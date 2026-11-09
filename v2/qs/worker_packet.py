"""Bounded current-collector packets, independent of process ownership."""

import hashlib
import json
import math
import time
from dataclasses import replace

from ..common import utils
from ..execution.budget import Budget, BudgetExhaustedError
from ..execution.resources import AllocationFailureError
from ..execution.work_budget import PollingBudget
from .factor_base import FactorBase
from .polynomial import Polynomial
from .relations import AtomicRelation, verify_atomic
from .residual_relations import ResidualAtom
from .sieve_collector import SieveCollector
from .sign_schedule import schedule_for
from .siqs import make_family

VERSION = "a7a-current-packet-v3"


def canonical(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    )


def digest(value):
    return hashlib.sha256(canonical(value).encode("ascii")).hexdigest()


def checked_int(value, name, low=0, high=2**63 - 1):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(name + " outside packet bounds")
    return value


def verification_options(config):
    return dict(
        residual_bound=config.residual_bound,
        large_prime_bound=getattr(config, "large_prime_bound", 0),
        large_product_bound=getattr(config, "large_product_bound", 0),
        small_prime_bound=getattr(config, "small_prime_bound", None),
        residual_arity=getattr(config, "residual_arity", 2),
        opaque_pairs=getattr(config, "opaque_pairs", False),
    )


def polynomial_at(
    base, primes, ordinal, config, budget, *, coefficient_cache=None
):
    """Reconstruct coefficients without consuming worker root tables."""
    if coefficient_cache is not None:
        return coefficient_cache.polynomial_at(
            base, primes, ordinal, config, budget
        )
    a = math.prod(primes)
    parity = config.polynomial_policy == "parity2" and base.n_prime % 4 == 1
    if parity:
        a *= 2
    roots = {
        e.prime: e.square_roots[0] for e in base.entries if e.prime in primes
    }
    schedule = schedule_for(
        len(primes),
        "swap_high"
        if config.sign_policy == "swap_high" and len(primes) >= 3
        else "control",
    )
    mask = schedule.mask(ordinal)
    b = a // 2 if parity else 0
    for index, prime in enumerate(primes):
        budget.consume(a.bit_length() + prime.bit_length() ** 2)
        q = a // prime
        term = q * pow(q, -1, prime) * roots[prime] % a
        b += -term if mask & (1 << index) else term
    b = (b + a // 2) % a - a // 2
    return Polynomial(base.n, base.multiplier, a, b)


def encode_atom(atom, ordinal):
    return [
        ordinal,
        atom.position,
        atom.sign,
        list(atom.exponents),
        atom.residual,
        list(atom.large_primes),
        getattr(atom, "opaque", None),
        isinstance(atom, ResidualAtom),
    ]


def decode_atom(record, polynomial, config, base, budget, *, verify=True):
    if not isinstance(record, list) or len(record) != 8:
        raise ValueError("invalid atom record")
    _, position, sign, exponents, residual, primes, opaque, extended = record
    if (
        not isinstance(exponents, list)
        or len(exponents) > len(base.entries)
        or not isinstance(primes, list)
        or len(primes) > 3
        or type(extended) is not bool
    ):
        raise ValueError("invalid atom dimensions")
    factors = tuple(tuple(pair) for pair in exponents)
    if extended:
        if residual != 1:
            raise ValueError("extended atom has a scalar residual")
        atom = ResidualAtom(
            polynomial,
            position,
            sign,
            factors,
            primes=tuple(primes),
            opaque=opaque,
        )
    else:
        if opaque is not None:
            raise ValueError("ordinary atom has an opaque residual")
        atom = AtomicRelation(
            polynomial,
            position,
            sign,
            factors,
            residual,
            large_primes=tuple(primes),
        )
    if verify:
        verify_atomic(
            atom, base, budget=budget, **verification_options(config)
        )
    return atom


def collect_packet(task, *, cancelled=None, publish_cpu=None):
    """Compute one leased prefix. A compatible service can call this directly.

    No relation/partial graph is reduced here. Recovery retains current roots,
    queues, score/hit metadata and private certified endpoint hints. Original
    atoms are exported; graph selection and retention belong to the parent.
    """
    config, settings = task["config"], task["settings"]
    started, cpu_started = time.monotonic(), time.process_time()
    # Root ownership includes final encoding after collection stops; the
    # private counter records only the reported mathematical prefix.
    raw_budget = Budget(
        seconds=task["seconds"],
        cpu_seconds=None,
        rss_bytes=None,
        cancelled=cancelled,
    )
    if settings.lease_budget:
        from .packet_lease import PacketLeaseBudget

        budget = PacketLeaseBudget(raw_budget, poll=publish_cpu)
    else:
        budget = PollingBudget(raw_budget, interval=64, poll=publish_cpu)
    rows, size, scanned, peak = [], 2, 0, 0
    if settings.packed_transport:
        from .worker_wire import PackedRows, legacy_length

        rows = task.get("_transport")
        if rows is None:
            rows = PackedRows(settings.packet_bytes)
        rows.reset()
    cursor = tuple(task["cursor"])
    divisor, reason = None, "complete"
    collector = None
    runtime = task.get("_runtime")
    family_reuse = collector_reuse = 0
    ordinal = cursor[0]
    try:
        budget.consume(0)
        record = task["base"]
        budget.consume(len(record[3]) * record[2].bit_length() ** 2)
        base = runtime.base if runtime is not None else FactorBase(*record)
        family = None
        if (
            runtime is not None
            and runtime.primes == tuple(task["primes"])
            and runtime.family is not None
            and runtime.family.next_index == ordinal
        ):
            family = runtime.family
            family.budget = budget
            family_reuse = 1
        if family is None:
            if runtime is not None:
                runtime.family = None
            family = make_family(
                base,
                tuple(task["primes"]),
                config=config,
                budget=budget,
            )
            family.next_index = ordinal
        if runtime is not None:
            runtime.family, runtime.primes = family, tuple(task["primes"])
        lo, hi = -config.half_width, config.half_width + 1
        cache_extra = (
            (runtime.bytes - base.workspace_bytes)
            if runtime is not None
            else 0
        )
        options = replace(config.collector)
        if task["split_grant"] is not None:
            options = replace(options, split_call_limit=task["split_grant"])
        if (
            runtime is not None
            and runtime.collector is not None
            and runtime.collector.config == options
        ):
            collector = runtime.collector
            collector.budget = budget
            runtime.reset_recovery(collector)
            collector_reuse = 1
        elif runtime is not None:
            # Incompatible capacity/lease settings cannot keep an old scratch
            # allocation alive while a replacement collector is constructed.
            runtime.collector = runtime.recovery = None

        def admit(atom, stats):
            nonlocal size
            if len(rows) >= settings.max_batch_atoms:
                return "packet_limit"
            row = encode_atom(atom, ordinal)
            length = (
                legacy_length(row)
                if settings.packed_transport
                else len(canonical(row))
            )
            if size + length + 1 > settings.packet_bytes:
                return "packet_limit"
            if settings.verification == "worker":
                if collector._checked is not None:
                    collector._checked.atom(atom, budget)
                else:
                    verify_atomic(
                        atom,
                        base,
                        budget=budget,
                        **verification_options(config.collector),
                    )
            budget.consume(length + len(atom.exponents) + 1)
            rows.append(row)
            size += length + 1
            stats["admitted_atoms"] += 1
            return None

        if collector is not None:
            collector._admit = admit
            if settings.grouped_lease_charging:
                from types import MethodType

                from .sieve_collector import SieveCollector as GenericCollector

                collector._sieve_direct = MethodType(
                    GenericCollector._sieve_direct, collector
                )
            collector._group_lease = (
                budget if settings.grouped_lease_charging else None
            )

        recovery_bound = False
        while ordinal < task["end"]:
            if (
                settings.whole_family
                and collector is not None
                and ordinal > task["cursor"][0]
                and ordinal % settings.polynomials_per_task == 0
            ):
                from .worker_state import reset_packet_hints

                # Whole-family scheduling must not extend packet-dependent
                # recovery hints or reset the packet's funded split counter.
                reset_packet_hints(collector)
            step = family.next()
            if collector is None:
                collector = SieveCollector(
                    step.polynomial,
                    base,
                    config=options,
                    budget=budget,
                    precomputed_roots=step.roots,
                )
                # Use the actual current collector, including DLP/extended
                # recovery. Only publication is replaced; no subclass bypass.
                collector._admit = admit
                if settings.grouped_lease_charging:
                    from types import MethodType

                    from .sieve_collector import (
                        SieveCollector as GenericCollector,
                    )

                    collector._sieve_direct = MethodType(
                        GenericCollector._sieve_direct, collector
                    )
                collector._group_lease = (
                    budget if settings.grouped_lease_charging else None
                )
                if runtime is not None:
                    runtime.collector = collector
                    runtime.recovery = getattr(
                        collector, "_divide_double", None
                    )
            else:
                collector.set_polynomial(step.polynomial, step.roots)
            if (
                not recovery_bound
                and collector._graph is not None
                and task["split_soft"]
            ):
                recover = (
                    runtime.recovery
                    if runtime is not None
                    else collector._divide_double
                )

                def leased_recovery(sieve, *args, **kwargs):
                    if sieve._split_calls >= task["split_grant"]:
                        budget.reason = "split_lease_limit"
                        raise BudgetExhaustedError(budget.reason)
                    return recover(sieve, *args, **kwargs)

                collector._divide_double = leased_recovery
                recovery_bound = True
            collector._set_plan_interval(lo, hi)
            position = cursor[1] if cursor[0] == ordinal else lo
            while position < hi:
                end = min(hi, position + config.batch_width)
                result = collector._collect_progress(position, end)
                scanned += result.stats["scanned"]
                position = result.next_position
                cursor = (ordinal, position)
                peak = max(
                    peak,
                    result.workspace_bytes
                    + family.workspace_bytes
                    - base.workspace_bytes,
                    result.workspace_bytes
                    + family.workspace_bytes
                    - base.workspace_bytes
                    + cache_extra,
                )
                if result.divisor is not None:
                    if not utils.valid_divisor(result.divisor, base.n):
                        raise ValueError("invalid packet divisor")
                    divisor = result.divisor
                    if task["fixed_work"]:
                        # The consumed divisor position has no unit atom.
                        continue
                if result.reason != "complete":
                    reason = result.reason
                    break
            if position != hi:
                break
            ordinal += 1
            cursor = (ordinal, lo)
        if ordinal == task["end"]:
            reason = "complete"
    except BudgetExhaustedError:
        reason = budget.reason
    except MemoryError as error:
        raise AllocationFailureError(
            "current SIQS packet", "qs_worker"
        ) from error
    publication_stop = None
    if settings.lease_budget:
        budget.reconcile()
        try:
            budget.consume(0)
        except BudgetExhaustedError:
            # The paid mathematical prefix remains publishable. Preserve the
            # stop at this boundary without relabeling a finished traversal.
            publication_stop = budget.reason
    if publish_cpu is not None:
        publish_cpu(raw_budget)
    packet = (
        rows.finish()
        if settings.packed_transport
        else canonical(rows).encode("ascii")
    )
    if len(packet) > settings.packet_bytes:
        raise ArithmeticError("packet reservation failed")
    report = dict(
        version=VERSION,
        kind="qs_collect",
        job=task["job"],
        id=task["id"],
        token=task["token"],
        attempt=task["attempt"],
        start=task["cursor"],
        cursor=list(cursor),
        reason=reason,
        divisor=divisor,
        work=budget.used,
        scanned=scanned,
        split_calls=0 if collector is None else collector._split_calls,
        packet=packet,
        sha256=hashlib.sha256(packet).hexdigest(),
        workspace_bytes=peak,
        cpu_seconds=time.process_time() - cpu_started,
        seconds=time.monotonic() - started,
        finished_at=time.monotonic(),
    )
    if settings.lease_budget:
        report.update(
            lease_policy="packet-parent-v1",
            publication_stop=publication_stop,
        )
    if settings.packed_transport:
        report["encoding"] = "packed-v1"
    if settings.grouped_lease_charging:
        report["grouped_marking_pairs"] = budget.marking_pairs
    if settings.reuse_state:
        from .worker_state import context_key

        report["context_key"] = (
            runtime.key
            if runtime is not None
            else context_key(
                dict(
                    config=config,
                    settings=settings,
                    base=task["base"],
                    job=task["job"],
                    execution=digest([task["job"], settings.record()]),
                )
            )
        )
        report["state_stats"] = dict(
            context_hits=int(runtime is not None),
            family_reuses=family_reuse,
            collector_reuses=collector_reuse,
        )
        if collector is not None and runtime is not None:
            # Rows are already encoded; do not keep the publication closure,
            # its old lease or its temporary row objects in the warm context.
            collector._admit = None
            collector._group_lease = None
    return report


def _process_cpu_reader(state):
    """Read one synchronized CPU snapshot against this run's fixed baseline."""
    # Pool ownership and quiescent finish prevent begin() from replacing the
    # baseline while a packet is active. Keep that immutable local copy; lock
    # the live array once per check rather than once for every scalar access.
    baseline = tuple(state["cpu_base"])
    counters = state["cpu"]
    values, lock = counters.get_obj(), counters.get_lock()

    def used():
        with lock:
            return sum(
                max(0, values[index] - prior)
                for index, prior in enumerate(baseline)
            )

    return used


def execute_packet(task, state=None):
    """Adapt the computational packet to the retained A7 executor."""
    from . import parallel

    if "context_key" in task:
        from .worker_state import expand_task

        task = expand_task(task, state)

    process = state is None
    slot = None
    if process:
        state, slot = parallel._PROCESS_STATE
    if getattr(task.get("settings"), "packed_transport", False):
        from .worker_wire import writer_for

        task = dict(
            task, _transport=writer_for(state, task["settings"].packet_bytes)
        )
    leased = getattr(task.get("settings"), "lease_budget", False)
    cpu_used = _process_cpu_reader(state) if process and not leased else None
    stop = state["stop"]
    if getattr(task.get("settings"), "private_stop", False):
        stop = stop.worker_event(slot)
    last_poll = 0.0

    def publish(budget):
        nonlocal last_poll
        now = time.monotonic()
        if process and now - last_poll >= 0.001:
            state["cpu"][slot] = time.process_time()
            state["rss"][slot] = parallel._rss()
            last_poll = now

    def cancelled():
        if stop.is_set():
            return True
        if leased:
            return False
        external = cpu_used() if process else 0
        limit = state["cpu_limit"].value
        return limit >= 0 and state["parent_cpu"].value + external >= limit

    try:
        return collect_packet(
            task,
            cancelled=cancelled,
            publish_cpu=None if leased else publish,
        )
    finally:
        if process:
            state["cpu"][slot] = time.process_time()
            state["rss"][slot] = parallel._rss()
