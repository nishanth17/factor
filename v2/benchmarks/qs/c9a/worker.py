"""One bounded PyPy worker, warmed and drift-checked in the same process."""

import gc
import hashlib
import platform
import statistics
import time
from contextlib import nullcontext
from dataclasses import asdict
from pathlib import Path

from .... import qs as candidate
from ....execution import budget as candidate_budget
from ...support.resources import PeakRSSReader
from .configuration import decode
from .control import control_modules
from .journal import Journal, RecorderExhaustedError, canonical, observe

VALIDATION_WORK = 10**7
WARM_N = 1000000007 * 1000000087
_RSS_READER = PeakRSSReader()


def runtime():
    if platform.python_implementation() != "PyPy" or tuple(
        platform.python_version_tuple()[:2]
    ) != ("3", "11"):
        raise RuntimeError("C9a requires PyPy Python 3.11")


def rss():
    return _RSS_READER()


def checked_split(n, divisor):
    if divisor is not None and not (
        type(divisor) is int and 1 < divisor < n and n % divisor == 0
    ):
        raise ArithmeticError("worker returned an improper divisor")
    return [divisor, n // divisor] if divisor is not None else []


def warm_once(qs, budgets, template, deadline):
    """Small validated same-path work; not a large-input feasibility sample."""
    values = dict(
        large_prime_bound=10000,
        large_product_bound=100000000,
        residual_bound=10000,
        block_width=128,
        division="bucket",
        score_policy="powers",
        memory_bytes=64 * 2**20,
        sparse_queue_bytes=2**20,
        power_plan_bytes=2**20,
    )
    selected = template["collector"]
    if "yield_policy" in selected:
        for key in (
            "preparation_policy",
            "cycle_policy",
            "cycle_search_edges",
            "duplicate_policy",
            "duplicate_bytes",
            "record_policy",
        ):
            values[key] = selected[key]
        values["small_prime_bound"] = (
            1000 if selected["small_prime_bound"] else None
        )
        if selected["yield_policy"] == "yield-v2":
            for key in (
                "recovery_policy",
                "split_policy",
                "batch_policy",
                "batch_size",
                "checked_bytes",
                "filter_policy",
                "extraction_policy",
            ):
                values[key] = selected[key]
            if selected.get("extension_policy"):
                from ....qs.extended_policy import (
                    BOUNDED_FIELDS,
                    EXTENSION_FIELDS,
                )

                values.update({key: selected[key] for key in EXTENSION_FIELDS})
                if selected["extension_policy"] == "extended-v2":
                    values.update(
                        {key: selected[key] for key in BOUNDED_FIELDS}
                    )
                    collector = qs.BoundedYieldSieveConfig(**values)
                else:
                    collector = qs.ExtendedYieldSieveConfig(**values)
            else:
                collector = qs.CombinedYieldSieveConfig(**values)
        else:
            collector = qs.YieldSieveConfig(**values)
    else:
        collector = qs.GrowingDoubleLargeSieveConfig(**values)
    arithmetic = {
        key: template[key]
        for key in (
            "polynomial_policy",
            "multiplier_policy",
            "root_update_policy",
        )
        if template.get(key, "standard") != "standard"
    }
    if arithmetic:
        arithmetic["assignment_policy"] = template["assignment_policy"]
        arithmetic["multiplier"] = 0
        if "parity_two" in selected["marking"].split("+"):
            collector = type(collector)(
                **dict(asdict(collector), marking="parity_two")
            )
    config = qs.SIQSConfig(
        base_bound=200,
        half_width=128,
        max_half_width=128,
        memory_bytes=64 * 2**20,
        collector=collector,
        sign_policy="swap_high",
        solve_policy="collection",
        **arithmetic,
    )
    budget = budgets.Budget(
        work_limit=10**9,
        seconds=max(0, deadline - time.monotonic()),
        cpu_seconds=5,
    )
    job = qs.SIQSJob(WARM_N, seed=7, config=config, budget=budget)
    try:
        return _validated_warm_job(qs, job, selected, values, budget)
    finally:
        # Warmup publishes no checkpoints. Explicitly retire its owned cold
        # files so repeated sentinel calls cannot accumulate orphaned stores.
        engine = job.engine
        arena = (
            None
            if engine is None
            else getattr(engine.collector, "_arena", None)
        )
        if arena is not None and getattr(arena, "spill", None) is not None:
            arena.spill.cleanup()


def _validated_warm_job(qs, job, selected, values, budget):
    result = job.run(max_blocks=2)
    checked_split(WARM_N, result.divisor)
    if result.reason not in ("paused", "factor_found"):
        raise RuntimeError("warmup did not complete its fixed validated work")
    atoms = (
        ()
        if job.engine is None
        else tuple(job.engine.collector._atoms.values())
    )
    for atom in atoms:
        options = dict(
            residual_bound=10000,
            large_prime_bound=10000,
            large_product_bound=100000000,
            budget=budget,
        )
        if "yield_policy" in selected:
            options["small_prime_bound"] = values["small_prime_bound"]
        if selected.get("extension_policy"):
            options.update(
                residual_arity=selected["residual_arity"],
                opaque_pairs=selected["opaque_pairs"],
            )
        qs.verify_atomic(atom, job.base, **options)
    signature = hashlib.sha256(
        canonical(
            dict(
                divisor=result.divisor,
                atoms=[
                    dict(
                        polynomial=asdict(a.polynomial),
                        position=a.position,
                        sign=a.sign,
                        exponents=tuple(a.exponents),
                        residual=a.residual,
                        large_primes=a.large_primes,
                        relation_id=a.relation_id,
                    )
                    for a in atoms
                ],
                position=result.next_position,
            )
        )
    ).hexdigest()
    return signature, budget.used


def sentinel(qs, budgets, template, deadline, expected=None):
    wall, cpu, work = [], [], 0
    for _ in range(9):
        start, cpu_start = time.monotonic(), time.process_time()
        for _ in range(4):
            signature, used = warm_once(qs, budgets, template, deadline)
            if expected is not None and signature != expected:
                raise ArithmeticError(
                    "warm/sentinel arithmetic signature changed"
                )
            expected = signature
            work += used
        wall.append(time.monotonic() - start)
        cpu.append(time.process_time() - cpu_start)
    return dict(wall=wall, cpu=cpu, work=work, signature=expected)


def drift(before, after):
    def spread(values):
        ordered = sorted(values)
        return (ordered[6] - ordered[2]) / statistics.median(ordered)

    ratios = {
        key: statistics.median(after[key]) / statistics.median(before[key])
        for key in ("wall", "cpu")
    }
    stable = all(
        spread(sample[key]) <= 0.15
        for sample in (before, after)
        for key in ("wall", "cpu")
    )
    return dict(
        ratios=ratios,
        stable=stable,
        valid=stable and all(0.85 <= r <= 1.15 for r in ratios.values()),
    )


def run(cell, *, capture_path=None):
    runtime()
    timing_policy = cell.get("timing_policy", "legacy")
    if timing_policy not in ("legacy", "retire-rewarm-v1"):
        raise ValueError("unknown worker timing policy")
    launch_wall, launch_cpu = time.monotonic(), time.process_time()
    qs, budgets = (
        control_modules()
        if cell["source"] == "immutable"
        else (candidate, candidate_budget)
    )
    warm_start, warm_cpu = time.monotonic(), time.process_time()
    warm_deadline = warm_start + cell["warm_seconds"] - cell["end_reserve"]
    warm_work, signature = 0, None
    while (
        min(time.monotonic() - warm_start, time.process_time() - warm_cpu) < 3
    ):
        if time.monotonic() >= warm_deadline:
            raise RuntimeError("validated warmup reservation exhausted")
        found, used = warm_once(
            qs, budgets, cell["configuration"], warm_deadline
        )
        if signature is not None and found != signature:
            raise ArithmeticError("warmup signature changed")
        signature, warm_work = found, warm_work + used
        if warm_work > cell["warm_work_limit"]:
            raise RuntimeError("warmup work reservation exhausted")
    before = sentinel(
        qs, budgets, cell["configuration"], warm_deadline, signature
    )
    warm_work += before["work"]
    warm_before_wall = time.monotonic() - warm_start
    warm_before_cpu = time.process_time() - warm_cpu
    if (
        warm_before_wall > cell["warm_seconds"] - cell["end_reserve"]
        or warm_before_cpu > cell["warm_seconds"] - cell["end_reserve"]
        or warm_work > cell["warm_work_limit"]
    ):
        raise RuntimeError("end sentinel reservation would be consumed")

    start, cpu_start = time.monotonic(), time.process_time()
    exceeded = False

    def cancelled():
        nonlocal exceeded
        exceeded = exceeded or rss() > cell["rss_bytes"]
        return exceeded

    budget = budgets.Budget(
        work_limit=cell["work_limit"] - VALIDATION_WORK,
        seconds=cell["seconds"],
        cpu_seconds=cell["seconds"],
        cancelled=cancelled,
    )
    writer = None
    job, result, reason, divisor, stats = None, None, None, None, {}
    try:
        config = decode(cell["configuration"], qs)
        if capture_path is not None:
            writer = Journal(
                capture_path,
                dict(
                    cell=cell,
                    format="c9a-atoms-v2"
                    if getattr(config.collector, "extension_policy", None)
                    else "c9a-atoms-v1",
                ),
                byte_limit=cell.get("observer_byte_limit", 256 * 2**20),
                work_limit=cell["observer_work_limit"],
            )
        with observe(qs, writer) if writer is not None else nullcontext():
            job = qs.SIQSJob(
                cell["n"], seed=cell["seed"], config=config, budget=budget
            )
            result = job.run(max_blocks=cell["max_blocks"])
        reason, divisor, stats = result.reason, result.divisor, result.stats
    except RecorderExhaustedError:
        reason = "recorder_limit"
    except MemoryError:
        reason = "native_memory_limit"
    finally:
        if writer is not None:
            writer.close(
                dict(reason=reason, divisor=divisor),
                complete=result is not None,
            )
    factors = checked_split(cell["n"], divisor)
    if factors and factors[0] * factors[1] != cell["n"]:
        raise ArithmeticError("worker reconstruction failed")
    seconds, cpu_seconds = (
        time.monotonic() - start,
        time.process_time() - cpu_start,
    )
    measured_work = budget.used + 4 * cell["n"].bit_length()
    observer_work = 0 if writer is None else writer.work
    if exceeded:
        reason = "process_rss_limit"
    actual_multiplier = None if job is None else job.multiplier
    actual_base_bound = (
        None if job is None or job.base is None else job.base.bound
    )
    next_position = None if result is None else result.next_position
    production_rss = rss()
    end_start, end_cpu = time.monotonic(), time.process_time()
    after = sentinel(
        qs,
        budgets,
        cell["configuration"],
        end_start + cell["end_reserve"],
        signature,
    )
    warm_work += after["work"]
    legacy_after = after
    retirement = None
    if timing_policy == "retire-rewarm-v1":
        # Preserve the immediate old screen, then measure the sentinel in a
        # separately declared, comparable post-retirement warm state. The
        # factoring interval above excludes all of this diagnostic overhead.
        cleanup_start, cleanup_cpu = time.monotonic(), time.process_time()
        engine = None if job is None else job.engine
        arena = (
            None
            if engine is None
            else getattr(engine.collector, "_arena", None)
        )
        if arena is not None and getattr(arena, "spill", None) is not None:
            arena.spill.cleanup()
        job = result = engine = arena = None
        gc.collect()
        retirement = dict(
            wall=time.monotonic() - cleanup_start,
            cpu=time.process_time() - cleanup_cpu,
        )
        rewarm_start, rewarm_cpu = time.monotonic(), time.process_time()
        end_deadline = end_start + cell["end_reserve"]
        while (
            min(
                time.monotonic() - rewarm_start,
                time.process_time() - rewarm_cpu,
            )
            < 3
        ):
            if (
                time.monotonic() >= end_deadline
                or warm_work >= cell["warm_work_limit"]
            ):
                raise RuntimeError("post-retirement warmup exhausted")
            found, used = warm_once(
                qs, budgets, cell["configuration"], end_deadline
            )
            if found != signature:
                raise ArithmeticError("post-retirement signature changed")
            warm_work += used
        retirement.update(
            rewarm_wall=time.monotonic() - rewarm_start,
            rewarm_cpu=time.process_time() - rewarm_cpu,
        )
        after = sentinel(
            qs, budgets, cell["configuration"], end_deadline, signature
        )
        warm_work += after["work"]
    if (
        time.monotonic() - end_start > cell["end_reserve"]
        or time.process_time() - end_cpu > cell["end_reserve"]
        or warm_work > cell["warm_work_limit"]
    ):
        raise RuntimeError("end sentinel reservation exhausted")
    screen = drift(before, after)
    return dict(
        fixture=cell["fixture"],
        arm=cell["arm"],
        n=cell["n"],
        seed=cell["seed"],
        digits=cell["digits"],
        split_found=divisor is not None,
        complete=False,
        divisor=divisor,
        factors=factors,
        remaining=[] if factors else [cell["n"]],
        reason=reason,
        seconds=seconds,
        cpu_seconds=cpu_seconds,
        measured_work=measured_work,
        warm_work=warm_work,
        observer_work=observer_work,
        work=measured_work + warm_work + observer_work,
        warm_wall=warm_before_wall + time.monotonic() - end_start,
        warm_cpu=warm_before_cpu + time.process_time() - end_cpu,
        process_wall=time.monotonic() - launch_wall,
        process_cpu=time.process_time() - launch_cpu,
        rss_bytes=rss(),
        production_rss_bytes=production_rss,
        rss_observer=_RSS_READER.identity,
        stats=stats,
        sentinel_before=before,
        sentinel_after=after,
        drift=screen,
        timing_policy=timing_policy,
        library_root=str(Path(qs.__file__).resolve().parents[1]),
        legacy_sentinel_after=legacy_after,
        legacy_drift=drift(before, legacy_after),
        retirement=retirement,
        timing_valid=screen["valid"] and reason != "recorder_limit",
        observation="capture" if writer else "clean",
        recorder=None
        if writer is None
        else dict(
            bytes=writer.bytes, work=writer.work, exhausted=writer.exhausted
        ),
        actual_multiplier=actual_multiplier,
        actual_base_bound=actual_base_bound,
        next_position=next_position,
    )
