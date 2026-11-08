"""Bounded portfolio with automatic ECM and checked owned-pool lifetimes."""

import hashlib
import json
import math
import random
import time
from dataclasses import asdict, dataclass, replace

from . import constants
from .common import arithmetic, prime_sieve, utils
from .common.arithmetic import isqrt
from .common.preprocessing import (
    fermat_step,
    integer_root,
    power_residue_possible,
    strip_twos,
)
from .ecm.programs import (
    PAIRED_VERSION,
    PROGRAM_VERSION,
    WHEEL_VERSION,
)
from .execution.allocation import (
    ALLOCATION_VERSION,
    ECMAllocation,
    HandoffRequiredError,
    PretestBudget,
    fallback_refusal,
)
from .execution.budget import Budget, BudgetExhaustedError
from .execution.portfolio_config import (
    _AUTO_SIQS as _AUTO_SIQS,
)
from .execution.portfolio_config import (
    PortfolioConfig as PortfolioConfig,
)
from .execution.portfolio_config import (
    _chain_bounds,
    _decode_siqs_config,
)
from .execution.portfolio_config import (
    _execution_context as _execution_context,
)
from .execution.resources import (
    decode_checkpoint,
    encode_checkpoint,
    removed_arguments,
)
from .execution.schedules import SieveContext
from .execution.stage_jobs import (
    advance_job,
    new_job,
    peek_prime,
    prime_cursor,
    take_prime,
)
from .factor import FactorizationResult, PrimeFactor
from .qs import SIQSConfig as SIQSConfig
from .qs import SIQSJob
from .qs.sss import SSSConfig, SSSJob

CHECKPOINT_VERSION = 10
OPTIONAL_CHAIN_CHECKPOINT_VERSION = 11
ALLOCATION_CHECKPOINT_VERSION = 12
AUTOMATIC_CHECKPOINT_VERSION = 13
SCHEDULE_VERSION = "half-open-prime-powers/ecm-even-baby-v1"


@dataclass(frozen=True)
class PortfolioRun:
    """Result, stop reason, consumed resources, trace, and resume payload."""

    result: FactorizationResult
    reason: str
    work_used: int
    wall_seconds: float
    cpu_seconds: float
    events: tuple
    dropped_events: int
    checkpoint: dict
    resource_stop: dict | None = None


def _canonical(value):
    """Encode exact decimal integers with stable keys for corruption checks."""
    return json.dumps(
        value,
        default=arithmetic.json_integer,
        sort_keys=True,
        separators=(",", ":"),
    )


def _tuples(value):
    """Restore Random's tuple-based state from JSON lists."""
    if isinstance(value, list):
        return tuple(_tuples(item) for item in value)
    return value


def _result(state):
    """Include current, pending, and exhausted cofactors in reconstruction."""
    remaining = list(state["remaining"])
    outstanding = list(state["pending"])
    if state["current"] is not None:
        outstanding.append([state["current"]["n"], state["current"]["mult"]])
    for n, multiplicity in outstanding:
        remaining.extend([n] * multiplicity)
    factors = tuple(
        PrimeFactor(int(value), record[0], utils.Primality(record[1]))
        for value, record in sorted(
            state["factors"].items(), key=lambda item: int(item[0])
        )
    )
    return FactorizationResult(
        state["original"],
        -1 if state["original"] < 0 else 1,
        factors,
        tuple(sorted(remaining)),
    )


def _add_factor(state, n, multiplicity, certainty):
    """Merge exact multiplicities without upgrading primality certainty."""
    key = str(n)
    if key in state["factors"]:
        state["factors"][key][0] += multiplicity
    else:
        state["factors"][key] = [multiplicity, certainty.value]


def _event(state, config, **record):
    """Keep a bounded event history; report omitted event counts explicitly."""
    if len(state["events"]) < config.trace_limit:
        state["events"].append(record)
    else:
        state["dropped_events"] += 1


def _split(state, divisor):
    """Replace a parent only after validating its exact divisor."""
    current = state["current"]
    n, multiplicity = current["n"], current["mult"]
    if not utils.valid_divisor(divisor, n):
        raise ValueError("portfolio split is not a proper divisor")
    state["pending"].extend(
        (
            [arithmetic.backend_for(n).integer(divisor), multiplicity],
            [arithmetic.divexact(n, divisor), multiplicity],
        )
    )
    state["current"] = None


def _classification_bases(n, policy):
    """Keep legacy random jobs under their original certainty/RNG contract."""
    if policy == utils.LEGACY_PRIMALITY_POLICY and (
        n >= utils.WORD_DETERMINISTIC_LIMIT
    ):
        return None
    bases = utils.deterministic_bases(n)
    return list(bases) if bases is not None else None


def _classify_step(current, config, budget, generator, policy):
    """Run at most one primality witness, preserving RNG position on pause."""
    n = current["n"]
    witness_state = current.get("prime_job")
    if witness_state is None:
        budget.consume(len(utils.SMALL_PRIMES))
        for prime in utils.SMALL_PRIMES:
            if n == prime:
                return utils.Primality.PROVEN.value
            if n % prime == 0:
                return utils.Primality.COMPOSITE.value

        if n < 41 * 41:
            return utils.Primality.PROVEN.value
        shifts = ((n - 1) & -(n - 1)).bit_length() - 1
        bases = _classification_bases(n, policy)

        current["prime_job"] = {
            "d": (n - 1) >> shifts,
            "s": shifts,
            "index": 0,
            "bases": bases,
            "tested": [],
        }
        return None

    bases = witness_state["bases"]
    index = witness_state["index"]
    rounds = len(bases) if bases else config.primality_rounds
    # Refusal precedes the random draw so resume sees the same witness.
    budget.consume(n.bit_length() + witness_state["s"])
    base = bases[index] if bases else generator.randrange(2, int(n) - 1)
    survives = utils._strong_probable_prime(
        n, base, witness_state["d"], witness_state["s"]
    )
    witness_state["index"] += 1
    if not survives:
        return utils.Primality.COMPOSITE.value
    witness_state["tested"].append(base)
    if witness_state["index"] == rounds:
        return (
            utils.Primality.PROVEN.value
            if bases
            else utils.Primality.PROBABLE.value
        )
    return None


def _handoff(state, config, budget, reason):
    """Commit handoff with spent work and partial-assignment evidence."""
    current = state["current"]
    job = current.get("job")
    if config.allocation is not None and current["stage"] == "trial":
        # Relation engines require power preprocessing even when optional
        # trial work cannot be admitted. Its exact roots use the shared ledger.
        current.update(stage="powers", pending_handoff=reason, job=None)
        return
    current.pop("pending_handoff", None)
    if config.allocation is not None:
        if job is not None:
            cursor = job.get("cursor")
            _event(
                state,
                config,
                stage=job["kind"],
                n=current["n"],
                seed=job["seed"],
                b1=job["b1"],
                b2=job["b2"],
                tier=current["tier"],
                attempt=current["attempt"],
                work=budget.used - job["start_work"],
                outcome="handoff",
                phase=job["phase"],
                cursor={key: cursor[key] for key in ("left", "next", "index")}
                if cursor is not None
                else None,
            )
        _event(
            state,
            config,
            stage="handoff",
            n=current["n"],
            outcome=reason,
            from_stage=current["stage"],
            work_used=budget.used,
            remaining_work=None,
        )
    current["job"] = None
    if config.siqs is not None:
        current.update(stage="siqs", siqs_start_work=budget.used)
    elif config.sss is not None:
        current.update(stage="sss", sss_start_work=budget.used)
    else:
        state["remaining"].extend([current["n"]] * current["mult"])
        state["current"] = None
        if config.allocation is not None and reason != "schedule_exhausted":
            state["stop_reason"] = "pretest_exhausted"


def _admit_fallback(current, config, budget):
    """Do not require a fresh full reservation for an already admitted job."""
    if config.allocation is None or current.get("fallback_admitted"):
        return
    reason = fallback_refusal(budget, config.allocation)
    if reason is not None:
        raise BudgetExhaustedError(reason)
    current["fallback_admitted"] = True


def _advance(
    state,
    config,
    budget,
    context,
    generator,
    siqs_runtime,
    programs,
    policy,
    ecm_executor=None,
):
    """Commit one portfolio transition or one resumable candidate action."""
    if state["current"] is None:
        if not state["pending"]:
            return
        budget.consume()
        n, multiplicity = state["pending"].pop()
        state["current"] = {
            "n": n,
            "mult": multiplicity,
            "stage": "twos",
            "attempt": 0,
            "tier": 0,
            "job": None,
            "fermat": 0,
            "power_index": 0,
        }

    current = state["current"]
    n = current["n"]
    if current["stage"] == "twos":
        budget.consume(n.bit_length())
        odd, exponent = strip_twos(n)
        if exponent:
            _add_factor(
                state, 2, exponent * current["mult"], utils.Primality.PROVEN
            )
        current.update(n=odd, stage="classify")
        if odd == 1:
            state["current"] = None
        return

    if current["stage"] == "classify":
        key = str(n)
        certainty = state["classifications"].get(key)
        if certainty is None:
            certainty = _classify_step(
                current, config, budget, generator, policy
            )
            if certainty is None:
                return
            state["classifications"][key] = certainty

        current.pop("prime_job", None)
        if certainty != utils.Primality.COMPOSITE.value:
            _add_factor(state, n, current["mult"], utils.Primality(certainty))
            state["current"] = None
        else:
            current.update(
                stage="trial", cursor=prime_cursor(3, config.trial_bound + 1)
            )
        return

    if current["stage"] == "trial":
        prime = peek_prime(current["cursor"], context, budget)
        if prime is None:
            current["stage"] = "powers"
            return
        cursor = current["cursor"]
        count = min(
            config.trial_chunk, len(cursor["values"]) - cursor["index"]
        )
        # Reserve a fixed chunk allowance, including an early successful
        # division. Resume therefore consumes exactly the same reservations.
        budget.consume(count * n.bit_length())
        root = isqrt(n)

        for _ in range(count):
            prime = cursor["values"][cursor["index"]]
            if prime > root:
                _add_factor(state, n, current["mult"], utils.Primality.PROVEN)
                state["current"] = None
                return
            exponent = 0
            while n % prime == 0:
                n = arithmetic.divexact(n, prime)
                exponent += 1
            take_prime(cursor)
            if exponent:
                _add_factor(
                    state,
                    prime,
                    exponent * current["mult"],
                    utils.Primality.PROVEN,
                )
                current.update(n=n, stage="classify")
                if n == 1:
                    state["current"] = None
                return

        return

    if current["stage"] == "powers":
        if "power_exponents" not in current:
            budget.consume(n.bit_length())
            current["power_exponents"] = prime_sieve.small_sieve(
                n.bit_length()
            )

        exponents = current["power_exponents"]
        position = current["power_index"]
        if position == len(exponents):
            current["stage"] = "fermat"
            return
        exponent = exponents[position]
        budget.consume(n.bit_length())
        if not power_residue_possible(n, exponent):
            current["power_index"] += 1
            return
        base = integer_root(n, exponent)
        current["power_index"] += 1
        if base**exponent == n:
            # Folding the exponent into multiplicity preserves every later
            # split without expanding repeated pending work.
            state["pending"].append([base, current["mult"] * exponent])
            state["current"] = None
        return

    if current["stage"] == "fermat":
        if current["fermat"] >= config.fermat_steps:
            current["stage"] = "rho"
            return
        budget.consume(n.bit_length())
        root = isqrt(n)
        candidate = root + int(root * root < n) + current["fermat"]
        divisor = fermat_step(n, candidate)
        current["fermat"] += 1
        if divisor is not None:
            _split(state, divisor)
        return

    if current["stage"] == "sss":
        _admit_fallback(current, config, budget)
        job = siqs_runtime.get("job")
        if job is None:
            if current.get("sss_checkpoint") is not None:
                job = SSSJob.from_checkpoint(
                    current["sss_checkpoint"], budget=budget, config=config.sss
                )
            else:
                if "sss_seed" not in current:
                    budget.consume()
                    current["sss_seed"] = generator.getrandbits(63)
                job = SSSJob(
                    n,
                    seed=current["sss_seed"],
                    config=config.sss,
                    budget=budget,
                )

            if job.n != n or job.seed != current["sss_seed"]:
                raise ValueError("SSS checkpoint differs from its parent")
            siqs_runtime["job"] = job

        result = job.run(batch_limit=1)
        if result.divisor is not None or result.reason in (
            "search_exhausted",
            "relation_limit",
            "atom_limit",
            "memory_limit",
        ):
            _event(
                state,
                config,
                stage=config.sss.mode,
                n=n,
                seed=job.seed,
                work=budget.used - current["sss_start_work"],
                outcome="factor" if result.divisor else result.reason,
                stats=result.stats,
            )
            siqs_runtime.clear()
            if result.divisor is not None:
                _split(state, result.divisor)
            else:
                state["remaining"].extend([n] * current["mult"])
                state["current"] = None
        elif result.reason != "paused":
            raise BudgetExhaustedError(result.reason)
        return

    if current["stage"] == "siqs":
        _admit_fallback(current, config, budget)
        job = siqs_runtime.get("job")
        if job is None:
            if current.get("siqs_checkpoint") is not None:
                job = SIQSJob.from_checkpoint(
                    current["siqs_checkpoint"],
                    budget=budget,
                    config=config.siqs,
                )
            else:
                if "siqs_seed" not in current:
                    budget.consume()
                    current["siqs_seed"] = generator.getrandbits(63)
                job = SIQSJob(
                    n,
                    seed=current["siqs_seed"],
                    config=config.siqs,
                    budget=budget,
                )

            if job.n != n or job.seed != current["siqs_seed"]:
                raise ValueError(
                    "SIQS checkpoint differs from its parent assignment"
                )
            siqs_runtime["job"] = job

        job.budget = budget
        result = job.run(max_blocks=1)
        if result.divisor is not None or job.finished_reason is not None:
            _event(
                state,
                config,
                stage="siqs",
                n=n,
                seed=job.seed,
                work=budget.used - current["siqs_start_work"],
                outcome="factor" if result.divisor else result.reason,
                stats=result.stats,
            )
            siqs_runtime.clear()
            if result.divisor is not None:
                _split(state, result.divisor)
            else:
                state["remaining"].extend([n] * current["mult"])
                state["current"] = None
        elif result.reason != "paused":
            raise BudgetExhaustedError(result.reason)
        return

    kind = current["stage"]
    if (
        kind == "ecm"
        and ecm_executor is not None
        and not current.get("ecm_serial", False)
    ):
        from .execution.worker_contract import LeaseRefusedError

        if current["tier"] >= len(config.ecm_tiers):
            _handoff(state, config, budget, "schedule_exhausted")
            return
        run = getattr(ecm_executor, "run_portfolio", ecm_executor.run)
        if hasattr(ecm_executor, "run_portfolio") and (
            current["job"] is not None or current["attempt"] or current["tier"]
        ):
            raise ValueError(
                "selector cannot replace existing serial ECM progress"
            )
        try:
            outcome = run(
                int(n),
                config=config,
                budget=budget,
                lineage=f"{state['original']}/{n}/{current['mult']}",
                checkpoint=current.get("worker_checkpoint"),
            )
        except LeaseRefusedError as error:
            raise BudgetExhaustedError(str(error)) from error
        if outcome is None:
            # Decide once per cofactor. Warming a pool or resuming later cannot
            # replace an incumbent curve prefix or move its fallback RNG draws.
            current["ecm_serial"] = True
            _event(
                state,
                config,
                stage="ecm_selection",
                n=n,
                outcome="legacy_serial",
                stats=ecm_executor.last_stats,
            )
        else:
            current["worker_checkpoint"] = outcome.checkpoint
            _event(
                state,
                config,
                stage="ecm_workers",
                n=n,
                outcome=outcome.reason,
                stats=outcome.stats,
            )
            if outcome.divisor is not None:
                _split(state, outcome.divisor)
            elif outcome.reason == "curves_exhausted":
                current.update(tier=len(config.ecm_tiers), attempt=0)
                current.pop("worker_checkpoint", None)
            else:
                raise BudgetExhaustedError(outcome.reason)
            return

    if kind == "rho":
        attempts, b1, b2 = config.rho_attempts, 0, 0
    elif kind == "pm1":
        attempts, b1, b2 = config.pm1_attempts, config.pm1_b1, config.pm1_b2
    elif current["tier"] < len(config.ecm_tiers):
        b1, b2, attempts = config.ecm_tiers[current["tier"]]
    else:
        _handoff(state, config, budget, "schedule_exhausted")
        return

    if current["attempt"] >= attempts:
        current["attempt"] = 0
        if kind == "rho":
            current["stage"] = "pm1"
        elif kind == "pm1":
            current["stage"] = "ecm"
        else:
            current["tier"] += 1
        return

    if current["job"] is None:
        budget.consume()
        seed = (
            current["attempt"] if kind == "pm1" else generator.getrandbits(63)
        )
        current["job"] = new_job(kind, n, seed, b1, b2)
        current["job"]["start_work"] = budget.used

    job = current["job"]
    started = time.perf_counter()
    job_context = (
        programs if kind == "ecm" and programs is not None else context
    )
    advance_job(job, budget, job_context, config)
    phase = job["phase"]
    timing_key = f"{kind}:{phase}"
    state["stage_seconds"][timing_key] = (
        state["stage_seconds"].get(timing_key, 0)
        + time.perf_counter()
        - started
    )
    if job["done"]:
        _event(
            state,
            config,
            stage=kind,
            n=n,
            seed=job["seed"],
            b1=b1,
            b2=b2,
            sigma=job.get("sigma"),
            base=job.get("base"),
            work=budget.used - job["start_work"],
            outcome="factor" if job["factor"] else "exhausted",
        )
        if job["factor"] is not None:
            _split(state, job["factor"])
        else:
            current.update(job=None, attempt=current["attempt"] + 1)


def _chain_job_supported(job, config):
    if job["kind"] != "ecm" or job["b1"] not in _chain_bounds(config):
        return False
    from .ecm.chains import supports_modulus

    return supports_modulus(job["n"])


def _chain_policy(config):
    if config.ecm_chain_mode == "off":
        return None
    from .ecm.chains import CHAIN_VERSION

    if config.ecm_chain_family != "auto":
        return _chain_identity(2000, config) + "/reuse8/chunk16/40-80digits-v1"
    return CHAIN_VERSION + "/reuse8/chunk16/40-80digits-v1"


def _chain_identity(bound, config):
    """Leave accepted default identities intact; pin explicit alternatives."""
    if config.ecm_chain_family == "auto":
        from .ecm.chains import identity

        return identity(bound, config.backend)
    from .ecm.chain_options import identity

    return identity(bound, config.backend, config.ecm_chain_family)


def _legacy_pm1(config):
    """Old snapshots retain their original chunk and gap execution."""
    return not config.pm1_attempts or (
        config.pm1_gap_mode == "cached" and config.pm1_chunk_size is None
    )


@encode_checkpoint
def _pack(state, config, budget, generator, policy):
    """Snapshot RNG, pending work, schedule identity, and resource use."""
    wall_used, cpu_used = budget.wall_used, budget.cpu_used
    configuration = asdict(config)
    if config.allocation is None:
        configuration.pop("allocation")
    if config.ecm_chain_family == "auto" or config.ecm_chain_mode == "off":
        configuration.pop("ecm_chain_family")
    version = 6
    if (
        _legacy_pm1(config)
        and config.ecm_chain_mode == "off"
        and config.backend == "python-int"
        and config.ecm_pair_distance is None
        and config.ecm_pair_wheel is None
    ):
        # Preserve mainline's native schemas. Version 5 was independently
        # used by P4.3 and P5.2; version 6 combines GMP and program identity.
        version = 5 if config.ecm_program_bytes else 4
        configuration.pop("backend")
        for name in ("siqs", "sss"):
            if configuration[name] is not None:
                configuration[name].pop("backend")
    if not config.ecm_program_bytes:
        configuration.pop("ecm_program_bytes")
    if config.ecm_pair_distance is None:
        configuration.pop("ecm_pair_distance")
    else:
        version = 7
    if config.ecm_pair_wheel is None:
        configuration.pop("ecm_pair_wheel")
    else:
        version = 8
    if _legacy_pm1(config) and config.ecm_chain_mode == "off":
        configuration.pop("pm1_gap_mode")
        configuration.pop("pm1_chunk_size")
    else:
        version = 9
    if config.ecm_chain_mode == "off":
        configuration.pop("ecm_chain_mode")
        configuration.pop("ecm_chain_bytes")
    else:
        version = CHECKPOINT_VERSION
        if config.ecm_chain_family != "auto":
            version = OPTIONAL_CHAIN_CHECKPOINT_VERSION
    if config.allocation is not None:
        version = ALLOCATION_CHECKPOINT_VERSION
        configuration = asdict(config)
    # Default scheduling keeps the legacy checkpoint identity and semantics.
    if config.siqs is not None and config.siqs.solve_policy == "window":
        configuration["siqs"].pop("solve_policy")
    payload = {
        "version": version,
        "primality": policy,
        "schedule": WHEEL_VERSION
        if config.ecm_pair_wheel is not None
        else PAIRED_VERSION
        if config.ecm_pair_distance is not None
        else PROGRAM_VERSION
        if config.ecm_program_bytes
        else SCHEDULE_VERSION,
        "chains": _chain_policy(config),
        "backend": arithmetic.get_backend(config.backend).identity,
        "config": configuration,
        "state": state,
        "rng": generator.getstate()
        if (state["pending"] or state["current"] is not None)
        else None,
        "work_used": budget.used,
        "wall_used": wall_used,
        "cpu_used": cpu_used,
        "limits": {
            "work": None,
            "seconds": budget.seconds,
            "cpu_seconds": budget.cpu_seconds,
            "remaining_work": None,
            "remaining_seconds": max(0, budget.seconds - wall_used)
            if budget.seconds is not None
            else None,
            "remaining_cpu_seconds": max(0, budget.cpu_seconds - cpu_used)
            if budget.cpu_seconds is not None
            else None,
        },
    }
    if config.allocation is not None:
        payload["allocation"] = ALLOCATION_VERSION
    if config.ecm_chain_mode == "off":
        payload.pop("chains")
    if "ecm_auto" in state:
        # Older readers must reject automatic state before they can ignore its
        # worker metadata. Retain the mathematical subschema for validation.
        payload["portfolio_version"] = payload["version"]
        payload["version"] = AUTOMATIC_CHECKPOINT_VERSION
    encoded = _canonical(payload)
    checkpoint = {
        "payload": json.loads(encoded),
        "sha256": hashlib.sha256(encoded.encode()).hexdigest(),
    }
    overhead = [budget.wall_used - wall_used, budget.cpu_used - cpu_used]
    # Charge encoding on resume too; repeated pauses cannot reset its cost.
    checkpoint["serialization_overhead"] = overhead
    checkpoint["overhead_sha256"] = hashlib.sha256(
        _canonical(overhead).encode()
    ).hexdigest()
    return checkpoint


def _verify_progress(current, config, policy):
    """Verify retained witnesses and complete prime buffers before reuse."""
    witness = current.get("prime_job")
    if witness:
        n = current["n"]
        expected_shifts = ((n - 1) & -(n - 1)).bit_length() - 1
        expected_bases = _classification_bases(n, policy)

        rounds = (
            len(expected_bases) if expected_bases else config.primality_rounds
        )
        utils.require_integer(witness["index"], "witness index", 0)
        if (
            current["stage"] != "classify"
            or n < 41 * 41
            or witness["s"] != expected_shifts
            or witness["d"] != (n - 1) >> expected_shifts
            or witness["bases"] != expected_bases
            or witness["index"] >= rounds
            or len(witness["tested"]) != witness["index"]
        ):
            raise ValueError("invalid primality progress metadata")

        for index, base in enumerate(witness["tested"]):
            utils.require_integer(base, "witness", 2)
            if expected_bases and base != expected_bases[index]:
                raise ValueError("invalid deterministic witness prefix")
            if not expected_bases and base >= n - 1:
                raise ValueError("random witness is outside its domain")
            if not utils._strong_probable_prime(
                n, base, witness["d"], witness["s"]
            ):
                raise ValueError("checkpoint retained a failed witness")

    cursors = [current.get("cursor")]
    if current["job"]:
        if current["job"]["kind"] != current["stage"]:
            raise ValueError("candidate kind disagrees with portfolio stage")
        job = current["job"]
        if job["kind"] == "pm1" and not _legacy_pm1(config):
            # Reject oversized new tables before backend rehydration copies
            # their contents. Full arithmetic verification remains charged.
            if (
                len(job.get("even_powers", [])) > 64
                or len(job.get("gap_powers", {})) > 64
            ):
                raise ValueError("p-1 gap storage exceeds its cap")
            if (job["b1"], job["b2"], job["seed"]) != (
                config.pm1_b1,
                config.pm1_b2,
                current["attempt"],
            ) or current["attempt"] >= config.pm1_attempts:
                raise ValueError("p-1 assignment disagrees with campaign")
        if config.ecm_program_bytes and job["kind"] == "ecm":
            tier = utils.require_integer(current["tier"], "ECM tier", 0)
            if tier >= len(config.ecm_tiers):
                raise ValueError(
                    "ECM program tier exceeds configured campaign"
                )
            b1, b2, curves = config.ecm_tiers[tier]
            attempt = utils.require_integer(
                current["attempt"], "ECM attempt", 0
            )
            if attempt >= curves or (job["b1"], job["b2"]) != (b1, b2):
                raise ValueError("ECM program bounds disagree with campaign")
            stage_one = job["phase"] in ("setup", "stage_one", "replay")
            if job["cursor"]["hi"] != (b1 + 1 if stage_one else b2 + 1):
                raise ValueError("ECM program cursor has the wrong endpoint")
        if "chain_identity" in job:
            if (
                job["kind"] != "ecm"
                or not _chain_job_supported(job, config)
                or job["chain_identity"] != _chain_identity(job["b1"], config)
            ):
                raise ValueError("incompatible ECM chain progress")
            utils.require_integer(job["chain_chunks"], "chain chunks", 1)
            utils.require_integer(job["chain_replays"], "chain replays", 0)
            if not job["chain_replays"] <= job["chain_chunks"] <= 303:
                raise ValueError("invalid ECM chain recovery counters")
        elif any(key in job for key in ("chain_chunks", "chain_replays")):
            raise ValueError("missing ECM chain progress identity")
        cursors.append(job.get("cursor"))
    if not any(cursor is not None for cursor in cursors):
        return
    verifier = SieveContext(config.max_hi, segment_size=config.segment_size)

    for cursor in cursors:
        if cursor is None:
            continue
        for name in ("left", "next", "hi", "index"):
            utils.require_integer(cursor[name], name, 0)
        # The odd-slot bound excludes 2, which can share a tiny first block
        # with every odd candidate. Exact buffer verification still follows.
        max_values = config.segment_size + int(
            cursor["left"] <= 2 < cursor["next"]
        )
        if not (
            cursor["left"] <= cursor["next"] <= cursor["hi"] <= config.max_hi
            and cursor["next"] - cursor["left"] <= 2 * config.segment_size
            and 0 <= cursor["index"] <= len(cursor["values"]) <= max_values
        ):
            raise ValueError("invalid buffered prime metadata")

        expected = list(verifier.primes(cursor["left"], cursor["next"]))
        if expected != cursor["values"]:
            raise ValueError("corrupt buffered prime values")

    if current["job"] and current["job"]["kind"] == "ecm":
        job = current["job"]
        if _chain_job_supported(job, config):
            if config.ecm_chain_family == "auto":
                from .ecm.chains import verify_progress

                verify_progress(job, config.backend, verifier)
            else:
                from .ecm.chain_options import verify_progress

                verify_progress(
                    job, config.backend, verifier, config.ecm_chain_family
                )

    if (
        current["job"]
        and current["job"]["kind"] == "ecm"
        and (
            config.ecm_pair_distance is not None
            or config.ecm_pair_wheel is not None
        )
    ):
        from .ecm.paired import verify_progress

        verify_progress(current["job"], config, verifier)


@decode_checkpoint
def _unpack(checkpoint, config):
    """Reject corrupt, incompatible, oversized, or inconsistent snapshots."""
    try:
        encoded = _canonical(checkpoint["payload"])
        payload = json.loads(encoded)
        overhead = checkpoint["serialization_overhead"]
        if len(overhead) != 2 or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0
            for value in overhead
        ):
            raise ValueError("invalid serialization resource metadata")

        if (
            hashlib.sha256(_canonical(overhead).encode()).hexdigest()
            != (checkpoint["overhead_sha256"])
        ):
            raise ValueError("serialization overhead checksum mismatch")

        if (
            hashlib.sha256(encoded.encode()).hexdigest()
            != checkpoint["sha256"]
        ):
            raise ValueError("checkpoint checksum mismatch")

        if payload["version"] == AUTOMATIC_CHECKPOINT_VERSION:
            base_version = payload.get("portfolio_version")
            if (
                type(base_version) is not int
                or not 2 <= base_version <= ALLOCATION_CHECKPOINT_VERSION
                or not isinstance(payload.get("state"), dict)
                or "ecm_auto" not in payload["state"]
                or "ecm_worker_identity" not in payload["state"]
            ):
                raise ValueError("invalid automatic checkpoint envelope")
            payload["version"] = base_version
        elif "portfolio_version" in payload or (
            isinstance(payload.get("state"), dict)
            and "ecm_auto" in payload["state"]
        ):
            raise ValueError("automatic checkpoint requires version 13")

        policy = payload.get("primality", utils.LEGACY_PRIMALITY_POLICY)
        if policy not in (
            utils.PRIMALITY_POLICY,
            utils.LEGACY_PRIMALITY_POLICY,
        ):
            raise ValueError("incompatible primality policy")

        expected_config = json.loads(_canonical(asdict(config)))
        if (
            config.siqs is not None
            and config.siqs.solve_policy == "window"
            and "solve_policy" not in (payload["config"].get("siqs") or {})
        ):
            expected_config["siqs"].pop("solve_policy")
        legacy = config.sss is None and (
            payload["version"] == 3
            or (payload["version"] == 2 and config.siqs is None)
        )
        if payload["version"] != ALLOCATION_CHECKPOINT_VERSION:
            if config.allocation is None:
                expected_config.pop("allocation")
            if (
                config.ecm_chain_family == "auto"
                or config.ecm_chain_mode == "off"
            ):
                expected_config.pop("ecm_chain_family")
            if config.ecm_chain_mode == "off":
                expected_config.pop("ecm_chain_mode")
                expected_config.pop("ecm_chain_bytes")
            if _legacy_pm1(config) and payload["version"] < 9:
                expected_config.pop("pm1_gap_mode")
                expected_config.pop("pm1_chunk_size")
            if config.ecm_pair_distance is None:
                expected_config.pop("ecm_pair_distance")
            if config.ecm_pair_wheel is None:
                expected_config.pop("ecm_pair_wheel")
            if not config.ecm_program_bytes and (
                "ecm_program_bytes" not in payload["config"]
            ):
                expected_config.pop("ecm_program_bytes")
            legacy = config.sss is None and (
                payload["version"] == 3
                or (payload["version"] == 2 and config.siqs is None)
            )
            if legacy and "sss" not in payload["config"]:
                expected_config.pop("sss")
            if payload["version"] == 2 and legacy:
                expected_config.pop("siqs")
            if payload["version"] < 6 and "backend" not in payload["config"]:
                expected_config.pop("backend")
                for name in ("siqs", "sss"):
                    if expected_config.get(name) is not None:
                        expected_config[name].pop("backend", None)

        if (
            type(payload["version"]) is not int
            or payload["version"]
            not in (
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                CHECKPOINT_VERSION,
                OPTIONAL_CHAIN_CHECKPOINT_VERSION,
                ALLOCATION_CHECKPOINT_VERSION,
            )
            or (payload["version"] < 9 and not _legacy_pm1(config))
            or (payload["version"] < 10 and config.ecm_chain_mode != "off")
            or (payload["version"] == 10 and config.ecm_chain_mode == "off")
            or (
                payload["version"] < 11
                and config.ecm_chain_mode != "off"
                and config.ecm_chain_family != "auto"
            )
            or (
                payload["version"] == 11
                and (
                    config.ecm_chain_mode == "off"
                    or config.ecm_chain_family == "auto"
                )
            )
            or (
                payload.get("allocation")
                != (
                    ALLOCATION_VERSION
                    if config.allocation is not None
                    else None
                )
            )
            or (
                (payload["version"] == ALLOCATION_CHECKPOINT_VERSION)
                != (config.allocation is not None)
            )
            or payload.get("chains") != _chain_policy(config)
            or (payload["version"] < 8 and config.ecm_pair_wheel is not None)
            or (
                payload["version"] < 7 and config.ecm_pair_distance is not None
            )
            or (payload["version"] in (2, 3) and not legacy)
            or (payload["version"] < 5 and config.ecm_program_bytes != 0)
            or (payload["version"] < 5 and config.backend != "python-int")
            or (
                payload["version"] == 5
                and "backend" in payload["config"]
                and config.ecm_program_bytes != 0
            )
            or payload["schedule"]
            != (
                WHEEL_VERSION
                if config.ecm_pair_wheel is not None
                else PAIRED_VERSION
                if config.ecm_pair_distance is not None
                else PROGRAM_VERSION
                if config.ecm_program_bytes
                else SCHEDULE_VERSION
            )
            or payload["backend"]
            != arithmetic.get_backend(config.backend).identity
            or payload["config"] != expected_config
        ):
            raise ValueError("incompatible checkpoint metadata")

        state = payload["state"]
        utils.require_integer(payload["work_used"], "work_used", 0)
        for name in ("wall_used", "cpu_used"):
            value = payload[name]
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError("invalid checkpoint resource metadata")
            if not math.isfinite(value) or value < 0:
                raise ValueError("invalid checkpoint resource metadata")

        payload["wall_used"] += overhead[0]
        payload["cpu_used"] += overhead[1]
        utils.require_integer(state["original"])
        if not state["original"] or abs(state["original"]).bit_length() > (
            config.max_input_bits
        ):
            raise ValueError("checkpoint input exceeds configured limit")
        if (
            max(
                len(state[name])
                for name in (
                    "factors",
                    "remaining",
                    "pending",
                    "classifications",
                )
            )
            > config.max_input_bits
        ):
            raise ValueError("checkpoint exceeds input-derived storage limit")

        for value, certainty in state["classifications"].items():
            n = int(value)
            utils.require_integer(n, "classified cofactor", 2)
            if n > abs(state["original"]):
                raise ValueError("classification belongs to another modulus")
            classification = utils.Primality(certainty)
            if classification is not utils.Primality.COMPOSITE:
                actual = utils.classify_prime(
                    n, tolerance=config.primality_rounds, rng=random.Random(0)
                )
                if actual is utils.Primality.COMPOSITE or (
                    classification is utils.Primality.PROVEN
                    and (
                        actual is not utils.Primality.PROVEN
                        or _classification_bases(n, policy) is None
                    )
                ):
                    raise ValueError("invalid cached primality evidence")

        for value, (multiplicity, certainty) in state["factors"].items():
            utils.require_integer(int(value), "factor", 2)
            utils.require_integer(multiplicity, "multiplicity", 1)
            if int(value) > abs(state["original"]) or multiplicity > (
                config.max_input_bits
            ):
                raise ValueError("checkpoint multiplicity exceeds input cap")
            if utils.Primality(certainty) is utils.Primality.COMPOSITE:
                raise ValueError("composite recorded as terminal factor")
            actual = utils.classify_prime(
                int(value),
                tolerance=config.primality_rounds,
                rng=random.Random(0),
            )
            if actual is utils.Primality.COMPOSITE or (
                certainty == utils.Primality.PROVEN.value
                and (
                    actual is not utils.Primality.PROVEN
                    or _classification_bases(int(value), policy) is None
                )
            ):
                raise ValueError(
                    "checkpoint contains invalid primality evidence"
                )

        for value, multiplicity in state["pending"]:
            utils.require_integer(value, "cofactor", 2)
            utils.require_integer(multiplicity, "multiplicity", 1)
            if (
                value > abs(state["original"])
                or multiplicity > config.max_input_bits
            ):
                raise ValueError("checkpoint multiplicity exceeds input cap")

        for value in state["remaining"]:
            utils.require_integer(value, "cofactor", 2)
            if value > abs(state["original"]):
                raise ValueError("checkpoint cofactor exceeds its input")
        current = state["current"]
        if current:
            utils.require_integer(current["n"], "cofactor", 2)
            utils.require_integer(current["mult"], "multiplicity", 1)
            if current["n"] > abs(state["original"]):
                raise ValueError("checkpoint cofactor exceeds its input")
            if current["mult"] > config.max_input_bits or current[
                "stage"
            ] not in (
                "twos",
                "classify",
                "trial",
                "powers",
                "fermat",
                "rho",
                "pm1",
                "ecm",
                "siqs",
                "sss",
            ):
                raise ValueError("invalid current-cofactor metadata")

        if config.allocation is not None:
            if state.get("stop_reason") not in (None, "pretest_exhausted"):
                raise ValueError("invalid allocation stop reason")
            if current:
                if "fallback_admitted" in current and (
                    type(current["fallback_admitted"]) is not bool
                    or current["stage"] not in ("siqs", "sss")
                ):
                    raise ValueError("invalid fallback admission metadata")
                if "pending_handoff" in current and (
                    current["stage"] not in ("powers", "fermat")
                    or current["pending_handoff"]
                    not in (
                        "pretest_work",
                        "pretest_wall",
                        "pretest_cpu",
                        "reserved_work",
                        "reserved_wall",
                        "reserved_cpu",
                    )
                ):
                    raise ValueError("invalid pending handoff metadata")

        if (
            current
            and "worker_checkpoint" in current
            and (
                not isinstance(state.get("ecm_worker_identity"), str)
                or current["stage"] != "ecm"
                or current["job"] is not None
            )
        ):
            raise ValueError("unbound ECM worker progress")

        if current and "ecm_serial" in current:
            if (
                current["ecm_serial"] is not True
                or not isinstance(state.get("ecm_worker_identity"), str)
                or "worker_checkpoint" in current
                or current["stage"] not in ("ecm", "siqs", "sss")
            ):
                raise ValueError("invalid incumbent ECM selection progress")

        if current and current["job"] and current["job"]["n"] != current["n"]:
            raise ValueError("candidate modulus disagrees with parent")
        if current:
            if current["stage"] == "siqs":
                if config.siqs is None or current["job"] is not None:
                    raise ValueError("incompatible SIQS dispatcher progress")
                utils.require_integer(
                    current["siqs_start_work"], "SIQS work start", 0
                )
            elif current["stage"] == "sss":
                if config.sss is None or current["job"] is not None:
                    raise ValueError("incompatible SSS dispatcher progress")
                utils.require_integer(
                    current["sss_start_work"], "SSS work start", 0
                )
            else:
                _verify_progress(current, config, policy)

        # Validate reconstruction after checking integer exponent bounds.
        _result(state)
        generator = random.Random()
        if payload["rng"] is not None:
            generator.setstate(_tuples(payload["rng"]))
    except (KeyError, TypeError, AttributeError, IndexError) as error:
        raise ValueError("malformed checkpoint") from error

    payload["primality"] = policy
    return payload, generator


def _promote_state(state, name):
    """Restore only arithmetic state; counters, cursors and RNG stay native."""
    backend = arithmetic.get_backend(name)
    state["pending"] = [
        [backend.integer(n), mult] for n, mult in state["pending"]
    ]
    state["remaining"] = [backend.integer(n) for n in state["remaining"]]
    current = state["current"]
    if current is not None:
        current["n"] = backend.integer(current["n"])
        if current.get("prime_job"):
            witness = current["prime_job"]
            witness["d"] = backend.integer(witness["d"])
        if current.get("job"):
            from .execution.stage_jobs import promote_job

            promote_job(current["job"], backend)


@removed_arguments(("ecm_memory_bytes",))
def factorize_bounded(
    n,
    *,
    seed=0,
    config=None,
    budget=None,
    checkpoint=None,
    stop_after_split=False,
    ecm_executor="auto",
):
    """Factor with automatic cost-gated ECM; None explicitly selects serial.

    Automatically created contexts are closed before returning or raising.
    Explicit contexts remain borrowed. Omitted budgets/memory receive finite
    parallel-capable grants; explicit allowances and legacy checkpoints win.
    """
    options = dict(
        seed=seed,
        config=config,
        budget=budget,
        checkpoint=checkpoint,
        stop_after_split=stop_after_split,
    )
    if isinstance(ecm_executor, str):
        if ecm_executor != "auto":
            raise ValueError("ECM executor must be auto, None or a context")
        from .ecm.auto import AutomaticECMOwner

        with AutomaticECMOwner(seed=seed) as owner:
            return _factorize_bounded(
                n, ecm_executor=owner, _owned_ecm=True, **options
            )
    return _factorize_bounded(n, ecm_executor=ecm_executor, **options)


def _factorize_bounded(
    n,
    *,
    seed=0,
    config=None,
    budget=None,
    checkpoint=None,
    stop_after_split=False,
    ecm_executor=None,
    _owned_ecm=False,
):
    """Factor under one budget, retaining exact state on pause/exhaustion.

    Resume with the same config and a Budget whose total allowances include
    consumed resources. Increasing an allowance explicitly grants extra work.
    Cancellation never removes a cofactor. Exhausted ECM schedules reach
    optional SIQS/SSS under the same allowance when configured.
    stop_after_split stops when a proper divisor is exposed, retaining children
    for a later full-factorization resume under the identical configuration.
    """
    utils.require_integer(n)
    if n == 0:
        raise ValueError("zero has no finite prime factorization")
    if type(stop_after_split) is not bool:
        raise TypeError("stop_after_split must be Boolean")
    if config is None:
        # A resume retains the saved p-1 executor; fresh calls use the accepted
        # default. Other configuration still has to match the snapshot.
        saved_config = {}
        if checkpoint is not None:
            try:
                saved_config = checkpoint["payload"]["config"]
                if type(saved_config) is not dict:
                    raise ValueError("malformed checkpoint config")
            except (KeyError, TypeError) as error:
                raise ValueError("malformed checkpoint") from error
        resume_options = {}
        if checkpoint is not None:
            resume_options = {
                "ecm_chain_mode": saved_config.get("ecm_chain_mode", "off"),
                "ecm_chain_family": saved_config.get(
                    "ecm_chain_family", "auto"
                ),
                "ecm_chain_bytes": saved_config.get("ecm_chain_bytes", 0),
                "ecm_program_bytes": saved_config.get("ecm_program_bytes", 0),
                "ecm_tiers": saved_config.get(
                    "ecm_tiers", ((constants.ECM_B1, constants.ECM_B2, 32),)
                ),
            }
        if (
            saved_config.get("allocation") is not None
            or saved_config.get("siqs") is not None
        ):
            saved_config = dict(saved_config)
            if saved_config.get("allocation") is not None:
                saved_config["allocation"] = ECMAllocation(
                    **saved_config["allocation"]
                )
            if saved_config.get("siqs") is not None:
                saved_config["siqs"] = _decode_siqs_config(
                    saved_config["siqs"]
                )
            if saved_config.get("sss") is not None:
                from .qs import SieveConfig

                values = dict(saved_config["sss"])
                values["collector"] = SieveConfig(**values["collector"])
                saved_config["sss"] = SSSConfig(**values)
            config = PortfolioConfig(**saved_config)
        else:
            config = PortfolioConfig(
                **resume_options,
                siqs=None if checkpoint is not None else _AUTO_SIQS,
                pm1_gap_mode=saved_config.get("pm1_gap_mode", "cached")
                if checkpoint is not None
                else "recurrence",
                pm1_chunk_size=saved_config.get("pm1_chunk_size")
                if checkpoint is not None
                else 64,
            )
    automatic_memory = config._automatic_memory
    if checkpoint is not None:
        try:
            _owned_ecm and ("ecm_auto" in checkpoint["payload"]["state"])
        except (KeyError, TypeError) as error:
            raise ValueError("malformed checkpoint") from error
    if budget is None:
        budget = Budget()
    if abs(n).bit_length() > config.max_input_bits:
        raise ValueError("input exceeds max_input_bits")
    if checkpoint is None:
        config = config.for_input(n)
    elif config._automatic_ecm_tiers:
        # An omitted tier setting means retain the checkpoint's exact plan,
        # including legacy fixed32; explicit tier conflicts still reject.
        try:
            tiers = checkpoint["payload"]["config"]["ecm_tiers"]
        except (KeyError, TypeError) as error:
            raise ValueError("malformed checkpoint") from error
        config = replace(config, ecm_tiers=tiers)

    if checkpoint is None:
        if budget.used or budget.prior_wall or budget.prior_cpu:
            raise ValueError("a fresh run requires an unused budget")
        generator = utils.resolve_rng(seed)
        policy = utils.PRIMALITY_POLICY
        state = {
            "original": n,
            "pending": [[abs(n), 1]] if abs(n) > 1 else [],
            "current": None,
            "factors": {},
            "remaining": [],
            "classifications": {},
            "events": [],
            "dropped_events": 0,
            "stage_seconds": {},
            "context_ready": False,
        }
    else:
        payload, generator = _unpack(checkpoint, config)
        policy = payload["primality"]
        state = payload["state"]
        if state["original"] != n:
            raise ValueError("checkpoint belongs to another input")
        budget.restore_charges(
            work=payload["work_used"],
            wall=payload["wall_used"],
            cpu=payload["cpu_used"],
        )

    if _owned_ecm:
        if checkpoint is not None and "ecm_auto" not in state:
            # Legacy serial/custom checkpoints never adopt a new stream.
            ecm_executor = None
        else:
            ecm_executor.configure(
                config,
                generator,
                state.get("ecm_auto"),
                automatic_memory=automatic_memory,
            )
            state["ecm_auto"] = ecm_executor.header

    # Opt-in search semantics never impersonate the incumbent RNG/checkpoint.
    saved_executor = state.get("ecm_worker_identity")
    supplied_executor = None if ecm_executor is None else ecm_executor.identity
    if checkpoint is not None and saved_executor != supplied_executor:
        raise ValueError(
            "ECM worker checkpoint requires its original executor/rule"
        )
    if ecm_executor is not None:
        state["ecm_worker_identity"] = supplied_executor

    _promote_state(state, config.backend)

    reason = "exhausted"
    siqs_runtime = {}
    optional_budget = (
        PretestBudget(
            budget, config.allocation, fallback=bool(config.siqs or config.sss)
        )
        if config.allocation is not None
        else budget
    )

    try:
        budget.consume(0)
        if checkpoint is not None and config.allocation is not None:
            # Verification has already reconstructed the bounded snapshot.
            # Charge its bytes on the shared ledger even if pretesting stops;
            # actual verification/decoding time is included by Budget's clock.
            budget.consume(len(_canonical(checkpoint["payload"]).encode()))
        context = programs = None
        if config.allocation is None:
            context, programs = _execution_context(state, config, budget)
        if checkpoint is not None and state["current"] is not None:
            try:
                job = state["current"].get("job")
                restore_budget = budget
                if config.allocation is not None and job is not None:
                    restore_budget = optional_budget
                if job is not None and _chain_job_supported(job, config):
                    # Bounded schedule/coordinate reconstruction was checked by
                    # unpack; a refused continuation must still pay that work.
                    restore_budget.consume(
                        2000 + job["n"].bit_length() + len(job["powers"])
                    )
                if (
                    job is not None
                    and job["kind"] == "pm1"
                    and not _legacy_pm1(config)
                ):
                    from .pm1.gaps import verify_powers

                    verify_powers(job, restore_budget, config.pm1_gap_mode)
            except HandoffRequiredError as error:
                _handoff(state, config, budget, str(error))
        while state["pending"] or state["current"] is not None:
            stage = (
                state["current"]["stage"] if state["current"] else "dispatch"
            )
            started = time.perf_counter()
            optional = stage in (
                "trial",
                "fermat",
                "rho",
                "pm1",
                "ecm",
            )
            ledger = budget
            if optional and config.allocation is not None:
                ledger = optional_budget
            try:
                current = state["current"]
                if (
                    current is not None
                    and current.get("pending_handoff")
                    and stage == "fermat"
                ):
                    _handoff(state, config, budget, current["pending_handoff"])
                    continue
                if ledger is not budget:
                    ledger.consume(0)
                if context is None and stage in ("trial", "pm1", "ecm"):
                    context, programs = _execution_context(
                        state, config, ledger
                    )
                _advance(
                    state,
                    config,
                    ledger,
                    context,
                    generator,
                    siqs_runtime,
                    programs,
                    policy,
                    ecm_executor,
                )
            except HandoffRequiredError as error:
                # Any retained chunk is already atomic and charged. Do not
                # retry a partial curve after the one-way fallback transition.
                _handoff(state, config, budget, str(error))
            state["stage_seconds"][stage] = (
                state["stage_seconds"].get(stage, 0)
                + time.perf_counter()
                - started
            )
            if stop_after_split:
                partial = _result(state)
                values = [factor.value for factor in partial.factors]
                values.extend(partial.remaining)
                if any(utils.valid_divisor(value, abs(n)) for value in values):
                    reason = "factor_found"
                    break
    except BudgetExhaustedError as error:
        reason = str(error)

    result = _result(state)
    if result.complete:
        reason = "complete"
    elif reason == "exhausted":
        reason = state.get("stop_reason", reason)
    # Live fallback state must reach the parent snapshot before packing.
    if siqs_runtime.get("job") is not None:
        stage = state["current"]["stage"]
        state["current"][stage + "_checkpoint"] = siqs_runtime[
            "job"
        ].checkpoint()

    if _owned_ecm and ecm_executor is not None:
        # Include owned-pool shutdown in the final clocks/checkpoint. The outer
        # owner also closes on exceptions before result construction finishes.
        ecm_executor.close()
    snapshot = _pack(state, config, budget, generator, policy)
    return PortfolioRun(
        result,
        reason,
        budget.used,
        budget.wall_used,
        budget.cpu_used,
        arithmetic.canonical(tuple(state["events"])),
        state["dropped_events"],
        snapshot,
        budget.stop,
    )
