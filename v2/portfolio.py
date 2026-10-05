"""Bounded serial portfolio with reconstructible, versioned checkpoints."""

import hashlib
import json
import math
import random
import time
from dataclasses import asdict, dataclass
from math import isqrt

from . import constants, prime_sieve, utils
from .budget import Budget, BudgetExhaustedError
from .factor import FactorizationResult, PrimeFactor
from .preprocessing import (
    fermat_step,
    integer_root,
    power_residue_possible,
    strip_twos,
)
from .qs import SIQSConfig, SIQSJob
from .qs.sss import SSSConfig, SSSJob
from .schedules import ScheduleCache, SieveContext
from .stage_jobs import (
    advance_job,
    new_job,
    peek_prime,
    prime_cursor,
    take_prime,
)

CHECKPOINT_VERSION = 4
SCHEDULE_VERSION = "half-open-prime-powers/ecm-even-baby-v1"


@dataclass(frozen=True)
class PortfolioConfig:
    """Finite allowances; ECM tiers are (inclusive B1, B2, curve count).

    Defaults retain the existing ECM bounds and serial execution. Alternative
    tiers/cutoffs are experimental until matched held-out evidence supports
    them. memory_bytes bounds conservative owned workspace estimates, not RSS.
    """

    trial_bound: int = constants.TRIAL_BOUND
    rho_attempts: int = 4
    rho_evaluations: int = 5000
    rho_batch: int = constants.RHO_BATCH_SIZE
    recovery_limit: int = constants.RHO_RECOVERY_LIMIT
    pm1_attempts: int = 1
    pm1_b1: int = constants.PM1_B1
    pm1_b2: int = constants.PM1_B2
    ecm_tiers: tuple = ((constants.ECM_B1, constants.ECM_B2, 32),)
    chunk_size: int = 16
    gcd_batch: int = constants.GCD_BATCH_SIZE
    primality_rounds: int = constants.PRIMALITY_ROUNDS
    fermat_steps: int = 0
    memory_bytes: int = 8_388_608
    segment_size: int = 1024
    trial_chunk: int = 64
    schedule_cache_bytes: int = 0
    max_input_bits: int = 4096
    trace_limit: int = 256
    rolling: bool = False
    siqs: SIQSConfig | None = None
    sss: SSSConfig | None = None

    def __post_init__(self):
        """Validate bounds and reserve storage before any allocation."""
        minima = {
            "trial_bound": 2,
            "rho_attempts": 0,
            "rho_evaluations": 0,
            "rho_batch": 1,
            "recovery_limit": 0,
            "pm1_attempts": 0,
            "pm1_b1": 2,
            "pm1_b2": self.pm1_b1,
            "chunk_size": 1,
            "gcd_batch": 1,
            "primality_rounds": 1,
            "fermat_steps": 0,
            "memory_bytes": 65_536,
            "segment_size": 1,
            "trial_chunk": 1,
            "schedule_cache_bytes": 0,
            "max_input_bits": 2,
            "trace_limit": 0,
        }
        for name, minimum in minima.items():
            utils.require_integer(getattr(self, name), name, minimum)
        if type(self.rolling) is not bool:
            raise TypeError("rolling must be Boolean")
        if 0 < self.schedule_cache_bytes < 4096:
            raise ValueError("schedule cache needs at least 4096 bytes")
        tiers = tuple(tuple(tier) for tier in self.ecm_tiers)
        if (
            len(tiers) > 64
            or max(
                self.chunk_size,
                self.gcd_batch,
                self.trial_chunk,
                self.primality_rounds,
            )
            > 256
        ):
            raise ValueError("limit tiers to 64 and arithmetic chunks to 256")
        for tier in tiers:
            if len(tier) != 3:
                raise ValueError("ECM tiers need B1, B2, and curves")
            b1, b2, curves = tier
            utils.require_integer(b1, "B1", 2)
            utils.require_integer(b2, "B2", b1)
            utils.require_integer(curves, "curves", 0)
        object.__setattr__(self, "ecm_tiers", tiers)
        if self.siqs is not None and self.sss is not None:
            raise ValueError("choose one relation fallback: siqs or sss")
        if self.sss is not None:
            if not isinstance(self.sss, SSSConfig):
                raise TypeError("sss must be an SSSConfig or None")
            if (
                self.sss.memory_bytes + self.workspace_reserve + 8192
                > self.memory_bytes
            ):
                raise MemoryError(
                    "portfolio/SSS coexistence exceeds memory cap"
                )
        if self.siqs is not None:
            if not isinstance(self.siqs, SIQSConfig):
                raise TypeError("siqs must be a SIQSConfig or None")
            if (
                self.siqs.memory_bytes + self.workspace_reserve + 8192
                > self.memory_bytes
            ):
                raise MemoryError(
                    "portfolio/SIQS coexistence exceeds memory cap"
                )
        if self.memory_bytes - self.workspace_reserve < 8192:
            raise MemoryError(
                "candidate/checkpoint reserve exceeds memory cap"
            )

    @property
    def max_hi(self):
        """Largest schedule endpoint, including optional exact-power roots."""
        return 1 + max(
            self.trial_bound,
            self.pm1_b2 if self.pm1_attempts else 2,
            self.max_input_bits,
            *(b2 for _, b2, curves in self.ecm_tiers if curves),
        )

    @property
    def workspace_reserve(self):
        """Conservative bound for cursor, point, result, and JSON storage.

        Reserve object overhead explicitly because PyPy cannot report useful
        sys.getsizeof values. A final serialized-size check also caps output.
        """
        coordinate_bytes = 128 + self.max_input_bits // 8
        distance = max(
            (
                min(isqrt(b2), (b1 - 1) // 2)
                for b1, b2, curves in self.ecm_tiers
                if curves
            ),
            default=0,
        )
        return (
            16_384
            + 8 * coordinate_bytes * (distance + self.gcd_batch)
            + 256 * self.segment_size
            + 256 * self.max_input_bits
            + 1024 * self.trace_limit
            + 256 * self.chunk_size
            + self.schedule_cache_bytes
        )


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


def _canonical(value):
    """Encode exact decimal integers with stable keys for corruption checks."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


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
        ([divisor, multiplicity], [n // divisor, multiplicity])
    )
    state["current"] = None


def _classify_step(current, config, budget, generator):
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
        deterministic = n < utils.DETERMINISTIC_LIMIT
        if not deterministic:
            bases = None
        elif n < 9_080_191:
            bases = [31, 73]
        elif n < 4_759_123_141:
            bases = [2, 7, 61]
        else:
            bases = list(utils.DETERMINISTIC_BASES)
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
    budget.consume(n.bit_length() + witness_state["s"])
    base = bases[index] if bases else generator.randrange(2, n - 1)
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


def _advance(state, config, budget, context, generator, siqs_runtime):
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
            certainty = _classify_step(current, config, budget, generator)
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
                n //= prime
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
    if kind == "rho":
        attempts, b1, b2 = config.rho_attempts, 0, 0
    elif kind == "pm1":
        attempts, b1, b2 = config.pm1_attempts, config.pm1_b1, config.pm1_b2
    elif current["tier"] < len(config.ecm_tiers):
        b1, b2, attempts = config.ecm_tiers[current["tier"]]
    else:
        if config.siqs is not None:
            current.update(stage="siqs", job=None, siqs_start_work=budget.used)
        elif config.sss is not None:
            current.update(stage="sss", job=None, sss_start_work=budget.used)
        else:
            state["remaining"].extend([n] * current["mult"])
            state["current"] = None
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
    advance_job(job, budget, context, config)
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


def _pack(state, config, budget, generator):
    """Snapshot RNG, pending work, schedule identity, and resource use."""
    wall_used, cpu_used = budget.wall_used, budget.cpu_used
    payload = {
        "version": CHECKPOINT_VERSION,
        "schedule": SCHEDULE_VERSION,
        "backend": "python-int",
        "config": asdict(config),
        "state": state,
        "rng": generator.getstate()
        if (state["pending"] or state["current"] is not None)
        else None,
        "work_used": budget.used,
        "wall_used": wall_used,
        "cpu_used": cpu_used,
        "limits": {
            "work": budget.work_limit,
            "seconds": budget.seconds,
            "cpu_seconds": budget.cpu_seconds,
            "remaining_work": budget.work_limit - budget.used,
            "remaining_seconds": max(0, budget.seconds - wall_used)
            if budget.seconds is not None
            else None,
            "remaining_cpu_seconds": max(0, budget.cpu_seconds - cpu_used)
            if budget.cpu_seconds is not None
            else None,
        },
    }
    encoded = _canonical(payload)
    if len(encoded.encode()) > config.memory_bytes // 2:
        raise MemoryError("serialized checkpoint exceeds output reserve")
    checkpoint = {
        "payload": json.loads(encoded),
        "sha256": hashlib.sha256(encoded.encode()).hexdigest(),
    }
    overhead = [budget.wall_used - wall_used, budget.cpu_used - cpu_used]
    checkpoint["serialization_overhead"] = overhead
    checkpoint["overhead_sha256"] = hashlib.sha256(
        _canonical(overhead).encode()
    ).hexdigest()
    return checkpoint


def _verify_progress(current, config):
    """Verify retained witnesses and complete prime buffers before reuse."""
    witness = current.get("prime_job")
    if witness:
        n = current["n"]
        expected_shifts = ((n - 1) & -(n - 1)).bit_length() - 1
        expected_bases = None
        if n < utils.DETERMINISTIC_LIMIT:
            if n < 9_080_191:
                expected_bases = [31, 73]
            elif n < 4_759_123_141:
                expected_bases = [2, 7, 61]
            else:
                expected_bases = list(utils.DETERMINISTIC_BASES)
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
        cursors.append(current["job"].get("cursor"))
    if not any(cursor is not None for cursor in cursors):
        return
    verifier = SieveContext(config.max_hi, segment_size=config.segment_size)
    for cursor in cursors:
        if cursor is None:
            continue
        for name in ("left", "next", "hi", "index"):
            utils.require_integer(cursor[name], name, 0)
        if not (
            cursor["left"] <= cursor["next"] <= cursor["hi"] <= config.max_hi
            and cursor["next"] - cursor["left"] <= 2 * config.segment_size
            and 0
            <= cursor["index"]
            <= len(cursor["values"])
            <= config.segment_size
        ):
            raise ValueError("invalid buffered prime metadata")
        expected = list(verifier.primes(cursor["left"], cursor["next"]))
        if expected != cursor["values"]:
            raise ValueError("corrupt buffered prime values")


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
        if len(encoded.encode()) > config.memory_bytes // 2:
            raise ValueError("checkpoint exceeds configured cap")
        if (
            hashlib.sha256(encoded.encode()).hexdigest()
            != checkpoint["sha256"]
        ):
            raise ValueError("checkpoint checksum mismatch")
        expected_config = json.loads(_canonical(asdict(config)))
        legacy = config.sss is None and (
            payload["version"] == 3
            or (payload["version"] == 2 and config.siqs is None)
        )
        if legacy and "sss" not in payload["config"]:
            expected_config.pop("sss")
        if payload["version"] == 2 and legacy:
            expected_config.pop("siqs")
        if (
            type(payload["version"]) is not int
            or payload["version"] not in (2, 3, CHECKPOINT_VERSION)
            or (payload["version"] in (2, 3) and not legacy)
            or payload["schedule"] != SCHEDULE_VERSION
            or payload["backend"] != "python-int"
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
                    and actual is not utils.Primality.PROVEN
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
                and actual is not utils.Primality.PROVEN
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
                _verify_progress(current, config)
        # Validate reconstruction after checking integer exponent bounds.
        _result(state)
        generator = random.Random()
        if payload["rng"] is not None:
            generator.setstate(_tuples(payload["rng"]))
    except (KeyError, TypeError, AttributeError, IndexError) as error:
        raise ValueError("malformed checkpoint") from error
    return payload, generator


def factorize_bounded(
    n,
    *,
    seed=0,
    config=None,
    budget=None,
    checkpoint=None,
    stop_after_split=False,
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
    config = PortfolioConfig() if config is None else config
    budget = Budget() if budget is None else budget
    if abs(n).bit_length() > config.max_input_bits:
        raise ValueError("input exceeds max_input_bits")
    if checkpoint is None:
        if budget.used or budget.prior_wall or budget.prior_cpu:
            raise ValueError("a fresh run requires an unused budget")
        generator = utils.resolve_rng(seed)
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
        state = payload["state"]
        if state["original"] != n:
            raise ValueError("checkpoint belongs to another input")
        if budget.work_limit < payload["work_used"]:
            raise ValueError("budget is below already consumed work")
        budget.used = payload["work_used"]
        budget.prior_wall = payload["wall_used"]
        budget.prior_cpu = payload["cpu_used"]
    reason = "exhausted"
    siqs_runtime = {}
    try:
        budget.consume(0)
        if not state.get("context_ready"):
            budget.consume((isqrt(config.max_hi - 1) + 1) // 2)
        context = SieveContext(
            config.max_hi,
            memory_bytes=config.memory_bytes
            - config.workspace_reserve
            - (
                (config.siqs or config.sss).memory_bytes
                if (config.siqs or config.sss)
                else 0
            ),
            segment_size=config.segment_size,
            rolling=config.rolling,
        )
        state["context_ready"] = True
        if config.schedule_cache_bytes:
            context = ScheduleCache(
                context, cache_bytes=config.schedule_cache_bytes
            )
        while state["pending"] or state["current"] is not None:
            stage = (
                state["current"]["stage"] if state["current"] else "dispatch"
            )
            started = time.perf_counter()
            _advance(state, config, budget, context, generator, siqs_runtime)
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
    if siqs_runtime.get("job") is not None:
        stage = state["current"]["stage"]
        state["current"][stage + "_checkpoint"] = siqs_runtime[
            "job"
        ].checkpoint()
    snapshot = _pack(state, config, budget, generator)
    return PortfolioRun(
        result,
        reason,
        budget.used,
        budget.wall_used,
        budget.cpu_used,
        tuple(state["events"]),
        state["dropped_events"],
        snapshot,
    )
