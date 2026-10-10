"""Bounded serial portfolio with reconstructible, versioned checkpoints."""

import hashlib
import json
import math
import random
import time
from dataclasses import asdict, dataclass, field

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
    ECMPrograms,
)
from .execution.budget import Budget, BudgetExhaustedError
from .execution.schedules import ScheduleCache, SieveContext
from .execution.stage_jobs import (
    advance_job,
    new_job,
    peek_prime,
    prime_cursor,
    take_prime,
)
from .factor import FactorizationResult, PrimeFactor
from .qs import SIQSConfig, SIQSJob
from .qs.sss import SSSConfig, SSSJob

CHECKPOINT_VERSION = 10
OPTIONAL_CHAIN_CHECKPOINT_VERSION = 11
SCHEDULE_VERSION = "half-open-prime-powers/ecm-even-baby-v1"


@dataclass(frozen=True)
class PortfolioConfig:
    """Finite allowances; ECM tiers are (inclusive B1, B2, curve count).

    Defaults retain the existing ECM bounds and serial execution. Alternative
    tiers/cutoffs are experimental until matched held-out evidence supports
    them. memory_bytes bounds conservative owned workspace estimates, not RSS.
    """

    backend: str = field(default="python-int", kw_only=True)
    pm1_gap_mode: str = field(default="recurrence", kw_only=True)
    pm1_chunk_size: int | None = field(default=64, kw_only=True)
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
    memory_bytes: int | None = None
    segment_size: int = 1024
    trial_chunk: int = 64
    schedule_cache_bytes: int = 0
    max_input_bits: int = 4096
    trace_limit: int = 256
    rolling: bool = False
    siqs: SIQSConfig | None = None
    sss: SSSConfig | None = None
    ecm_program_bytes: int = 0
    ecm_chain_mode: str = field(default="auto", kw_only=True)
    ecm_chain_bytes: int = field(default=0, kw_only=True)
    ecm_chain_family: str = field(default="auto", kw_only=True)
    ecm_pair_distance: int | None = None
    ecm_pair_wheel: int | None = None

    def __post_init__(self):
        """Validate bounds and reserve storage before any allocation."""
        automatic_memory = self.memory_bytes is None
        if automatic_memory:
            object.__setattr__(self, "memory_bytes", 8_388_608)
        if self.pm1_gap_mode not in ("cached", "recurrence"):
            raise ValueError("p-1 gap mode must be cached or recurrence")
        if self.pm1_chunk_size is not None:
            utils.require_integer(self.pm1_chunk_size, "pm1_chunk_size", 1)
            if self.pm1_chunk_size > 256:
                raise ValueError("limit p-1 chunks to 256 primes")
        arithmetic.get_backend(self.backend)
        if self.ecm_chain_mode not in ("auto", "off", "reuse"):
            raise ValueError("ECM chain mode must be auto, off or reuse")
        if self.ecm_chain_family not in ("auto", "lucas", "cf"):
            raise ValueError("ECM chain family must be auto, lucas or cf")
        utils.require_integer(self.ecm_chain_bytes, "ecm_chain_bytes", 0)
        if self.ecm_chain_mode == "reuse":
            from .ecm.chains import MIN_MEMORY_BYTES

            if not self.ecm_program_bytes:
                raise ValueError(
                    "ECM chains require bounded reusable programs"
                )
            if self.ecm_chain_bytes < MIN_MEMORY_BYTES:
                raise MemoryError(
                    "ECM chain preparation exceeds configured cap"
                )
        elif self.ecm_chain_bytes and self.ecm_chain_mode == "off":
            raise ValueError("ECM chain memory requires reuse mode")
        for fallback in (self.siqs, self.sss):
            if fallback is not None and fallback.backend != self.backend:
                raise ValueError("portfolio and fallback backends must match")

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
            "ecm_program_bytes": 0,
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

        if self.ecm_chain_mode == "auto":
            from .ecm.chain_options import default_options
            from .ecm.chains import MIN_MEMORY_BYTES

            # Store the resolved policy, never an auto decision, so resume
            # retains its executor even if future defaults change again.
            for name, value in default_options(
                self, tiers, automatic_memory=automatic_memory
            ).items():
                object.__setattr__(self, name, value)
            if self.ecm_chain_mode == "reuse":
                utils.require_integer(
                    self.ecm_chain_bytes, "ecm_chain_bytes", MIN_MEMORY_BYTES
                )
            elif self.ecm_chain_bytes:
                raise ValueError("ECM chain memory requires reuse mode")

        if self.ecm_program_bytes:
            if self.ecm_program_bytes < 4096 + 256 * self.segment_size:
                raise MemoryError("ECM program scratch exceeds configured cap")
            if any(b2 >= 2**64 for _, b2, curves in tiers if curves):
                raise ValueError("ECM programs require B2 below 2**64")

        if self.ecm_pair_distance is not None:
            distance = utils.require_integer(
                self.ecm_pair_distance, "ecm_pair_distance", 0
            )
            if not self.ecm_program_bytes:
                raise ValueError("paired ECM requires a reusable program cap")
            if distance and (distance < 2 or distance % 2):
                raise ValueError("paired D must be zero or positive and even")
            if distance and any(
                2 * distance >= b1 - (b1 % 2 == 0)
                for b1, _, curves in tiers
                if curves
            ):
                raise ValueError("paired D requires positive initialization")
            if any(
                b2 + distance >= 2**64 for _, b2, curves in tiers if curves
            ):
                raise ValueError("paired centers must fit packed words")

        if self.ecm_pair_wheel is not None:
            wheel = utils.require_integer(
                self.ecm_pair_wheel, "ecm_pair_wheel", 2
            )
            if wheel % 2 or wheel > 2 * self.segment_size:
                raise ValueError("even wheel must fit one prime segment")
            if (
                self.ecm_pair_distance is not None
                or not self.ecm_program_bytes
            ):
                raise ValueError(
                    "wheel requires programs and excludes legacy D"
                )
            if any(
                b2 + wheel // 2 >= 2**64 for _, b2, curves in tiers if curves
            ):
                raise ValueError("wheel centers must fit packed words")

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
        paired_reserve = 0
        if (
            self.ecm_pair_distance is not None
            or self.ecm_pair_wheel is not None
        ):
            distance = (
                self.ecm_pair_wheel // 2 + 6
                if self.ecm_pair_wheel
                else self.ecm_pair_distance // 2 + 3
            )
            # Construction dict/sort/packing, decoded certificates, retained
            # replay records and their JSON copies coexist with point tables.
            paired_reserve = 8192 + 2048 * (
                self.segment_size + self.gcd_batch + 1
            )
        return (
            16_384
            + 8 * coordinate_bytes * (distance + self.gcd_batch)
            + 256 * self.segment_size
            + 256 * self.max_input_bits
            + 1024 * self.trace_limit
            + 256 * max(self.chunk_size, self.pm1_chunk_size or 0)
            + (
                8 * coordinate_bytes * 64
                if self.pm1_attempts and self.pm1_gap_mode == "recurrence"
                else 0
            )
            + self.schedule_cache_bytes
            + self.ecm_program_bytes
            + self.ecm_chain_bytes
            + paired_reserve
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


def _advance(
    state, config, budget, context, generator, siqs_runtime, programs, policy
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


def _chain_bounds(config):
    """Keep fresh/unsupported schedules on the measured B4 ladder."""
    if config.ecm_chain_mode == "off" or config.chunk_size != 16:
        return ()
    return tuple(
        b1 for b1, _, curves in config.ecm_tiers if b1 == 2000 and curves >= 8
    )


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


def _pack(state, config, budget, generator, policy):
    """Snapshot RNG, pending work, schedule identity, and resource use."""
    wall_used, cpu_used = budget.wall_used, budget.cpu_used
    configuration = asdict(config)
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
    if config.ecm_chain_mode == "off":
        payload.pop("chains")
    encoded = _canonical(payload)
    if len(encoded.encode()) > config.memory_bytes // 2:
        raise MemoryError("serialized checkpoint exceeds output reserve")
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

        policy = payload.get("primality", utils.LEGACY_PRIMALITY_POLICY)
        if policy not in (
            utils.PRIMALITY_POLICY,
            utils.LEGACY_PRIMALITY_POLICY,
        ):
            raise ValueError("incompatible primality policy")

        expected_config = json.loads(_canonical(asdict(config)))
        if config.ecm_chain_family == "auto" or config.ecm_chain_mode == "off":
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
                "memory_bytes": saved_config.get("memory_bytes", 8_388_608),
            }
        config = PortfolioConfig(
            **resume_options,
            pm1_gap_mode=saved_config.get("pm1_gap_mode", "cached")
            if checkpoint is not None
            else "recurrence",
            pm1_chunk_size=saved_config.get("pm1_chunk_size")
            if checkpoint is not None
            else 64,
        )
    budget = Budget() if budget is None else budget
    if abs(n).bit_length() > config.max_input_bits:
        raise ValueError("input exceeds max_input_bits")
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
        if budget.work_limit < payload["work_used"]:
            raise ValueError("budget is below already consumed work")
        budget.used = payload["work_used"]
        budget.prior_wall = payload["wall_used"]
        budget.prior_cpu = payload["cpu_used"]

    _promote_state(state, config.backend)

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
        programs = (
            ECMPrograms(context, memory_bytes=config.ecm_program_bytes)
            if config.ecm_program_bytes
            else None
        )
        if config.ecm_chain_mode != "off":
            if config.ecm_chain_family == "auto":
                from .ecm.chains import ChainPlans

                programs.chains = ChainPlans(
                    config.ecm_chain_bytes,
                    config.backend,
                    _chain_bounds(config),
                )
            else:
                from .ecm.chain_options import ChainPlans

                programs.chains = ChainPlans(
                    config.ecm_chain_bytes,
                    config.backend,
                    _chain_bounds(config),
                    config.ecm_chain_family,
                )
        if checkpoint is not None and state["current"] is not None:
            job = state["current"].get("job")
            if job is not None and _chain_job_supported(job, config):
                # Bounded schedule/coordinate reconstruction was checked by
                # unpack; a refused continuation must still pay that work.
                budget.consume(
                    2000 + job["n"].bit_length() + len(job["powers"])
                )
            if (
                job is not None
                and job["kind"] == "pm1"
                and not _legacy_pm1(config)
            ):
                from .pm1.gaps import verify_powers

                verify_powers(job, budget, config.pm1_gap_mode)
        while state["pending"] or state["current"] is not None:
            stage = (
                state["current"]["stage"] if state["current"] else "dispatch"
            )
            started = time.perf_counter()
            _advance(
                state,
                config,
                budget,
                context,
                generator,
                siqs_runtime,
                programs,
                policy,
            )
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
    # Live fallback state must reach the parent snapshot before packing.
    if siqs_runtime.get("job") is not None:
        stage = state["current"]["stage"]
        state["current"][stage + "_checkpoint"] = siqs_runtime[
            "job"
        ].checkpoint()

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
    )
