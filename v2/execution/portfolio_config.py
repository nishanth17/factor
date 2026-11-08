"""Portfolio configuration and ECM metadata without eager relation engines."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

from .. import constants
from ..common import arithmetic, utils
from ..common.arithmetic import isqrt
from ..ecm.programs import ECMPrograms
from .allocation import ECMAllocation
from .resources import removed_arguments
from .schedules import ScheduleCache, SieveContext

if TYPE_CHECKING:
    from ..qs import SIQSConfig
    from ..qs.sss import SSSConfig

_AUTO_SIQS = object()


def _decode_siqs_config(values):
    """Rebuild a concrete SIQS child from an asdict-style mapping."""
    from ..qs import (
        CombinedYieldSieveConfig,
        CountedSIQSConfig,
        DoubleLargeSieveConfig,
        GrowingDoubleLargeSieveConfig,
        SieveConfig,
        SIQSConfig,
        YieldSieveConfig,
    )
    from ..qs.extended_policy import (
        BoundedYieldSieveConfig,
        ExtendedYieldSieveConfig,
    )

    decoded = dict(values)
    collector = decoded["collector"]
    if collector.get("extension_policy") == "extended-v2":
        collector_type = BoundedYieldSieveConfig
    elif collector.get("extension_policy") == "extended-v1":
        collector_type = ExtendedYieldSieveConfig
    elif collector.get("yield_policy") == "yield-v2":
        collector_type = CombinedYieldSieveConfig
    elif "yield_policy" in collector or "small_prime_bound" in collector:
        collector_type = YieldSieveConfig
    elif collector.get("storage_policy") == "grow":
        collector_type = GrowingDoubleLargeSieveConfig
    elif "large_prime_bound" in collector:
        collector_type = DoubleLargeSieveConfig
    else:
        collector_type = SieveConfig
    decoded["collector"] = collector_type(**collector)
    config_type = CountedSIQSConfig if "base_count" in decoded else SIQSConfig
    return config_type(**decoded)


@removed_arguments(("memory_bytes",))
@dataclass(frozen=True)
class PortfolioConfig:
    """Finite allowances; ECM tiers are (inclusive B1, B2, curve count).

    Omitted ECM tiers resolve once from the original input via for_input().
    Explicit tiers and saved numerical schedules retain their identity.
    memory_bytes bounds conservative owned workspace estimates, not RSS.
    """

    allocation: ECMAllocation | None = field(default=None, kw_only=True)
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
    ecm_tiers: tuple | None = None
    chunk_size: int = 16
    gcd_batch: int = constants.GCD_BATCH_SIZE
    primality_rounds: int = constants.PRIMALITY_ROUNDS
    fermat_steps: int = 0
    segment_size: int = 1024
    trial_chunk: int = 64
    schedule_cache_bytes: int = 0
    max_input_bits: int = 4096
    trace_limit: int = 256
    rolling: bool = False
    siqs: SIQSConfig | None = _AUTO_SIQS
    sss: SSSConfig | None = None
    ecm_program_bytes: int = 0
    ecm_chain_mode: str = field(default="auto", kw_only=True)
    ecm_chain_bytes: int = field(default=0, kw_only=True)
    ecm_chain_family: str = field(default="auto", kw_only=True)
    ecm_pair_distance: int | None = None
    ecm_pair_wheel: int | None = None

    def __post_init__(self):
        """Validate bounds and reserve storage before any allocation."""
        if isinstance(self.siqs, dict):
            object.__setattr__(self, "siqs", _decode_siqs_config(self.siqs))

        # Routing intent is transient; only concrete numerical tiers enter
        # asdict/checkpoints. Reconstructing a saved config never reselects.
        object.__setattr__(
            self, "_automatic_ecm_tiers", self.ecm_tiers is None
        )
        if self.ecm_tiers is None:
            object.__setattr__(
                self,
                "ecm_tiers",
                ((constants.ECM_B1, constants.ECM_B2, 32),),
            )

        automatic_siqs = self.siqs is _AUTO_SIQS
        object.__setattr__(
            self,
            "_automatic_siqs",
            automatic_siqs and self._automatic_ecm_tiers and self.sss is None,
        )
        if automatic_siqs:
            object.__setattr__(self, "siqs", None)
        object.__setattr__(self, "_automatic_memory", True)

        if self.allocation is not None:
            if not isinstance(self.allocation, ECMAllocation):
                raise TypeError("allocation must be an ECMAllocation or None")
            if not (self.siqs or self.sss) and any(
                (
                    self.allocation.fallback_work,
                    self.allocation.fallback_seconds,
                    self.allocation.fallback_cpu_seconds,
                )
            ):
                raise ValueError(
                    "fallback reservations require a relation engine"
                )
        automatic_memory = True
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
            if not self.ecm_program_bytes:
                raise ValueError(
                    "ECM chains require bounded reusable programs"
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
            from ..ecm.chain_options import default_options

            # Store the resolved policy, never an auto decision, so resume
            # retains its executor even if future defaults change again.
            for name, value in default_options(
                self, tiers, automatic_memory=automatic_memory
            ).items():
                object.__setattr__(self, name, value)
            if self.ecm_chain_mode == "reuse":
                utils.require_integer(
                    self.ecm_chain_bytes, "ecm_chain_bytes", 0
                )
            elif self.ecm_chain_bytes:
                raise ValueError("ECM chain memory requires reuse mode")

        if self.ecm_program_bytes:
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
            from ..qs.sss import SSSConfig

            if not isinstance(self.sss, SSSConfig):
                raise TypeError("sss must be an SSSConfig or None")

        if self.siqs is not None:
            from ..qs import SIQSConfig

            if not isinstance(self.siqs, SIQSConfig):
                raise TypeError("siqs must be a SIQSConfig or None")

    def for_input(self, n):
        """Resolve the same size-based numerical and cache recipes once."""
        utils.require_integer(n)
        if not self._automatic_ecm_tiers:
            return self
        from ..ecm.defaults import resolve_defaults
        from ..qs.defaults import default_siqs_config, uses_siqs_defaults

        selected = resolve_defaults(self, n)
        if self._automatic_siqs and uses_siqs_defaults(
            n, backend=self.backend
        ):
            return replace(
                selected, siqs=default_siqs_config(n, backend=self.backend)
            )
        return selected

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


def _chain_bounds(config):
    """Keep fresh/unsupported schedules on the measured B4 ladder."""
    if config.ecm_chain_mode == "off" or config.chunk_size != 16:
        return ()
    return tuple(
        b1 for b1, _, curves in config.ecm_tiers if b1 == 2000 and curves >= 8
    )


def _execution_context(state, config, budget):
    """Prepare finite schedules only when the selected policy will use them."""
    if config.allocation is not None or not state.get("context_ready"):
        budget.consume((isqrt(config.max_hi - 1) + 1) // 2)
    context = SieveContext(
        config.max_hi, segment_size=config.segment_size, rolling=config.rolling
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
            from ..ecm.chains import ChainPlans

            programs.chains = ChainPlans(
                config.ecm_chain_bytes,
                config.backend,
                _chain_bounds(config),
            )
        else:
            from ..ecm.chain_options import ChainPlans

            programs.chains = ChainPlans(
                config.ecm_chain_bytes,
                config.backend,
                _chain_bounds(config),
                config.ecm_chain_family,
            )
    return context, programs


# Keep the public import/pickle identity while sharing one config definition.
PortfolioConfig.__module__ = "v2.portfolio"
