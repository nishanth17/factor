"""Explicit finite policies for the expanded useful-yield experiment."""

from dataclasses import dataclass, field

from ..common import utils
from .sieve_collector import CombinedYieldSieveConfig


@dataclass(frozen=True)
class ExtendedYieldSieveConfig(CombinedYieldSieveConfig):
    """Versioned extensions; production and earlier combined policies persist.

    Alternatives are eligible operations, not unconditional extra passes.
    Disk and threaded recovery require explicit resources from the caller.
    """

    extension_policy: str = field(default="extended-v1", kw_only=True)
    lazy_cycles: bool = field(default=True, kw_only=True)
    indexed_recovery: bool = field(default=True, kw_only=True)
    refinement_width: int = field(default=128, kw_only=True)
    outcome_feedback: bool = field(default=True, kw_only=True)
    incremental_bytes: int = field(default=8 * 2**20, kw_only=True)
    split_dispatch: bool = field(default=True, kw_only=True)
    split_cache_entries: int = field(default=256, kw_only=True)
    threshold_epochs: int = field(default=4, kw_only=True)
    threshold_exposure: int = field(default=1024, kw_only=True)
    smooth_batch_size: int = field(default=32, kw_only=True)
    opaque_pairs: bool = field(default=True, kw_only=True)
    opaque_conversions: int = field(default=4, kw_only=True)
    residual_arity: int = field(default=3, kw_only=True)
    cycle_reduction_passes: int = field(default=2, kw_only=True)
    cycle_reduction_choices: int = field(default=8, kw_only=True)
    recompute_prime_limit: int = field(default=13, kw_only=True)
    merge_weight: int = field(default=4, kw_only=True)
    merge_fill: int = field(default=0, kw_only=True)
    sparse_incidence: bool = field(default=True, kw_only=True)
    solver_panel_bits: int = field(default=4, kw_only=True)
    recovery_workers: int = field(default=0, kw_only=True)
    spill_directory: str | None = field(default=None, kw_only=True)
    spill_bytes: int = field(default=0, kw_only=True)
    spill_cache_bytes: int = field(default=2**20, kw_only=True)

    def __post_init__(self):
        # Validate before the inherited product-envelope exponentiation.
        if type(self.residual_arity) is not int or self.residual_arity not in (
            2,
            3,
        ):
            raise ValueError("invalid expanded residual arity")
        super().__post_init__()
        if self.extension_policy != "extended-v1":
            raise ValueError("unknown extension policy")
        for name in (
            "lazy_cycles",
            "indexed_recovery",
            "outcome_feedback",
            "split_dispatch",
            "opaque_pairs",
            "sparse_incidence",
        ):
            if type(getattr(self, name)) is not bool:
                raise TypeError(name + " must be Boolean")
        limits = {
            "split_cache_entries": (0, 4096),
            "incremental_bytes": (0, 64 * 2**20),
            "refinement_width": (0, 4096),
            "threshold_epochs": (0, 16),
            "threshold_exposure": (32, 2**20),
            "smooth_batch_size": (0, 256),
            "residual_arity": (2, 3),
            "opaque_conversions": (0, 16),
            "cycle_reduction_passes": (0, 4),
            "cycle_reduction_choices": (1, 32),
            "recompute_prime_limit": (0, 257),
            "merge_weight": (2, 8),
            "merge_fill": (0, 4096),
            "solver_panel_bits": (0, 8),
            "recovery_workers": (0, 1),
            "spill_bytes": (0, 4 * 2**30),
            "spill_cache_bytes": (65536, 64 * 2**20),
        }
        for name, (minimum, maximum) in limits.items():
            value = getattr(self, name)
            utils.require_integer(value, name, minimum)
            if value > maximum:
                raise ValueError(name + " exceeds its finite cap")
        if self.record_policy != "indexed-v2":
            raise ValueError("expanded storage requires indexed records")
        if self.cycle_policy != "short-v2":
            raise ValueError(
                "expanded incidence uses its own bounded reduction controls"
            )
        if self.recovery_policy != "checked-v2":
            raise ValueError("expanded arithmetic requires checked recovery")
        if bool(self.spill_bytes) != (self.spill_directory is not None):
            raise ValueError("spill requires a directory and finite disk cap")
        if self.spill_directory is not None and (
            not isinstance(self.spill_directory, str)
            or not self.spill_directory.startswith("/")
        ):
            raise ValueError("spill directory must be an absolute path")


EXTENSION_FIELDS = frozenset(
    ExtendedYieldSieveConfig.__dataclass_fields__
) - frozenset(CombinedYieldSieveConfig.__dataclass_fields__)
