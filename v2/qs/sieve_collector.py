"""Bounded conservative score sieving and single-large-prime matching."""

from array import array
from dataclasses import dataclass, field

from ..common import utils
from ..common.arithmetic import gcd, isqrt
from ..execution.budget import BudgetExhaustedError
from ..execution.resources import (
    accounted_job,
    removed_arguments,
    resolve_budget,
)
from ..execution.work_budget import PollingBudget
from . import collector_kernels as kernels
from .families import verify_polynomial_roots
from .hit_words import MAX_EXTENDED_HIT_BYTES, MAX_HIT_BYTES, WordHits
from .polynomial import (
    MAX_COEFFICIENT_BITS,
    MAX_POSITION_BITS,
    checked_position,
    polynomial_roots,
)
from .power_reuse import (
    MAX_EXTENDED_TRANSFER_BYTES,
    MAX_TRANSFER_BYTES,
    RootTransfers,
)
from .power_sieve import HenselInverses, prime_power_roots
from .relations import (
    AtomicRelation,
    _combination_workspace,
    checked_residual_bound,
    combine_relations,
    combined_storage_reserve,
    verify_atomic,
)
from .score import log_bounds
from .storage import (
    RelationPublicationError,
    require_published,
    storage_limit,
)

MAX_BLOCK_WIDTH = 4096
MAX_WINDOW_WIDTH = 1_000_000
MAX_STORED_ATOMS = 65536
MAX_METADATA_CHUNK = 1024
MAX_SHARED_INVERSE_BYTES = 8 * 2**20
MAX_WIDE_INVERSE_BYTES = 32 * 2**20


@removed_arguments(
    ("max_atoms", "max_relations", "memory_bytes"),
    positional=(10, "max_relations"),
)
@dataclass(frozen=True)
class SieveConfig:
    """Finite collector controls; threshold_extra > 0 is intentionally lossy.

    Scores upper-bound factor-base log contributions, including all prime
    powers. Zero extra guarantees candidate coverage, not admission: exact
    division and the verifier decide admission. bytearray scores saturate at
    255; thresholds saturate too, allowing extra false candidates, never
    missed admissible values. Cutoffs omit marking only, never exact recovery.
    """

    block_width: int = 256
    metadata_chunk: int = 64
    score_backend: str = "list"
    marking: str = "sparse"
    division: str = "roots"
    small_prime_cutoff: int = 0
    threshold_extra: int = 0
    score_policy: str = "adaptive"
    sparse_queue_bytes: int = field(default=0, kw_only=True)
    power_plan_bytes: int = field(default=0, kw_only=True)
    residual_bound: int = 500
    preparation_cache_bytes: int = field(default=2 * 2**20, kw_only=True)
    target_rows: int | None = field(default=1024, kw_only=True)
    max_partials: int = 256

    def __post_init__(self):
        """Validate bounded choices before allocating buffers or roots."""
        limits = {
            "block_width": (1, MAX_BLOCK_WIDTH),
            "metadata_chunk": (1, MAX_METADATA_CHUNK),
            "small_prime_cutoff": (0, 100000),
            "threshold_extra": (0, 65536),
            "max_partials": (0, MAX_STORED_ATOMS),
            "power_plan_bytes": (0, 16 * 1024 * 1024),
            "sparse_queue_bytes": (0, 32 * 1024 * 1024),
        }
        for name, (minimum, maximum) in limits.items():
            value = getattr(self, name)
            if getattr(self, "storage_policy", None) == "grow" and name in (
                "max_partials",
            ):
                storage_limit(value)
                continue
            utils.require_integer(value, name, minimum)
            if value > maximum:
                raise ValueError(f"{name} exceeds the collector limit")

        checked_residual_bound(self.residual_bound)
        if self.score_backend not in ("list", "bytearray", "array"):
            raise ValueError("unknown score backend")
        if not self.marking:
            raise ValueError("unknown marking policy")
        object.__setattr__(self, "_kernels", kernels.components(self.marking))
        if self.division not in ("full", "roots", "bucket", "resieve"):
            raise ValueError("unknown division policy")
        if self.score_policy not in (
            "adaptive",
            "conservative",
            "candidate",
            "powers",
            "fixed",
        ):
            raise ValueError("unknown score policy")

    @property
    def partial_limit(self):
        return storage_limit(self.max_partials)


@removed_arguments(
    ("memory_bytes", "max_atoms", "max_relations"),
    positional=(10, "max_relations"),
)
@dataclass(frozen=True)
class DoubleLargeSieveConfig(SieveConfig):
    """Opt-in two-prime bounds and a finite cumulative splitting quota.

    candidate_bound=0 uses the full admissible product domain. A smaller
    explicit allowance intentionally loses DLP candidates while retaining
    conservative SLP coverage. Default SieveConfig serialization is unchanged.
    """

    large_prime_bound: int = field(kw_only=True)
    large_product_bound: int = field(kw_only=True)
    candidate_bound: int = field(default=0, kw_only=True)
    split_call_limit: int = field(default=131072, kw_only=True)

    def __post_init__(self):
        super().__post_init__()
        checked_residual_bound(self.large_prime_bound)
        utils.require_integer(self.large_prime_bound, "large_prime_bound", 2)
        utils.require_integer(
            self.large_product_bound, "large_product_bound", 4
        )
        if self.large_product_bound > self.large_prime_bound ** getattr(
            self, "residual_arity", 2
        ):
            raise ValueError("product bound exceeds the endpoint envelope")
        utils.require_integer(self.candidate_bound, "candidate_bound", 0)
        if self.candidate_bound and not (
            self.residual_bound
            <= self.candidate_bound
            <= max(self.residual_bound, self.large_product_bound)
        ):
            raise ValueError("candidate bound must preserve SLP coverage")
        if getattr(self, "storage_policy", None) == "grow":
            storage_limit(self.split_call_limit)
        else:
            utils.require_integer(self.split_call_limit, "split_call_limit", 0)
            if self.split_call_limit > 1048576:
                raise ValueError("split-call limit exceeds the collector cap")

    @property
    def split_limit(self):
        return storage_limit(self.split_call_limit)


@removed_arguments(
    ("max_atoms", "max_relations", "memory_bytes"),
    positional=(10, "max_relations"),
)
@dataclass(frozen=True)
class GrowingDoubleLargeSieveConfig(DoubleLargeSieveConfig):
    """Retain partial history with allocator-grown containers by default.

    None removes a storage or cumulative split-count quota. Integer limits
    remain explicit opt-ins; work/time and per-residual effort are unchanged.
    The policy marker distinguishes this execution from legacy capped DLP.
    """

    target_rows: int | None = field(default=None, kw_only=True)
    storage_policy: str = field(default="grow", kw_only=True)
    max_partials: int | None = None
    split_call_limit: int | None = field(default=None, kw_only=True)

    def __post_init__(self):
        if self.storage_policy != "grow":
            raise ValueError("unknown growing-storage policy")
        super().__post_init__()


@removed_arguments(
    ("memory_bytes", "max_atoms", "max_relations"),
    positional=(10, "max_relations"),
)
@dataclass(frozen=True)
class YieldSieveConfig(GrowingDoubleLargeSieveConfig):
    """Independent useful-yield policies with checked original provenance.

    For proven endpoints p <= q, require p <= small_prime_bound and
    q <= large_prime_bound, including p == q. The product ceiling remains
    binding. This changes admission only; it never widens reporting itself.
    """

    small_prime_bound: int | None = field(default=None, kw_only=True)
    yield_policy: str = field(default="yield-v1", kw_only=True)
    preparation_policy: str = field(default="fixed", kw_only=True)
    cycle_policy: str = field(default="forest", kw_only=True)
    cycle_search_edges: int = field(default=4096, kw_only=True)
    duplicate_policy: str = field(default="identity", kw_only=True)
    duplicate_bytes: int = field(default=16 * 2**20, kw_only=True)
    record_policy: str = field(default="objects", kw_only=True)

    def __post_init__(self):
        super().__post_init__()
        if self.yield_policy not in ("yield-v1", "yield-v2"):
            raise ValueError("unknown useful-yield policy")
        if self.yield_policy == "yield-v1" and any(
            getattr(self, name).endswith("v2")
            for name in (
                "preparation_policy",
                "cycle_policy",
                "duplicate_policy",
                "record_policy",
            )
        ):
            raise ValueError("combined policies require yield-v2")
        if self.yield_policy == "yield-v2" and not isinstance(
            self, CombinedYieldSieveConfig
        ):
            raise ValueError("yield-v2 requires the combined configuration")
        if self.small_prime_bound is not None:
            utils.require_integer(
                self.small_prime_bound, "small_prime_bound", 2
            )
            if self.small_prime_bound > self.large_prime_bound:
                raise ValueError("smaller endpoint bound exceeds larger bound")
        for name, choices in (
            ("preparation_policy", ("fixed", "target-v1", "active-v2")),
            ("cycle_policy", ("forest", "short-v1", "short-v2")),
            ("duplicate_policy", ("identity", "pair-v1", "congruence-v2")),
            ("record_policy", ("objects", "packed-v1", "indexed-v2")),
        ):
            if getattr(self, name) not in choices:
                raise ValueError("unknown " + name)
        for name, maximum in (
            ("cycle_search_edges", 65536),
            ("duplicate_bytes", 64 * 2**20),
        ):
            utils.require_integer(getattr(self, name), name, 0)
            if getattr(self, name) > maximum:
                raise ValueError(name + " exceeds its finite cap")


@removed_arguments(
    ("memory_bytes", "max_atoms", "max_relations"),
    positional=(10, "max_relations"),
)
@dataclass(frozen=True)
class CombinedYieldSieveConfig(YieldSieveConfig):
    """Versioned combined challenger; every auxiliary path has a finite cap."""

    yield_policy: str = field(default="yield-v2", kw_only=True)
    preparation_policy: str = field(default="active-v2", kw_only=True)
    cycle_policy: str = field(default="short-v2", kw_only=True)
    duplicate_policy: str = field(default="congruence-v2", kw_only=True)
    record_policy: str = field(default="indexed-v2", kw_only=True)
    recovery_policy: str = field(default="checked-v2", kw_only=True)
    split_policy: str = field(default="portfolio-v2", kw_only=True)
    batch_policy: str = field(default="shared-v2", kw_only=True)
    batch_size: int = field(default=64, kw_only=True)
    checked_bytes: int = field(default=8 * 2**20, kw_only=True)
    filter_policy: str = field(default="light-v2", kw_only=True)
    extraction_policy: str = field(default="joint-v2", kw_only=True)

    def __post_init__(self):
        super().__post_init__()
        if self.yield_policy != "yield-v2":
            raise ValueError("combined policies require yield-v2")
        for name, choices in (
            ("recovery_policy", ("reference", "checked-v2")),
            (
                "split_policy",
                (
                    "reference",
                    "rho-v2",
                    "portfolio-v2",
                    "ecm-v2",
                    "rho-first-v3",
                ),
            ),
            ("batch_policy", ("off", "known-v2", "shared-v2")),
            ("filter_policy", ("reference", "fill-v2", "light-v2")),
            ("extraction_policy", ("reference", "joint-v2")),
        ):
            if getattr(self, name) not in choices:
                raise ValueError("unknown " + name)
        for name, low, high in (
            ("batch_size", 2, 256),
            ("checked_bytes", 0, 64 * 2**20),
        ):
            utils.require_integer(getattr(self, name), name, low)
            if getattr(self, name) > high:
                raise ValueError(name + " exceeds its finite cap")


@dataclass(frozen=True)
class SieveResult:
    """Cumulative checked store and per-call diagnostics for [lo,hi).

    atoms includes full/combined provenance and unmatched partials. FIFO
    eviction touches unmatched atoms only. next_position is the first
    uncommitted position; resuming on the same collector retains matches.
    Returned snapshots own references; release old snapshots before resuming
    if their memory is to be excluded from the collector's reservation.
    """

    atoms: tuple
    full_relations: tuple
    combined_relations: tuple
    partial_ids: tuple
    divisor: int | None
    next_position: int
    reason: str
    stats: dict
    workspace_bytes: int
    resource_stop: dict | None = None


@dataclass(frozen=True)
class _CollectionProgress:
    """Internal status without references to the cumulative relation store."""

    divisor: int | None
    next_position: int
    reason: str
    stats: dict
    workspace_bytes: int


class SieveCollector:
    """Cache one polynomial's roots/logs and reuse bounded working blocks.

    Constructor setup can raise ValueError, MemoryError or budget exhaustion.
    Collection refusals return a checked prefix. Single-thread use only;
    state is in-memory, not the serialized SIQS checkpoint of P3.4.
    """

    @property
    def budget(self):
        return self._budget

    @budget.setter
    def budget(self, value):
        self._budget = value
        arena = getattr(self, "_arena", None)
        if arena is not None and arena.spill is not None:
            # Reads of old cold records belong to the resumed allowance too.
            arena.spill.budget = value

    def __init__(
        self,
        polynomial,
        factor_base,
        *,
        config=None,
        budget=None,
        precomputed_roots=None,
    ):
        """Reserve metadata/buffers and validate A's factor-base support."""
        self.config = config if config is not None else SieveConfig()
        if (
            isinstance(self.config, CombinedYieldSieveConfig)
            and type(factor_base.n) is not int
        ):
            raise ValueError(
                "combined useful-yield policies require native integers"
            )
        double = isinstance(self.config, DoubleLargeSieveConfig)
        if double and type(self) is not SieveCollector:
            raise ValueError("DLP supports the serial SieveCollector only")
        self._graph = None
        self._short_cycles = None
        if double:
            # Keep ordinary SLP and owned historic loader dependencies intact.
            from .large_prime_collector import (
                COUNTERS,
                admit_large_relation,
                divide_large_residual,
            )
            from .large_primes import LargePrimeForest, graph_reserve

            if isinstance(self.config, GrowingDoubleLargeSieveConfig):
                self._graph = LargePrimeForest()
                self._graph_reserve = graph_reserve
            else:
                self._graph = LargePrimeForest()
                self._graph_reserve = graph_reserve
            self._divide_double = divide_large_residual
            self._admit_double = admit_large_relation
            self._double_counters = COUNTERS
        self._split_calls = 0
        self._candidate_bound = self.config.residual_bound
        if double:
            self._candidate_bound = self.config.candidate_bound or max(
                self.config.residual_bound, self.config.large_product_bound
            )
        self.polynomial, self.factor_base = polynomial, factor_base
        self.budget = resolve_budget(budget)
        if (polynomial.n, polynomial.multiplier) != (
            factor_base.n,
            factor_base.multiplier,
        ):
            raise ValueError("polynomial and factor-base identity mismatch")

        count = len(factor_base.entries)
        # Include roots, Python object overhead, score integers, hit bitsets,
        # candidate temporaries, and two bounded translation slices. Account
        # for combination scratch before retaining any atom.
        self._workspace = factor_base.workspace_bytes + 32768 + count * 512
        self._workspace += self.config.block_width * (128 + (count + 7) // 8)
        self._workspace += 128 * (polynomial.n_prime.bit_length() + 16384)
        self._workspace += self.config.power_plan_bytes
        self._workspace += kernels.workspace(
            self.config.marking, self.config.block_width, count
        )
        if double:
            self._workspace += self._graph_reserve(0)
        if not getattr(self.config, "extension_policy", None) and getattr(
            self.config, "cycle_policy", "forest"
        ) in (
            "short-v1",
            "short-v2",
        ):
            from .short_cycles import ShortCycleIndex

            self._short_cycles = ShortCycleIndex(
                self.config.cycle_search_edges,
                strategy=self.config.cycle_policy,
            )
            self._workspace += self._short_cycles.scratch_reserve
            self._workspace += self._short_cycles.edge_reserve(0)
        if self.config.score_policy in ("powers", "fixed"):
            # Two capped root lists coexist while lifting; coefficients and
            # positions have finite bit bounds even far from zero.
            self._workspace += 128 * (
                128 + 4 * (MAX_COEFFICIENT_BITS + 2 * MAX_POSITION_BITS + 4)
            )
            self._workspace += 512 * count

        if precomputed_roots is not None:
            verify_polynomial_roots(
                polynomial, factor_base, precomputed_roots, budget=self.budget
            )
        remaining, a_exponents, roots, logs = (
            polynomial.supported_a,
            [],
            [],
            [],
        )

        for index, entry in enumerate(factor_base.entries):
            self.budget.consume(polynomial.a.bit_length() + 1)
            exponent = 0
            while remaining % entry.prime == 0:
                remaining //= entry.prime
                exponent += 1
            a_exponents.append(exponent)
            roots.append(
                precomputed_roots[index]
                if precomputed_roots is not None
                else polynomial_roots(
                    polynomial, factor_base, entry, budget=self.budget
                )
            )
            # ceil(log2 p) uses exact bit lengths, with 2 exactly one bit.
            logs.append(
                log_bounds(entry.prime)[1]
                if self.config.score_policy == "fixed"
                else (entry.prime - 1).bit_length()
            )

        if remaining != 1:
            raise ValueError("A must factor completely over the factor base")
        self._roots, self._logs = tuple(roots), tuple(logs)
        self._a_exponents = tuple(a_exponents)
        self._a_support = tuple(
            i for i, value in enumerate(a_exponents) if value
        )
        self._a_factor_value = polynomial.supported_a
        self._a_factor_primes = factor_base.primes
        self._skipped = tuple(
            i
            for i, roots in enumerate(self._roots)
            if roots.prime < self.config.small_prime_cutoff
        )
        width = self.config.block_width
        if self.config.score_backend == "bytearray":
            self._scores = bytearray(width)
        elif self.config.score_backend == "array":
            self._scores = array("I", [0]) * width
        else:
            self._scores = [0] * width

        self._hits = [0] * width
        self._group_lease = None
        self._resieved = {}
        self._scratch_bytes = 0
        self._scratch_peak_bytes = 0
        self._omitted_allowance = 0
        self._lift_inverses = {}
        from .power_sieve_native import prime_power_roots as native_power_roots

        self._power_roots = (
            prime_power_roots
            if "pruned_transfers" in self.config._kernels
            else native_power_roots
        )
        if "presieve_powers" in self.config._kernels:
            self._presieve = kernels.PowerPresieve()
        elif "presieve" in self.config._kernels:
            self._presieve = kernels.Presieve()
        else:
            self._presieve = None
        self._power_plans = {}
        self._plan_refused = set()
        self._plan_bytes = 0
        self._plan_window = None
        self._atoms, self._pending = {}, {}
        self._full, self._combined, self._rows = [], [], []
        self._atom_bytes = {}
        self._checked = self._arena = self._endpoints = None
        if isinstance(self.config, CombinedYieldSieveConfig):
            from .cofactor_policy import KnownEndpoints
            from .indexed_relations import RecordArena
            from .yield_collector import divide_yield_residual
            from .yield_support import CheckedStore, residual_floor

            self._workspace += self.config.checked_bytes + 65536
            if self.config.duplicate_policy == "congruence-v2":
                self._workspace += self.config.duplicate_bytes // 2
            # Lookahead norms, exponents and product/remainder temporaries
            # coexist with scoring buffers and retained graph/row ownership.
            self._batch_reserve = (
                self.config.batch_size
                * (
                    8192
                    + 256 * len(factor_base.entries)
                    + 16 * (MAX_COEFFICIENT_BITS + 2 * MAX_POSITION_BITS + 4)
                )
                if self.config.batch_policy == "shared-v2"
                else 0
            )
            if (
                getattr(self.config, "extension_policy", None) == "extended-v2"
                and self.config.batch_policy == "shared-v2"
            ):
                from .yield_memory import base_batch_reserve

                self._batch_base = base_batch_reserve(
                    factor_base, self.config.batch_size, self.budget
                )
                self._batch_reserve = self._batch_base
            self._workspace += self._batch_reserve
            self._checked = CheckedStore(factor_base, self.config)
            self._endpoints = KnownEndpoints()
            if self.config.record_policy == "indexed-v2":
                self._arena = RecordArena(
                    factor_base, config=self.config, budget=self.budget
                )
            self._divide_double = divide_yield_residual
            self._residual_floor = (
                residual_floor(
                    factor_base,
                    self.budget,
                )
                if self.config.recovery_policy == "checked-v2"
                else 2
            )
            self._incoming_key = None
            if getattr(self.config, "extension_policy", None):
                from .cofactor_policy import SplitMemo
                from .incremental_structure import IncrementalStructure
                from .residual_kernel import ResidualKernel
                from .yield_feedback import ReportingController

                self._workspace += (
                    4096 * self.config.split_cache_entries + 2**21
                )
                self._workspace += (
                    self.config.spill_cache_bytes
                    if (self.config.spill_directory is not None)
                    else 0
                )
                self._workspace += self.config.incremental_bytes
                if self.config.extension_policy == "extended-v2":
                    from .opaque_staging import OpaqueStaging

                    self._opaque_staging = OpaqueStaging(
                        self.config, factor_base
                    )
                    self._workspace += self.config.opaque_cache_bytes
                    self._slab_bytes = 0
                self._workspace += 16 * 2**20 * self.config.recovery_workers
                self._split_memo = SplitMemo(self.config.split_cache_entries)
                self._conversion_attempts = 0
                self._conversion_events = []
                self._conversion_bytes = 0
                self._incremental = IncrementalStructure(
                    self.config.incremental_bytes
                )
                self._graph = ResidualKernel()
                self._short_cycles = None
                self._reporting = ReportingController(
                    self._candidate_bound,
                    max(
                        self._candidate_bound, self.config.large_product_bound
                    ),
                    self.config.threshold_exposure,
                    self.config.threshold_epochs,
                )

        # Optional derived helpers compete for the same owned allowance.
        # Their refusal preserves the original exact arithmetic/shift paths.
        self._shared_lift_inverses = None
        self._reuse_workspace = 0
        inverse_bytes = 256 + 512 * count
        wide_inverses = "wide_inverses" in self.config._kernels
        inverse_limit = (
            MAX_WIDE_INVERSE_BYTES
            if wide_inverses
            else MAX_SHARED_INVERSE_BYTES
        )
        if (
            self.config.score_policy in ("powers", "fixed")
            and inverse_bytes <= inverse_limit
        ):
            self._workspace += inverse_bytes
            self._reuse_workspace += inverse_bytes
            self._shared_lift_inverses = HenselInverses(count)
        self._hit_masks = ()
        if self.config.division in ("bucket", "resieve"):
            cap = 4 * 2**20
            masks, reserved = [], 256
            for index in range(count):
                size = 120 + 8 * ((index + 63) // 64)
                if reserved + size > cap:
                    break
                self.budget.consume(1 + (index + 63) // 64)
                masks.append(1 << index)
                reserved += size
            if masks:
                self._hit_masks = tuple(masks)
                # Include the temporary list during construction as well.
                self._workspace += reserved
                self._reuse_workspace += reserved
        self._power_transfers = None
        wide_transfers = "wide_transfers" in self.config._kernels
        transfer_limit = (
            MAX_EXTENDED_TRANSFER_BYTES
            if wide_transfers
            else MAX_TRANSFER_BYTES
        )
        # Multiple powers of one prime own separate cached pairs. A larger
        # experiment therefore funds more than one entry per base member;
        # the default capacity and refusal behavior remain unchanged.
        transfer_bytes = min(
            transfer_limit, 65536 + (8192 if wide_transfers else 4096) * count
        )
        if self.config.score_policy in ("powers", "fixed"):
            # Leave the mandatory old/new polynomial metadata coexistence
            # funded; an optional cache must not prevent the next switch.
            self._workspace += transfer_bytes
            self._reuse_workspace += transfer_bytes
            self._power_transfers = RootTransfers(
                transfer_bytes,
                capacity_limit=transfer_limit,
                pruned_replay="pruned_transfers" in self.config._kernels,
            )
            self._power_transfers.bind(polynomial)
        self._word_hits = None
        hit_bytes = WordHits.reservation(count, width)
        hit_limit = (
            MAX_EXTENDED_HIT_BYTES
            if "wide_hits" in self.config._kernels
            else MAX_HIT_BYTES
        )
        if (
            count > 60
            and self.config.score_backend == "list"
            and self.config.division == "bucket"
            and self.config.score_policy in ("powers", "fixed")
            and hit_bytes <= hit_limit
        ):
            # Existing helpers keep priority. Construction and lazy reads
            # have explicit extra work; scarce memory keeps the list path.
            word_type = WordHits
            if "byte_stamps" in self.config._kernels:
                from .byte_hit_words import ByteStampWordHits

                word_type = ByteStampWordHits
            self._word_hits = word_type(
                count, self._hits, self.budget, byte_limit=hit_limit
            )
            self._hits = self._word_hits
            self._workspace += hit_bytes
            self._reuse_workspace += hit_bytes
        self._sparse_split = None
        if self.config.sparse_queue_bytes:
            target = polynomial.n_prime
            self.budget.consume(count + target.bit_length() ** 2 + 1)
            self._sparse_split = next(
                (
                    i
                    for i, entry in enumerate(factor_base.entries)
                    if entry.prime >= self.config.block_width
                ),
                count,
            )
            self._sparse_singular = any(
                target % entry.prime == 0
                for entry in factor_base.entries[self._sparse_split :]
            )
            self._sparse_square = isqrt(target) ** 2 == target

        # Bind optional implementations once. Default block/root traversals
        # must not repeatedly test disabled experimental mechanisms.
        from types import MethodType

        from . import scalar_kernels
        from .sparse_queue import sieve as optional_sparse_sieve

        if "presieve_powers" not in self.config._kernels:
            self._sieve_direct = MethodType(scalar_kernels.direct_scalar, self)
        if "scan_fused" not in self.config._kernels:
            self._collect = MethodType(scalar_kernels.collect_scalar, self)
        self._sparse_sieve = (
            scalar_kernels.sieve_native
            if "packed_queue" not in self.config._kernels
            else optional_sparse_sieve
        )

    def _recovery_mask(self, index):
        """Reuse a funded immutable prefix; larger indices stay exact."""
        masks = self._hit_masks
        return masks[index] if index < len(masks) else 1 << index

    @property
    def matrix_relations(self):
        """Return an immutable prefix in mixed full/matched admission order."""
        return tuple(self._rows)

    @property
    def partial_ids(self):
        """Unowned atoms only; emitted cycle provenance stays pinned."""
        if getattr(self.config, "extension_policy", None):
            return tuple(
                self._arena.atoms[i].relation_id for i in self._graph.unowned
            )
        return (
            tuple(self._graph.unowned)
            if self._graph is not None
            else tuple(self._pending.values())
        )

    def set_polynomial(self, polynomial, roots, *, retain_relations=True):
        """Switch metadata atomically; retain verified cross-family rows.

        The target/base are fixed. Reserve coexistence with old metadata before
        rebuilding; budget or memory refusal leaves the published state intact.
        """
        require_published(self)
        if type(retain_relations) is not bool:
            raise TypeError("retain_relations must be Boolean")
        extra = 512 * len(self.factor_base.entries) + 32768
        verify_polynomial_roots(
            polynomial, self.factor_base, roots, budget=self.budget
        )
        supported_a = polynomial.supported_a
        primes = self.factor_base.primes
        if (
            supported_a == self._a_factor_value
            and primes == self._a_factor_primes
        ):
            # Gray steps change B but retain A. The immutable exponent tuple
            # is already checked against exactly this supported A and base.
            self.budget.consume(1)
            exponents, support = self._a_exponents, self._a_support
        else:
            remaining, exponents = supported_a, []
            for entry in self.factor_base.entries:
                self.budget.consume(polynomial.a.bit_length() + 1)
                exponent = 0
                while remaining % entry.prime == 0:
                    remaining //= entry.prime
                    exponent += 1
                exponents.append(exponent)

            if remaining != 1:
                raise ValueError(
                    "A must factor completely over the factor base"
                )
            exponents = tuple(exponents)
            support = tuple(i for i, value in enumerate(exponents) if value)
        self._switch_peak = max(
            getattr(self, "_switch_peak", 0), self._workspace + extra
        )
        self.budget.publication()
        try:
            if not retain_relations:
                if self._checked is not None:
                    from .indexed_relations import RecordArena

                    self._checked.clear()
                    if self._arena is not None:
                        # Old public snapshots keep their original arena alive.
                        spill = self._arena.spill
                        if (
                            getattr(self.config, "extension_policy", None)
                            == "extended-v2"
                        ):
                            self._workspace -= self._slab_bytes
                            self._slab_bytes = 0
                        self._arena = RecordArena(
                            self.factor_base,
                            config=self.config,
                            budget=self.budget,
                        )
                        self._arena.spill = spill
                    if getattr(self.config, "extension_policy", None):
                        from .incremental_structure import IncrementalStructure

                        self._incremental = IncrementalStructure(
                            self.config.incremental_bytes
                        )
                        self._conversion_events = []
                        self._workspace -= self._conversion_bytes
                        self._conversion_bytes = 0
                        if self.config.extension_policy == "extended-v2":
                            from .opaque_staging import OpaqueStaging

                            self._opaque_staging = OpaqueStaging(
                                self.config, self.factor_base
                            )
                cache = getattr(self, "_preparation_cache", None)
                if cache is not None:
                    cache.clear()
                self._workspace -= sum(self._atom_bytes.values())
                if self._short_cycles is not None:
                    self._workspace -= 2048 * len(self._short_cycles.edges)
                    self._short_cycles.edges.clear()
                    self._short_cycles.neighbors.clear()
                    self._short_cycles.weights.clear()
                if self._graph is not None:
                    if not getattr(self.config, "extension_policy", None):
                        self._workspace -= self._graph_reserve(
                            len(self._graph.edges)
                        ) - self._graph_reserve(0)
                    self._graph = type(self._graph)()
                self._atoms.clear()
                self._pending.clear()
                self._atom_bytes.clear()
                self._full.clear()
                self._combined.clear()
                self._rows.clear()

            if self.config.sparse_queue_bytes:
                from .sparse_queue import clear

                clear(self)
            self.polynomial, self._roots = polynomial, roots
            if self._power_transfers is not None:
                self._power_transfers.bind(polynomial)
            self._a_exponents, self._a_support = exponents, support
            self._a_factor_value, self._a_factor_primes = supported_a, primes
            self._lift_inverses.clear()
            self._power_plans.clear()
            self._plan_refused.clear()
            self._plan_bytes = 0
            self._plan_window = None
            self._resieved.clear()
        except MemoryError as error:
            self._publication_failed = True
            raise RelationPublicationError(
                "polynomial metadata publication", "sieve_collector"
            ) from error

    def _power_marks(self, index, maximum, lo, hi, stats):
        """Reuse a charged, capped plan bound to this polynomial/window.

        The whole collect interval bounds valuations; each replay intersects
        marks with its working block. A refusal falls back to block streaming.
        Construction scratch retains the existing two-root-list reservation;
        each retained tuple/list entry is checked before allocation.
        """
        roots, weight = self._roots[index], self._logs[index]
        cap = self.config.power_plan_bytes
        if cap and index not in self._plan_refused:
            if index not in self._power_plans:
                window_lo, window_hi, window_maximum = self._plan_window
                plan, reserve = [], 256

                for mark in self._power_roots(
                    self.polynomial,
                    roots,
                    window_maximum,
                    weight,
                    self.budget,
                    lo=window_lo,
                    hi=window_hi,
                    inverses=self._lift_inverses,
                    shared_inverses=self._shared_lift_inverses,
                    root_transfers=self._power_transfers,
                ):
                    modulus, residues, _ = mark
                    words = (modulus.bit_length() + 7) // 8
                    size = 384 + (4 + len(residues)) * (96 + 8 * words)
                    if self._plan_bytes + reserve + size > cap:
                        self._plan_refused.add(index)
                        stats["plan_refusals"] += 1
                        break
                    self.budget.consume(len(residues) + words + 1)
                    reserve += size
                    plan.append(mark)
                else:
                    if self._plan_bytes + reserve <= cap:
                        self.budget.consume(len(plan) + 1)
                        # Only a complete lift sequence may enter the cache;
                        # a truncated prefix could miss higher valuations.
                        self._power_plans[index] = tuple(plan)
                        self._plan_bytes += reserve
                        stats["plan_builds"] += 1
                    else:
                        self._plan_refused.add(index)
                        stats["plan_refusals"] += 1

            if index in self._power_plans:
                plan = self._power_plans[index]
                self.budget.consume(len(plan) + 1)
                stats["plan_replays"] += 1
                return plan

        return self._power_roots(
            self.polynomial,
            roots,
            maximum,
            weight,
            self.budget,
            lo=lo,
            hi=hi,
            inverses=self._lift_inverses,
            shared_inverses=self._shared_lift_inverses,
            root_transfers=self._power_transfers,
        )

    def _set_plan_interval(self, lo, hi):
        """Bind reusable plans to a verified whole polynomial interval."""
        bits = max(abs(lo).bit_length(), abs(hi).bit_length())
        coefficient_bits = max(
            self.polynomial.a.bit_length(),
            abs(self.polynomial.b).bit_length(),
            abs(self.polynomial.c).bit_length(),
        )
        self.budget.consume(4 * (coefficient_bits + 2 * bits + 3))
        key = lo, hi, self._bounds(lo, hi)[1]
        # Lift depth and hit coverage depend on the whole interval. Reusing
        # a narrower window's plan could drop newly eligible prime powers.
        if key != self._plan_window:
            self._power_plans.clear()
            self._plan_refused.clear()
            self._plan_bytes = 0
            self._plan_window = key

    def _lower_score(self, value):
        if self.config.score_policy == "fixed" and value:
            return log_bounds(value)[0]
        return max(0, value.bit_length() - 1)

    def _residual_score(self):
        if self.config.score_policy == "fixed":
            return log_bounds(self._candidate_bound)[1]
        return (self._candidate_bound - 1).bit_length()

    def _clip_threshold(self, threshold):
        if self.config.score_backend == "bytearray":
            return min(255, threshold)
        if self.config.score_backend == "array":
            return min(2**32 - 1, threshold)
        return threshold

    def _bounds(self, lo, hi):
        """Bound |F| with endpoints and the two integer vertex neighbors.

        A positive quadratic's extrema on integer positions occur here.
        A sign crossing gets lower bound zero, conservatively selecting all
        positions. This evaluates only O(1) big integers per working block.
        """
        vertex = -self.polynomial.b // self.polynomial.a
        positions = {lo, hi - 1}
        for position in (vertex, vertex + 1):
            if lo <= position < hi:
                positions.add(position)
        values = [self.polynomial.value(position) for position in positions]
        minimum, maximum = min(values), max(values)
        lower = 0 if minimum <= 0 <= maximum else min(map(abs, values))
        return lower, max(map(abs, values))

    def _weights(self, maximum):
        """Upper-bound all p-power log contributions using exact powers.

        If p divides F, its valuation is at most floor(log_p(max |F|)).
        Multiplying that bound by ceil(log2 p) cannot underestimate any
        smooth contribution. No Hensel depth cutoff can lose high powers.
        """
        weights = []

        for entry, log in zip(self.factor_base.entries, self._logs):
            self.budget.consume(maximum.bit_length() + 1)
            exponent, power = 0, entry.prime
            while power <= maximum:
                exponent += 1
                power *= entry.prime
            weights.append(exponent * log)

        return weights

    def _add(self, indices, weight, stats):
        """Update scores; byte translations allocate at most block width."""
        if self.config.score_backend == "bytearray":
            if isinstance(indices, range) and indices:
                table = bytes(min(255, value + weight) for value in range(256))
                start, stop, step = indices.start, indices.stop, indices.step
                segment = self._scores[start:stop:step]
                self._scores[start:stop:step] = segment.translate(table)
                stats["slice_bytes"] += 2 * len(segment) + 256
                stats["slice_allocations"] += 3
            else:
                for index in indices:
                    self._scores[index] = min(
                        255, self._scores[index] + weight
                    )
        elif self.config.score_backend == "array":
            for index in indices:
                self._scores[index] = min(
                    2**32 - 1, self._scores[index] + weight
                )
        else:
            for index in indices:
                self._scores[index] += weight

    def _sieve(self, lo, hi, stats):
        self._block_norm_bits = None
        if (
            getattr(self.config, "extension_policy", None) == "extended-v2"
            and self.config.batch_policy == "shared-v2"
        ):
            from .yield_memory import reserve_batch

            reserve_batch(self, lo, hi)
        self._piecewise_thresholds = None
        if self.config.sparse_queue_bytes:
            threshold = self._sparse_sieve(self, lo, hi, stats)
        else:
            threshold = self._sieve_direct(lo, hi, stats)
        if self._word_hits is not None and self._hits is self._word_hits:
            self._word_hits.seal()
        width = getattr(self.config, "refinement_width", 0)
        if (
            width
            and self.config.score_policy in ("candidate", "powers")
            and not self._skipped
        ):
            if (
                self._block_norm_bits is not None
                and self._block_norm_bits[0] == self._block_norm_bits[1]
            ):
                # Every subwindow has the same norm bit length, so the
                # refined threshold cannot exceed this block's threshold.
                return threshold
            refined = []
            for start in range(lo, hi, width):
                self.budget.consume(
                    8 * self.factor_base.n_prime.bit_length() + 1
                )
                lower, _ = self._bounds(start, min(start + width, hi))
                needed = max(0, lower.bit_length() - 1)
                needed -= (
                    self._candidate_bound - 1
                ).bit_length() + self._omitted_allowance
                needed = max(0, needed + self.config.threshold_extra)
                refined.append(max(threshold, self._clip_threshold(needed)))
            self._piecewise_thresholds = tuple(refined)
        return threshold

    def _piecewise_passes(self, offset):
        bounds = getattr(self, "_piecewise_thresholds", None)
        if bounds is None:
            return True
        self.budget.consume(1)
        return (
            self._scores[offset]
            >= bounds[offset // self.config.refinement_width]
        )

    def _sieve_direct(self, lo, hi, stats):
        """Mark complete root hits in bounded prime chunks, without F(x)."""
        width = hi - lo
        self._omitted_allowance = 0
        self.budget.consume(width)
        two_levels, initial_score, initial_hit = 0, 0, 0
        powers = self.config.score_policy in ("powers", "fixed")
        parity_two = "parity_two" in self.config._kernels
        if parity_two:
            lower, maximum = self._bounds(lo, hi)
            two_levels, initial_score, initial_hit = (
                kernels.parity_two_initial(self, maximum)
            )
        period_initializes = (
            powers
            and self._presieve is not None
            and self.config.score_backend == "list"
        )
        words = self._word_hits if self._hits is self._word_hits else None
        grouped_budget = self._group_lease
        if words is not None:
            words.begin(initial_hit, self.budget)
        # Both periodic implementations replace the full active slices,
        # including empty-prime periods. Avoid clearing those slices first.
        if not period_initializes:
            for index in range(width):
                self._scores[index] = initial_score
                if words is None:
                    self._hits[index] = initial_hit
        if not parity_two:
            lower, maximum = self._bounds(lo, hi)
        self._block_norm_bits = lower.bit_length(), maximum.bit_length()
        if two_levels:
            stats["parity_two_blocks"] = stats.get("parity_two_blocks", 0) + 1
            stats["parity_two_seeded_hits"] = (
                stats.get("parity_two_seeded_hits", 0) + two_levels * width
            )
        lower_threshold = max(0, lower.bit_length() - 1)
        lower_threshold -= (self._candidate_bound - 1).bit_length()
        if (
            self.config.score_policy == "adaptive"
            and lower_threshold <= 0
            and self.config.threshold_extra == 0
            and self.config.division != "bucket"
        ):
            # Every score is nonnegative. A zero lower threshold selects
            # every position, so weight calculation/marking cannot affect
            # coverage. Root division or candidate-only resieving recovers
            # exact exponents later; buckets still need their hit metadata.
            stats["blocks"] += 1
            stats["threshold_min"] = 0
            stats["skipped_score_blocks"] = (
                stats.get("skipped_score_blocks", 0) + 1
            )
            return 0

        weights = self._logs if powers else self._weights(maximum)
        omitted = 0
        buckets = (
            self.config.marking == "bucket" or self.config.division == "bucket"
        )
        chunk = self.config.metadata_chunk
        unroll = "unroll4" in self.config._kernels
        power_period = "presieve_powers" in self.config._kernels
        density = (
            "density" in self.config._kernels or unroll
        ) and self.config.score_backend == "list"
        presieved = (
            self._presieve.apply(self, lo, hi, maximum, buckets)
            if powers and self._presieve is not None
            else frozenset()
        )

        for start in range(0, len(self._roots), chunk):
            for index in range(start, min(start + chunk, len(self._roots))):
                roots, weight = self._roots[index], weights[index]
                prime = roots.prime
                recovery_mask = self._recovery_mask(index) if buckets else 0
                skipped = prime < self.config.small_prime_cutoff
                if skipped:
                    # Power scores recover skipped valuations exactly at the
                    # candidate. A block-wide allowance remains conservative,
                    # but cannot use bit_length * log as the refined threshold.
                    omitted += (
                        weight * maximum.bit_length() if powers else weight
                    )

                if powers and not skipped:
                    first_level = True
                    marks = (
                        self._power_marks(index, maximum, lo, hi, stats)
                        if self.config.power_plan_bytes
                        else self._power_roots(
                            self.polynomial,
                            roots,
                            maximum,
                            weight,
                            self.budget,
                            lo=lo,
                            hi=hi,
                            inverses=self._lift_inverses,
                            shared_inverses=self._shared_lift_inverses,
                            root_transfers=self._power_transfers,
                        )
                    )

                    for modulus, lifted, power_weight in marks:
                        if (
                            index == 0
                            and two_levels
                            and power_weight == weight
                            and modulus in (2, 4)
                            and modulus <= (1 << two_levels)
                        ):
                            # Universal normal levels were fused into buffer
                            # initialization. Keep semantic hit counters; never
                            # skip a later conservative singular overmark.
                            self.budget.consume(1)
                            if first_level:
                                stats["root_hits"] += width
                            stats["power_hits"] += width
                            first_level = False
                            continue
                        if (
                            grouped_budget is not None
                            and density
                            and first_level
                            and buckets
                            and words is not None
                            and len(lifted) == 2
                            and index not in presieved
                            and not power_period
                        ):
                            first_hits = range(
                                (lifted[0] - lo) % modulus, width, modulus
                            )
                            second_hits = range(
                                (lifted[1] - lo) % modulus, width, modulus
                            )
                            first_count = len(first_hits)
                            second_count = len(second_hits)
                            if grouped_budget._precharge_marking_pair(
                                first_count + 1, second_count + 1
                            ):
                                # These are the same two logical actions,
                                # already inside a checked polling interval.
                                # Preserve their original marking order and
                                # keep scalar refusal behavior on boundaries.
                                stats["root_hits"] += (
                                    first_count + second_count
                                )
                                stats["power_hits"] += (
                                    first_count + second_count
                                )
                                words.add(
                                    index,
                                    first_hits,
                                    self._scores,
                                    power_weight,
                                )
                                words.add(
                                    index,
                                    second_hits,
                                    self._scores,
                                    power_weight,
                                )
                                first_level = False
                                continue
                        for root in lifted:
                            if (first_level and index in presieved) or (
                                power_period
                                and (index, modulus) in presieved
                                and power_weight == weight
                            ):
                                hits = range(
                                    (root - lo) % modulus, width, modulus
                                )
                                self.budget.consume(len(hits) + 1)
                                if first_level:
                                    stats["root_hits"] += len(hits)
                                stats["power_hits"] += len(hits)
                                continue
                            if first_level and buckets and words is not None:
                                hits = range(
                                    (root - lo) % modulus, width, modulus
                                )
                                self.budget.consume(len(hits) + 1)
                                stats["root_hits"] += len(hits)
                                stats["power_hits"] += len(hits)
                                words.add(
                                    index, hits, self._scores, power_weight
                                )
                                continue
                            if density:
                                # Sparse powers have at most one hit. Keep
                                # the scalar reservation, including empty
                                # ranges, before touching either buffer.
                                offset = (root - lo) % modulus
                                if modulus >= width:
                                    count = int(offset < width)
                                    self.budget.consume(count + 1)
                                    if first_level:
                                        stats["root_hits"] += count
                                    stats["power_hits"] += count
                                    if count:
                                        if first_level and buckets:
                                            self._hits[offset] |= recovery_mask
                                        self._scores[offset] += power_weight
                                    continue

                                hits = range(offset, width, modulus)
                                self.budget.consume(len(hits) + 1)
                                stats["power_hits"] += len(hits)
                                if first_level:
                                    stats["root_hits"] += len(hits)
                                if unroll:
                                    kernels.add_unrolled(
                                        self._scores,
                                        self._hits,
                                        offset,
                                        width,
                                        modulus,
                                        power_weight,
                                        recovery_mask
                                        if first_level and buckets
                                        else 0,
                                    )
                                    continue
                                if first_level and buckets:
                                    # The same certified first-level hits
                                    # supply both recovery bits and scores.
                                    bit = recovery_mask
                                    for hit in hits:
                                        self._hits[hit] |= bit
                                        self._scores[hit] += power_weight
                                else:
                                    for hit in hits:
                                        self._scores[hit] += power_weight
                                continue

                            hits = range((root - lo) % modulus, width, modulus)
                            self.budget.consume(len(hits) + 1)
                            if first_level:
                                stats["root_hits"] += len(hits)

                            stats["power_hits"] += len(hits)
                            if (
                                first_level
                                and buckets
                                and self.config.score_backend == "list"
                            ):
                                # Both writes use the same certified range;
                                # avoid traversing all first-level hits twice.
                                scores, hit_masks = self._scores, self._hits
                                for hit in hits:
                                    hit_masks[hit] |= recovery_mask
                                    scores[hit] += power_weight
                            else:
                                if first_level and buckets:
                                    for hit in hits:
                                        self._hits[hit] |= recovery_mask
                                self._add(hits, power_weight, stats)

                        first_level = False

                    continue

                residues = (0,) if roots.all_positions else roots.roots
                step = 1 if roots.all_positions else prime

                for root in residues:
                    offset = (root - lo) % step
                    hits = range(offset, width, step)
                    self.budget.consume(len(hits) + 1)
                    stats["root_hits"] += len(hits)
                    if buckets:
                        if words is not None:
                            words.add(index, hits)
                        else:
                            for hit in hits:
                                self._hits[hit] |= recovery_mask
                    if skipped or powers or self.config.marking == "bucket":
                        continue
                    if (
                        self.config.marking == "sparse"
                        or bool(self.config._kernels)
                    ) and step >= width:
                        hits = (offset,) if offset < width else ()
                    self._add(hits, weight, stats)

        if not powers and self.config.marking == "bucket":
            for offset in range(width):
                bits = self._hits[offset]

                while bits:
                    bit = bits & -bits
                    index = bit.bit_length() - 1
                    if (
                        self._roots[index].prime
                        >= self.config.small_prime_cutoff
                    ):
                        self._add((offset,), weights[index], stats)

                    bits ^= bit

        # log2 |F| - log2 residual <= log2 factor-base part. Lower-bound
        # the first term and upper-bound the second; omitted primes get a
        # proved universal allowance. Positive extra deliberately loses it.
        threshold = (
            self._lower_score(lower)
            if self.config.score_policy == "fixed"
            else max(0, lower.bit_length() - 1)
        )
        self._omitted_allowance = omitted
        threshold -= (
            self._residual_score()
            if self.config.score_policy == "fixed"
            else (self._candidate_bound - 1).bit_length()
        ) + omitted
        threshold = max(0, threshold + self.config.threshold_extra)
        if self.config.score_backend != "list":
            threshold = self._clip_threshold(threshold)
        stats["blocks"] += 1
        stats["threshold_min"] = min(stats["threshold_min"], threshold)
        stats["threshold_max"] = max(stats["threshold_max"], threshold)
        return threshold

    def _candidate_passes(self, value, offset):
        """Refine a safe threshold using an already evaluated candidate norm.

        floor(log2 |F(x)|) replaces the block-wide lower bound. Scores still
        upper-bound all prime-power contributions; omitted primes retain
        their full allowance. Thus zero-extra refinement cannot reject an
        admissible norm. Clip both sides for byte scores as in block scoring.
        """
        if self.config.score_policy not in ("candidate", "powers", "fixed"):
            return True
        if self.config.score_policy == "fixed":
            threshold = self._lower_score(abs(value)) - self._residual_score()
        else:
            threshold = abs(value).bit_length() - 1
            threshold -= (self._candidate_bound - 1).bit_length()

        if self.config.score_policy in ("powers", "fixed") and self._skipped:
            remaining, allowance = abs(value), 0

            for index in self._skipped:
                self.budget.consume(remaining.bit_length() + 1)
                prime, log = self._roots[index].prime, self._logs[index]
                while remaining and remaining % prime == 0:
                    remaining //= prime
                    allowance += log

            threshold -= allowance
        else:
            threshold -= self._omitted_allowance

        threshold = max(0, threshold + self.config.threshold_extra)
        if self.config.score_backend != "list":
            threshold = self._clip_threshold(threshold)
        return self._scores[offset] >= threshold

    def _resieve(self, lo, hi, threshold, stats):
        """Recover exact candidate exponents in a separate root-hit pass.

        Only selected candidate values are formed. Iterate prime/root hits
        across this block, divide every valuation to completion, and retain
        A's exponents separately. Refusal changes scratch only; no position
        or checked store is committed before the full recovery pass returns.
        """
        self._resieved.clear()
        _, maximum = self._bounds(lo, hi)
        # Every recovered prime divides supported_A * F(x). The product of
        # the smallest base primes bounds the number of distinct entries,
        # including A support, without reserving a dense base per position.
        # At F(x)=0 the recovered map contains A support only.
        bits = (
            self.polynomial.supported_a.bit_length()
            + max(1, maximum).bit_length()
        )
        product, support = 1, 0

        for entry in self.factor_base.entries:
            self.budget.consume(
                product.bit_length() + entry.prime.bit_length()
            )
            product *= entry.prime
            if product.bit_length() > bits:
                break
            support += 1

        scratch = (hi - lo) * (1024 + 128 * support + 4 * maximum.bit_length())
        self._scratch_bytes = scratch
        self._scratch_peak_bytes = max(self._scratch_peak_bytes, scratch)

        values, remaining, recovered = {}, {}, {}

        for offset in range(hi - lo):
            if self._scores[offset] >= threshold:
                self.budget.consume(self.polynomial.n_prime.bit_length() + 1)
                value = self.polynomial.value(lo + offset)
                if not self._candidate_passes(value, offset):
                    continue
                values[offset] = value
                remaining[offset] = abs(value)
                recovered[offset] = {
                    (
                        i + 1
                        if getattr(self.config, "indexed_recovery", False)
                        else self._roots[i].prime
                    ): self._a_exponents[i]
                    for i in self._a_support
                }

        for index, (roots, a_exponent) in enumerate(
            zip(self._roots, self._a_exponents)
        ):
            self.budget.consume(1)
            prime = roots.prime
            residues = (0,) if roots.all_positions else roots.roots
            step = 1 if roots.all_positions else prime

            for root in residues:
                for offset in range((root - lo) % step, hi - lo, step):
                    if offset not in values or remaining[offset] == 0:
                        continue

                    self.budget.consume(remaining[offset].bit_length() + 1)

                    exponent = 0
                    while remaining[offset] % prime == 0:
                        remaining[offset] //= prime
                        exponent += 1
                    if exponent:
                        key = (
                            index + 1
                            if getattr(self.config, "indexed_recovery", False)
                            else prime
                        )
                        recovered[offset][key] = a_exponent + exponent
                    stats["division_primes"] += 1
                    stats["division_steps"] += exponent

        self._resieved = {
            offset: (
                values[offset],
                remaining[offset],
                sorted(recovered[offset].items()),
            )
            for offset in values
        }
        return True

    def _divide(self, position, offset, stats):
        """Recover every exponent exactly; roots/buckets certify hit coverage.

        Root exclusion is exact, independent of scores and omitted primes.
        All valuations at every included prime are divided to completion.
        Full factor-base division remains an explicit reference control.
        """
        if self.config.division == "resieve":
            if offset not in self._resieved:
                self.budget.consume()
                stats["candidates"] += 1
                stats["refined_rejections"] = (
                    stats.get("refined_rejections", 0) + 1
                )
                return None, None

            value, remaining, exponents = self._resieved[offset]
        else:
            value = self.polynomial.value(position)
            remaining, exponents = abs(value), []
        smooth_residual = getattr(self, "_smooth_hints", {}).get(position, 1)
        if smooth_residual != 1 and self.config.division != "resieve":
            remaining, remainder = divmod(remaining, smooth_residual)
            if remainder:
                raise ArithmeticError(
                    "batch residual does not divide its norm"
                )

        self.budget.consume(abs(value).bit_length() + 2)
        stats["candidates"] += 1
        if value == 0:
            stats["zeros"] += 1
            divisor = gcd(
                self.polynomial.u_value(position), self.factor_base.n
            )
            return None, divisor if utils.valid_divisor(
                divisor, self.factor_base.n
            ) else None

        if not self._candidate_passes(value, offset):
            stats["refined_rejections"] = (
                stats.get("refined_rejections", 0) + 1
            )
            return None, None

        if self.config.division == "resieve":
            indices = ()
        elif self.config.division == "bucket":
            hit_bits = (
                self._word_hits.mask(offset, self.budget)
                if self._word_hits is not None
                and self._hits is self._word_hits
                else self._hits[offset]
            )
            bits = hit_bits
            for index in self._a_support:
                bits |= self._recovery_mask(index)
            indices = []
            while bits:
                bit = bits & -bits
                indices.append(bit.bit_length() - 1)
                bits ^= bit
        else:
            indices = range(len(self._roots))

        self.budget.consume(len(indices) * (abs(value).bit_length() + 1))
        roots_to_divide = ((index, self._roots[index]) for index in indices)

        for index, roots in roots_to_divide:
            prime = roots.prime
            if self.config.division == "full":
                hit = True
            elif self.config.division == "bucket":
                hit = bool(hit_bits & self._recovery_mask(index))
            else:
                hit = roots.all_positions or position % prime in roots.roots

            exponent = self._a_exponents[index]
            if hit:
                stats["division_primes"] += 1
                while remaining % prime == 0:
                    remaining //= prime
                    exponent += 1
                    stats["division_steps"] += 1

            if exponent:
                exponents.append(
                    (
                        index + 1
                        if getattr(self.config, "indexed_recovery", False)
                        else prime,
                        exponent,
                    )
                )

        if getattr(self.config, "indexed_recovery", False):
            from .indexed_relations import IndexedFactors

            exponents = IndexedFactors(exponents, self.factor_base.primes)
        if self.config.division != "resieve":
            remaining *= smooth_residual

        if self._graph is not None:
            return self._divide_double(
                self, position, value, remaining, exponents, stats
            )
        if remaining > self.config.residual_bound:
            return None, None
        self.budget.consume(remaining.bit_length() ** 2)
        if remaining != 1 and utils.classify_prime(remaining) != (
            utils.Primality.PROVEN
        ):
            stats["composite_residuals"] += 1
            return None, None

        atom = AtomicRelation(
            self.polynomial,
            position,
            -1 if value < 0 else 1,
            tuple(exponents),
            remaining,
        )
        verify_atomic(
            atom,
            self.factor_base,
            residual_bound=self.config.residual_bound,
            budget=self.budget,
        )
        if remaining != 1:
            self.budget.consume(self.factor_base.n.bit_length())
            divisor = gcd(remaining, self.factor_base.n)
            if utils.valid_divisor(divisor, self.factor_base.n):
                return atom, divisor
            if divisor != 1:
                raise ValueError("nonunit residual has no proper split")

        return atom, None

    def _admit(self, atom, stats):
        """Commit only after verification/combination and all cap checks.

        FIFO eviction cannot reach pinned atoms. Refusal leaves store state
        untouched, so the same position can be retried with a new Budget.
        Expanded admission may first commit a separately recorded opaque
        conversion; refusal then preserves that conversion and retries only
        the uncommitted incoming atom.
        """
        require_published(self)
        if self._checked is not None:
            from .yield_collector import prepare_admission

            if not prepare_admission(self, atom, stats):
                return None
        if getattr(self.config, "extension_policy", None):
            from .extended_collector import admit

            return admit(self, atom, stats)
        if self._graph is not None and atom.residual_product != 1:
            return self._admit_double(self, atom, stats)
        if atom.relation_id in self._atoms:
            stats["duplicates"] += 1
            return None
        partner = (
            self._pending.get(atom.residual) if atom.residual != 1 else None
        )
        accepted = atom.residual == 1 or partner is not None
        evicted = None
        if not accepted:
            if self.config.max_partials == 0:
                stats["dropped_partials"] += 1
                return None
            if (
                self.config.max_partials is not None
                and len(self._pending) >= self.config.max_partials
            ):
                evicted = next(iter(self._pending.values()))

        released = self._atom_bytes[evicted] if evicted else 0

        reserve = 4096 + 256 * len(atom.exponents)
        reserve += 16 * (
            abs(atom.position).bit_length()
            + self.polynomial.n_prime.bit_length()
        )
        combination_atoms, scratch = None, 0
        if partner is not None:
            combination_atoms = (self._atoms[partner], atom)
            reserve += combined_storage_reserve(
                sum(len(a.exponents) for a in combination_atoms),
                self.factor_base.n.bit_length(),
            )
            scratch = (
                _combination_workspace(
                    combination_atoms,
                    self.factor_base,
                )
                - self.factor_base.workspace_bytes
            )

        combination = None
        if partner is not None:
            self._scratch_peak_bytes = max(
                self._scratch_peak_bytes, self._scratch_bytes + scratch
            )
            combination = combine_relations(
                combination_atoms,
                self.factor_base,
                budget=self.budget,
            )
            if combination.divisor is not None:
                raise ArithmeticError(
                    "previously checked unit residual changed"
                )

        compact_slot = None
        if self._arena is not None:
            atom, _, compact_slot = self._arena.plan(atom, None, self.budget)
        elif getattr(self.config, "record_policy", "objects") == "packed-v1":
            from .packed_relations import pack_relation

            atom = pack_relation(atom, self.budget)
        self.budget.consume(len(atom.exponents) + 1)
        # A failed insertion must not expose a partially published smooth
        # row (or legacy SLP match) as a resumable memory-limit result.
        self.budget.publication()
        self._publication_failed = True
        try:
            evicted_atoms = (self._atoms[evicted],) if evicted else ()
            if evicted:
                del self._pending[self._atoms[evicted].residual]
                del self._atoms[evicted]
                del self._atom_bytes[evicted]
                stats["evictions"] += 1

            self._atoms[atom.relation_id] = atom
            if self._checked is not None:
                from .yield_collector import publish_auxiliary

                publish_auxiliary(
                    self, atom, None, compact_slot, evicted_atoms
                )
            self._atom_bytes[atom.relation_id] = reserve
            self._workspace += reserve - released
            stats["admitted_atoms"] += 1
            if atom.residual == 1:
                self._full.append(atom)
                self._rows.append(atom)
            elif partner is not None:
                del self._pending[atom.residual]
                self._combined.append(combination.relation)
                self._rows.append(combination.relation)
                stats["matches"] += 1
            else:
                self._pending[atom.residual] = atom.relation_id

        except MemoryError as error:
            raise RelationPublicationError(
                "allocation failed during relation publication"
            ) from error
        else:
            self._publication_failed = False

        return None

    @accounted_job
    def collect(self, lo, hi):
        """Collect a bounded window using reusable half-open working blocks.

        A refused block or position returns its first uncommitted x. Scoring
        scratch may change on refusal; accepted store state never does.
        Eviction and thresholds may lose yield only as documented in config.
        """
        return self._collect(lo, hi, snapshot=True)

    def _collect_progress(self, lo, hi):
        """Collect for an internal caller that already owns the live store."""
        return self._collect(lo, hi, snapshot=False)

    def _collect(self, lo, hi, *, snapshot):
        """Share validation, traversal and publication for both results."""
        require_published(self)
        checked_position(lo)
        checked_position(hi)
        if not 0 <= hi - lo <= MAX_WINDOW_WIDTH:
            raise ValueError("window must be ordered and at most 1000000 wide")
        stats = dict.fromkeys(
            (
                "blocks",
                "scanned",
                "candidates",
                "root_hits",
                "power_hits",
                "division_primes",
                "division_steps",
                "zeros",
                "composite_residuals",
                "admitted_atoms",
                "matches",
                "duplicates",
                "evictions",
                "dropped_partials",
                "slice_bytes",
                "slice_allocations",
                "threshold_max",
                "plan_builds",
                "plan_replays",
                "plan_refusals",
            ),
            0,
        )
        if self._graph is not None:
            stats.update(dict.fromkeys(self._double_counters, 0))
        stats["threshold_min"] = 2**63
        position, reason, divisor = lo, "complete", None
        span_budget = (
            self.budget if isinstance(self.budget, PollingBudget) else None
        )
        reporting = (
            self._reporting
            if getattr(self.config, "extension_policy", None)
            else None
        )
        kernel_before = kernels.counters(self)
        batcher = None
        if (
            self._checked is not None
            and self.config.batch_policy == "shared-v2"
        ):
            from .cofactor_policy import CofactorBatch

            batcher = CofactorBatch(self.config.batch_size)

        try:
            if (
                (
                    self.config.power_plan_bytes
                    or self.config.sparse_queue_bytes
                )
                and lo < hi
                and (
                    self._plan_window is None
                    or lo < self._plan_window[0]
                    or hi > self._plan_window[1]
                )
            ):
                self._set_plan_interval(lo, hi)

            while position < hi:
                block_lo = position
                block_hi = min(hi, block_lo + self.config.block_width)
                if getattr(self.config, "extension_policy", None):
                    # Resume completes the same physical block before making
                    # a reporting decision; a pause cannot create a new epoch.
                    block_hi = min(
                        hi,
                        (position // self.config.block_width + 1)
                        * self.config.block_width,
                    )
                    if position % self.config.block_width == 0:
                        self._candidate_bound = self._reporting.decide()
                if batcher is not None:
                    batcher.records.clear()
                    batcher.hints.clear()
                threshold = self._sieve(block_lo, block_hi, stats)
                if self.config.division == "resieve":
                    self._resieve(block_lo, block_hi, threshold, stats)
                candidates = None
                if "scan_fused" in self.config._kernels:
                    self.budget.consume(block_hi - block_lo + 1)
                    candidates = kernels.word_candidates(
                        self._scores,
                        threshold,
                        block_hi - block_lo,
                        self.budget,
                    )
                    next_candidate = next(candidates, -1)
                while position < block_hi:
                    offset = position - block_lo
                    selected = (
                        self._scores[offset] >= threshold
                        if candidates is None
                        else offset == next_candidate
                    )
                    if not selected and candidates is None and span_budget:
                        # Remaining is already bounded by64. Cache the
                        # ledger reference, never its mutable limit or usage.
                        amount = span_budget._remaining
                        if amount < 1:
                            amount = 1
                        stop = min(block_hi - block_lo, offset + amount)
                        end = offset + 1
                        while end < stop and self._scores[end] < threshold:
                            end += 1
                        amount = end - offset
                        recorder_span = getattr(
                            self.budget, "recorder_span", None
                        )
                        if recorder_span is not None:
                            amount = recorder_span(amount)
                        # Publish only affordable certified rejects. The
                        # adapter counts these logical actions individually,
                        # so batching cannot enlarge its 64-action poll gap.
                        span_budget.consume_span(amount)
                        position += amount
                        stats["scanned"] += amount
                        if reporting is not None:
                            reporting.record(amount, 0, 0)
                        continue
                    if not selected and candidates is not None:
                        end = (
                            block_hi
                            if next_candidate < 0
                            else block_lo + next_candidate
                        )
                        amount = min(
                            kernels.SCAN_CHUNK,
                            end - position,
                        )
                        recorder_span = getattr(
                            self.budget, "recorder_span", None
                        )
                        if recorder_span is not None:
                            amount = recorder_span(amount)
                        # These positions are certified rejects. Reserve one
                        # unit each before committing this bounded span; a
                        # refused poll leaves its first position untouched.
                        self.budget.consume(amount)
                        position += amount
                        stats["scanned"] += amount
                        if getattr(self.config, "extension_policy", None):
                            self._reporting.record(amount, 0, 0)
                        continue
                    if selected:
                        # Candidate division reserves the position unit and
                        # recovery work in one atomic check. The former two
                        # adjacent polls did no intervening published work.
                        if not self._piecewise_passes(offset):
                            self.budget.consume()
                            atom, divisor = None, None
                        elif batcher is None:
                            atom, divisor = self._divide(
                                position, offset, stats
                            )
                        else:
                            atom, divisor = batcher.recover(
                                self,
                                position,
                                offset,
                                block_hi,
                                threshold,
                                stats,
                            )
                        if divisor is not None:
                            position += 1
                            stats["scanned"] += 1
                            if getattr(self.config, "extension_policy", None):
                                self._reporting.record(1, 1, 0)
                            reason = "factor_found"
                            break

                        if atom is not None:
                            refusal = self._admit(atom, stats)
                            if refusal:
                                reason = refusal
                                break
                    else:
                        self.budget.consume()

                    position += 1
                    stats["scanned"] += 1
                    if getattr(self.config, "extension_policy", None):
                        self._reporting.record(1, int(selected), 0)
                    # Searching another group can poll and refuse. Publish
                    # this candidate's cursor before doing that derived work.
                    if selected and candidates is not None:
                        next_candidate = next(candidates, -1)

                if reason != "complete":
                    break
                if (
                    getattr(self.config, "extension_policy", None)
                    and position % self.config.block_width == 0
                ):
                    self._candidate_bound = self._reporting.decide()
        except BudgetExhaustedError as error:
            reason = self.budget.reason or str(error)

        if stats["blocks"] == 0:
            stats["threshold_min"] = 0
        self._resieved.clear()
        self._scratch_bytes = 0
        stats["plan_bytes"] = self._plan_bytes
        stats.update(
            {
                key: value - kernel_before.get(key, 0)
                for key, value in kernels.counters(self).items()
            }
        )
        if not snapshot:
            # QSJob consumes status only. Expanding cumulative records here
            # copies a growing history after every bounded collection batch.
            return _CollectionProgress(
                divisor,
                position,
                reason,
                stats,
                self._workspace + self._scratch_peak_bytes,
            )
        return SieveResult(
            tuple(self._atoms.values()),
            tuple(self._full),
            tuple(self._combined),
            self.partial_ids,
            divisor,
            position,
            reason,
            stats,
            self._workspace + self._scratch_peak_bytes,
            resource_stop=getattr(self.budget, "stop", None),
        )
