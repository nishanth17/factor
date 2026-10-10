"""Bounded multi-polynomial QS, shared verified stores and finite recovery."""

import time
from dataclasses import dataclass, field, replace
from functools import partial

from ..common import arithmetic, utils
from ..execution.budget import Budget, BudgetExhaustedError
from .assignment_stream import AssignmentStream
from .external_square import (
    CoefficientDivisorError,
    CoefficientExhaustedError,
    external_square,
)
from .factor_base import (
    MAX_FACTOR_BASE_BOUND,
    build_factor_base,
    checked_target,
)
from .families import PolynomialFamily, family_assignments
from .multiplier import select_multiplier
from .pipeline import QSJob, QSResult
from .polynomial import (
    a_target,
    mpqs_polynomial,
    polynomial_roots,
    qs_polynomial,
)
from .sieve_collector import SieveCollector, SieveConfig


@dataclass(frozen=True)
class SIQSConfig:
    """Finite experimental parameters; growth retains the same factor base."""

    backend: str = field(default="python-int", kw_only=True)
    mode: str = "siqs"
    base_bound: int = 1000
    multiplier: int = 1
    half_width: int = 512
    factor_count: int = 3
    family_count: int = 16
    pool_size: int = 16
    assignment_policy: str = "reference"
    external_coefficients: bool = False
    coefficient_trials: int = 4096
    polynomials_per_family: int = 0
    diverse: bool = True
    shared_relations: bool = True
    growth_steps: int = 0
    max_half_width: int = 8192
    max_stalled: int = 16
    max_trivial: int = 128
    weight_two: bool = True
    row_excess: int = 2
    batch_width: int = 256
    filter_row_growth: int = 1
    tested_dependencies: bool = False
    memory_bytes: int = 32 * 1024 * 1024
    checkpoint_bytes: int = 1024 * 1024
    collector: SieveConfig = field(
        default_factory=lambda: SieveConfig(
            residual_bound=10000,
            max_atoms=4096,
            max_relations=2048,
            max_partials=512,
        )
    )

    def __post_init__(self):
        """Bound schedules and reserve checkpoint/registry storage."""
        if self.backend not in ("python-int", "gmpy2-mpz"):
            raise ValueError("unknown arithmetic backend")

        limits = {
            "base_bound": (3, MAX_FACTOR_BASE_BOUND),
            "multiplier": (0, 1000000),
            "half_width": (1, 499999),
            "factor_count": (
                1,
                8 if self.assignment_policy == "reference" else 32,
            ),
            "family_count": (1, 2**31 - 1 if self.streaming else 64),
            "pool_size": (self.factor_count, 4096 if self.streaming else 128),
            "coefficient_trials": (1, 65536),
            "polynomials_per_family": (0, 2**31),
            "growth_steps": (0, 4),
            "max_half_width": (self.half_width, 499999),
            "max_stalled": (1, 2**31 - 1),
            "max_trivial": (1, 2**31 - 1),
            "row_excess": (0, 4096),
            "batch_width": (1, 4096),
            "filter_row_growth": (1, 4096),
            "checkpoint_bytes": (
                4096,
                (64 if self.streaming else 16) * 1024 * 1024,
            ),
            "memory_bytes": (4 * 1024 * 1024, 4 * 1024 * 1024 * 1024),
        }
        for name, (low, high) in limits.items():
            value = utils.require_integer(getattr(self, name), name, low)
            if value > high:
                raise ValueError(f"{name} exceeds the SIQS limit")
        for name in (
            "diverse",
            "shared_relations",
            "weight_two",
            "external_coefficients",
            "tested_dependencies",
        ):
            if type(getattr(self, name)) is not bool:
                raise TypeError(f"{name} must be Boolean")

        if self.assignment_policy not in ("reference", "nearest", "flyer"):
            raise ValueError("unknown assignment policy")
        if self.assignment_policy != "reference" and self.mode != "siqs":
            raise ValueError("streaming A policies require SIQS")
        if self.external_coefficients and self.mode != "mpqs":
            raise ValueError("external square coefficients require MPQS")
        if self.streaming and (self.growth_steps or not self.diverse):
            raise ValueError(
                "streaming jobs use a fixed width and distinct assignments"
            )
        if self.polynomials_per_family and not self.streaming:
            raise ValueError("a Gray quota requires a streaming job")
        if self.mode not in ("qs", "mpqs", "siqs"):
            raise ValueError("unknown polynomial mode")
        if not isinstance(self.collector, SieveConfig):
            raise TypeError("collector must be a SieveConfig")
        if self.metadata_reserve + 3 * 1024 * 1024 > self.memory_bytes:
            raise MemoryError(
                "SIQS metadata/checkpoint reservation exceeds cap"
            )

    @property
    def streaming(self):
        """Whether assignments have a constant-storage, extendable prefix."""
        return (
            self.assignment_policy != "reference" or self.external_coefficients
        )

    @property
    def gray_limit(self):
        """Finite per-family quota, independent of total assignment count."""
        size = 1 << (self.factor_count - 1) if self.mode == "siqs" else 1
        return min(size, self.polynomials_per_family or size)

    @property
    def polynomial_limit(self):
        """Finite bound including every permitted width epoch."""
        count = self.family_count * self.gray_limit
        return count * (self.growth_steps + 1)

    @property
    def metadata_reserve(self):
        """Registry, assignments and coexistence with checkpoint encoding."""
        registry = 0 if self.streaming else 512 * self.polynomial_limit
        return 65536 + 8 * self.checkpoint_bytes + registry


class SIQSJob:
    """One setup-to-extraction allowance, with a shared single-prime store.

    The base and multiplier stay fixed during bounded width recovery; there
    is no spill or implicit base remapping. Checkpoints rebuild charged root,
    matrix and elimination caches from compact progress, never reset resources.
    """

    def __init__(self, n, *, seed=7, config=None, budget=None):
        self.config = config if config is not None else SIQSConfig()
        checked_target(n, max(1, self.config.multiplier))
        utils.require_integer(seed, "seed", 0)
        if seed >= 2**64:
            raise ValueError("seed exceeds 64 bits")
        self.n = arithmetic.get_backend(self.config.backend).integer(n)
        self.seed = seed
        self.budget = budget if budget is not None else Budget()
        self.base = self.family = self.engine = self.pending_step = None
        self.assignments = None
        self.coefficient_cursor = 0
        self.multiplier = None
        self.family_index = self.epoch = self.stalled = 0
        self.half_width = self.config.half_width
        self.active = False
        self.seen = set()
        self.finished_reason = self.divisor = None
        self.poly_start_rows = 0
        self.stats = dict(
            polynomials=0,
            duplicate_polynomials=0,
            distinct_polynomials=0,
            recoveries=0,
            blocks=0,
            root_reconstructions=0,
            matrix_reconstructions=0,
        )
        self.stats["stage_seconds"] = {}

    def _timed(self, stage, call):
        started = time.perf_counter()

        try:
            return call()
        finally:
            timings = self.stats["stage_seconds"]
            timings[stage] = (
                timings.get(stage, 0) + time.perf_counter() - started
            )

    def _setup(self):
        if self.multiplier is None:
            if self.config.multiplier:
                self.multiplier = self.config.multiplier
            else:
                choice = self._timed(
                    "multiplier",
                    lambda: select_multiplier(self.n, budget=self.budget),
                )
                self.multiplier, self.divisor = (
                    choice.multiplier,
                    choice.divisor,
                )
                self.stats["multiplier_scores"] = choice.scores
                if self.divisor is not None:
                    return

        if self.base is None:
            setup = self._timed(
                "base",
                lambda: build_factor_base(
                    self.n,
                    multiplier=self.multiplier,
                    bound=self.config.base_bound,
                    budget=self.budget,
                    memory_bytes=self.config.memory_bytes
                    - self.config.metadata_reserve,
                ),
            )
            self.base, self.divisor = setup.factor_base, setup.divisor
            if self.divisor is not None:
                return

        if "capacity" not in self.stats:
            from .capacity import capacity_report

            self.stats["capacity"] = capacity_report(self.base, self.config)
        if self.assignments is None:
            if self.config.mode == "siqs" and self.config.streaming:
                self.assignments = self._timed(
                    "assignments",
                    lambda: AssignmentStream(
                        self.base,
                        self.half_width,
                        factor_count=self.config.factor_count,
                        family_count=self.config.family_count,
                        pool_size=self.config.pool_size,
                        seed=self.seed,
                        policy=self.config.assignment_policy,
                        budget=self.budget,
                        memory_bytes=self.config.memory_bytes
                        - self.config.metadata_reserve,
                    ),
                )
            elif self.config.external_coefficients:
                self.assignments = range(self.config.family_count)
            elif self.config.mode == "siqs":
                count = self.config.family_count if self.config.diverse else 1
                self.assignments = self._timed(
                    "assignments",
                    lambda: family_assignments(
                        self.base,
                        self.half_width,
                        factor_count=self.config.factor_count,
                        family_count=count,
                        pool_size=self.config.pool_size,
                        seed=(self.seed + self.epoch) % 2**64,
                        budget=self.budget,
                        memory_bytes=self.config.memory_bytes
                        - self.config.metadata_reserve,
                    ),
                )
            elif self.config.mode == "mpqs":
                target = a_target(self.base.n_prime, self.half_width)
                self.budget.consume(len(self.base.entries))
                primes = sorted(
                    (
                        p
                        for p in self.base.primes
                        if p != 2 and self.base.n_prime % p
                    ),
                    key=lambda p: (abs(p * p - target), p),
                )
                self.assignments = tuple(
                    (p,) for p in primes[: self.config.family_count]
                )
            else:
                self.assignments = ((),)

            self.stats["assignment_count"] = len(self.assignments)
            if isinstance(self.assignments, AssignmentStream):
                self.stats["assignment_space"] = self.assignments.total

    def _next_step(self):
        if isinstance(self.assignments, AssignmentStream):
            self.assignments.budget = self.budget
        if self.config.mode != "siqs":
            if self.family_index >= len(self.assignments):
                return None
            if self.config.external_coefficients:
                polynomial, cursor, certainty = external_square(
                    self.base,
                    self.half_width,
                    budget=self.budget,
                    cursor=self.coefficient_cursor,
                    trials=self.config.coefficient_trials,
                )
                roots = tuple(
                    polynomial_roots(
                        polynomial, self.base, entry, budget=self.budget
                    )
                    for entry in self.base.entries
                )
                self.coefficient_cursor = cursor
                self.stats["coefficient_certainty"] = certainty
                return polynomial, roots

            polynomial = (
                qs_polynomial(self.base)
                if self.config.mode == "qs"
                else mpqs_polynomial(
                    self.base,
                    self.half_width,
                    budget=self.budget,
                    prime=self.assignments[self.family_index][0],
                )
            )
            roots = tuple(
                polynomial_roots(
                    polynomial, self.base, entry, budget=self.budget
                )
                for entry in self.base.entries
            )
            return polynomial, roots

        while self.family_index < len(self.assignments):
            if self.family is None:
                self.family = PolynomialFamily(
                    self.base,
                    self.assignments[self.family_index],
                    budget=self.budget,
                    memory_bytes=self.config.memory_bytes
                    - self.config.metadata_reserve,
                )

            self.family.budget = self.budget
            step = (
                self.family.next()
                if self.family.next_index < self.config.gray_limit
                else None
            )
            if step is not None:
                return step.polynomial, step.roots
            self.family_index += 1
            self.family = None

        return None

    def _activate(self):
        if self.pending_step is None:
            self.pending_step = self._timed("roots", self._next_step)
        if self.pending_step is None:
            return False
        polynomial, roots = self.pending_step
        key = (
            polynomial.a,
            min(polynomial.b % polynomial.a, -polynomial.b % polynomial.a),
            self.half_width,
        )
        if key in self.seen:
            self.stats["duplicate_polynomials"] += 1
            self.pending_step = None
            if self.config.mode != "siqs":
                self.family_index += 1
            return True

        extra = (
            self.config.metadata_reserve + 65536 + 640 * len(self.base.entries)
        )
        extra += 16 * (
            max(2, self.config.factor_count) * self.base.bound.bit_length()
            + self.n.bit_length()
        )
        config = replace(
            self.config.collector,
            memory_bytes=self.config.memory_bytes - extra,
        )
        if self.engine is None:
            self.engine = QSJob(
                polynomial,
                self.base,
                -self.half_width,
                self.half_width + 1,
                config=config,
                budget=self.budget,
                weight_two=self.config.weight_two,
                row_excess=self.config.row_excess,
                batch_width=self.config.batch_width,
                filter_row_growth=self.config.filter_row_growth,
                tested_dependencies=self.config.tested_dependencies,
                collector_class=partial(
                    SieveCollector, precomputed_roots=roots
                ),
            )
        else:
            self.engine.collector.budget = self.budget
            self.engine.collector.set_polynomial(
                polynomial,
                roots,
                retain_relations=self.config.shared_relations,
            )
            self.engine.lo = self.engine.next_position = -self.half_width
            self.engine.hi = self.half_width + 1
            self.engine.final_solve_done = False
            if self.config.collector.power_plan_bytes:
                self.engine.collector._set_plan_interval(
                    -self.half_width, self.half_width + 1
                )
            if not self.config.shared_relations:
                self.engine.last_count = self.engine.last_solved_count = -1

        if not self.config.streaming:
            self.seen.add(key)
        self.stats["polynomials"] += 1
        self.stats["a_minimum_used"] = min(
            self.stats.get("a_minimum_used", polynomial.a), polynomial.a
        )
        self.stats["a_maximum_used"] = max(
            self.stats.get("a_maximum_used", polynomial.a), polynomial.a
        )
        self.stats["distinct_polynomials"] = (
            self.stats["polynomials"]
            if self.config.streaming
            else len({(a, b) for a, b, _ in self.seen})
        )
        self.poly_start_rows = len(self.engine.collector._full) + len(
            self.engine.collector._combined
        )
        self.pending_step = None
        self.active = True
        return True

    def _recover(self, reason):
        if (
            self.epoch < self.config.growth_steps
            and self.half_width < self.config.max_half_width
        ):
            self.budget.consume(len(self.base.entries) + 1)
            self.epoch += 1
            self.half_width = min(
                self.config.max_half_width, 2 * self.half_width
            )
            self.assignments = self.family = self.pending_step = None
            self.family_index = self.stalled = 0
            self.stats["recoveries"] += 1
            self.active = False
        else:
            self.finished_reason = reason

    def run(self, *, max_blocks=None):
        """Run or pause at block boundaries; all exhausted results retain n."""
        if max_blocks is not None:
            utils.require_integer(max_blocks, "max_blocks", 0)
        blocks = 0
        reason = self.finished_reason or "paused"

        try:
            while self.finished_reason is None and self.divisor is None:
                if max_blocks is not None and blocks >= max_blocks:
                    break
                self.budget.consume(0)
                self._setup()
                if self.divisor is not None:
                    break
                if not self.active:
                    if not self._activate():
                        exhausted = "families_exhausted"
                        if (
                            isinstance(self.assignments, AssignmentStream)
                            and len(self.assignments)
                            < self.config.family_count
                        ):
                            exhausted = "assignment_space_exhausted"

                        self._recover(exhausted)

                    continue

                self.engine.budget = self.budget
                result = self._timed(
                    "collection_to_extraction",
                    lambda: self.engine.run(batch_limit=1),
                )
                blocks += 1
                self.stats["blocks"] += 1
                self.divisor = result.divisor
                reason = result.reason
                if reason == "window_exhausted":
                    rows = len(self.engine.collector._full) + len(
                        self.engine.collector._combined
                    )
                    self.stalled = (
                        self.stalled + 1 if rows == self.poly_start_rows else 0
                    )
                    self.active = False
                    if self.config.mode != "siqs":
                        self.family_index += 1
                    if self.stalled >= self.config.max_stalled:
                        self._recover("stalled_yield")
                    elif (
                        self.engine.stats["trivial_dependencies"]
                        >= self.config.max_trivial
                    ):
                        self._recover("trivial_dependency_limit")
                elif reason in (
                    "relation_limit",
                    "atom_limit",
                    "memory_limit",
                ):
                    self.finished_reason = reason
                elif reason not in ("paused", "factor_found"):
                    break
        except CoefficientDivisorError as found:
            self.divisor = found.divisor
        except CoefficientExhaustedError:
            self.finished_reason = reason = "coefficient_limit"
        except BudgetExhaustedError:
            reason = self.budget.reason
        except MemoryError:
            self.finished_reason = reason = "memory_limit"

        if self.divisor is not None:
            if not utils.valid_divisor(self.divisor, self.n):
                raise ArithmeticError("SIQS returned an invalid split")
            reason = "factor_found"
        elif self.finished_reason is not None:
            reason = self.finished_reason
        elif (
            max_blocks is not None
            and blocks >= max_blocks
            and reason in ("paused", "window_exhausted")
        ):
            reason = "paused"

        if self.engine is not None:
            self.stats["engine"] = dict(self.engine.stats)
            self.stats["relations"] = len(self.engine.collector._full) + len(
                self.engine.collector._combined
            )
            self.stats["partials"] = len(self.engine.collector.partial_ids)
            graph = self.engine.collector._graph
            if graph is not None:
                self.stats["large_prime_graph"] = dict(
                    forest_edges=len(graph.edges),
                    forest_vertices=len(graph.adjacency),
                    components=len(graph.sizes),
                    unowned_edges=len(graph.unowned),
                    owned_forest_edges=len(graph.edges) - len(graph.unowned),
                    split_calls=self.engine.collector._split_calls,
                )
            self.stats["workspace_bytes"] = (
                self.config.memory_bytes
                - self.engine.config.memory_bytes
                + max(
                    self.engine._workspace_peak,
                    self.engine.collector._workspace,
                    getattr(self.engine.collector, "_switch_peak", 0),
                )
            )

        self.stats["work_used"] = self.budget.used
        return QSResult(
            self.divisor,
            self.n // self.divisor if self.divisor else self.n,
            reason,
            self.engine.next_position if self.engine else 0,
            dict(self.stats),
        )

    def checkpoint(self):
        """Return a compact full-job envelope with checked-store identity."""
        from .checkpoint import pack_job

        return pack_job(self)

    @classmethod
    def from_checkpoint(
        cls, checkpoint, *, budget, config=None, allow_extension=False
    ):
        """Restore state under an allowance retaining prior resources."""
        from .checkpoint import restore_job

        try:
            return restore_job(
                checkpoint,
                budget=budget,
                config=config,
                allow_extension=allow_extension,
            )
        except (
            KeyError,
            TypeError,
            IndexError,
            AttributeError,
            RecursionError,
        ) as error:
            raise ValueError("malformed SIQS checkpoint") from error
