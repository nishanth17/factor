"""Bounded fixed-polynomial QS collection and verified extraction."""

import time
from dataclasses import dataclass

from .. import utils
from ..budget import Budget, BudgetExhaustedError
from ..work_budget import PollingBudget
from .extraction import (
    DependencyExtractor,
    _prepare_relations,
    _VerificationCache,
)
from .linear_algebra import DependencySolver, filter_matrix
from .polynomial import checked_position
from .sieve_collector import SieveCollector, SieveConfig


@dataclass(frozen=True)
class QSResult:
    """Proper divisor or unfinished cofactor, with cumulative costs."""

    divisor: int | None
    cofactor: int
    reason: str
    next_position: int
    stats: dict


class QSJob:
    """In-memory resumable fixed-polynomial job, independent of dispatch.

    Collection batches are half-open. Singleton-filtered usable-row excess
    controls solve attempts, with zero-row kernels tried immediately. Trivial
    dependencies permit further collection. Window/store/budget caps bound
    every retry; SIQS families and serialized checkpoints belong to P3.4.
    """

    def __init__(
        self,
        polynomial,
        factor_base,
        lo,
        hi,
        *,
        config=None,
        budget=None,
        weight_two=False,
        pivot="highest",
        row_excess=2,
        batch_width=256,
        filter_row_growth=1,
        tested_dependencies=False,
        collector_class=SieveCollector,
    ):
        """Build collection state; propagate setup failures."""
        utils.require_integer(batch_width, "batch_width", 1)
        if batch_width > 4096:
            raise ValueError("batch_width exceeds 4096")
        utils.require_integer(row_excess, "row_excess", 0)
        if row_excess > 4096:
            raise ValueError("row_excess exceeds 4096")
        utils.require_integer(filter_row_growth, "filter_row_growth", 1)
        if filter_row_growth > 4096:
            raise ValueError("filter_row_growth exceeds 4096")
        if type(tested_dependencies) is not bool:
            raise TypeError("tested_dependencies must be Boolean")
        checked_position(lo)
        checked_position(hi)
        if not 0 <= hi - lo <= 1_000_000:
            raise ValueError("window must be ordered and at most 1000000 wide")
        self.budget = budget if budget is not None else Budget()
        self.config = config if config is not None else SieveConfig()
        self.collector = collector_class(
            polynomial,
            factor_base,
            config=self.config,
            budget=self.budget,
        )
        if getattr(self.config, "power_plan_bytes", 0):
            self.collector._set_plan_interval(lo, hi)
        # Only retained rows are cached. The collector owns the reservation,
        # including while matrix scratch coexists or a polynomial is switched.
        cache_bytes = min(2 * 2**20, 2048 * self.config.max_relations)
        self.collector._preparation_cache = None
        if self.collector._workspace + cache_bytes < self.config.memory_bytes:
            self.collector._preparation_cache = _VerificationCache(cache_bytes)
            self.collector._workspace += cache_bytes
        self.lo, self.hi, self.next_position = lo, hi, lo
        self.weight_two, self.pivot = weight_two, pivot
        self.row_excess, self.batch_width = row_excess, batch_width
        self.filter_row_growth = filter_row_growth
        self.tested_dependencies = tested_dependencies
        self.solver, self.extractor, self.prepared = None, None, None
        self.last_count = -1
        self.last_solved_count = -1
        self.final_solve_done = False
        self.divisor = None
        self.storage_reason = None
        self.storage_solve_done = False
        self._workspace_peak = self.collector._workspace
        self.stats = {
            "filter_calls": 0,
            "solve_calls": 0,
            "trivial_dependencies": 0,
            "collected_positions": 0,
            "stage_seconds": {},
            "collector_counts": {},
        }

    def _timed(self, stage, call):
        """Retain inclusive stage costs, including refused private work."""
        started = time.perf_counter()
        try:
            return call()
        finally:
            times = self.stats["stage_seconds"]
            times[stage] = times.get(stage, 0) + time.perf_counter() - started

    def _solve(self, final):
        """Prepare changed stores; preserve pending solve/extraction state."""
        relations = self.collector.matrix_relations
        count = len(relations)
        if self.solver is None and self.extractor is None:
            if count == self.last_count and (
                not final or count == self.last_solved_count
            ):
                return None
            if (
                not final
                and self.last_count >= 0
                and count - self.last_count < self.filter_row_growth
            ):
                return None

            store = dict(self.collector._atoms)
            prepared = self._timed(
                "preparation",
                lambda: _prepare_relations(
                    relations,
                    self.collector.factor_base,
                    store,
                    budget=self.budget,
                    memory_bytes=self.config.memory_bytes,
                    retained_workspace_bytes=self.collector._workspace,
                    verification_cache=self.collector._preparation_cache,
                ),
            )
            # Collection buffers and preparation scratch coexist. Only
            # their common immutable base and pinned atoms are shared.
            live_workspace = (
                self.collector._workspace
                + prepared.workspace_bytes
                - prepared.shared_workspace_bytes
            )
            self._workspace_peak = max(self._workspace_peak, live_workspace)
            remaining = self.config.memory_bytes - live_workspace
            matrix = self._timed(
                "filtering",
                lambda: filter_matrix(
                    prepared.rows,
                    weight_two=self.weight_two,
                    budget=self.budget,
                    memory_bytes=max(0, remaining),
                ),
            )
            self._workspace_peak = max(
                self._workspace_peak,
                live_workspace + matrix.workspace_bytes,
            )
            self.stats["filter_calls"] += 1
            self.stats["last_filter"] = matrix.stats
            excess = len(matrix.rows) - matrix.stats["output_columns"]
            if (
                not final
                and not matrix.zero_dependencies
                and (excess < self.row_excess)
            ):
                # Defer only a growing store; final/storage-cap attempts must
                # still test kernels when no further collection is possible.
                self.last_count = count
                return None

            self.prepared = prepared
            self.solver = DependencySolver(
                matrix,
                pivot=self.pivot,
                budget=self.budget,
            )
            self.stats["solve_calls"] += 1

        if self.extractor is None:
            self.solver.budget = self.budget
            self.budget.consume(0)
            dependencies = self._timed("elimination", self.solver.run)
            self.extractor = DependencyExtractor(
                self.prepared,
                dependencies,
                budget=self.budget,
                tested_cache=(
                    self.collector._preparation_cache
                    if self.tested_dependencies
                    else None
                ),
            )

        self.extractor.budget = self.budget
        self.budget.consume(0)
        divisor = self._timed("extraction", self.extractor.run)
        skipped = getattr(self.extractor, "cache_skips", 0)
        self.stats["trivial_dependencies"] += skipped + sum(
            trial.divisor is None for trial in self.extractor.trials
        )
        self.stats["tested_dependency_skips"] = (
            self.stats.get("tested_dependency_skips", 0) + skipped
        )
        self.stats["last_solver"] = {
            "xors": self.solver.xors,
            "peak_pivot_nonzeros": self.solver.peak_nonzeros,
            "dependencies": len(self.solver.dependencies),
        }
        self.last_count = count
        self.last_solved_count = count
        self.solver, self.extractor, self.prepared = None, None, None
        return divisor

    def _finish_storage(self):
        """Try retained rows once; budget refusal preserves pending work.

        A full store cannot collect further, but its checked rows may already
        yield a factor. A matrix memory refusal preserves an explicit stop
        without repeating preparation on every unchanged run.
        """
        if not self.storage_solve_done:
            try:
                self.divisor = self._solve(True)
            except MemoryError:
                self.storage_reason = "memory_limit"
            self.storage_solve_done = True

        return self.storage_reason

    def run(self, *, batch_limit=None):
        """Amortize external polls while retaining the caller's work ledger."""
        original = self.budget
        if not isinstance(original, PollingBudget):
            self.budget = PollingBudget(original)
        try:
            return self._run(batch_limit=batch_limit)
        finally:
            self.budget = original
            self.collector.budget = original

    def _run(self, *, batch_limit=None):
        """Return factor/exhaustion state; resume after extending one budget.

        Retain consumed work and prior active wall/CPU when replacing budgets.
        MemoryError/invalid provenance propagate; the checked store survives.
        Storage caps trigger a final bounded extraction attempt. Successful
        and unchanged storage-exhausted results return without repeating work.
        """
        if batch_limit is not None:
            utils.require_integer(batch_limit, "batch_limit", 1)
        batches = 0
        reason = "window_exhausted"

        try:
            while self.divisor is None:
                if self.storage_reason is not None:
                    reason = self._finish_storage()
                    break
                final = self.next_position >= self.hi
                if final and self.final_solve_done:
                    break
                if self.solver is not None or self.extractor is not None:
                    self.divisor = self._solve(final)
                    if self.divisor is not None:
                        break
                if not final:
                    self.collector.budget = self.budget
                    stop = min(self.hi, self.next_position + self.batch_width)
                    result = self._timed(
                        "collection",
                        lambda: self.collector.collect(
                            self.next_position, stop
                        ),
                    )
                    counts = self.stats["collector_counts"]
                    for key, value in result.stats.items():
                        if not key.startswith("threshold"):
                            counts[key] = counts.get(key, 0) + value
                    self._workspace_peak = max(
                        self._workspace_peak,
                        result.workspace_bytes,
                    )
                    self.stats["collected_positions"] += result.stats[
                        "scanned"
                    ]
                    self.next_position = result.next_position
                    if result.divisor is not None:
                        self.divisor = result.divisor
                        break
                    batches += 1
                    collection_reason = result.reason
                    # Release caller snapshots before any filtering, including
                    # the final attempt at a storage cap. Otherwise evicted
                    # atoms remain pinned outside the owned reservation.
                    del result
                    if collection_reason != "complete":
                        reason = collection_reason
                        if reason in (
                            "relation_limit",
                            "atom_limit",
                            "memory_limit",
                        ):
                            self.storage_reason = reason
                            reason = self._finish_storage()

                        break

                final = self.next_position >= self.hi
                self.divisor = self._solve(final)
                if final:
                    self.final_solve_done = True
                elif batch_limit is not None and batches >= batch_limit:
                    reason = "paused"
                    break

            self.budget.consume(0)
        except BudgetExhaustedError:
            reason = self.budget.reason

        if self.divisor is not None:
            if not utils.valid_divisor(
                self.divisor, self.collector.factor_base.n
            ):
                raise ArithmeticError("extraction returned an invalid divisor")
            reason = "factor_found"

        self.stats["work_used"] = self.budget.used
        self.stats["workspace_bytes"] = self._workspace_peak
        return QSResult(
            self.divisor,
            (
                self.collector.factor_base.n // self.divisor
                if self.divisor
                else self.collector.factor_base.n
            ),
            reason,
            self.next_position,
            dict(self.stats),
        )
