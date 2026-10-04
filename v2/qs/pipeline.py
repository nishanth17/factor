"""Bounded fixed-polynomial QS collection and verified extraction."""

from dataclasses import dataclass

from .. import utils
from ..budget import Budget, BudgetExhaustedError
from .extraction import DependencyExtractor, prepare_relations
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
        collector_class=SieveCollector,
    ):
        """Build collection state; propagate setup failures."""
        utils.require_integer(batch_width, "batch_width", 1)
        if batch_width > 4096:
            raise ValueError("batch_width exceeds 4096")
        utils.require_integer(row_excess, "row_excess", 0)
        if row_excess > 4096:
            raise ValueError("row_excess exceeds 4096")
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
        self.lo, self.hi, self.next_position = lo, hi, lo
        self.weight_two, self.pivot = weight_two, pivot
        self.row_excess, self.batch_width = row_excess, batch_width
        self.solver, self.extractor, self.prepared = None, None, None
        self.last_count = -1
        self.last_solved_count = -1
        self.final_solve_done = False
        self.divisor = None
        self._workspace_peak = self.collector._workspace
        self.stats = {
            "filter_calls": 0,
            "solve_calls": 0,
            "trivial_dependencies": 0,
            "collected_positions": 0,
        }

    def _solve(self, final):
        """Prepare changed stores; preserve pending solve/extraction state."""
        relations = tuple(self.collector._full + self.collector._combined)
        count = len(relations)
        if self.solver is None and self.extractor is None:
            if count == self.last_count and (
                not final or count == self.last_solved_count
            ):
                return None
            store = dict(self.collector._atoms)
            prepared = prepare_relations(
                relations,
                self.collector.factor_base,
                store,
                budget=self.budget,
                memory_bytes=self.config.memory_bytes,
            )
            remaining = self.config.memory_bytes - max(
                prepared.workspace_bytes,
                self.collector._workspace,
            )
            matrix = filter_matrix(
                prepared.rows,
                weight_two=self.weight_two,
                budget=self.budget,
                memory_bytes=max(0, remaining),
            )
            self._workspace_peak = max(
                self._workspace_peak,
                max(prepared.workspace_bytes, self.collector._workspace)
                + matrix.workspace_bytes,
            )
            self.stats["filter_calls"] += 1
            self.stats["last_filter"] = matrix.stats
            excess = len(matrix.rows) - matrix.stats["output_columns"]
            if (
                not final
                and not matrix.zero_dependencies
                and (excess < self.row_excess)
            ):
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
            dependencies = self.solver.run()
            self.extractor = DependencyExtractor(
                self.prepared,
                dependencies,
                budget=self.budget,
            )
        self.extractor.budget = self.budget
        divisor = self.extractor.run()
        self.stats["trivial_dependencies"] += sum(
            trial.divisor is None for trial in self.extractor.trials
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

    def run(self):
        """Return factor/exhaustion state; resume after extending one budget.

        Retain consumed work and prior active wall/CPU when replacing budgets.
        MemoryError/invalid provenance propagate; the checked store survives.
        A successful result returns again without repeating work.
        """
        reason = "window_exhausted"
        try:
            while self.divisor is None:
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
                    result = self.collector.collect(self.next_position, stop)
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
                    if result.reason != "complete":
                        reason = result.reason
                        break
                    # Drop caller snapshots before filtering or collecting
                    # another batch. Otherwise FIFO-evicted partial atoms
                    # stay alive outside the collector's owned reservation.
                    del result
                final = self.next_position >= self.hi
                self.divisor = self._solve(final)
                if final:
                    self.final_solve_done = True
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
