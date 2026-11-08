"""Recorder work frontiers, separate from production resource limits.

Historical runners use algorithmic counters as deterministic diagnostic endpoints.
These controls belong to the recorder, not to Factor's production resource API.
"""

from ...common import utils
from ...execution.budget import Budget as ResourceBudget
from ...execution.budget import BudgetExhaustedError


class RecordingBudget(ResourceBudget):
    """Inject a stop before a recorder-selected action frontier."""

    frontier_reason = "work_limit"

    def __init__(self, *, work_limit=2_000_000, **kwargs):
        self.work_limit = utils.require_integer(work_limit, "recorder work", 0)
        if kwargs.get("used", 0) > work_limit:
            raise ValueError("recorder frontier is below consumed work")
        wall_start = kwargs.pop("_wall_start", None)
        cpu_start = kwargs.pop("_cpu_start", None)
        super().__init__(**kwargs)
        if wall_start is not None:
            self._wall_start = wall_start
        if cpu_start is not None:
            self._cpu_start = cpu_start

    @property
    def used(self):
        return self._recorded_work

    @used.setter
    def used(self, value):
        if value > self.work_limit:
            self.reason = self.frontier_reason
            raise BudgetExhaustedError(self.reason)
        self._recorded_work = value
