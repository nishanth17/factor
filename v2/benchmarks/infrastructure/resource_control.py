"""Recorder work frontiers, separate from production resource limits.

Historical runners use algorithmic counters as diagnostic endpoints.
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


RETIRED_CONFIG_FIELDS = frozenset(
    (
        "memory_bytes",
        "checkpoint_bytes",
        "max_atoms",
        "max_relations",
        "max_rows",
        "max_candidates",
        "max_tree_bits",
        "max_tree_nodes",
        "max_bits",
        "max_nodes",
        "edge_limit",
        "assignment_work",
        "assignment_cpu",
        "startup_work",
        "worker_memory",
        "parent_memory",
        "parent_memory_bytes",
        "worker_memory_bytes",
        "frame_bytes",
        "max_runs",
        "residue_bytes",
        "spill_bytes",
        "ecm_memory_bytes",
    )
)


def adapt_resource_configuration(factory, *args, **values):
    """Translate historical harness metadata at the explicit resource boundary.

    Production constructors still reject every removed argument. Immutable
    historical engines keep their original accepted fields. The adapter
    translates mixed row controls to their concrete retained cache/target
    fields, and keeps cache bytes on genuine retention objects.
    """
    import inspect

    parameters = inspect.signature(factory).parameters
    fields = getattr(factory, "__dataclass_fields__", {})
    values = dict(values)
    if "preparation_cache_bytes" in fields and "max_relations" in values:
        rows = values["max_relations"]
        values.setdefault(
            "preparation_cache_bytes",
            2 * 2**20 if rows is None else min(2 * 2**20, 2048 * rows),
        )
        values.setdefault("target_rows", rows)
    if "preparation_cache_bytes" not in parameters:
        values.pop("preparation_cache_bytes", None)
        if "target_rows" in values and "max_relations" in parameters:
            values["max_relations"] = values.pop("target_rows")
        else:
            values.pop("target_rows", None)
    if "recipe" in values and "recipe" not in parameters:
        values["memory_bytes"] = (
            2**30 if values.pop("recipe") == "calibrated" else 128 * 2**20
        )
    for name in RETIRED_CONFIG_FIELDS:
        if name in values and name not in parameters:
            values.pop(name)
    return factory(*args, **values)


def recording_budget(factory, **kwargs):
    """Keep a historical runner frontier at the recorder boundary."""
    if factory is ResourceBudget:
        return RecordingBudget(**kwargs)
    return factory(**kwargs)
