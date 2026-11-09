"""Recorder work frontiers, separate from production resource limits.

Historical runners use algorithmic counters as diagnostic endpoints.
These controls belong to the recorder, not to Factor's production resource API.
"""

import json
from dataclasses import fields, is_dataclass
from dataclasses import replace as dataclass_replace
from pathlib import Path

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

    def recorder_span(self, amount):
        """Keep coalesced diagnostic spans at the recorder frontier."""
        return max(1, min(amount, self.work_limit - self.used))

    @property
    def used(self):
        return self._recorded_work

    @used.setter
    def used(self, value):
        if value > self.work_limit:
            self._check()
            self.reason = self.frontier_reason
            raise BudgetExhaustedError(self.reason)
        self._recorded_work = value

    def restore_charges(self, *, work, wall, cpu):
        if work > self.work_limit:
            raise ValueError("recorder frontier is below saved work")
        return super().restore_charges(work=work, wall=wall, cpu=cpu)


RETIRED_CONFIG_FIELDS = frozenset(
    (
        "assignment_seconds",
        "parent_work_reserve",
        "parent_cpu_reserve",
        "shutdown_work",
        "shutdown_cpu",
        "helper_cpu",
        "runtime_memory_bytes",
        "coexistence_bytes",
        "context_limit",
        "receipt_limit",
        "max_attempts",
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
    historical = {
        key: value
        for key, value in values.items()
        if key in RETIRED_CONFIG_FIELDS
    }
    recipe = values.get("recipe")
    if (
        getattr(factory, "__name__", "")
        in ("default_siqs_config", "legacy_siqs_config")
        and "recipe" in parameters
    ):
        legacy_memory = historical.get("memory_bytes")
        recipe = recipe or (
            "calibrated" if legacy_memory in (None, 2**30) else "reference"
        )
        values["recipe"] = recipe
        historical["memory_bytes"] = (
            (2**30 if recipe == "calibrated" else 128 * 2**20)
            if legacy_memory is None
            else legacy_memory
        )
        historical.setdefault("checkpoint_bytes", 2**20)
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
    result = factory(*args, **values)
    if is_dataclass(result):
        if recipe is not None:
            historical.setdefault(
                "memory_bytes",
                2**30 if recipe == "calibrated" else 128 * 2**20,
            )
        _attach_metadata(result, historical)
    return result


def recording_budget(factory, **kwargs):
    """Keep a historical runner frontier at the recorder boundary."""
    if factory is ResourceBudget:
        return RecordingBudget(**kwargs)
    return factory(**kwargs)


# Explicit alias for isolated historical runtime loaders.
Budget = RecordingBudget


def historical_asdict(value):
    """Serialize recorder metadata in its historical, source-pinned schema.

    Retired grants are metadata only here. Production configuration and its
    checkpoint serializers continue to expose only the new resource API.
    """
    if is_dataclass(value):
        result = {
            f.name: historical_asdict(getattr(value, f.name))
            for f in fields(value)
        }
        defaults_path = (
            Path(__file__).resolve().parents[1]
            / "inputs/controls/h3_historical_resources.json"
        )
        defaults = json.loads(defaults_path.read_text())["defaults"]
        metadata = defaults.get(type(value).__name__, {})
        metadata = dict(
            metadata, **getattr(value, "_historical_resources", {})
        )
        if (
            type(value).__name__ in ("SIQSConfig", "CountedSIQSConfig")
            and getattr(value.collector, "storage_policy", None) != "grow"
        ):
            if metadata.get("memory_bytes") is None:
                metadata["memory_bytes"] = 32 * 2**20
            if metadata.get("checkpoint_bytes") is None:
                metadata["checkpoint_bytes"] = 2**20
        metadata = {
            key: item for key, item in metadata.items() if key not in result
        }
        if "preparation_cache_bytes" in result:
            result.pop("preparation_cache_bytes")
            metadata["max_relations"] = result.pop("target_rows")
        result.update(metadata)
        return result
    if isinstance(value, dict):
        return {key: historical_asdict(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return type(value)(historical_asdict(item) for item in value)
    return value


def historical_replace(value, **changes):
    """Retain recorder-only grant metadata across literal policy changes."""
    metadata = dict(getattr(value, "_historical_resources", {}))
    accepted = {f.name for f in fields(value)}
    if "target_rows" in accepted and "max_relations" in changes:
        rows = changes["max_relations"]
        changes.setdefault("target_rows", rows)
        changes.setdefault(
            "preparation_cache_bytes",
            2 * 2**20 if rows is None else min(2 * 2**20, 2048 * rows),
        )
    for key in tuple(changes):
        if key in RETIRED_CONFIG_FIELDS and key not in accepted:
            metadata[key] = changes.pop(key)
    result = dataclass_replace(value, **changes)
    _attach_metadata(result, metadata)
    return result


def _attach_metadata(result, supplied):
    defaults_path = (
        Path(__file__).resolve().parents[1]
        / "inputs/controls/h3_historical_resources.json"
    )
    defaults = json.loads(defaults_path.read_text())["defaults"]
    metadata = dict(defaults.get(type(result).__name__, {}), **supplied)
    object.__setattr__(result, "_historical_resources", metadata)
    actual = {f.name for f in fields(result)}
    for key, value in metadata.items():
        if key not in actual:
            object.__setattr__(result, key, value)


def decode_historical_siqs(values):
    """Bridge a historical recorder recipe to current numerical settings."""
    from v2 import qs
    from v2.benchmarks.qs.c9a.configuration import decode

    return decode(values, qs)
