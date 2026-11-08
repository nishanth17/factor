"""Deterministic interruption seams for arithmetic/checkpoint tests.

The old work-frontier fixtures inject deterministic stops into the production
owner. This preserves precise interrupted prefixes without reintroducing a
production work quota. Resource contract tests import the real Budget directly.
"""

from v2.benchmarks.infrastructure.resource_control import RecordingBudget
from v2.ecm.worker_state import ActionBudget as ResourceActionBudget
from v2.ecm.worker_state import (
    ProspectiveActionBudget as ResourceProspectiveBudget,
)


class InterruptBudget(RecordingBudget):
    frontier_reason = "work_limit"

    def _sample(self, **kwargs):
        # Historical arithmetic fixtures inject action-frontier stops. Root
        # RSS/live CPU sampling has dedicated provider tests in H3.
        return None

    def track_process(self, process):
        from multiprocessing.process import BaseProcess

        if isinstance(process, BaseProcess):
            return super().track_process(process)

    def acknowledge(self, pid, cpu_total):
        if any(key[0] == pid for key in self._children):
            return super().acknowledge(pid, cpu_total)


class ActionBudget(ResourceActionBudget, RecordingBudget):
    """Action polling with a test-only interruption frontier."""


class ProspectiveActionBudget(ResourceProspectiveBudget, RecordingBudget):
    """Prospective polling with a test-only interruption frontier."""


def recipe_values(value):
    """Normalize only removed grants in a historical recipe comparison."""
    from dataclasses import asdict, is_dataclass

    if is_dataclass(value):
        value = asdict(value)
    if isinstance(value, dict):
        value = dict(value)
        if "max_relations" in value:
            rows = value["max_relations"]
            value.setdefault(
                "preparation_cache_bytes",
                2 * 2**20 if rows is None else min(2 * 2**20, 2048 * rows),
            )
            value.setdefault("target_rows", rows)
        return {
            key: recipe_values(item)
            for key, item in value.items()
            if key
            not in (
                "memory_bytes",
                "checkpoint_bytes",
                "max_atoms",
                "max_relations",
            )
        }
    if isinstance(value, (tuple, list)):
        return type(value)(recipe_values(item) for item in value)
    return value
