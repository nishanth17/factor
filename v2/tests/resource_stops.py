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
