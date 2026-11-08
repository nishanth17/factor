"""Deterministic interruption seams for arithmetic/checkpoint tests.

The old work-frontier fixtures now inject deterministic stops into the new production
owner. This preserves precise interrupted prefixes without reintroducing a
production work quota. Resource contract tests import the real Budget directly.
"""

from v2.benchmarks.infrastructure.resource_control import RecordingBudget


class InterruptBudget(RecordingBudget):
    frontier_reason = "work_limit"
