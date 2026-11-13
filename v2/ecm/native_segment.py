"""Bounded native work with real root observations and bulk settlement."""

import os

from ..execution.budget import BudgetExhaustedError
from ..execution.ecm_workers import _base_budget
from ..execution.resources import MacResources
from .native_resources import native_bulk_allowed, poll_native_budget

UINT64_MAX = 2**64 - 1
REASONS = {
    0: "not_entered",
    1: "completed",
    2: "metadata",
    3: "quantum",
    4: "children",
    5: "rss",
    6: "observer_error",
    7: "cancelled",
    8: "frontier",
}


class NativeSegment:
    """One coordinator-owned policy/result buffer; never a worker callback."""

    def __init__(self, owner):
        self.ffi, self.lib = owner.ffi, owner.lib
        self.blocks = owner.blocks
        self.policy = self.ffi.new("b4_segment_policy *")
        self.info = self.ffi.new("b4_segment_info *")
        self.returns = {}
        self.observations = self.waves = self.refills = 0
        self._recorded = True

    def close(self):
        self.ffi.release(self.policy)
        self.ffi.release(self.info)

    def prepare(self, budget):
        if not native_bulk_allowed(budget):
            return False
        base = _base_budget(budget)
        pid = os.getpid()
        if set(base.roots) != {pid} or base.roots[pid] is None:
            return False

        def remaining(limit, used):
            if limit is None:
                return UINT64_MAX
            return min(UINT64_MAX, max(0, int((limit - used) * 10**9)))

        self.policy.pid = pid
        self.policy.birth = base.roots[pid]
        self.policy.wall_ns = remaining(base.seconds, base.wall_used)
        self.policy.cpu_ns = remaining(base.cpu_seconds, base.cpu_used)
        self.policy.rss_limited = int(base.rss_bytes is not None)
        self.policy.rss = (
            0 if base.rss_bytes is None else min(UINT64_MAX, base.rss_bytes)
        )
        return True

    def settle(self, packets, budget):
        work = [int(self.info.work[i]) for i in range(len(packets))]
        deltas = [
            value - packet.charged for value, packet in zip(work, packets)
        ]
        metadata_work = int(self.info.metadata_work)
        metadata_delta = metadata_work - self.blocks.charged_metadata
        if any(delta < 0 for delta in deltas) or metadata_delta < 0:
            raise ValueError("native segment incurred work moved backwards")
        # Validate every peer before mutating any charged cursor. Recording is
        # unconditional and precedes cancellation, observation or publication.
        for packet, value in zip(packets, work):
            packet.charged = value
        self.blocks.charged_metadata = metadata_work
        _base_budget(budget).record_incurred_work(sum(deltas) + metadata_delta)

    def begin(self, packets):
        # A signal can arrive before C entry. Seed the result with already
        # settled cursors, so the unconditional finally records zero new work
        # in that case. C replaces the snapshot only after it takes ownership.
        size = self.ffi.sizeof("b4_segment_info")
        self.ffi.memmove(self.info, b"\x00" * size, size)
        for index, packet in enumerate(packets):
            self.info.work[index] = packet.charged
        self.info.metadata_work = self.blocks.charged_metadata
        self._recorded = False

    def record(self):
        if self._recorded:
            return
        self._recorded = True
        reason = REASONS.get(int(self.info.reason), "invalid")
        self.returns[reason] = self.returns.get(reason, 0) + 1
        self.observations += int(self.info.observations)
        self.waves += int(self.info.waves)
        self.refills += int(self.info.refills)

    def check(self, budget, owner):
        if self.info.reason == 6:
            MacResources._proc_error(
                "native segment observation",
                self.policy.pid,
                int(self.info.error),
            )
        poll_native_budget(budget, owner)
        if self.info.reason == 5:
            base = _base_budget(budget)
            # An actual in-segment RSS refusal survives a later drop in RSS.
            base.observed_rss = int(self.info.rss)
            base.reason = "rss_limit"
            base.stop = dict(
                reason=base.reason,
                limit="rss",
                observed=base.observed_rss,
                limit_value=base.rss_bytes,
                scope=base.scope,
            )
            raise BudgetExhaustedError(base.reason, base.stop)


def segment_for(owner, budget):
    if (
        not getattr(owner, "coarse_segments", True)
        or getattr(owner, "atomic_actions", False)
        or not native_bulk_allowed(budget)
    ):
        return None
    segment = getattr(owner, "_segment", None)
    if segment is None:
        segment = owner._segment = NativeSegment(owner)
    return segment if segment.prepare(budget) else None
