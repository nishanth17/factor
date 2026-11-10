"""The persistent pthread scheduler driving independent whole native curves."""

from ..execution.budget import BudgetExhaustedError
from .native_curve import BlockStore, NativeCurvePacket
from .native_resources import close_resource_checks, poll_native_budget
from .native_team import (
    _PROCESS_TEAM_LOCK,
    SCHEDULES,
    NativeTeam,
    TeamRefusedError,
)


class NativeCurveTeam(NativeTeam):
    """Borrow private contexts for one drained generation; refill in Python."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.blocks = BlockStore()
        self.bindings = {}
        self.phase_clocks = False
        self.atomic_actions = False

    def _grow(self, count, n, budget):
        while len(self.packets) < count:
            binding = (
                None
                if not self.packets
                else (self.ffi, self.lib, self.packets[0].manifest)
            )
            packet = NativeCurvePacket(
                _variant="pthread",
                _binding=binding,
                _programs=self.programs,
                _blocks=self.blocks,
            )
            self.packets.append(packet)
            try:
                packet.available(n, budget)
                if self.ffi is None:
                    self.ffi, self.lib = packet.ffi, packet.lib
                    self.contexts = self.ffi.new("b4_context *[32]")
                    self.statuses = self.ffi.new("int[32]")
                    self.actual = self.ffi.new("unsigned *")
                    self.owner = self.lib.b4_pthread_create()
                    if self.owner == self.ffi.NULL:
                        self.owner = None
                        raise MemoryError(
                            "whole native pthread owner creation failed"
                        )
            except BaseException:
                # The packet is already reachable by close(), even when a
                # partial preparation fails. No worker has borrowed it yet.
                raise

    def advance_jobs(self, jobs, budget, context, config):
        if self.closed or not self.lock.acquire(blocking=False):
            raise RuntimeError("native curve team closed or active")
        if not _PROCESS_TEAM_LOCK.acquire(blocking=False):
            self.lock.release()
            raise RuntimeError("another native team is active in this process")
        borrowed = []
        try:
            if not jobs or len(jobs) > self.capacity:
                raise ValueError("whole native curve queue exceeds capacity")
            poll_native_budget(budget, self)
            self._grow(self.capacity, jobs[0][1]["n"], budget)
            if self.lib.b4_pthread_reset(self.owner):
                raise RuntimeError("whole native pthread owner is active")
            positions = {position for position, _ in jobs}
            for position in list(self.bindings):
                if position not in positions:
                    self.bindings.pop(position).discard(sync=False)
            available = [
                p for p in self.packets if p not in self.bindings.values()
            ]
            for position, job in jobs:
                if position not in self.bindings:
                    packet = available.pop(0)
                    packet.atomic_actions, packet.phase_clocks = (
                        self.atomic_actions,
                        self.phase_clocks,
                    )
                    packet.bind_job(job, budget, context, config)
                    self.bindings[position] = packet
                elif self.bindings[position].job is not job:
                    self.bindings[position].bind_job(
                        job, budget, context, config
                    )
            prepared = [self.bindings[position] for position, _ in jobs]
            for index, packet in enumerate(prepared):
                if not packet.call_lock.acquire(blocking=False):
                    raise RuntimeError(
                        "whole native curve is already borrowed"
                    )
                borrowed.append(packet)
                self.contexts[index] = packet.context
            try:
                status = self.lib.b4_pthread_step(
                    self.owner,
                    self.contexts,
                    len(prepared),
                    self.workers,
                    SCHEDULES[self.schedule],
                    self.chunk_size,
                    5,
                    self.packets[0].manifest.get("slice_actions", 128),
                    self.statuses,
                    self.actual,
                )
                self.slices += 1
            finally:
                # All workers have drained before foreign return, including
                # signal reentry. Charge all peers before a cancellable poll.
                for packet in prepared:
                    packet.settle(budget)
            if status == -3:
                raise BudgetExhaustedError("cancelled")
            if status == -2 or int(self.actual[0]) != self.workers:
                raise TeamRefusedError("whole native pthread startup refused")
            if status == -4:
                raise RuntimeError("whole native pthread worker failed")
            if status or any(
                self.statuses[index] not in (0, 1, 2)
                for index in range(len(prepared))
            ):
                raise ValueError(
                    "whole native pthread generation returned invalid status"
                )
            poll_native_budget(budget, self)
            # Validate every ready divisor before any stable-order publication.
            ready = []
            for index, ((position, job), packet) in enumerate(
                zip(jobs, prepared)
            ):
                if self.statuses[index] == 2:
                    packet.feed(budget)
                elif self.statuses[index] == 1:
                    packet.sync(
                        allow_result=True, compact=not self.atomic_actions
                    )
                    ready.append((position, job))
                elif self.atomic_actions:
                    packet.sync()
            return ready
        finally:
            for packet in reversed(borrowed):
                packet.call_lock.release()
            _PROCESS_TEAM_LOCK.release()
            self.lock.release()

    def clear_jobs(self):
        for packet in self.packets:
            packet.discard(sync=False)
        self.bindings.clear()

    def _before_close(self):
        # Parked workers have no borrowed contexts between generations. Stop
        # them, release curve metadata references, free shared blocks, and only
        # then allow the first packet to unload the DSO.
        close_resource_checks(self)
        self.clear_jobs()
        self.blocks.close()
