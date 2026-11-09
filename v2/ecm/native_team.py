"""Reusable private contexts and one persistent C-only pthread team."""

import threading

from ..execution.budget import BudgetExhaustedError
from ..execution.resources import removed_arguments
from .native_packet import RESERVE, NativePacket, ProgramStore, poll_budget

SCHEDULES = {"dynamic": 0, "static": 1}
MAX_TASKS = 32
TEAM_OVERHEAD = 64 * 2**20
_PROCESS_TEAM_LOCK = threading.Lock()


class TeamRefusedError(RuntimeError):
    """The runtime refused the funded team size before native arithmetic."""


@removed_arguments(("memory_bytes",))
class NativeTeam:
    """No Python workers/callbacks inside a parallel region; one active owner.

    Requests are (kind, scalar-or-entries, point, n, a24). Mathematical work
    must already be reserved. Transport/setup are additionally charged here.
    publish runs in the Python owner only after the region barrier and checked
    export. A true return stops further slices and discards unfinished state.
    """

    def __init__(
        self, *, workers=2, schedule="dynamic", capacity=8, chunk_size=1
    ):
        if type(workers) is not int or workers not in (1, 2, 4):
            raise ValueError("native team supports 1, 2 or 4 workers")
        if schedule not in SCHEDULES:
            raise ValueError("unknown native team schedule")
        if type(capacity) is not int or not 1 <= capacity <= MAX_TASKS:
            raise ValueError("native team capacity must be 1..32")
        if type(chunk_size) is not int or not 1 <= chunk_size <= MAX_TASKS:
            raise ValueError("native team chunk size must be 1..32")
        self.chunk_size = chunk_size
        required = TEAM_OVERHEAD + capacity * RESERVE
        self.workers, self.schedule, self.capacity = (
            workers,
            schedule,
            capacity,
        )
        self.memory_bytes = required
        self.packets = []
        self.programs = ProgramStore()
        self.ffi = self.lib = self.contexts = self.statuses = self.actual = (
            None
        )
        self.lock = threading.Lock()
        self.closed = False
        self.owner = None
        self.slices = self.discarded = self.completed = 0
        self.lifecycle = (0, 0, 0, 0)

    def _grow(self, count, n, budget):
        while len(self.packets) < count:
            binding = (
                None
                if not self.packets
                else (self.ffi, self.lib, self.packets[0].manifest)
            )
            packet = NativePacket(
                _variant="pthread", _binding=binding, _programs=self.programs
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
                        raise MemoryError("pthread owner creation failed")
            except BaseException:
                # The private packet is registered before any foreign setup,
                # so interruption cannot orphan its native context or lock.
                if packet._lock.locked():
                    packet.end()
                packet.close()
                self.packets.pop()
                if not self.packets:
                    if self.ffi is not None:
                        if self.owner is not None:
                            self.lib.b4_pthread_destroy(self.owner)
                            self.owner = None
                        for owner in (
                            self.actual,
                            self.statuses,
                            self.contexts,
                        ):
                            if owner is not None:
                                self.ffi.release(owner)
                    self.contexts = self.statuses = self.actual = None
                    self.ffi = self.lib = None
                raise

    def execute(self, requests, budget, publish):
        if self.closed or not self.lock.acquire(blocking=False):
            raise RuntimeError("native team closed or active")
        if not _PROCESS_TEAM_LOCK.acquire(blocking=False):
            self.lock.release()
            raise RuntimeError(
                "another native team is already active in this process"
            )
        prepared = []
        try:
            if not 0 <= len(requests) <= self.capacity:
                raise ValueError("native team queue exceeds reservation")
            if not callable(publish) or any(
                not isinstance(request, (tuple, list)) or len(request) != 5
                for request in requests
            ):
                raise ValueError("invalid team request or publication owner")
            poll_budget(budget)
            if not requests:
                return
            self._grow(len(requests), requests[0][3], budget)
            if self.lib.b4_pthread_reset(self.owner):
                raise RuntimeError("pthread owner is active")
            for packet, request in zip(self.packets, requests):
                kind, payload, point, n, a24 = request
                if kind not in ("ladder", "records"):
                    raise ValueError("unknown native chunk kind")
                method = (
                    packet.prepare_ladder
                    if kind == "ladder"
                    else packet.prepare_records
                )
                prepared.append(packet)
                method(payload, point, n, a24, budget)
            pending = list(range(len(prepared)))
            while pending:
                poll_budget(budget)
                for slot, index in enumerate(pending):
                    self.contexts[slot] = prepared[index].context
                status = self.lib.b4_pthread_step(
                    self.owner,
                    self.contexts,
                    len(pending),
                    self.workers,
                    SCHEDULES[self.schedule],
                    self.chunk_size,
                    5,
                    128,
                    self.statuses,
                    self.actual,
                )
                # The generation has drained, including when CFFI raises a
                # signal on reentry. Parked workers no longer borrow contexts.
                self.slices += 1
                if status == -3:
                    raise BudgetExhaustedError("cancelled")
                if status == -4:
                    raise RuntimeError("pthread worker failed")
                if status == -2 or int(self.actual[0]) != self.workers:
                    raise TeamRefusedError("pthread team size/startup refused")
                if status or any(
                    self.statuses[i] not in (0, 1) for i in range(len(pending))
                ):
                    raise ValueError("native team returned an invalid status")
                poll_budget(budget)
                ready = []
                for slot, index in enumerate(pending):
                    if not self.statuses[slot]:
                        continue
                    kind, payload, _, n, _ = requests[index]
                    result = prepared[index].checked(
                        prepared[index]._export(),
                        n,
                        len(payload) if kind == "records" else None,
                    )
                    ready.append((index, result))
                # Validate the complete returned set before publishing any
                # result. Stable queue order resolves simultaneous successes.
                finished = {index for index, _ in ready}
                for offset, (index, result) in enumerate(ready):
                    if publish(index, result):
                        self.completed += offset + 1
                        self.discarded += len(pending) - offset - 1
                        return
                self.completed += len(ready)
                pending = [index for index in pending if index not in finished]
        finally:
            for packet in reversed(self.packets):
                if packet._lock.locked():
                    packet.end()
            _PROCESS_TEAM_LOCK.release()
            self.lock.release()

    def close(self):
        if not self.lock.acquire(blocking=False):
            raise RuntimeError("cannot close an active native team")
        try:
            if self.owner is not None:
                if self.lib.b4_pthread_shutdown(self.owner):
                    raise RuntimeError("cannot stop an active pthread team")
                self.lifecycle = tuple(
                    int(self.lib.b4_pthread_counter(self.owner, kind))
                    for kind in range(4)
                )
                self.lib.b4_pthread_destroy(self.owner)
                self.owner = None
            # Borrowing packets share the first packet's DSO handle. Destroy
            # them before closing that handle, exactly once, in its owner.
            for packet in reversed(self.packets):
                packet.close()
            self.packets.clear()
            if self.ffi is not None:
                self.programs.close(self.ffi)
                for owner in (self.actual, self.statuses, self.contexts):
                    if owner is not None:
                        self.ffi.release(owner)
            self.contexts = self.statuses = self.actual = None
            self.closed = True
        finally:
            self.lock.release()

    def cancel(self):
        """Cooperatively stop an active C dispatch, with no Python callback."""
        # close is forbidden during execute. Callers must keep this owner live
        # until their cancel call returns; close/cancel cannot race each other.
        if self.closed:
            raise RuntimeError("native team is closed")
        if self.owner is not None:
            self.lib.b4_pthread_cancel(self.owner)
