"""One C-only OpenMP team for bounded, pre-funded ECM continuation chunks."""

import threading

from .native_packet import RESERVE, NativePacket, ProgramStore, poll_budget

SCHEDULES = {"dynamic": 0, "static": 1, "guided": 2}
MAX_TASKS = 32
TEAM_OVERHEAD = 64 * 2**20
_PROCESS_TEAM_LOCK = threading.Lock()


class TeamRefusedError(RuntimeError):
    """The runtime refused the funded team size before native arithmetic."""


class NativeTeam:
    """No Python workers/callbacks inside a parallel region; one active owner.

    Requests are (kind, scalar-or-entries, point, n, a24). Mathematical work
    must already be reserved. Transport/setup are additionally charged here.
    publish runs in the Python owner only after the region barrier and checked
    export. A true return stops further slices and discards unfinished state.
    """

    def __init__(
        self, *, workers=2, schedule="dynamic", capacity=8, memory_bytes=None
    ):
        if type(workers) is not int or workers not in (1, 2, 4):
            raise ValueError("native team supports 1, 2 or 4 workers")
        if schedule not in SCHEDULES:
            raise ValueError("unknown native team schedule")
        if type(capacity) is not int or not 1 <= capacity <= MAX_TASKS:
            raise ValueError("native team capacity must be 1..32")
        required = TEAM_OVERHEAD + capacity * RESERVE
        if (
            memory_bytes is None
            or type(memory_bytes) is not int
            or memory_bytes < required
        ):
            raise MemoryError("native team memory reservation refused")
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
        self.slices = self.discarded = self.completed = 0

    def _grow(self, count, n, budget):
        while len(self.packets) < count:
            binding = (
                None
                if not self.packets
                else (self.ffi, self.lib, self.packets[0].manifest)
            )
            packet = NativePacket(
                _variant="openmp", _binding=binding, _programs=self.programs
            )
            self.packets.append(packet)
            try:
                packet.available(n, budget)
                if self.ffi is None:
                    self.ffi, self.lib = packet.ffi, packet.lib
                    self.contexts = self.ffi.new("b4_context *[32]")
                    self.statuses = self.ffi.new("int[32]")
                    self.actual = self.ffi.new("unsigned *")
            except BaseException:
                # The private packet is registered before any foreign setup,
                # so interruption cannot orphan its native context or lock.
                if packet._lock.locked():
                    packet.end()
                packet.close()
                self.packets.pop()
                if not self.packets:
                    if self.ffi is not None:
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
            if not 1 <= len(requests) <= self.capacity:
                raise ValueError("native team queue exceeds reservation")
            if not callable(publish) or any(
                not isinstance(request, (tuple, list)) or len(request) != 5
                for request in requests
            ):
                raise ValueError("invalid team request or publication owner")
            self._grow(len(requests), requests[0][3], budget)
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
                status = self.lib.b4_team_step(
                    self.contexts,
                    len(pending),
                    self.workers,
                    SCHEDULES[self.schedule],
                    5,
                    128,
                    self.statuses,
                    self.actual,
                )
                # No region survives this return, including a signal raised
                # by CFFI on reentry. Borrowed data can now be released.
                self.slices += 1
                if status == -2 or int(self.actual[0]) != self.workers:
                    raise TeamRefusedError("OpenMP team size/startup refused")
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
            for packet in self.packets:
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
