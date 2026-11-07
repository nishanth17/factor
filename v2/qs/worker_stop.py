"""Opt-in current-packet stop flags without a shared worker-reader lock."""

import multiprocessing
import os


class SingleReaderFlag:
    """One worker reads; only the parent sets and clears while quiescent."""

    def __init__(self, context):
        self.flag = context.Semaphore(0)
        self.raised = False
        self.owner_pid = os.getpid()

    def is_set(self):
        # Only this worker can temporarily consume the token. The parent
        # never removes it until all packets and snapshot barriers drain.
        if self.flag.acquire(False):
            self.flag.release()
            return True
        return False

    def set(self):
        if os.getpid() != self.owner_pid:
            raise ValueError("only the owner parent may set a stop flag")
        if not self.raised:
            self.flag.release()
            self.raised = True

    def clear(self):
        if os.getpid() != self.owner_pid:
            raise ValueError("only the owner parent may clear a stop flag")
        while self.flag.acquire(False):
            pass
        self.raised = False


class WorkerStops:
    """Parent broadcast; each worker reads its independently locked Event."""

    def __init__(self, events):
        self.events = events

    def set(self):
        for event in self.events:
            event.set()

    def clear(self):
        for event in self.events:
            event.clear()

    def is_set(self):
        return any(event.is_set() for event in self.events)

    def worker_event(self, slot):
        return self.events[slot]


class SemaphoreStops(WorkerStops):
    """Keep setup checks on the same single-reader flag as packet checks."""

    def is_set(self):
        from .parallel import _PROCESS_STATE

        if _PROCESS_STATE is not None and _PROCESS_STATE[0]["stop"] is self:
            return self.events[_PROCESS_STATE[1]].is_set()
        # Parent observation must not become a second semaphore reader.
        return any(event.raised for event in self.events)


def configure_stop(pool, private, backend="event"):
    """Pin before spawn; never replace handles already owned by a worker."""
    if backend not in ("event", "semaphore") or (
        backend == "semaphore" and not private
    ):
        raise ValueError("invalid private stop backend")
    current = (
        getattr(pool, "_qs_private_stop", False),
        getattr(pool, "_qs_stop_backend", "event"),
    )
    if (private, backend) == current:
        return
    if pool.closed or pool.lock.locked():
        raise ValueError("cannot change active or closed pool stop flags")
    if pool.mode != "process" or pool._owned_processes:
        raise ValueError("stop policy requires a fresh process pool")
    context = multiprocessing.get_context("spawn")
    if backend == "semaphore":
        pool.state["stop"] = SemaphoreStops(
            [SingleReaderFlag(context) for _ in range(pool.workers)]
        )
    else:
        pool.state["stop"] = (
            WorkerStops([context.Event() for _ in range(pool.workers)])
            if private
            else context.Event()
        )
    pool._qs_private_stop = private
    pool._qs_stop_backend = backend
