"""All-worker bounded recovery and sampled stop-poll diagnostics only."""

import json
import os
import time
from functools import wraps
from pathlib import Path

_observer = None


class PacketCounters:
    """Aggregate packets, retaining only eight longest recovery calls."""

    def __init__(self, interval=256):
        if interval < 1 or interval & (interval - 1):
            raise ValueError("poll sample interval must be a power of two")
        self.interval = interval
        self.active = False

    def begin(self):
        self.values = {k: [0, 0, 0] for k in ("recovery", "split")}
        self.depth = {k: 0 for k in self.values}
        self.longest = {k: [] for k in self.values}
        self.poll_calls = self.poll_samples = self.poll_wall = (
            self.poll_cpu
        ) = 0
        self.poll_histogram = [0] * 64
        self.poll_max = 0
        self.active = True

    def wrap(self, name, function):
        @wraps(function)
        def measured(*args, **kwargs):
            if not self.active or self.depth[name]:
                return function(*args, **kwargs)
            self.depth[name] += 1
            wall, cpu = time.monotonic_ns(), time.process_time_ns()
            try:
                return function(*args, **kwargs)
            finally:
                cpu = time.process_time_ns() - cpu
                wall = time.monotonic_ns() - wall
                self.depth[name] -= 1
                row = self.values[name]
                row[0] += 1
                row[1] += wall
                row[2] += cpu
                longest = self.longest[name]
                if len(longest) < 8 or wall > longest[-1][0]:
                    longest.append([wall, cpu])
                    longest.sort(reverse=True)
                    del longest[8:]

        return measured

    def stop_poll(self, function):
        if not self.active:
            return function()
        self.poll_calls += 1
        if self.poll_calls & (self.interval - 1):
            return function()
        wall, cpu = time.monotonic_ns(), time.process_time_ns()
        try:
            return function()
        finally:
            cpu = time.process_time_ns() - cpu
            wall = time.monotonic_ns() - wall
            self.poll_samples += 1
            self.poll_wall += wall
            self.poll_cpu += cpu
            self.poll_histogram[min(63, wall.bit_length())] += 1
            self.poll_max = max(self.poll_max, wall)

    def finish(self):
        self.active = False
        return dict(
            recovery=self.values["recovery"],
            split=self.values["split"],
            longest=self.longest,
            poll_calls=self.poll_calls,
            poll_samples=self.poll_samples,
            poll_wall_ns=self.poll_wall,
            poll_cpu_ns=self.poll_cpu,
            poll_histogram=self.poll_histogram,
            poll_max_ns=self.poll_max,
            poll_interval=self.interval,
        )


class Recorder:
    """One bounded append-only journal per owned process."""

    def __init__(self, role):
        self.folder = Path(os.environ["FACTOR_A7A_PACKET_DIR"])
        self.path = self.folder / (str(os.getpid()) + ".jsonl")
        self.stream = self.path.open("x")
        self.rows = self.size = 0
        self.write("ready", role=role, pid=os.getpid(), ppid=os.getppid())

    def write(self, kind, **values):
        row = json.dumps(dict(kind=kind, at=time.monotonic(), **values)) + "\n"
        if self.rows >= 16384 or self.size + len(row) > 4 * 2**20:
            (self.folder / (str(os.getpid()) + ".failure.json")).write_text(
                json.dumps(dict(error="packet recorder allowance exhausted"))
            )
            raise RuntimeError("packet recorder allowance exhausted")
        self.stream.write(row)
        self.stream.flush()
        self.rows += 1
        self.size += len(row)


def install_worker(state):
    global _observer
    from ....qs import cofactor_policy, worker_packet, yield_collector

    recorder, counters = Recorder("worker"), PacketCounters()
    _observer = counters
    # Lazy imports in expanded recovery see these wrappers. Depth prevents
    # split_residual -> split_cofactor from counting the same span twice.
    for name in ("split_cofactor", "split_residual"):
        original = getattr(cofactor_policy, name)
        wrapper = counters.wrap("split", original)
        setattr(cofactor_policy, name, wrapper)
        if getattr(yield_collector, name, None) is original:
            setattr(yield_collector, name, wrapper)
    yield_collector.divide_yield_residual = counters.wrap(
        "recovery", yield_collector.divide_yield_residual
    )
    stop = state["stop"]
    if hasattr(stop, "worker_event"):
        from ....qs.parallel import _PROCESS_STATE

        stop = stop.worker_event(_PROCESS_STATE[1])
    original_stop = stop.is_set
    stop.is_set = lambda: counters.stop_poll(original_stop)
    original_packet = worker_packet.collect_packet

    def packet(task, **kwargs):
        counters.begin()
        wall, cpu = time.monotonic(), time.process_time()
        report = None
        try:
            report = original_packet(task, **kwargs)
            return report
        finally:
            recorder.write(
                "packet",
                job=task["job"],
                index=task["id"],
                token=task["token"],
                start=wall,
                wall=time.monotonic() - wall,
                cpu=time.process_time() - cpu,
                reason=None if report is None else report["reason"],
                report_cpu=None if report is None else report["cpu_seconds"],
                counters=counters.finish(),
            )

    worker_packet.collect_packet = packet


def initialize_worker(state):
    from ....qs.parallel import _initialize_process

    _initialize_process(state)
    install_worker(state)


def validate_capture(folder, workers):
    if list(Path(folder).glob("*.failure.json")):
        raise RuntimeError("packet recorder exhausted its allowance")
    roles, packets = [], 0
    for path in Path(folder).glob("*.jsonl"):
        with path.open() as stream:
            first = json.loads(next(stream))
            if first["kind"] != "ready":
                raise RuntimeError("packet recorder missing owner identity")
            roles.append(first["role"])
            for line in stream:
                row = json.loads(line)
                if row["kind"] == "packet":
                    packets += 1
    if (
        roles.count("worker") != workers
        or roles.count("parent") != 1
        or not packets
    ):
        raise RuntimeError("packet recorder lacks all owned workers")


def main():
    from ....qs import parallel
    from ....qs.current_parallel import CurrentParallelSIQSJob
    from .cell import main as cell_main

    recorder = Recorder("parent")
    parallel._initialize_process = initialize_worker
    for name in ("_receive", "_merge"):
        original = getattr(CurrentParallelSIQSJob, name)

        def bind(function, category):
            def observed(owner, item, *args, **kwargs):
                result = function(owner, item, *args, **kwargs)
                recorder.write(
                    category,
                    job=owner.identity,
                    index=item["id"] if category == "receive" else item,
                    frontier=owner.next_assignment,
                    pending=len(owner.pending),
                )
                return result

            return observed

        setattr(CurrentParallelSIQSJob, name, bind(original, name[1:]))
    cell_main()


if __name__ == "__main__":
    main()
