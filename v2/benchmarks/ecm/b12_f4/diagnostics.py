"""Opt-in observation of owned ECM processes; never imported by production."""

import gc
import json
import os
import resource
import threading
import time
from pathlib import Path


class Capture:
    """Bounded process samples and optional RPython heap dumps."""

    def __init__(self, request, role):
        self.request = request
        self.mode = request["observer"]
        self.role = role
        self.directory = Path(request["capture"]) / str(os.getpid())
        self.directory.mkdir(parents=True, exist_ok=False)
        self.marker = Path(request["capture"]) / "measuring"
        self.profile = None
        self.phase = "warmup/idle"
        self.measured = False
        self.start_wall = time.monotonic()
        self.start_cpu = time.process_time()
        self.phase_wall = self.start_wall
        self.phase_cpu = self.start_cpu
        self.phases = {}
        self.jobs = 0
        self.heaps = []
        self.sampling = False
        self.gc_counts = dict(minor=0, major_steps=0, major=0)
        self.gc_seconds = dict(minor=0.0, major_steps=0.0)
        self.dump_seconds = 0.0
        self.gc_by_phase = {}
        self.before = self.memory()
        if self.mode in ("counters", "heap"):
            gc.hooks.on_gc_minor = self.minor
            gc.hooks.on_gc_collect_step = self.major_step
            gc.hooks.on_gc_collect = self.major
        if self.mode == "cpu" and request["mode"] == "cold":
            self.start_profile()

    def memory(self):
        stats = gc.get_stats()
        # Numeric fields follow this pinned PyPy build's GcStats.
        raw = stats._s
        names = (
            "total_gc_memory",
            "total_allocated_memory",
            "peak_memory",
            "total_arena_memory",
            "total_rawmalloced_memory",
            "nursery_size",
            "jit_backend_used",
            "jit_backend_allocated",
            "total_gc_time",
        )
        return {name: getattr(raw, name) for name in names}

    def start_profile(self):
        if self.sampling:
            return
        import _vmprof

        self.profile = (self.directory / "cpu.vmprof").open("xb")
        # CPU timer at 2 ms; PyPy supplies symbols. Extra tracing stays off.
        _vmprof.enable(self.profile.fileno(), 0.002, 0, 0, 0, 0)
        self.sampling = True

    def stop_profile(self):
        if self.sampling:
            import _vmprof

            _vmprof.disable()
            self.profile.close()
            self.profile = None
            self.sampling = False

    def minor(self, stats):
        self.gc_counts["minor"] += stats.count
        self.gc_seconds["minor"] += stats.duration
        row = self.gc_by_phase.setdefault(self.phase, [0, 0.0, 0, 0.0])
        row[0] += stats.count
        row[1] += stats.duration

    def major_step(self, stats):
        self.gc_counts["major_steps"] += stats.count
        self.gc_seconds["major_steps"] += stats.duration
        row = self.gc_by_phase.setdefault(self.phase, [0, 0.0, 0, 0.0])
        row[2] += stats.count
        row[3] += stats.duration

    def major(self, stats):
        self.gc_counts["major"] += stats.count

    def transition(self, phase):
        now, cpu = time.monotonic(), time.process_time()
        previous = self.phases.setdefault(self.phase, [0.0, 0.0, 0])
        previous[0] += now - self.phase_wall
        previous[1] += cpu - self.phase_cpu
        previous[2] += 1
        label = ("measured/" if self.measured else "warmup/") + phase
        self.phase, self.phase_wall, self.phase_cpu = label, now, cpu
        if (
            self.mode == "heap"
            and phase == "stage_two"
            and self.marker.exists()
            and not any(h["label"] == "stage_two_live" for h in self.heaps)
        ):
            self.dump("stage_two_live")

    def begin_job(self):
        self.measured = self.marker.exists()
        self.jobs += 1
        if self.jobs > 2048:
            raise RuntimeError("recorder job allowance exhausted")
        if self.mode == "cpu":
            if not self.sampling and self.marker.exists():
                self.start_profile()
        else:
            self.transition("prepare")

    def end_job(self):
        if self.mode in ("counters", "heap"):
            self.transition("idle")
        if (
            self.mode == "heap"
            and any(h["label"] == "stage_two_live" for h in self.heaps)
            and not any(h["label"] == "after_curve" for h in self.heaps)
        ):
            self.dump("after_curve")

    def dump(self, label):
        if len(self.heaps) >= 4:
            raise RuntimeError("recorder heap count exhausted")
        path = self.directory / (label + ".heap")
        if path.exists():
            raise RuntimeError("heap snapshot cannot overwrite evidence")
        before, started, cpu = (
            self.memory(),
            time.monotonic(),
            time.process_time(),
        )
        # Heap walking forces collection and perturbs execution timing.
        gc.dump_rpy_heap(str(path))
        elapsed = time.monotonic() - started
        self.dump_seconds += elapsed
        self.heaps.append(
            dict(
                label=label,
                file=str(path),
                bytes=path.stat().st_size,
                seconds=elapsed,
                cpu=time.process_time() - cpu,
                before=before,
                after=self.memory(),
            )
        )

    def finish(self):
        self.stop_profile()
        if self.mode in ("counters", "heap"):
            self.transition("finished")
            gc.hooks.reset()
        result = dict(
            pid=os.getpid(),
            role=self.role,
            observer=self.mode,
            elapsed=time.monotonic() - self.start_wall,
            cpu=time.process_time() - self.start_cpu,
            jobs=self.jobs,
            phases=self.phases,
            gc_counts=self.gc_counts,
            gc_seconds=self.gc_seconds,
            gc_by_phase=self.gc_by_phase,
            memory_before=self.before,
            memory_after=self.memory(),
            heaps=self.heaps,
            dump_seconds=self.dump_seconds,
            hwm_rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        )
        with (self.directory / "process.json").open("x") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)


def install_engine_hooks(capture):
    """Observe stages while preserving mathematical state and results."""
    from ....execution.ecm_worker import CurveEngine
    from ....execution.ecm_worker_state import ActionBudget

    execute, enter = CurveEngine.execute, ActionBudget.enter_phase

    def observed_execute(engine, *args, **kwargs):
        capture.begin_job()
        try:
            return execute(engine, *args, **kwargs)
        finally:
            capture.end_job()

    def observed_enter(budget, phase):
        if phase != budget.phase:
            capture.transition(phase)
        return enter(budget, phase)

    CurveEngine.execute = observed_execute
    if capture.mode in ("counters", "heap"):
        ActionBudget.enter_phase = observed_enter


def captured_worker(*args):
    """Spawn target observing this worker and its own import costs."""
    request = json.loads(Path(os.environ["FACTOR_ECM_DIAGNOSTIC"]).read_text())
    capture = Capture(request, "worker")
    try:
        from ....execution.ecm_worker import worker_main

        install_engine_hooks(capture)
        worker_main(*args)
    finally:
        capture.finish()


class RSSSampler:
    """Sample known owned PIDs via macOS libproc; never enumerate processes."""

    def __init__(self, executor):
        import ctypes

        class TaskInfo(ctypes.Structure):
            _fields_ = [
                (name, ctypes.c_uint64)
                for name in (
                    "virtual_size",
                    "resident_size",
                    "total_user",
                    "total_system",
                    "threads_user",
                    "threads_system",
                )
            ] + [
                (name, ctypes.c_int32)
                for name in (
                    "policy",
                    "faults",
                    "pageins",
                    "cow_faults",
                    "messages_sent",
                    "messages_received",
                    "syscalls_mach",
                    "syscalls_unix",
                    "csw",
                    "threadnum",
                    "numrunning",
                    "priority",
                )
            ]

        self.ctypes, self.task_type = ctypes, TaskInfo
        self.libproc = ctypes.CDLL("/usr/lib/libproc.dylib")
        self.libproc.proc_pidinfo.argtypes = [
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_uint64,
            ctypes.c_void_p,
            ctypes.c_int,
        ]
        self.libproc.proc_pidinfo.restype = ctypes.c_int
        self.executor = executor
        self.stop_event = threading.Event()
        self.samples, self.error, self.over_limit = [], None, False
        self.started = time.monotonic()
        self.thread = threading.Thread(target=self.run, daemon=True)

    def run(self):
        try:
            while not self.stop_event.is_set():
                if len(self.samples) >= 1000:
                    raise RuntimeError(
                        "RSS recorder sample allowance exhausted"
                    )
                pids = [os.getpid()] + [
                    slot["process"].pid
                    for slot in self.executor.slots
                    if slot["process"].pid is not None
                ]
                values = {}
                for pid in pids:
                    task = self.task_type()
                    size = self.ctypes.sizeof(task)
                    copied = self.libproc.proc_pidinfo(
                        pid, 4, 0, self.ctypes.byref(task), size
                    )
                    if copied == size:
                        values[pid] = task.resident_size
                total = sum(values.values())
                self.samples.append(
                    dict(
                        elapsed=time.monotonic() - self.started,
                        rss=values,
                        total=total,
                        complete=len(values) == len(pids),
                    )
                )
                self.over_limit |= total > 8 * 2**30
                self.stop_event.wait(0.05)
        except BaseException as error:
            self.error = repr(error)
            self.over_limit = True

    def close(self):
        self.stop_event.set()
        self.thread.join(timeout=1)
        if self.thread.is_alive() or self.error:
            raise RuntimeError("RSS recorder failure: " + str(self.error))
