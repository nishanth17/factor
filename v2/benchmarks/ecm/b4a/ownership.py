"""macOS owned-group watchdog, isolated from the arithmetic interpreter.

The SDK's public libproc ABI supplies process CPU (including native threads),
live RSS and reaped-child CPU. Unsupported/denied observation is a harness
failure. Process-wide counters already include threads; never add thread CPU.
"""

import ctypes
import json
import os
import select
import signal
import sys
import time
from pathlib import Path


class Usage(ctypes.Structure):
    _fields_ = [("uuid", ctypes.c_ubyte * 16)] + [
        (name, ctypes.c_uint64)
        for name in (
            "user",
            "system",
            "idle",
            "interrupts",
            "pageins",
            "wired",
            "rss",
            "footprint",
            "start",
            "exit",
            "child_user",
            "child_system",
            "child_idle",
            "child_interrupts",
            "child_pageins",
            "child_elapsed",
            "read",
            "written",
        )
    ]


class TaskInfo(ctypes.Structure):
    """Public proc_taskinfo layout from the installed macOS SDK."""

    _fields_ = [
        (name, ctypes.c_uint64)
        for name in (
            "virtual",
            "resident",
            "user",
            "system",
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
            "sent",
            "received",
            "mach_calls",
            "unix_calls",
            "switches",
            "thread_count",
            "running",
            "priority",
        )
    ]


class GroupObserver:
    """Bounded PID inventory; no ps subprocesses or interpreter callbacks."""

    def __init__(self):
        if sys.platform != "darwin":
            raise RuntimeError("B4a ownership observer requires macOS")
        self.lib = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
        self.lib.proc_listpids.argtypes = [
            ctypes.c_uint32,
            ctypes.c_uint32,
            ctypes.c_void_p,
            ctypes.c_int,
        ]
        self.lib.proc_listpids.restype = ctypes.c_int
        self.lib.proc_pid_rusage.argtypes = [
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_void_p,
        ]
        self.lib.proc_pid_rusage.restype = ctypes.c_int
        self.lib.proc_pidinfo.argtypes = [
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_uint64,
            ctypes.c_void_p,
            ctypes.c_int,
        ]
        self.lib.proc_pidinfo.restype = ctypes.c_int

    def sample(self, group):
        pids = (ctypes.c_int * 64)()
        ctypes.set_errno(0)
        size = self.lib.proc_listpids(2, group, pids, ctypes.sizeof(pids))
        if size < 0 or size >= ctypes.sizeof(pids):
            raise RuntimeError("owned PID recorder refused/exhausted")
        if size == 0:
            try:
                os.killpg(group, 0)
            except ProcessLookupError:
                return dict(pids=[], cpu=0, rss=0, threads=0)
            except PermissionError as error:
                raise RuntimeError("owned PID observation denied") from error
            raise RuntimeError(
                "empty PID observation for an existing owned group"
            )
        inventory = [
            pid
            for pid in pids[: size // ctypes.sizeof(ctypes.c_int)]
            if pid > 0
        ]
        if not inventory:
            raise RuntimeError("invalid nonempty PID observation")
        records, threads = [], 0
        for pid in inventory:
            usage = Usage()
            if self.lib.proc_pid_rusage(pid, 2, ctypes.byref(usage)):
                if (
                    ctypes.get_errno() == 3
                ):  # ESRCH: process exited during sample.
                    continue
                raise RuntimeError("owned process CPU/RSS observation denied")
            info = TaskInfo()
            size = self.lib.proc_pidinfo(
                pid, 4, 0, ctypes.byref(info), ctypes.sizeof(info)
            )
            if size != ctypes.sizeof(info):
                if ctypes.get_errno() == 3:
                    continue
                raise RuntimeError("owned thread inventory denied/truncated")
            if not 0 <= info.thread_count <= 4096:
                raise RuntimeError("invalid owned thread count")
            threads += info.thread_count
            records.append((pid, usage))
        # Reaped descendants are counted by their still-live ancestor. A live
        # child has not yet been added to its parent's reaped-child counters.
        cpu = (
            sum(
                u.user + u.system + u.child_user + u.child_system
                for _, u in records
            )
            / 1e9
        )
        return dict(
            pids=inventory,
            cpu=cpu,
            rss=sum(u.rss for _, u in records),
            threads=threads,
        )


def terminate(group):
    try:
        os.killpg(group, signal.SIGTERM)
    except ProcessLookupError:
        return


def guard(request, pipe, receipt):
    """Parent EOF or any limit stops the entire owned group, then drains it."""
    observer = GroupObserver()
    started = time.monotonic()
    reason, peak_rss, peak_cpu, peak_processes = None, 0, 0, 0
    empty, peak_threads = False, 0
    try:
        if request.get("ready"):
            Path(request["ready"]).write_text("ready\n")
        while True:
            sample = observer.sample(request["pgid"])
            if not sample["pids"]:
                empty = True
                break
            peak_rss = max(peak_rss, sample["rss"])
            peak_cpu = max(peak_cpu, sample["cpu"])
            peak_processes = max(peak_processes, len(sample["pids"]))
            peak_threads = max(peak_threads, sample["threads"])
            output = sum(
                Path(path).stat().st_size
                for path in request["outputs"]
                if Path(path).exists()
            )
            if time.monotonic() >= request["deadline"]:
                reason = "wall_limit"
            elif peak_cpu >= request["cpu"]:
                reason = "cpu_limit"
            elif peak_rss > request["rss"]:
                reason = "rss_limit"
            elif peak_processes > request["processes"]:
                reason = "process_limit"
            elif peak_threads > request.get("threads", 32):
                reason = "thread_limit"
            elif output > request["output"]:
                reason = "recorder_output_limit"
            elif select.select([pipe], [], [], 0.05)[0] and not os.read(
                pipe, 1
            ):
                reason = "controller_closed"
            if reason:
                break
    except BaseException as error:
        reason = "recorder_failure: " + repr(error)
    finally:
        if not empty:
            terminate(request["pgid"])
            deadline = time.monotonic() + request["shutdown"]
            while time.monotonic() < deadline:
                try:
                    if not observer.sample(request["pgid"])["pids"]:
                        empty = True
                        break
                except RuntimeError:
                    pass
                time.sleep(0.02)
            if not empty:
                try:
                    os.killpg(request["pgid"], signal.SIGKILL)
                except ProcessLookupError:
                    pass
                # The parent must additionally verify group disappearance
                # after reaping the command. A sent signal is not that proof.
        Path(receipt).write_text(
            json.dumps(
                dict(
                    reason=reason,
                    group_empty=empty,
                    sampled_peak_rss=peak_rss,
                    sampled_cpu=peak_cpu,
                    sampled_peak_processes=peak_processes,
                    sampled_peak_threads=peak_threads,
                    wall=time.monotonic() - started,
                    cpu_scope=(
                        "process CPU incl native threads + reaped "
                        "descendants; live snapshots"
                    ),
                ),
                indent=2,
            )
            + "\n"
        )
