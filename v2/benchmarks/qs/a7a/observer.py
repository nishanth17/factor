"""Optional two-Hz RSS samples of this call's owned processes only."""

import ctypes
import errno
import os
import subprocess
import sys
import threading
import time


class _TaskInfo(ctypes.Structure):
    # Darwin sys/proc_info.h: six uint64_t followed by twelve int32_t.
    # This declaration reads the existing OS ABI; it allocates no task data.
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
            "context_switches",
            "thread_count",
            "running_threads",
            "priority",
        )
    ]


class DarwinRSS:
    """Read resident bytes for supplied owned PIDs through libproc."""

    def __init__(self):
        if sys.platform != "darwin" or ctypes.sizeof(_TaskInfo) != 96:
            raise RuntimeError("unsupported Darwin task-info ABI")
        self.library = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
        self.query = self.library.proc_pidinfo
        self.query.argtypes = [
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_uint64,
            ctypes.c_void_p,
            ctypes.c_int,
        ]
        self.query.restype = ctypes.c_int

    def __call__(self, pids):
        rss, missing = {}, []
        for pid in pids:
            info = _TaskInfo()
            ctypes.set_errno(0)
            size = self.query(
                pid, 4, 0, ctypes.pointer(info), ctypes.sizeof(info)
            )
            if size != ctypes.sizeof(info):
                code = ctypes.get_errno()
                if not size and code == errno.ESRCH:
                    missing.append(pid)
                    continue
                raise OSError(
                    code or errno.EIO,
                    f"proc_pidinfo({pid}) returned {size} task-info bytes",
                )
            rss[pid] = int(info.resident_size)
        return rss, missing


def validate_observation(result):
    """A failed recorder must not masquerade as production cancellation."""
    if result.get("observer_failure") is not None:
        raise RuntimeError(
            "RSS recorder failure: " + result["observer_failure"]
        )
    if result.get("observer_exhausted", False):
        raise RuntimeError("diagnostic owned RSS exceeded its allowance")
    if result["cell"].get("observe", False) and not result.get("rss_samples"):
        raise RuntimeError("RSS recorder emitted no samples")


class OwnedRSS:
    """Bounded outer samples; child HWM sums are never called aggregate RSS."""

    def __init__(self, pool, *, limit, samples=16000):
        self.pool, self.limit, self.maximum = pool, limit, samples
        self.rows, self.error, self.exceeded = [], None, False
        self.reader = DarwinRSS() if sys.platform == "darwin" else None
        self.backend = "darwin-libproc" if self.reader is not None else "ps"
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        while not self.stop.is_set():
            processes = []
            if self.pool is not None and self.pool.mode == "process":
                mapping = self.pool.executor._processes or {}
                processes = list(mapping.copy().values())
            pids = [os.getpid()] + [p.pid for p in processes if p.pid]
            try:
                from multiprocessing import resource_tracker

                helper = resource_tracker._resource_tracker._pid
                if helper is not None:
                    pids.append(helper)
                if self.reader is not None:
                    rss, missing = self.reader(pids)
                else:
                    result = subprocess.run(
                        [
                            "ps",
                            "-o",
                            "pid=,rss=",
                            "-p",
                            ",".join(map(str, pids)),
                        ],
                        capture_output=True,
                        text=True,
                        timeout=1,
                    )
                    values = [
                        line.split() for line in result.stdout.splitlines()
                    ]
                    rss = {
                        int(pid): 1024 * int(value) for pid, value in values
                    }
                    missing = sorted(set(pids) - rss.keys())
                if os.getpid() not in rss:
                    raise RuntimeError("RSS recorder lost its parent process")
                total = sum(rss.values())
                self.rows.append(
                    dict(
                        time=time.monotonic(),
                        bytes=total,
                        processes=rss,
                        helper=helper,
                        missing=missing,
                        backend=self.backend,
                    )
                )
                if total > self.limit:
                    self.exceeded = True
                    return
                if len(self.rows) >= self.maximum:
                    self.error = "RSS recorder sample allowance exhausted"
                    return
            except Exception as error:
                self.error = str(error)
                return
            self.stop.wait(0.5)

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *args):
        self.stop.set()
        self.thread.join(timeout=2)
