"""Bounded, read-only Mach observations of this owner's idle native threads.

The SDK's THREAD_BASIC_INFO reports time_value seconds and microseconds.
Opaque pthread handles are translated only while their owning team is alive;
ports and ranks do not identify physical CPU cores.
"""

import ctypes
import time


class _TimeValue(ctypes.Structure):
    _fields_ = [("seconds", ctypes.c_int), ("microseconds", ctypes.c_int)]

    def value(self):
        if self.seconds < 0 or not 0 <= self.microseconds < 1_000_000:
            raise ValueError("invalid Mach thread CPU time")
        return self.seconds + self.microseconds / 1_000_000


class _ThreadInfo(ctypes.Structure):
    _fields_ = [
        ("user", _TimeValue),
        ("system", _TimeValue),
        ("cpu_usage", ctypes.c_int),
        ("policy", ctypes.c_int),
        ("run_state", ctypes.c_int),
        ("flags", ctypes.c_int),
        ("suspend_count", ctypes.c_int),
        ("sleep_time", ctypes.c_int),
    ]


class NativeThreadProbe:
    """Observe a live drained pthread team; never inspect freed handles."""

    def __init__(self, team):
        if team.closed or team.variant != "pthread":
            raise ValueError("thread observation requires a live pthread team")
        self.team = team
        activity = team.worker_activity()
        if not activity["actual"]:
            raise ValueError("thread observation requires started workers")
        self.identities = [row["identity"] for row in activity["workers"]]
        self.system = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
        self.system.pthread_mach_thread_np.argtypes = [ctypes.c_void_p]
        self.system.pthread_mach_thread_np.restype = ctypes.c_uint
        self.system.thread_info.argtypes = [
            ctypes.c_uint,
            ctypes.c_int,
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_uint),
        ]
        self.system.thread_info.restype = ctypes.c_int
        self.ports = [
            self.system.pthread_mach_thread_np(token)
            for token in self.identities
        ]
        if not all(self.ports) or len(set(self.ports)) != len(self.ports):
            raise ValueError("native worker ports are missing or repeated")
        if ctypes.sizeof(_ThreadInfo) != 40:
            raise ValueError("unexpected SDK thread-basic layout")

    def snapshot(self, label):
        if self.team.closed:
            raise RuntimeError("native owner closed before thread observation")
        activity = self.team.worker_activity()
        if self.identities != [row["identity"] for row in activity["workers"]]:
            raise RuntimeError("native compute-worker identities changed")
        rows = []
        for token, port in zip(self.identities, self.ports):
            info = _ThreadInfo()
            count = ctypes.c_uint(ctypes.sizeof(info) // 4)
            error = self.system.thread_info(
                port, 3, ctypes.byref(info), ctypes.byref(count)
            )
            if error or count.value != ctypes.sizeof(info) // 4:
                raise OSError(error, "Mach thread_info failed")
            rows.append(
                dict(
                    identity=token,
                    port=port,
                    cpu=info.user.value() + info.system.value(),
                    run_state=info.run_state,
                )
            )
        return dict(
            label=label,
            at=time.monotonic(),
            threads=rows,
            activity=activity,
        )
