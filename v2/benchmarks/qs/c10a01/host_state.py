"""Bounded host-state receipts and a native-code execution-rate canary."""

import ctypes
import gc
import hashlib
import os
import statistics
import subprocess
import sys
import time

CANARY_ITERATIONS = 200_000


def parse_power(source, profiles):
    """Extract active power source and its configured low-power setting."""
    active = None
    for line in source.splitlines():
        if line.startswith("Now drawing from '") and line.endswith("'"):
            active = line.split("'", 2)[1]
            break
    modes, section = {}, None
    for line in profiles.splitlines():
        if line.endswith("Power:") and not line.startswith(" "):
            section = line[:-1]
        fields = line.split()
        if len(fields) == 2 and fields[0] == "lowpowermode":
            modes[section] = int(fields[1])
    if active not in modes or modes[active] not in (0, 1):
        raise RuntimeError("active power mode is unavailable")
    return dict(source=active, low_power=modes[active])


def power_state():
    if sys.platform != "darwin":
        raise RuntimeError("this host-state adapter requires macOS")
    outputs = [
        subprocess.check_output(["pmset", "-g", option], text=True, timeout=3)
        for option in ("batt", "custom")
    ]
    if sum(len(value.encode()) for value in outputs) > 16384:
        raise RuntimeError("recorder_exhaustion: host state")
    return parse_power(*outputs)


def runtime_state():
    """Read requested thread QoS and GC settings without changing either."""
    library = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
    library.qos_class_self.argtypes = []
    library.qos_class_self.restype = ctypes.c_uint
    return dict(
        qos=library.qos_class_self(),
        nursery=gc.get_stats(memory_pressure=False).nursery_size,
        environment={
            name: value
            for name, value in os.environ.items()
            if name.startswith("PYPY_GC_") or name == "PYPYJIT"
        },
    )


def native_canary():
    """Use a fixed OpenSSL workload, outside PyPy's compiled Python loops."""
    samples, reference = [], None
    for _ in range(9):
        wall, cpu = time.monotonic(), time.process_time()
        value = hashlib.pbkdf2_hmac(
            "sha256", b"c10a01-host-canary", b"fixed-salt", CANARY_ITERATIONS
        ).hex()
        samples.append(
            dict(wall=time.monotonic() - wall, cpu=time.process_time() - cpu)
        )
        if reference is None:
            reference = value
        elif value != reference:
            raise ArithmeticError("native canary output changed")
    return dict(
        digest=reference,
        iterations=CANARY_ITERATIONS,
        samples=samples,
        medians={
            key: statistics.median(s[key] for s in samples)
            for key in ("wall", "cpu")
        },
        # This is a separate diagnostic charge, not SIQS Budget.used.
        work=9 * CANARY_ITERATIONS,
    )


def matches(before, after):
    """A source/mode change invalidates clean timing without losing factors."""
    return (
        before["power"] == after["power"]
        and before["runtime"] == after["runtime"]
    )
