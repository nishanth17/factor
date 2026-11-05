"""Process RSS high-water reads without changing observation cadence."""

import os
import platform
import resource
import struct
import sys
import threading


def _darwin_reader():
    """Use the stable 64-bit Darwin getrusage ABI when CFFI is available."""
    from cffi import FFI

    ffi = FFI()
    ffi.cdef(
        """
        struct timeval { long tv_sec; int tv_usec; };
        struct rusage {
            struct timeval ru_utime; struct timeval ru_stime;
            long ru_maxrss; long ru_ixrss; long ru_idrss; long ru_isrss;
            long ru_minflt; long ru_majflt; long ru_nswap; long ru_inblock;
            long ru_oublock; long ru_msgsnd; long ru_msgrcv; long ru_nsignals;
            long ru_nvcsw; long ru_nivcsw;
        };
        int getrusage(int, struct rusage *);
        """
    )
    if (
        ffi.sizeof("struct rusage") != 144
        or ffi.offsetof("struct rusage", "ru_maxrss") != 32
    ):
        raise ValueError("unsupported Darwin rusage layout")
    library = ffi.dlopen("/usr/lib/libSystem.B.dylib")
    buffers = threading.local()

    def read():
        buffer = getattr(buffers, "usage", None)
        if buffer is None:
            buffer = buffers.usage = ffi.new("struct rusage *")
        if library.getrusage(resource.RUSAGE_SELF, buffer):
            error = ffi.errno
            raise OSError(error, os.strerror(error))
        return int(buffer.ru_maxrss)

    # Confirm the units/offset against the ordinary reader before admission.
    # The counter is monotone, so intervening allocation cannot invalidate
    # this bracket as long as the native reader observes that same counter.
    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    value = read()
    after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if not before <= value <= after:
        raise ValueError("native RSS reader disagrees with resource")
    return read


class PeakRSSReader:
    """Read the same OS counter every call, with a portable quiet fallback.

    The native reader avoids materializing unused rusage fields and keeps a
    separate buffer per thread. It neither caches RSS nor delays a limit check.
    """

    def __init__(self, *, native=True):
        self.identity = "resource-getrusage"
        self._native = None
        self._scale = 1 if sys.platform == "darwin" else 1024
        if (
            native
            and sys.platform == "darwin"
            and platform.machine() in ("arm64", "x86_64")
            and struct.calcsize("P") == 8
        ):
            try:
                self._native = _darwin_reader()
            except (ImportError, OSError, ValueError):
                pass
            else:
                self.identity = "darwin-getrusage-cffi"

    def __call__(self):
        if self._native is not None:
            return self._native()
        value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return int(value * self._scale)
