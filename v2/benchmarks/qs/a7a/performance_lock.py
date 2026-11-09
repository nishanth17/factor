"""Open the shared lock, retaining a supervising owner's descriptor."""

import os
from pathlib import Path

LOCK_PATH = Path("/private/tmp/factor-performance.lock")


def open_lock():
    """Return a private handle; callers still acquire flock normally."""
    inherited = os.environ.get("FACTOR_A7A_PERFORMANCE_LOCK_FD")
    if inherited is None:
        return LOCK_PATH.open("r+")

    stream = os.fdopen(os.dup(int(inherited)), "r+")
    actual, expected = os.fstat(stream.fileno()), LOCK_PATH.stat()
    if (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino):
        stream.close()
        raise ValueError("inherited performance lock is a different file")
    return stream
