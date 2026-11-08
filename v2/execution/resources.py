"""Top-level observed resources and atomic, optional checkpoint persistence."""

import contextvars
import ctypes
import errno
import inspect
import json
import math
import os
import sys
import tempfile
from dataclasses import dataclass
from functools import wraps

RESOURCE_VERSION = "observed-resources/saved-history-v2"
DEFAULT_RSS = object()


class RemovedResourceLimitError(TypeError):
    """A former internal refusal argument was supplied."""

    def __init__(
        self, argument, replacement="Budget(rss_bytes=...) at the root"
    ):
        self.argument, self.replacement = argument, replacement
        super().__init__(f"{argument} was removed; use {replacement}")


class UnsupportedResourceVersionError(ValueError):
    """The checkpoint predates the current accounting contract."""


def require_resource_version(payload):
    if payload.get("resource_version") != RESOURCE_VERSION:
        raise UnsupportedResourceVersionError(
            "unsupported-resource-version: "
            + repr(payload.get("resource_version"))
        )


class ResourceProviderError(RuntimeError):
    """An OS resource observation failed; the guard remains enabled."""


@dataclass(frozen=True)
class ProcessSample:
    pid: int
    birth: int
    rss: int
    cpu: float


class _RUsage(ctypes.Structure):
    _fields_ = [("uuid", ctypes.c_ubyte * 16)] + [
        (name, ctypes.c_uint64)
        for name in (
            "user",
            "system",
            "idle",
            "interrupt",
            "pageins",
            "wired",
            "resident",
            "physical",
            "start",
            "exit",
        )
    ]


class MacResources:
    """Current libproc RSS/CPU for the parent and its live descendants."""

    def __init__(self):
        self.lib = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
        self.system = ctypes.CDLL(None, use_errno=True)
        self.lib.proc_listchildpids.argtypes = [
            ctypes.c_int,
            ctypes.c_void_p,
            ctypes.c_int,
        ]
        self.lib.proc_listchildpids.restype = ctypes.c_int
        self.lib.proc_pid_rusage.argtypes = [
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_void_p,
        ]
        self.lib.proc_pid_rusage.restype = ctypes.c_int
        self.system.sysctlbyname.argtypes = [
            ctypes.c_char_p,
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_size_t),
            ctypes.c_void_p,
            ctypes.c_size_t,
        ]
        self.system.sysctlbyname.restype = ctypes.c_int

    def physical_memory(self):
        value = ctypes.c_uint64()
        size = ctypes.c_size_t(ctypes.sizeof(value))
        if self.system.sysctlbyname(
            b"hw.memsize", ctypes.byref(value), ctypes.byref(size), None, 0
        ):
            raise ResourceProviderError("sysctlbyname(hw.memsize) failed")
        return value.value

    def sample(self, roots):
        pending, seen, result = list(roots), set(), []
        while pending:
            pid = pending.pop()
            if pid in seen:
                continue
            seen.add(pid)
            usage = _RUsage()
            if self.lib.proc_pid_rusage(pid, 0, ctypes.byref(usage)):
                error = ctypes.get_errno()
                if error == errno.ESRCH:
                    continue
                raise ResourceProviderError(f"proc_pid_rusage({pid}): {error}")
            if isinstance(roots, dict):
                expected_birth = roots.get(pid)
                if (
                    expected_birth is not None
                    and expected_birth != usage.start
                ):
                    continue
            result.append(
                ProcessSample(
                    pid,
                    usage.start,
                    usage.resident,
                    (usage.user + usage.system) / 1_000_000_000,
                )
            )
            # Size can change between calls. Retry with the returned size;
            # this grows an observation buffer, never a production quota.
            size = self.lib.proc_listchildpids(pid, None, 0)
            if size < 0:
                raise ResourceProviderError(
                    f"proc_listchildpids({pid}) failed"
                )
            while size:
                children = (ctypes.c_int * (size + 1))()
                count = self.lib.proc_listchildpids(
                    pid, children, ctypes.sizeof(children)
                )
                if count < 0:
                    raise ResourceProviderError(
                        f"proc_listchildpids({pid}) failed"
                    )
                if count < len(children):
                    pending.extend(children[:count])
                    break
                size = count + 1
        return result


def default_provider():
    if sys.platform != "darwin":
        raise ResourceProviderError(
            "RSS provider requires macOS; supply a resource provider "
            "or explicitly disable rss_bytes"
        )
    return MacResources()


def validate_seconds(value, name):
    if value is not None and (
        type(value) not in (int, float)
        or not math.isfinite(value)
        or value < 0
    ):
        raise ValueError(name + " must be finite and nonnegative")


def atomic_save(path, payload, *, checkpoint_bytes=None, events=None):
    """Save published state atomically; oversize leaves the old file valid.

    The caller supplies a snapshot built from published state. Encoding and
    failed saves belong to its active accounting interval. A skipped save
    preserves only the old checkpoint's saved charges (saved-history model).
    """
    if checkpoint_bytes is not None and (
        type(checkpoint_bytes) is not int or checkpoint_bytes < 0
    ):
        raise ValueError("checkpoint_bytes must be nonnegative or None")
    encoded = json.dumps(
        payload, separators=(",", ":"), allow_nan=False
    ).encode()
    if checkpoint_bytes is not None and len(encoded) > checkpoint_bytes:
        event = dict(
            reason="checkpoint_size_skip",
            limit="checkpoint_bytes",
            requested=len(encoded),
            limit_value=checkpoint_bytes,
            scope="checkpoint_save",
        )
        if events is not None:
            events.append(event)
        return event
    path = os.path.abspath(path)
    directory = os.path.dirname(path)
    descriptor, temporary = tempfile.mkstemp(prefix=".factor-", dir=directory)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(encoded)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        directory_fd = os.open(directory, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return None


def removed_arguments(names, replacement="Budget(rss_bytes=...) at the root"):
    """Reject retired keyword arguments before calling the revised API."""

    def decorate(target):
        original = target.__init__ if isinstance(target, type) else target

        parameters = tuple(inspect.signature(original).parameters.values())
        positional_count = sum(
            parameter.kind
            in (
                inspect.Parameter.POSITIONAL_ONLY,
                inspect.Parameter.POSITIONAL_OR_KEYWORD,
            )
            for parameter in parameters
        )
        variadic = any(
            parameter.kind == inspect.Parameter.VAR_POSITIONAL
            for parameter in parameters
        )

        @wraps(original)
        def checked(*args, **kwargs):
            if not variadic and len(args) > positional_count:
                raise RemovedResourceLimitError(names[0], replacement)
            for name in names:
                if name in kwargs:
                    raise RemovedResourceLimitError(name, replacement)
            return original(*args, **kwargs)

        if isinstance(target, type):
            target.__init__ = checked
            return target
        return checked

    return decorate


class AllocationFailureError(RuntimeError):
    """Allocation interrupted an operation; a poisoned owner cannot resume."""

    reason = "allocation_failure"

    def __init__(self, operation, owner=None):
        self.operation, self.owner = operation, owner
        self.payload = dict(
            reason=self.reason, operation=operation, owner=owner
        )
        super().__init__("allocation_failure: " + operation)


class UnsupportedEncodingError(ValueError):
    """A retained arithmetic representation cannot encode this object."""


def encode_checkpoint(function):
    """Bind every new snapshot to the saved-history resource schema."""

    @wraps(function)
    def encode(*args, **kwargs):
        result = function(*args, **kwargs)
        result["resource_version"] = RESOURCE_VERSION
        return result

    return encode


def decode_checkpoint(function):
    """Reject pre-resource-model checkpoints before any replay or decoding."""
    signature = inspect.signature(function)

    @wraps(function)
    def decode(*args, **kwargs):
        bound = signature.bind(*args, **kwargs)
        checkpoint = bound.arguments["checkpoint"]
        if not isinstance(checkpoint, dict):
            raise ValueError("checkpoint must be an object")
        require_resource_version(checkpoint)
        bound.arguments["checkpoint"] = {
            key: value
            for key, value in checkpoint.items()
            if key != "resource_version"
        }
        return function(*bound.args, **bound.kwargs)

    return decode


_ACTIVE_OWNER = contextvars.ContextVar("factor_resource_owner", default=None)


def active_owner():
    """Return the borrowed root during a public factoring call, if any."""
    return _ACTIVE_OWNER.get()


def accounted_entry(function):
    """Own validation, restore, computation, encoding and final cleanup."""

    @wraps(function)
    def run(*args, **kwargs):
        from .budget import Budget

        previous = _ACTIVE_OWNER.get()
        owner = previous or kwargs.get("budget") or Budget()
        own_interval = previous is None
        if own_interval:
            owner.start()
        kwargs["budget"] = owner
        token = _ACTIVE_OWNER.set(owner)
        try:
            return function(*args, **kwargs)
        finally:
            _ACTIVE_OWNER.reset(token)
            if own_interval:
                owner.close()

    return run


def accounted_job(function):
    """Exclude idle intervals between direct retained-job invocations."""

    @wraps(function)
    def run(self, *args, **kwargs):
        owner = getattr(self.budget, "budget", self.budget)
        previous = _ACTIVE_OWNER.get()
        own_interval = previous is None
        if own_interval:
            owner.start()
        token = _ACTIVE_OWNER.set(owner)
        try:
            return function(self, *args, **kwargs)
        finally:
            _ACTIVE_OWNER.reset(token)
            if own_interval:
                owner.close()

    return run
