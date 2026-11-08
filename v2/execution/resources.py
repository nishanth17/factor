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
    if not isinstance(payload, dict):
        raise UnsupportedResourceVersionError(
            "unsupported-resource-version: checkpoint must be an object"
        )
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


class _SigInfo(ctypes.Structure):
    # macOS arm64 siginfo_t from the host SDK's sys/signal.h.
    _fields_ = [
        (name, ctypes.c_int)
        for name in ("signo", "error", "code", "pid", "uid", "status")
    ] + [
        ("address", ctypes.c_void_p),
        ("value", ctypes.c_void_p),
        ("band", ctypes.c_long),
        ("reserved", ctypes.c_ulong * 7),
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
        self.system.waitid.argtypes = [
            ctypes.c_int,
            ctypes.c_uint,
            ctypes.POINTER(_SigInfo),
            ctypes.c_int,
        ]
        self.system.waitid.restype = ctypes.c_int

    def peek_exit(self, pid, flag):
        """Read final self CPU before the sole reaper consumes the zombie.

        macOS wait4 rolls reaped descendant CPU into its returned rusage.
        libproc V0 retains self CPU separately, including for zombies.
        WNOWAIT leaves reaping and status ownership with the original Popen.
        """
        info = _SigInfo()
        # P_PID=1, WEXITED=4, WNOWAIT=32 are macOS ABI constants.
        if self.system.waitid(1, pid, ctypes.byref(info), 4 | 32 | flag):
            error = ctypes.get_errno()
            raise OSError(error, os.strerror(error))
        if info.pid == 0:
            return None
        usage = _RUsage()
        if self.lib.proc_pid_rusage(pid, 0, ctypes.byref(usage)):
            raise ResourceProviderError("final proc_pid_rusage failed")
        return ProcessSample(
            pid, usage.start, 0, (usage.user + usage.system) / 1_000_000_000
        )

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
                    0 if usage.exit else usage.resident,
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
    if value is None:
        return
    try:
        valid = (
            type(value) in (int, float) and math.isfinite(value) and value >= 0
        )
    except OverflowError:
        valid = False
    if not valid:
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
    try:
        encoded = json.dumps(
            payload, separators=(",", ":"), allow_nan=False
        ).encode()
    except MemoryError as error:
        raise AllocationFailureError(
            "checkpoint encoding", "checkpoint_save"
        ) from error
    except (TypeError, ValueError) as error:
        raise UnsupportedEncodingError(
            "checkpoint encoding: " + str(error)
        ) from error
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

    def __reduce__(self):
        return type(self), (self.operation, self.owner), self.__dict__


class UnsupportedEncodingError(ValueError):
    """A retained arithmetic representation cannot encode this object."""

    reason = "unsupported_encoding"


def record_checkpoint(owner, result):
    """Remember the complete recorded prefix, including encoding overhead."""
    checkpoint = getattr(result, "checkpoint", result)
    if not isinstance(checkpoint, dict):
        return
    payload = checkpoint.get("payload", {})
    header = checkpoint.get("header", {})
    payload = payload if isinstance(payload, dict) else {}
    header = header if isinstance(header, dict) else {}
    records = (
        checkpoint.get("resources"),
        payload.get("resources"),
        payload,
        checkpoint.get("budget"),
        header.get("budget"),
        checkpoint,
    )
    for values in records:
        if not isinstance(values, dict):
            continue
        if {"work_used", "wall_used", "cpu_used"} <= values.keys():
            charges = dict(
                work=values["work_used"],
                wall=values["wall_used"],
                cpu=values["cpu_used"],
            )
        elif {"work", "wall", "cpu"} <= values.keys():
            charges = {key: values[key] for key in ("work", "wall", "cpu")}
        else:
            continue
        overhead = checkpoint.get("serialization_overhead")
        if isinstance(overhead, dict):
            charges["wall"] += overhead["wall"]
            charges["cpu"] += overhead["cpu"]
        elif isinstance(overhead, (list, tuple)):
            charges["wall"] += overhead[0]
            charges["cpu"] += overhead[1]
        owner.record_charges(**charges)
        return


def encode_checkpoint(function):
    """Encode published state within its owner's active accounting scope."""

    @wraps(function)
    def encode(*args, **kwargs):
        previous = active_owner()
        supplied = getattr(args[0], "budget", None) if args else None
        owner = previous or (root_budget(supplied) if supplied else None)
        owned = owner is not None and previous is None
        if owned:
            owner.start()
        token = _ACTIVE_OWNER.set(owner)
        try:
            result = function(*args, **kwargs)
            result["resource_version"] = RESOURCE_VERSION
            if owner is not None:
                record_checkpoint(owner, result)
            return result
        except MemoryError as error:
            raise AllocationFailureError(
                function.__qualname__, "checkpoint_encoding"
            ) from error
        finally:
            _ACTIVE_OWNER.reset(token)
            if owned:
                owner.close()

    return encode


def decode_checkpoint(function):
    """Restore within one owner and reject old resource schemas first."""
    signature = inspect.signature(function)

    @wraps(function)
    def decode(*args, **kwargs):
        bound = signature.bind(*args, **kwargs)
        previous = active_owner()
        supplied = (
            resolve_budget(bound.arguments.get("budget"))
            if "budget" in signature.parameters
            else None
        )
        owner = previous or (root_budget(supplied) if supplied else None)
        owned = owner is not None and previous is None
        if owned:
            owner.start()
        token = _ACTIVE_OWNER.set(owner)
        try:
            checkpoint = bound.arguments["checkpoint"]
            require_resource_version(checkpoint)
            bound.arguments["checkpoint"] = {
                key: value
                for key, value in checkpoint.items()
                if key != "resource_version"
            }
            if supplied is not None:
                bound.arguments["budget"] = supplied
            return function(*bound.args, **bound.kwargs)
        except MemoryError as error:
            raise AllocationFailureError(
                function.__qualname__, "checkpoint_restore"
            ) from error
        finally:
            _ACTIVE_OWNER.reset(token)
            if owned:
                owner.close()

    return decode


_ACTIVE_OWNER = contextvars.ContextVar("factor_resource_owner", default=None)


def active_owner():
    """Return the borrowed root during a public factoring call, if any."""
    return _ACTIVE_OWNER.get()


def root_budget(budget):
    """Find the root while retaining stage/polling views at call sites."""
    while hasattr(budget, "budget"):
        budget = budget.budget
    return budget


def resolve_budget(budget=None):
    """Use an explicit view/counter, otherwise borrow the active owner."""
    if budget is not None:
        return budget
    previous = active_owner()
    if previous is not None:
        return previous
    from .budget import Budget

    return Budget()


def accounted_entry(function):
    """Own validation, restore, computation, encoding and final cleanup."""

    @wraps(function)
    def run(*args, **kwargs):
        previous = _ACTIVE_OWNER.get()
        supplied = kwargs.get("budget")
        if previous is not None and (
            supplied is None or root_budget(supplied) is not previous
        ):
            supplied = previous
        supplied = resolve_budget(supplied)
        owner = previous or root_budget(supplied)
        own_interval = previous is None
        if own_interval:
            owner.start()
        kwargs["budget"] = supplied
        token = _ACTIVE_OWNER.set(owner)
        try:
            result = function(*args, **kwargs)
            record_checkpoint(owner, result)
            return result
        except MemoryError as error:
            raise AllocationFailureError(
                function.__qualname__, "unpublished_scratch"
            ) from error
        finally:
            _ACTIVE_OWNER.reset(token)
            if own_interval:
                owner.close()

    return run


def accounted_job(function):
    """Exclude idle intervals between direct retained-job invocations."""

    @wraps(function)
    def run(self, *args, **kwargs):
        owner = root_budget(self.budget)
        previous = _ACTIVE_OWNER.get()
        own_interval = previous is None
        if own_interval:
            owner.start()
        token = _ACTIVE_OWNER.set(owner)
        try:
            return function(self, *args, **kwargs)
        except MemoryError as error:
            raise AllocationFailureError(
                function.__qualname__, type(self).__name__
            ) from error
        finally:
            _ACTIVE_OWNER.reset(token)
            if own_interval:
                owner.close()

    return run
