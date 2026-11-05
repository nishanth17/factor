"""Owned, lazy automatic ECM selection for ordinary portfolio calls."""

import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

AUTO_VERSION = "b12-f4-auto-v1"
POLICY_PATH = Path(__file__).with_name("ecm_worker_policy.json")
POLICY_SHA256 = (
    "1e01666ddbe89b7250c3ccdce9e5f0444b26794bed9d4226cb4c236b0dd4f637"
)
MEASURED_PROFILE_SHA256 = (
    "16d52a6e55270d3cfea467f2a5ca41dd00ad1d8d052f2263861656a99b4b00cf"
)
DEFAULT_MEMORY_BYTES = 8 * 2**30
DEFAULT_WORK = 10**13
DEFAULT_CPU_SECONDS = 4096
# The smallest measured pool has four workers plus one spawn helper. Each
# needs ten startup and four shutdown CPU seconds. Metadata and transport
# additionally need finite work before any useful curve can be dispatched.
MIN_STARTUP_CPU = 70
MIN_STARTUP_WORK = 8_000_000
SELECTOR_OWNERSHIP = 4 * 2**20 + 2 * 256 * 2**10 + 65536 + 3 * 16 * 2**10


def _digest(value):
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _size_supported(n):
    # Every packaged anchor is in this interval, and the cost model never
    # extrapolates. Rejecting outside it needs neither imports nor JSON.
    return 10**29 <= n < 10**100


def _runtime_supported():
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        return False
    if sys.platform != "darwin":
        return False
    import platform

    return platform.machine() == "arm64"


def selected_context(
    *,
    seed=0,
    prefer_warm=False,
    parent_memory=1040 * 2**20,
    memory_bytes=DEFAULT_MEMORY_BYTES,
):
    """Return a caller-owned context using the packaged production policy."""
    from .ecm_selection import SelectedECMContext

    return SelectedECMContext(
        POLICY_PATH,
        POLICY_SHA256,
        seed=seed,
        prefer_warm=prefer_warm,
        parent_memory=parent_memory,
        memory_bytes=memory_bytes,
    )


class AutomaticECMOwner:
    """Own a cold-intent pool without runtime imports on easy exits."""

    def __init__(self, *, seed=0, memory_bytes=None):
        if memory_bytes is not None and (
            type(memory_bytes) is not int or not 0 <= memory_bytes < 2**63
        ):
            raise ValueError("ECM aggregate memory must be a finite allowance")
        self.requested_seed = seed
        self.requested_memory = memory_bytes
        self.context = None
        self.closed = False
        self.last_stats = {}
        self.disabled_reason = None

    def configure(self, config, generator, saved, *, automatic_memory):
        """Bind a policy/seed/grant without advancing the portfolio RNG."""
        implementation = hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest()
        if saved is not None:
            if (
                type(saved) is not dict
                or set(saved)
                != {
                    "version",
                    "profile",
                    "seed",
                    "memory_bytes",
                    "implementation",
                }
                or saved["implementation"] != implementation
                or saved["version"] != AUTO_VERSION
                or saved["profile"] != POLICY_SHA256
                or type(saved["seed"]) is not int
                or not 0 <= saved["seed"] < 2**256
                or type(saved["memory_bytes"]) is not int
                or not 0 <= saved["memory_bytes"] <= DEFAULT_MEMORY_BYTES
            ):
                raise ValueError("incompatible automatic ECM checkpoint")
            memory = saved["memory_bytes"]
            if (
                self.requested_memory is not None
                and min(self.requested_memory, DEFAULT_MEMORY_BYTES) != memory
            ):
                raise ValueError("automatic ECM resume cannot change memory")
            seed = saved["seed"]
        else:
            # Retain the measured integer-seed mapping. Other legacy Random
            # seed types use getstate without advancing the parent generator.
            seed = self.requested_seed
            if type(seed) is not int or not 0 <= seed < 2**256:
                seed = int(_digest([AUTO_VERSION, generator.getstate()]), 16)
            memory = self.requested_memory
            if memory is None:
                memory = (
                    DEFAULT_MEMORY_BYTES
                    if automatic_memory
                    else config.memory_bytes
                )
            memory = min(memory, DEFAULT_MEMORY_BYTES)
        self.memory_bytes = memory
        self.seed = seed
        self.parent_memory = config.memory_bytes
        self.settings = SimpleNamespace(
            parent_memory=config.memory_bytes + SELECTOR_OWNERSHIP
        )
        self.header = dict(
            version=AUTO_VERSION,
            profile=POLICY_SHA256,
            seed=seed,
            memory_bytes=memory,
            implementation=implementation,
        )
        self.identity = _digest([self.header, self.settings.parent_memory])

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        if not self.closed:
            if self.context is not None:
                self.context.close()
            self.closed = True

    def _serial(self, reason):
        self.last_stats = dict(
            engaged=False,
            legacy_serial=True,
            engagement=dict(engage=False, state="cold", reason=reason),
        )
        return None

    def run_portfolio(self, n, *, config, budget, checkpoint=None, **kwargs):
        if self.closed:
            raise RuntimeError("automatic ECM owner is closed")
        budget.consume(0)
        base = getattr(budget, "budget", budget)
        if checkpoint is None and self.context is None:
            reason = self.disabled_reason
            if reason is None and not _size_supported(n):
                reason = "unmeasured_selection_cost"
            if reason is None and config.backend != "python-int":
                reason = "unsupported_automatic_backend"
            if reason is None and not _runtime_supported():
                reason = "unsupported_automatic_runtime"
            if reason is None and (
                base.seconds is None
                or base.cpu_seconds is None
                or base.cpu_seconds - base.cpu_used <= MIN_STARTUP_CPU
                or base.work_limit - base.used <= MIN_STARTUP_WORK
            ):
                reason = "insufficient_automatic_startup_allowance"
            if reason is None and (
                self.settings.parent_memory > 4 * 2**30
                or self.memory_bytes
                < self.settings.parent_memory + 5 * 64 * 2**20
            ):
                reason = "worker_ownership_refused"
            if reason is not None:
                return self._serial(reason)
        if self.context is None:
            import multiprocessing

            if multiprocessing.parent_process() is not None:
                if checkpoint is not None:
                    raise ValueError("automatic ECM cannot resume in a worker")
                return self._serial("nested_process")
            self.context = selected_context(
                seed=self.seed,
                parent_memory=self.parent_memory,
                memory_bytes=self.memory_bytes,
            )
            try:
                self.context._load(budget)
            except (OSError, ValueError):
                if checkpoint is not None:
                    raise
                # A missing/stale built-in estimate cannot start a different
                # search. Preserve incumbent serial behavior for this call.
                self.context.close()
                self.context = None
                self.disabled_reason = "unavailable_automatic_policy"
                return self._serial(self.disabled_reason)
        result = self.context.run_portfolio(
            n, config=config, budget=budget, checkpoint=checkpoint, **kwargs
        )
        self.last_stats = self.context.last_stats
        return result

    run = run_portfolio
