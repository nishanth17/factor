"""Owned, lazy automatic ECM selection for ordinary portfolio calls."""

import hashlib
import json
import sys

from ..execution.resources import removed_arguments

AUTO_VERSION = "b12-f4-auto-observed-v3"
POLICY_PATH = None
POLICY_SHA256 = (
    "38b28f423811ceb3242dcc989f8e14f70759c96fa9944cd2fa33bbb846327032"
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


def cold_serial_default(n, config):
    """Decline measured cheap cold plans without loading a cost catalogue.

    This is a numerical-plan guard, not a factor-size prediction. Larger
    explicit quotas/bounds, paired arithmetic and allocation policies retain
    ordinary selection. Warm and resumed owners do not use this shortcut.
    """
    if (
        config.backend != "python-int"
        or config.allocation is not None
        or config.ecm_pair_distance is not None
        or config.ecm_pair_wheel is not None
    ):
        return False
    tiers = tuple(map(tuple, config.ecm_tiers))
    if tiers == ((2000, 50000, 64),):
        return 10**29 <= n < 10**30
    if tiers == ((2000, 147396, 32),):
        return 10**49 <= n < 10**50 or 10**69 <= n < 10**70
    return False


def _runtime_supported():
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        return False
    if sys.platform != "darwin":
        return False
    import platform

    return platform.machine() == "arm64"


@removed_arguments(("parent_memory", "memory_bytes"))
def selected_context(*, seed=0, prefer_warm=False):
    """Return a retained context using the unchanged packaged cost policy."""
    from .selection import SelectedECMContext

    return SelectedECMContext(
        POLICY_PATH, POLICY_SHA256, seed=seed, prefer_warm=prefer_warm
    )


class AutomaticECMOwner:
    """Own a cold-intent pool without runtime imports on easy exits."""

    @removed_arguments(("memory_bytes",))
    def __init__(self, *, seed=0):
        self.requested_seed = seed
        self.context = None
        self.closed = False
        self.last_stats = {}
        self.disabled_reason = None
        self.native_context = None
        self.native_admission = None
        self.native_disabled_reason = None
        self.retained_intent = False

    def configure(
        self,
        config,
        generator,
        saved,
        *,
        automatic_memory=True,
        retained_intent=False,
    ):
        """Bind numerical policy/seed without advancing the parent RNG."""
        if saved is None:
            # Negative optional-build decisions last for one call. A retained
            # owner can discover a newly supplied build on its next call.
            self.native_disabled_reason = None
        self.whole_native = config.ecm_kernel == "auto-full"
        self.retained_intent = retained_intent
        version = "b4a-auto-full-v1" if self.whole_native else AUTO_VERSION
        implementation = (
            "ecm-auto-full-v1" if self.whole_native else "ecm-auto-observed-v3"
        )
        if saved is not None:
            if (
                set(saved) != {"version", "profile", "seed", "implementation"}
                or saved["version"] != version
                or saved["profile"] != POLICY_SHA256
                or saved["implementation"] != implementation
                or type(saved["seed"]) is not int
                or not 0 <= saved["seed"] < 2**256
            ):
                raise ValueError("incompatible automatic ECM checkpoint")
            seed = saved["seed"]
        else:
            seed = self.requested_seed
            if type(seed) is not int or not 0 <= seed < 2**256:
                seed = int(_digest([AUTO_VERSION, generator.getstate()]), 16)
        self.seed = seed
        self.header = dict(
            version=version,
            profile=POLICY_SHA256,
            seed=seed,
            implementation=implementation,
        )
        self.identity = _digest(self.header)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        if not self.closed:
            if self.context is not None:
                self.context.close()
            if self.native_context is not None:
                self.native_context.close()
            self.closed = True

    def _serial(self, reason):
        self.last_stats = dict(
            engaged=False,
            legacy_serial=True,
            engagement=dict(engage=False, state="cold", reason=reason),
        )
        return None

    def run_portfolio(self, n, *, config, budget, checkpoint=None, **kwargs):
        if config.ecm_kernel == "auto-full":
            from .native_auto import run_portfolio

            return run_portfolio(
                self,
                n,
                config=config,
                budget=budget,
                checkpoint=checkpoint,
                **kwargs,
            )
        return self._run_incumbent(
            n, config=config, budget=budget, checkpoint=checkpoint, **kwargs
        )

    def _run_incumbent(self, n, *, config, budget, checkpoint=None, **kwargs):
        if self.closed:
            raise RuntimeError("automatic ECM owner is closed")
        budget.consume(0)
        if checkpoint is None and self.context is None:
            reason = self.disabled_reason
            if reason is None and not _size_supported(n):
                reason = "unmeasured_selection_cost"
            if reason is None and config.backend != "python-int":
                reason = "unsupported_automatic_backend"
            if reason is None and cold_serial_default(n, config):
                reason = "cold_serial_default"
            if reason is None and not _runtime_supported():
                reason = "unsupported_automatic_runtime"
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
