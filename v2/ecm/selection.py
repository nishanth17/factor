"""Lazy opt-in ECM configuration selection from a bounded measured profile."""

import hashlib
import json
import threading
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path

from ..execution.resources import accounted_entry, removed_arguments
from ..execution.worker_contract import (
    SEARCH_VERSION,
    VERSION,
    CostRow,
    canonical,
    digest,
    finite,
    integer,
)

PROFILE_VERSION = "b12-f4-selection-v3"
MAX_DIGIT_GAP = 20
PROFILE_BYTES = 256 * 2**10
MODEL_BYTES = 4 * 2**20
HEADER_BYTES = 16 * 2**10
PROFILE_SCRATCH = 2 * PROFILE_BYTES + 65536
COMPATIBILITY = "ecm-python-policy-v1"


def _checksum(value):
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(char not in "0123456789abcdef" for char in value)
    ):
        raise ValueError("invalid SHA256")
    return value


def _key(digits, options):
    integer(digits, "decimal digits", 2, 1234)
    return digits, _checksum(options)


@dataclass(frozen=True)
class SelectionCosts:
    """Only observable regime keys and measured costs enter decisions."""

    digits: int
    options: str
    seconds_per_curve: tuple
    cold_seconds: float | None
    warm_seconds: float | None
    pool_seconds: float | None

    def __post_init__(self):
        _key(self.digits, self.options)
        checked = CostRow(
            seconds_per_curve=self.seconds_per_curve,
            cold_seconds=self.cold_seconds,
            warm_seconds=self.warm_seconds,
            pool_seconds=self.pool_seconds,
        )
        for name in (
            "seconds_per_curve",
            "cold_seconds",
            "warm_seconds",
            "pool_seconds",
        ):
            object.__setattr__(self, name, getattr(checked, name))


def _path_costs(digits, options, costs):
    """Do not fit through an accepted PRAC/ladder eligibility discontinuity."""
    chains = (
        options["ecm_chain_mode"] != "off"
        and options["ecm_chain_bytes"] > 0
        and options["chunk_size"] == 16
        and any(
            b1 == 2000 and count >= 8 for b1, _, count in options["ecm_tiers"]
        )
    )
    if not chains:
        return costs
    from .chains import MAX_MODULUS, MIN_MODULUS

    lower, upper = len(str(MIN_MODULUS)), len(str(MAX_MODULUS)) - 1
    supported = lower <= digits <= upper
    return tuple(
        row for row in costs if (lower <= row.digits <= upper) == supported
    )


def _bracket(digits, options, costs):
    """Bound interpolation by actual measured sizes of the identical plan."""
    matching = [row for row in costs if row.options == options]
    lower = max(
        (row for row in matching if row.digits < digits),
        key=lambda row: row.digits,
        default=None,
    )
    upper = min(
        (row for row in matching if row.digits > digits),
        key=lambda row: row.digits,
        default=None,
    )
    if lower is None or upper is None:
        return None
    if upper.digits - lower.digits > MAX_DIGIT_GAP:
        return None
    return lower, upper


def _size_cost(digits, options, costs, rate_costs):
    """Interpolate serial rates; crossovers remain bound to the pool bundle."""
    matching = [row for row in costs if row.options == options]
    exact = next((row for row in matching if row.digits == digits), None)
    if exact is not None:
        return exact, dict(size_model="measured", anchor_digits=[digits])
    measured_rate = next(
        (
            row
            for row in rate_costs
            if row.digits == digits and row.options == options
        ),
        None,
    )
    bounds = _bracket(digits, options, rate_costs)
    if not matching or (measured_rate is None and bounds is None):
        return None, {}
    if measured_rate is not None:
        rates = measured_rate.seconds_per_curve
        provenance = dict(
            size_model="measured-rate-v1", anchor_digits=[digits]
        )
    else:
        lower, upper = bounds
        if not lower.seconds_per_curve or len(lower.seconds_per_curve) != len(
            upper.seconds_per_curve
        ):
            return None, {}
        weight = (digits - lower.digits) / (upper.digits - lower.digits)
        rates = tuple(
            a + weight * (b - a)
            for a, b in zip(lower.seconds_per_curve, upper.seconds_per_curve)
        )
        provenance = dict(
            size_model="bracket-linear-v1",
            anchor_digits=[lower.digits, upper.digits],
        )
    # A retained context never borrows another bundle's crossover. Its nearest
    # compatible anchor supplies overhead; this remains an unconfirmed model.
    anchor = min(
        matching, key=lambda row: (abs(row.digits - digits), row.digits)
    )
    cost = replace(anchor, digits=digits, seconds_per_curve=rates)
    return cost, dict(provenance, threshold_anchor=anchor.digits)


@dataclass(frozen=True)
class SelectionRule:
    """Share serial rates across bundles, never pool-specific crossovers."""

    profile: str
    execution: str
    costs: tuple = ()
    rate_costs: tuple = ()
    workers: int = 1

    def decision(self, manifest, **kwargs):
        digits, options = len(str(manifest.n)), manifest.options
        cost, provenance = _size_cost(
            digits,
            profile_options_digest(options, self.costs + self.rate_costs),
            _path_costs(digits, options, self.costs),
            _path_costs(digits, options, self.rate_costs),
        )
        ready = kwargs.get("ready", False)
        started = kwargs.get("pool_started", False)
        state = "warm" if ready else "pool" if started else "cold"
        if cost is None or (state == "pool" and cost.pool_seconds is None):
            return dict(
                engage=False,
                state=state,
                estimate=None,
                threshold=None,
                reason="unmeasured_selection_cost",
                **provenance,
            )
        # Recompute the minimum useful dispatch wave for an interpolated size.
        # Existing exact anchors retain their literal calibrated thresholds.
        if provenance["size_model"] != "measured":
            counts = [tier[2] for tier in manifest.tiers]
            total = sum(counts)
            if total and len(counts) == len(cost.seconds_per_curve):
                wave = (
                    self.workers
                    * sum(
                        count * rate
                        for count, rate in zip(counts, cost.seconds_per_curve)
                    )
                    / total
                )
                thresholds = {}
                for name in ("cold_seconds", "warm_seconds", "pool_seconds"):
                    value = getattr(cost, name)
                    thresholds[name] = (
                        None if value is None else max(value, wave)
                    )
                cost = replace(cost, **thresholds)
        decision = CostRow(
            schedule=manifest.schedule,
            cold_seconds=cost.cold_seconds,
            warm_seconds=cost.warm_seconds,
            pool_seconds=cost.pool_seconds,
            seconds_per_curve=cost.seconds_per_curve,
            execution=self.execution,
        ).decision(manifest, **kwargs)
        return dict(decision, **provenance)


@dataclass(frozen=True)
class ParentOwnership:
    parent_memory: int


def _saved_budget(budget, saved):
    """Restore cumulative outer charges before reserving a new header."""
    base = getattr(budget, "budget", budget)
    if set(saved) != {"work", "cpu", "wall"}:
        raise ValueError("invalid selector budget")
    integer(saved["work"], "saved selector work")
    for key in ("cpu", "wall"):
        finite(saved[key], "saved selector " + key)
    base.restore_charges(
        work=saved["work"], cpu=saved["cpu"], wall=saved["wall"]
    )


class SelectedECMContext:
    """Quiet lazy adapter; a retained pool keeps its first selected settings.

    Supply a version-compatible measured profile and its expected SHA256.
    Direct component calls retain SEARCHv1 serial replay. Portfolio declines
    return to the incumbent arithmetic/RNG path without starting a search.
    Closing is explicit, including when factorization raises; use a context
    manager.
    """

    @removed_arguments(("parent_memory", "memory_bytes"))
    def __init__(
        self,
        profile,
        profile_sha256,
        *,
        seed=0,
        prefer_warm=False,
    ):
        integer(seed, "method seed", 0, 2**256 - 1)
        if type(prefer_warm) is not bool:
            raise ValueError("pool preference must be Boolean")
        self.profile = None if profile is None else Path(profile)
        self.profile_sha256 = _checksum(profile_sha256)
        self.seed, self.prefer_warm = seed, prefer_warm
        self.identity = digest(
            [PROFILE_VERSION, VERSION, self.profile_sha256, seed, prefer_warm]
        )
        self.inner = None
        self.catalogue = None
        self.selected = None
        self.origin = None
        self.mode = None
        self.closed = False
        self.lock = threading.Lock()
        self.last_stats, self.lifetime_stats = {}, {}
        self.selector_work = 0
        self.selector_cpu = 0.0

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    @accounted_entry
    def close(self, *, budget=None):
        with self.lock:
            if self.inner is not None:
                self.inner.close()
            self.catalogue = None
            self.closed = True

    def _load(self, budget):
        from ..execution.worker_contract import WorkerSettings

        if self.catalogue is not None:
            return
        budget.consume(0)
        if self.profile is None:
            from .auto import POLICY_SHA256
            from .worker_policy import policy_data

            if self.profile_sha256 != POLICY_SHA256:
                raise ValueError("bundled policy identity mismatch")

            budget.consume(1)
            self.selector_work += 1
            data = policy_data()
        else:
            size = self.profile.stat().st_size
            budget.consume(size)
            self.selector_work += size
            with self.profile.open("rb") as stream:
                raw = stream.read()
            if (
                len(raw) != size
                or hashlib.sha256(raw).hexdigest() != self.profile_sha256
            ):
                raise ValueError("profile bytes changed")
            data = json.loads(raw)
        if set(data) != {
            "version",
            "status",
            "executor_version",
            "search_version",
            "compatibility",
            "serial_settings",
            "rows",
        } or (
            data["version"] != PROFILE_VERSION
            or data["executor_version"] not in (VERSION, "b12-f4-v7")
            or data["search_version"] != SEARCH_VERSION
            or data["status"]
            not in ("training", "frozen-for-confirmation", "confirmed")
            or data["compatibility"] != COMPATIBILITY
        ):
            raise ValueError("unsupported selector profile")
        rows = data["rows"]
        if not isinstance(rows, list) or len(rows) > 256:
            raise ValueError("selector row capacity exceeded")
        keys, declared = set(), {}
        serial = WorkerSettings(**profile_settings(data["serial_settings"]))
        if serial.workers != 1:
            raise ValueError("serial fallback must own one context")
        for row in rows:
            if set(row) != {"digits", "options", "tiers", "cold", "warm"}:
                raise ValueError("invalid observable selection row")
            tiers = row["tiers"]
            if not isinstance(tiers, list) or not 1 <= len(tiers) <= 32:
                raise ValueError("bounded declared tiers required")
            for tier in tiers:
                if not isinstance(tier, list) or len(tier) != 3:
                    raise ValueError("invalid declared tier")
                b1, b2, count = tier
                integer(b1, "B1", 2, 2**32 - 1)
                integer(b2, "B2", b1, 2**32 - 1)
                integer(count, "curve quota", 0, 4096)
            if sum(tier[2] for tier in tiers) > 4096:
                raise ValueError("selector curve quota exceeds contract")
            key = _key(row["digits"], row["options"])
            if key in keys:
                raise ValueError("duplicate selector regime")
            keys.add(key)
            plan = tuple(map(tuple, tiers))
            if declared.setdefault(row["options"], plan) != plan:
                raise ValueError(
                    "identical option keys declare different tiers"
                )
            for mode in ("cold", "warm"):
                candidate = row[mode]
                if set(candidate) != {"settings", "costs"}:
                    raise ValueError("invalid selection candidate")
                candidate["settings"] = asdict(
                    WorkerSettings(**profile_settings(candidate["settings"]))
                )
                costs = SelectionCosts(**candidate["costs"])
                if (costs.digits, costs.options) != key:
                    raise ValueError("costs belong to another regime")
        self.catalogue = data

    def _select(self, n, config, saved, origin, *, portfolio=False):
        # The selector coexists with the caller's parent workspace. Refuse the
        # combined reservation before constructing bounded worker settings;
        # portfolio callers can then retain their incumbent serial search.
        from ..execution.worker_contract import WorkerSettings
        from .workers import arithmetic_options

        data = self.catalogue
        mode = "warm" if self.prefer_warm else "cold"
        mathematical = arithmetic_options(config, WorkerSettings())
        options = profile_options_digest(mathematical, data["rows"])
        key = len(str(n)), options
        if origin is None:
            origin = next(
                (
                    index
                    for index, entry in enumerate(data["rows"])
                    if (entry["digits"], entry["options"]) == key
                ),
                -1,
            )
            if origin == -1:
                rates = tuple(
                    SelectionCosts(**entry[mode]["costs"])
                    for entry in data["rows"]
                )
                bounds = _bracket(
                    key[0], key[1], _path_costs(key[0], mathematical, rates)
                )
                if bounds is not None:
                    anchor = min(
                        bounds,
                        key=lambda cost: (
                            abs(cost.digits - key[0]),
                            cost.digits,
                        ),
                    )
                    origin = next(
                        i
                        for i, entry in enumerate(data["rows"])
                        if (entry["digits"], entry["options"])
                        == (anchor.digits, anchor.options)
                    )
        integer(origin, "selected origin", -1, len(data["rows"]) - 1)
        row = None if origin == -1 else data["rows"][origin]
        if portfolio and row is None:
            return None
        candidate = (
            data["serial_settings"] if row is None else (row[mode]["settings"])
        )
        tiers = config.ecm_tiers if row is None else row["tiers"]
        if (
            saved is None
            and row is not None
            and tuple(map(tuple, tiers))
            != (tuple(map(tuple, config.ecm_tiers)))
        ):
            raise ValueError("profile declared tiers differ from mathematics")
        sum(tier[2] for tier in tiers)
        settings = WorkerSettings(**profile_settings(candidate))
        if saved is not None and asdict(settings) != saved:
            raise ValueError("saved selector settings differ from its policy")
        if self.inner is not None:
            if settings != self.inner.settings or origin != self.origin:
                raise ValueError("retained selector cannot change workers")
            return self.inner.settings, self.inner.rule, self.origin
        normalized = asdict(WorkerSettings(**profile_settings(candidate)))
        normalized.pop("parent_memory", None)
        normalized.pop("checkpoint_bytes", None)
        costs = []
        for entry in data["rows"]:
            proposed = profile_settings(entry[mode]["settings"])
            proposed.pop("parent_memory", None)
            proposed.pop("checkpoint_bytes", None)
            if proposed == normalized:
                costs.append(SelectionCosts(**entry[mode]["costs"]))
        rule = SelectionRule(
            self.profile_sha256,
            digest([VERSION, "cost-policy", asdict(settings)]),
            tuple(costs),
            tuple(
                SelectionCosts(**entry[mode]["costs"])
                for entry in data["rows"]
            ),
            settings.workers,
        )
        return settings, rule, origin

    def _portfolio_decision(self, n, config, budget, lineage, settings, rule):
        from ..execution.worker_contract import Manifest
        from .workers import arithmetic_options

        options = arithmetic_options(config, settings)
        if not any(
            row.options == profile_options_digest(options, rule.costs)
            for row in rule.costs
        ):
            return dict(
                engage=False,
                state="pool" if self.inner is not None else "cold",
                estimate=None,
                threshold=None,
                reason="unmeasured_selection_cost",
                funded=False,
            )
        manifest = Manifest(n, self.seed, lineage, options)
        started = self.inner is not None and self.inner.pool_started()
        ready = self.inner is not None and self.inner.ready(manifest)
        decision = rule.decision(
            manifest,
            ready=ready,
            pool_started=started,
            execution=rule.execution,
        )
        budget.consume(0)
        funded = True
        return dict(decision, funded=funded)

    def _stats(self, started, stats):
        self.selector_cpu += max(
            0,
            time.process_time() - started - stats.get("parent_cpu", 0),
        )
        self.last_stats = dict(
            stats,
            selector_profile=self.profile_sha256,
            selector_work=self.selector_work,
            selector_cpu=self.selector_cpu,
        )
        self.lifetime_stats = (
            {} if self.inner is None else dict(self.inner.lifetime_stats)
        )
        self.lifetime_stats["charged_work"] = (
            self.lifetime_stats.get("charged_work", 0) + self.selector_work
        )
        return self.last_stats

    @accounted_entry
    def run_portfolio(self, n, **kwargs):
        """Return None to retain incumbent curves and fallback RNG draws."""
        return self._invoke(n, portfolio=True, **kwargs)

    @accounted_entry
    def run(
        self,
        n,
        *,
        config,
        budget,
        lineage="root",
        checkpoint=None,
        fixed_work=False,
    ):
        return self._invoke(
            n,
            config=config,
            budget=budget,
            lineage=lineage,
            checkpoint=checkpoint,
            fixed_work=fixed_work,
            portfolio=False,
        )

    def _invoke(
        self,
        n,
        *,
        config,
        budget,
        lineage="root",
        checkpoint=None,
        fixed_work=False,
        portfolio,
    ):
        if not self.lock.acquire(blocking=False):
            raise RuntimeError("selector context already active")
        try:
            return self._run(
                n, config, budget, lineage, checkpoint, fixed_work, portfolio
            )
        finally:
            self.lock.release()

    def _run(
        self, n, config, budget, lineage, checkpoint, fixed_work, portfolio
    ):
        if self.closed:
            raise RuntimeError("selector context is closed")
        actual_identity = digest(
            [
                PROFILE_VERSION,
                VERSION,
                self.profile_sha256,
                self.seed,
                self.prefer_warm,
            ]
        )
        if actual_identity != self.identity:
            raise ValueError("selector configuration mutated")
        mode = "portfolio" if portfolio else "component"
        if self.mode is not None and self.mode != mode:
            raise ValueError(
                "cannot mix component and portfolio context modes"
            )
        self.mode = mode
        integer(n, "cofactor", 2, 2**4096 - 1)
        started = time.process_time()
        budget.consume(0)
        saved, inner_checkpoint, origin = None, None, None
        if checkpoint is not None:
            if (
                not isinstance(checkpoint, dict)
                or len(checkpoint) != 5
                or set(checkpoint)
                != {
                    "header",
                    "sha256",
                    "inner",
                    "serialization_overhead",
                    "overhead_sha256",
                }
            ):
                raise ValueError("invalid selector checkpoint fields")
            header = checkpoint["header"]
            overhead = checkpoint["serialization_overhead"]
            if (
                not isinstance(header, dict)
                or len(header) != 7
                or set(header)
                != {
                    "version",
                    "mode",
                    "identity",
                    "origin",
                    "settings",
                    "budget",
                    "inner",
                }
                or header["mode"] not in ("portfolio", "component")
                or not isinstance(header["settings"], dict)
                or len(header["settings"]) > 64
                or not isinstance(header["budget"], dict)
                or set(header["budget"]) != {"work", "cpu", "wall"}
                or any(
                    not isinstance(key, str) or len(key) > 64
                    for key in header["settings"]
                )
                or not isinstance(overhead, dict)
                or len(overhead) != 2
                or set(overhead) != {"cpu", "wall"}
                or any(
                    type(value) not in (int, float, bool, str)
                    or isinstance(value, str)
                    and len(value) > 64
                    or isinstance(value, int)
                    and value.bit_length() > 4096
                    for value in header["settings"].values()
                )
            ):
                raise ValueError("unbounded selector header")
            integer(header["origin"], "selected origin", -1, 255)
            for name, value in overhead.items():
                finite(value, "selector serialization " + name)
            _checksum(header["identity"])
            _checksum(header["inner"])
            if header["version"] != PROFILE_VERSION:
                raise ValueError("unsupported selector checkpoint")
            for name, value in header["budget"].items():
                finite(value, "selector budget " + name)
            if (
                len(canonical(header)) > HEADER_BYTES
                or checkpoint["sha256"] != digest(header)
                or checkpoint["overhead_sha256"] != digest(overhead)
                or set(overhead) != {"cpu", "wall"}
                or header["version"] != PROFILE_VERSION
                or header["identity"] != self.identity
                or header["mode"] != mode
                or set(header)
                != {
                    "version",
                    "mode",
                    "identity",
                    "origin",
                    "settings",
                    "budget",
                    "inner",
                }
            ):
                raise ValueError("corrupt or changed selector checkpoint")
            saved = header["settings"]
            origin = integer(header["origin"], "selected origin", -1, 255)
            inner_checkpoint = checkpoint["inner"]
            if inner_checkpoint["sha256"] != header["inner"]:
                raise ValueError("selector inner checkpoint changed")
            charges = dict(header["budget"])
            for key in ("cpu", "wall"):
                finite(overhead[key], "selector serialization " + key)
                charges[key] += overhead[key]
            _saved_budget(budget, charges)
        if (
            portfolio
            and checkpoint is None
            and self.inner is None
            and not self.prefer_warm
        ):
            from .auto import POLICY_PATH, POLICY_SHA256, cold_serial_default

            # Only the packaged cold policy has this predeclared serial rule.
            # Custom catalogues and checked resumes must still validate their
            # own compatibility/policy identity before selecting any search.
            if (
                self.profile is POLICY_PATH
                and self.profile_sha256 == POLICY_SHA256
                and cold_serial_default(n, config)
            ):
                self._stats(
                    started,
                    dict(
                        engaged=False,
                        legacy_serial=True,
                        cofactor_digits=len(str(n)),
                        engagement=dict(
                            engage=False,
                            state="cold",
                            reason="cold_serial_default",
                        ),
                    ),
                )
                return None
        self._load(budget)
        budget.consume(HEADER_BYTES)
        self.selector_work += HEADER_BYTES
        if self.inner is None or saved is not None:
            try:
                plan = self._select(
                    n,
                    config,
                    saved,
                    origin,
                    portfolio=portfolio,
                )
                refusal = "unmeasured_selection_cost"
            except MemoryError:
                if not portfolio or checkpoint is not None:
                    raise
                plan, refusal = None, "worker_ownership_refused"
            if plan is None:
                self._stats(
                    started,
                    dict(
                        engaged=False,
                        legacy_serial=True,
                        cofactor_digits=len(str(n)),
                        engagement=dict(
                            engage=False,
                            state="cold",
                            estimate=None,
                            threshold=None,
                            reason=refusal,
                            funded=False,
                        ),
                    ),
                )
                return None
            settings, rule, origin = plan
        else:
            settings, rule, origin = (
                self.inner.settings,
                self.inner.rule,
                self.origin,
            )
        if portfolio and checkpoint is None:
            decision = self._portfolio_decision(
                n, config, budget, lineage, settings, rule
            )
            if not decision["engage"] or not decision["funded"]:
                self._stats(
                    started,
                    dict(
                        engaged=False,
                        legacy_serial=True,
                        engagement=decision,
                        cofactor_digits=len(str(n)),
                        configured_workers=settings.workers,
                    ),
                )
                return None
        if self.inner is None:
            from .workers import ECMExecutor

            # Cost selection already ran. A later resource/cancellation refusal
            # must pause this search, never silently execute SEARCHv1 serial.
            self.inner = ECMExecutor(
                seed=self.seed,
                settings=settings,
                rule=rule,
                force_workers=portfolio,
            )
            self.selected, self.origin = asdict(settings), origin
        result = self.inner.run(
            n,
            config=config,
            budget=budget,
            lineage=lineage,
            checkpoint=inner_checkpoint,
            fixed_work=fixed_work,
        )
        base = getattr(budget, "budget", budget)
        header = dict(
            version=PROFILE_VERSION,
            mode=mode,
            identity=self.identity,
            settings=self.selected,
            origin=self.origin,
            budget=dict(
                work=base.used, cpu=base.cpu_used, wall=base.wall_used
            ),
            inner=result.checkpoint["sha256"],
        )
        wall, cpu = time.monotonic(), time.process_time()
        if len(canonical(header)) > HEADER_BYTES:
            raise ValueError("selector checkpoint header exceeds ownership")
        checkpoint = dict(
            header=header, sha256=digest(header), inner=result.checkpoint
        )
        overhead = dict(
            wall=time.monotonic() - wall, cpu=time.process_time() - cpu
        )
        checkpoint.update(
            serialization_overhead=overhead, overhead_sha256=digest(overhead)
        )
        stats = self._stats(started, result.stats)
        return replace(result, checkpoint=checkpoint, stats=stats)


REMOVED_PROFILE_FIELDS = frozenset(
    (
        "assignment_work",
        "assignment_cpu",
        "startup_work",
        "worker_memory",
        "parent_memory",
        "memory_bytes",
        "frame_bytes",
        "checkpoint_bytes",
        "max_runs",
        "residue_bytes",
    )
)


def profile_settings(values):
    """Adapt the pinned historical cost model, never a user configuration."""
    return {
        key: value
        for key, value in values.items()
        if key not in REMOVED_PROFILE_FIELDS
    }


def profile_options_digest(options, rows):
    """Match historical recipe identities without imposing a byte grant."""
    current = digest(options)
    identities = {
        row["options"] if isinstance(row, dict) else row.options
        for row in rows
    }
    if current in identities:
        return current
    historical = digest(dict(options, memory_bytes=64 * 2**20))
    return historical if historical in identities else current
