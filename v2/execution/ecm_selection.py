"""Lazy opt-in ECM configuration selection from a bounded measured profile."""

import hashlib
import json
import threading
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path

from .worker_contract import (
    SEARCH_VERSION,
    VERSION,
    CostRow,
    canonical,
    digest,
    finite,
    integer,
)

PROFILE_VERSION = "b12-f4-selection-v1"
PROFILE_BYTES = 256 * 2**10
MODEL_BYTES = 4 * 2**20
HEADER_BYTES = 16 * 2**10
SOURCE_BYTES = 2 * 2**20
PROFILE_SCRATCH = 2 * PROFILE_BYTES + 65536
SOURCE_ROOT = Path(__file__).resolve().parents[1]
# Whole-call costs include shared imports, including the QS fallback.
SOURCE_NAMES = frozenset(
    (
        "__init__.py",
        "common/__init__.py",
        "common/arithmetic.py",
        "common/preprocessing.py",
        "common/prime_sieve.py",
        "common/utils.py",
        "constants.py",
        "ecm/__init__.py",
        "ecm/chain_options.py",
        "ecm/chain_records.py",
        "ecm/chains.py",
        "ecm/core.py",
        "ecm/paired.py",
        "ecm/prac.py",
        "ecm/programs.py",
        "ecm/wheel.py",
        "execution/__init__.py",
        "execution/allocation.py",
        "execution/budget.py",
        "execution/ecm_defaults.py",
        "execution/ecm_presets.py",
        "execution/ecm_selection.py",
        "execution/ecm_stage_two.py",
        "execution/ecm_worker.py",
        "execution/ecm_worker_state.py",
        "execution/ecm_workers.py",
        "execution/portfolio_config.py",
        "execution/schedules.py",
        "execution/stage_jobs.py",
        "execution/work_budget.py",
        "execution/worker_contract.py",
        "factor.py",
        "pm1/__init__.py",
        "pm1/bounded.py",
        "pm1/core.py",
        "pm1/gaps.py",
        "pm1/tuning.py",
        "portfolio.py",
        "qs/__init__.py",
        "qs/assignment_stream.py",
        "qs/capacity.py",
        "qs/checkpoint.py",
        "qs/cofactor_policy.py",
        "qs/cold_factors.py",
        "qs/collector_kernels.py",
        "qs/coordinate_family.py",
        "qs/count_policy.py",
        "qs/defaults.py",
        "qs/equal_parity.py",
        "qs/extended_checkpoint.py",
        "qs/extended_collector.py",
        "qs/extended_policy.py",
        "qs/external_square.py",
        "qs/extraction.py",
        "qs/factor_base.py",
        "qs/factor_witness.py",
        "qs/families.py",
        "qs/hit_words.py",
        "qs/incremental_structure.py",
        "qs/indexed_relations.py",
        "qs/large_prime_checkpoint.py",
        "qs/large_prime_collector.py",
        "qs/large_primes.py",
        "qs/lazy_relations.py",
        "qs/linear_algebra.py",
        "qs/multiplier.py",
        "qs/norm_multiplier.py",
        "qs/opaque_conversion.py",
        "qs/opaque_staging.py",
        "qs/packed_relations.py",
        "qs/parallel.py",
        "qs/pipeline.py",
        "qs/polynomial.py",
        "qs/power_reuse.py",
        "qs/power_sieve.py",
        "qs/product_window.py",
        "qs/recovery_batch.py",
        "qs/reference_collector.py",
        "qs/relations.py",
        "qs/residual_kernel.py",
        "qs/residual_relations.py",
        "qs/root_updates.py",
        "qs/score.py",
        "qs/short_cycles.py",
        "qs/sieve_collector.py",
        "qs/sign_schedule.py",
        "qs/siqs.py",
        "qs/smooth_batch.py",
        "qs/sparse_queue.py",
        "qs/sss.py",
        "qs/sss_checkpoint.py",
        "qs/storage.py",
        "qs/structured_matrix.py",
        "qs/yield_collector.py",
        "qs/yield_feedback.py",
        "qs/yield_memory.py",
        "qs/yield_support.py",
        "rho/__init__.py",
        "rho/brent.py",
    )
)


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


@dataclass(frozen=True)
class SelectionRule:
    """A fixed context may use only rows calibrated for its own settings."""

    profile: str
    execution: str
    costs: tuple = ()

    def decision(self, manifest, **kwargs):
        key = len(str(manifest.n)), digest(manifest.options)
        cost = next(
            (row for row in self.costs if (row.digits, row.options) == key),
            None,
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
            )
        return CostRow(
            schedule=manifest.schedule,
            cold_seconds=cost.cold_seconds,
            warm_seconds=cost.warm_seconds,
            pool_seconds=cost.pool_seconds,
            seconds_per_curve=cost.seconds_per_curve,
            execution=self.execution,
        ).decision(manifest, **kwargs)


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
    if (
        saved["work"] > base.work_limit
        or base.cpu_seconds is None
        or saved["cpu"] > base.cpu_seconds
        or base.seconds is None
        or saved["wall"] > base.seconds
    ):
        raise ValueError("resume grant below selector charges")
    if base.used == base.prior_cpu == base.prior_wall == 0:
        base.used = saved["work"]
        base.prior_cpu = saved["cpu"]
        base.prior_wall = saved["wall"]
    elif (
        base.used < saved["work"]
        or base.cpu_used < saved["cpu"]
        or base.wall_used < saved["wall"]
    ):
        raise ValueError("resume would reset selector charges")


class SelectedECMContext:
    """Quiet lazy adapter; a retained pool keeps its first selected settings.

    Supply a source-compatible measured profile and its expected SHA256.
    Serial choices use the explicit v1 ECM stream, not legacy portfolio ECM.
    Unknown mathematical regimes remain serial. Closing is explicit, including
    when factorization raises; use this context as a context manager.
    """

    def __init__(
        self,
        profile,
        profile_sha256,
        *,
        seed=0,
        prefer_warm=False,
        parent_memory=1040 * 2**20,
    ):
        integer(seed, "method seed", 0, 2**256 - 1)
        integer(parent_memory, "parent memory", 2**20, 4 * 2**30)
        if type(prefer_warm) is not bool:
            raise ValueError("pool preference must be Boolean")
        self.profile = Path(profile)
        self.profile_sha256 = _checksum(profile_sha256)
        self.seed, self.prefer_warm = seed, prefer_warm
        self.settings = ParentOwnership(
            parent_memory + MODEL_BYTES + PROFILE_SCRATCH + 3 * HEADER_BYTES
        )
        self.identity = digest(
            [
                PROFILE_VERSION,
                VERSION,
                self.profile_sha256,
                seed,
                prefer_warm,
                self.settings.parent_memory,
            ]
        )
        self.inner = None
        self.catalogue = None
        self.selected = None
        self.origin = None
        self.closed = False
        self.lock = threading.Lock()
        self.last_stats, self.lifetime_stats = {}, {}
        self.selector_work = 0
        self.selector_cpu = 0.0

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        with self.lock:
            if self.inner is not None:
                self.inner.close()
            self.catalogue = None
            self.closed = True

    def _load(self, budget):
        from .worker_contract import WorkerSettings

        if self.catalogue is not None:
            return
        budget.consume(0)
        size = self.profile.stat().st_size
        if not 0 < size <= PROFILE_BYTES:
            raise ValueError("profile exceeds bounded ownership")
        budget.consume(size)
        self.selector_work += size
        with self.profile.open("rb") as stream:
            raw = stream.read(PROFILE_BYTES + 1)
        if len(raw) != size or hashlib.sha256(raw).hexdigest() != (
            self.profile_sha256
        ):
            raise ValueError("profile bytes changed")
        data = json.loads(raw)
        if set(data) != {
            "version",
            "status",
            "executor_version",
            "search_version",
            "source_sha256",
            "serial_settings",
            "rows",
        } or (
            data["version"] != PROFILE_VERSION
            or data["executor_version"] != VERSION
            or data["search_version"] != SEARCH_VERSION
            or data["status"]
            not in ("training", "frozen-for-confirmation", "confirmed")
            or set(data["source_sha256"]) != SOURCE_NAMES
        ):
            raise ValueError("unsupported selector profile")
        total = 0
        for name, expected in data["source_sha256"].items():
            _checksum(expected)
            path = SOURCE_ROOT / name
            size = path.stat().st_size
            total += size
            if total > SOURCE_BYTES:
                raise ValueError("source compatibility input exceeds cap")
            budget.consume(size)
            self.selector_work += size
            actual, consumed = hashlib.sha256(), 0
            with path.open("rb") as stream:
                while True:
                    chunk = stream.read(min(65536, size - consumed + 1))
                    if not chunk:
                        break
                    consumed += len(chunk)
                    if consumed > size:
                        raise ValueError(
                            "source grew during compatibility check"
                        )
                    actual.update(chunk)
            if consumed != size or actual.hexdigest() != expected:
                raise ValueError("stale selector source: " + name)
        rows = data["rows"]
        if not isinstance(rows, list) or len(rows) > 256:
            raise ValueError("selector row capacity exceeded")
        keys = set()
        serial = WorkerSettings(**data["serial_settings"])
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
            for mode in ("cold", "warm"):
                candidate = row[mode]
                if set(candidate) != {"settings", "costs"}:
                    raise ValueError("invalid selection candidate")
                candidate["settings"] = asdict(
                    WorkerSettings(**candidate["settings"])
                )
                costs = SelectionCosts(**candidate["costs"])
                if (costs.digits, costs.options) != key:
                    raise ValueError("costs belong to another regime")
        self.catalogue = data

    def _select(self, n, config, saved, origin):
        from .ecm_workers import ECMExecutor, arithmetic_options
        from .worker_contract import WorkerSettings

        data = self.catalogue
        mode = "warm" if self.prefer_warm else "cold"
        options = digest(
            arithmetic_options(
                config, WorkerSettings(worker_memory=64 * 2**20)
            )
        )
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
        integer(origin, "selected origin", -1, len(data["rows"]) - 1)
        row = None if origin == -1 else data["rows"][origin]
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
        curves = sum(tier[2] for tier in tiers)
        settings = WorkerSettings(**candidate)
        if settings.worker_memory != 64 * 2**20:
            raise ValueError("selector mathematical memory changed")
        settings = replace(settings, parent_memory=self.settings.parent_memory)
        required = settings.checkpoint_bound(curves)
        settings = replace(
            settings,
            checkpoint_bytes=max(
                2**20, (required + 2**20 - 1) // 2**20 * 2**20
            ),
        )
        if settings.checkpoint_bound(sum(t[2] for t in config.ecm_tiers)) > (
            settings.checkpoint_bytes
        ):
            raise MemoryError("retained selector checkpoint capacity refused")
        if saved is not None and asdict(settings) != saved:
            raise ValueError("saved selector settings differ from its policy")
        if self.inner is not None:
            if settings != self.inner.settings or origin != self.origin:
                raise ValueError("retained selector cannot change workers")
            return
        normalized = asdict(WorkerSettings(**candidate))
        normalized.pop("parent_memory", None)
        normalized.pop("checkpoint_bytes", None)
        costs = []
        for entry in data["rows"]:
            proposed = dict(entry[mode]["settings"])
            proposed.pop("parent_memory", None)
            proposed.pop("checkpoint_bytes", None)
            if proposed == normalized:
                costs.append(SelectionCosts(**entry[mode]["costs"]))
        rule = SelectionRule(
            self.profile_sha256,
            digest([VERSION, "cost-policy", asdict(settings)]),
            tuple(costs),
        )
        self.inner = ECMExecutor(seed=self.seed, settings=settings, rule=rule)
        self.selected = asdict(settings)
        self.origin = origin

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
        if not self.lock.acquire(blocking=False):
            raise RuntimeError("selector context already active")
        try:
            return self._run(
                n, config, budget, lineage, checkpoint, fixed_work
            )
        finally:
            self.lock.release()

    def _run(self, n, config, budget, lineage, checkpoint, fixed_work):
        if self.closed:
            raise RuntimeError("selector context is closed")
        actual_identity = digest(
            [
                PROFILE_VERSION,
                VERSION,
                self.profile_sha256,
                self.seed,
                self.prefer_warm,
                self.settings.parent_memory,
            ]
        )
        if actual_identity != self.identity:
            raise ValueError("selector configuration mutated")
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
                or len(header) != 6
                or set(header)
                != {
                    "version",
                    "identity",
                    "origin",
                    "settings",
                    "budget",
                    "inner",
                }
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
                or set(header)
                != {
                    "version",
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
        self._load(budget)
        budget.consume(HEADER_BYTES)
        self.selector_work += HEADER_BYTES
        if self.inner is None or saved is not None:
            self._select(n, config, saved, origin)
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
        self.selector_cpu += max(
            0,
            time.process_time() - started - result.stats.get("parent_cpu", 0),
        )
        self.last_stats = dict(
            result.stats,
            selector_profile=self.profile_sha256,
            selector_work=self.selector_work,
            selector_cpu=self.selector_cpu,
        )
        self.lifetime_stats = dict(self.inner.lifetime_stats)
        self.lifetime_stats["charged_work"] = (
            self.lifetime_stats.get("charged_work", 0) + self.selector_work
        )
        return replace(result, checkpoint=checkpoint, stats=self.last_stats)
