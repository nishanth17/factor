"""Pure B12 identities, framing and leases; no factorer or process imports."""

import hashlib
import json
import math
import random
import struct
from dataclasses import asdict, dataclass

from .resources import (
    RemovedResourceLimitError,
    decode_checkpoint,
    encode_checkpoint,
    removed_arguments,
)

SEARCH_VERSION = "b12-f4-v1"
VERSION = "b12-f4-observed-v9"
RESIDUE_VERSION = "ecm-stage-boundary-v1"
RESULT_BYTES = 4096


def integer(value, name, low=0, high=None):
    """Reject booleans and unbounded protocol integers."""
    if (
        type(value) is not int
        or value < low
        or (high is not None and value > high)
    ):
        raise ValueError(name + " outside finite bounds")
    return value


def finite(value, name, low=0):
    if (
        type(value) not in (int, float)
        or not math.isfinite(value)
        or value < low
    ):
        raise ValueError(name + " must be finite and nonnegative")
    return value


def canonical(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    )


def digest(value):
    return hashlib.sha256(canonical(value).encode("ascii")).hexdigest()


@removed_arguments(
    (
        "assignment_cpu",
        "assignment_work",
        "checkpoint_bytes",
        "frame_bytes",
        "max_runs",
        "memory_bytes",
        "parent_memory",
        "residue_bytes",
        "startup_work",
        "worker_memory",
    ),
    positional=(1, "assignment_work"),
)
@dataclass(frozen=True)
class WorkerSettings:
    """Finite pool ownership; acceleration choices remain explicitly opt-in."""

    workers: int = 1
    startup_seconds: float = 10.0
    drain_seconds: float = 1.0
    kill_seconds: float = 1.0
    max_attempts: int = 3
    assignment_policy: str = "dynamic"
    packet_curves: int = 1
    poll_every: int = 1
    metadata_slots: int = 1
    decoded_cache_bytes: int = 0
    packed_cache_bytes: int = 0
    manifest_cache_bytes: int = 0
    worker_imports: str = "portfolio"
    stage_two_kernel: str = "legacy"
    stage_two_cache_bytes: int = 0
    budget_mode: str = "legacy"
    bind_context: bool = False
    packet_policy: str = "fixed"
    stage_policy: str = "whole"
    residue_capacity: int = 16
    allow_wide_workers: bool = False
    parent_poll_seconds: float = 0.01

    def __post_init__(self):
        integer(self.workers, "workers", 1, 16)
        if type(self.allow_wide_workers) is not bool:
            raise ValueError("wider-worker opt-in must be Boolean")
        if self.workers not in (1, 2, 4) and not self.allow_wide_workers:
            raise ValueError("workers outside 1/2/4 require explicit opt-in")
        for name, low, high in (
            ("max_attempts", 1, 3),
            ("packet_curves", 1, 4),
            ("poll_every", 1, 64),
            ("metadata_slots", 1, 2),
            ("residue_capacity", 1, 64),
        ):
            integer(getattr(self, name), name, low, high)
        for name in ("startup_seconds", "drain_seconds", "kill_seconds"):
            finite(getattr(self, name), name, 0.001)
        finite(self.parent_poll_seconds, "parent poll", 0.001)
        if self.parent_poll_seconds > 0.05:
            raise ValueError("parent poll exceeds cancellation bound")
        integer(self.decoded_cache_bytes, "decoded cache", 0, 64 * 2**20)
        integer(self.packed_cache_bytes, "packed cache", 0, 16 * 2**20)
        integer(self.manifest_cache_bytes, "manifest cache", 0, 16 * 2**20)
        integer(self.stage_two_cache_bytes, "stage-two cache", 0, 64 * 2**20)
        if self.budget_mode not in ("legacy", "prospective"):
            raise ValueError("unknown worker budget mode")
        if type(self.bind_context) is not bool:
            raise ValueError("context binding must be Boolean")
        if self.packet_policy not in ("fixed", "tail_one"):
            raise ValueError("unknown packet policy")
        if self.worker_imports not in ("portfolio", "ecm"):
            raise ValueError("unknown worker import path")
        if self.stage_two_kernel not in ("legacy", "reduced"):
            raise ValueError("unknown stage-two kernel")
        if self.assignment_policy not in ("dynamic", "static"):
            raise ValueError("unknown curve assignment policy")
        if self.stage_policy not in ("whole", "staged"):
            raise ValueError("unknown stage policy")

    @property
    def owned_bytes(self):
        """Diagnostic retained-cache estimate; current RSS is observed."""
        return self.manifest_cache_bytes + (
            self.workers + 1
        ) * self.metadata_slots * (
            self.decoded_cache_bytes
            + self.packed_cache_bytes
            + self.stage_two_cache_bytes
        )

    def checkpoint_bound(self, curves):
        """Diagnostic estimate; saving is never required for dispatch."""
        return 65536 + curves * (RESULT_BYTES + 256 + self.max_attempts * 2048)


@dataclass(frozen=True)
class CostRow:
    """One source/config/size-bound row; bad estimates disable that row."""

    schedule: str = ""
    cold_seconds: float | None = None
    warm_seconds: float | None = None
    seconds_per_curve: tuple = ()
    pool_seconds: float | None = None
    seconds_per_stage_two: tuple = ()
    execution: str = ""

    def __post_init__(self):
        for name in ("cold_seconds", "warm_seconds", "pool_seconds"):
            try:
                finite(getattr(self, name), name, 0.001)
            except (ValueError, TypeError):
                object.__setattr__(self, name, None)
        for name in ("seconds_per_curve", "seconds_per_stage_two"):
            try:
                costs = tuple(getattr(self, name))
                if len(costs) > 32:
                    raise ValueError("too many cost rows")
                for cost in costs:
                    finite(cost, "curve cost", 0.000001)
            except (ValueError, TypeError):
                costs = ()
            object.__setattr__(self, name, costs)
        for name in ("schedule", "execution"):
            if (
                not isinstance(getattr(self, name), str)
                or len(getattr(self, name)) > 64
            ):
                object.__setattr__(self, name, "")

    def decision(
        self,
        manifest,
        *,
        ready,
        pool_started=False,
        remaining=None,
        remaining_stage_two=None,
        execution=None,
    ):
        state = "warm" if ready else "pool" if pool_started else "cold"
        threshold = (
            self.warm_seconds
            if ready
            else (
                self.pool_seconds
                if pool_started and self.pool_seconds is not None
                else self.cold_seconds
            )
        )
        result = dict(
            engage=False,
            state=state,
            estimate=None,
            threshold=threshold,
            reason="invalid_cost",
        )
        try:
            finite(threshold, "threshold", 0.001)
            if self.schedule != manifest.schedule or (
                execution is not None and self.execution != execution
            ):
                return result
            if len(self.seconds_per_curve) != len(manifest.tiers):
                return result
            counts = (
                tuple(t[2] for t in manifest.tiers)
                if remaining is None
                else tuple(remaining)
            )
            continuations = (
                (0,) * len(counts)
                if remaining_stage_two is None
                else tuple(remaining_stage_two)
            )
            if len(counts) != len(manifest.tiers) or len(continuations) != len(
                counts
            ):
                return result
            if any(continuations) and len(self.seconds_per_stage_two) != len(
                counts
            ):
                return result
            estimate = 0.0
            for i, (count, continued, tier) in enumerate(
                zip(counts, continuations, manifest.tiers)
            ):
                integer(count, "remaining curves", 0, tier[2])
                integer(
                    continued, "remaining continuations", 0, tier[2] - count
                )
                estimate += self.seconds_per_curve[i] * count
                if continued:
                    estimate += self.seconds_per_stage_two[i] * continued
            finite(estimate, "total estimate")
        except (ValueError, TypeError, OverflowError):
            return result
        return dict(
            result,
            engage=estimate >= threshold,
            estimate=estimate,
            reason="threshold",
        )


@dataclass(frozen=True)
class EngagementRule(CostRow):
    """Bounded cost table; no factor labels or future outcomes are inputs."""

    cost_table: tuple = ()

    def __post_init__(self):
        super().__post_init__()
        try:
            rows = tuple(self.cost_table)
            if len(rows) > 32:
                raise ValueError("cost table capacity exceeded")
            rows = tuple(
                row
                if isinstance(row, CostRow)
                and not isinstance(row, EngagementRule)
                else CostRow(**row)
                for row in rows
            )
            keys = [row.schedule for row in rows]
            if len(set(keys)) != len(keys) or (
                self.schedule and self.schedule in keys
            ):
                raise ValueError("duplicate cost identities")
        except (ValueError, TypeError):
            rows = ()
        object.__setattr__(self, "cost_table", rows)

    def decision(self, manifest, **kwargs):
        if self.schedule == manifest.schedule:
            return super().decision(manifest, **kwargs)
        for row in self.cost_table:
            if row.schedule == manifest.schedule:
                return row.decision(manifest, **kwargs)
        return dict(
            engage=False,
            state="unknown",
            estimate=None,
            threshold=None,
            reason="missing_cost",
        )

    def engages(self, manifest, **kwargs):
        return self.decision(manifest, **kwargs)["engage"]


def metadata_identity(options):
    """Only immutable schedule data: no point, modulus or modulus size."""
    return digest([VERSION, "immutable-ecm-metadata", options])


class Manifest:
    """Immutable mathematical search; no worker count enters curve identity."""

    def __init__(self, n, seed, lineage, options):
        integer(n, "cofactor", 2, 2**4096 - 1)
        integer(seed, "method seed", 0, 2**256 - 1)
        if not isinstance(lineage, str):
            raise ValueError("bounded cofactor lineage required")
        # Canonical copying prevents callers mutating the frozen schedule.
        self._options = canonical(options)
        tiers = options["ecm_tiers"]
        if not isinstance(tiers, (list, tuple)) or not 1 <= len(tiers) <= 32:
            raise ValueError("finite ECM tiers required")
        for b1, b2, curves in tiers:
            integer(b1, "B1", 2, 2**32 - 1)
            integer(b2, "B2", b1, 2**32 - 1)
            integer(curves, "curves", 0, 4096)
        if sum(t[2] for t in tiers) > 4096:
            raise ValueError("curve manifest exceeds 4096 assignments")
        if options["backend"] != "python-int":
            raise ValueError(
                "initial worker executor supports python-int only"
            )
        self.n, self.seed, self.lineage = n, seed, lineage
        self.tiers = tuple(tuple(t) for t in tiers)
        self.schedule = digest([VERSION, "cost", options, n.bit_length()])
        self.metadata = metadata_identity(options)
        # The former 64 MiB worker grant entered the established curve hash.
        # Retain that hash input solely as a compatibility identity; no byte
        # allowance is constructed or enforced from it. Removing a resource
        # argument must not silently change sigmas for the same root seed.
        search_options = (
            dict(options, memory_bytes=64 * 2**20)
            if "gcd_batch" in options
            else options
        )
        self.identity = digest(
            [SEARCH_VERSION, n, seed, lineage, search_options]
        )
        self._assignments = {}
        self.jobs = tuple(
            (tier, index)
            for tier, (_, _, curves) in enumerate(self.tiers)
            for index in range(curves)
        )

    @property
    def options(self):
        return json.loads(self._options)

    def assignment(self, position):
        integer(position, "position", 0, len(self.jobs) - 1)
        if position in self._assignments:
            return dict(self._assignments[position])
        tier, index = self.jobs[position]
        b1, b2, _ = self.tiers[tier]
        identity = digest([self.identity, tier, index])
        seed = int(digest([SEARCH_VERSION, "curve", identity]), 16)
        result = dict(
            id=identity,
            position=position,
            tier=tier,
            curve=index,
            n=self.n,
            seed=seed,
            b1=b1,
            b2=b2,
            sigma=random.Random(seed).randrange(6, 2**63),
        )
        self._assignments[position] = result
        return dict(result)


class LeaseRefusedError(Exception):
    """A lifecycle watchdog or numerical retry boundary was reached."""


class Ledger:
    """Assignment identities and observed reports; no prepaid allowance.

    The root Budget observes process CPU, including failed/killed workers.
    Reports here count algorithmic work and support checked replay; they do
    not reserve, refund or independently charge CPU into the root owner.
    """

    def __init__(self, *old, **removed):
        if old or removed:
            name = next(iter(removed)) if removed else "work"
            raise RemovedResourceLimitError(name, "the shared root Budget")
        self.work = 0
        self.cpu = self.parent_cpu = 0.0
        self.active = {}
        self.next_token = self.lost = 0
        self.receipts = []

    def charge(self, work=0, *, parent_cpu=None):
        integer(work, "controller work")
        if parent_cpu is not None:
            finite(parent_cpu, "parent CPU")
            if parent_cpu < self.parent_cpu:
                raise ValueError("parent CPU counter decreased")
            self.parent_cpu = parent_cpu
        self.work += work

    @removed_arguments(("work", "cpu"))
    def reserve(self, identity):
        """Issue a token before dispatch; the token carries no resources."""
        self.next_token += 1
        token = self.next_token
        self.active[token] = identity
        self.receipts.append(dict(token=token, id=identity, status="issued"))
        return token

    def settle(self, token, identity, work, cpu):
        integer(work, "reported work")
        finite(cpu, "reported CPU")
        if self.active[token] != identity:
            raise ValueError("worker identity changed")
        del self.active[token]
        self.receipts[token - 1].update(status="reported", work=work, cpu=cpu)
        self.work += work
        self.cpu += cpu

    def lose(self, token):
        self.active.pop(token)
        self.receipts[token - 1].update(status="lost", work=0, cpu=0.0)
        self.lost += 1


@removed_arguments(("limit",))
def checked_residue(residue, assignment, metadata):
    """Validate a canonical boundary before any stage-two tables exist."""
    if not isinstance(residue, dict):
        raise ValueError("oversized/noncanonical residue")
    required = {
        "version",
        "id",
        "metadata",
        "n",
        "b1",
        "b2",
        "sigma",
        "value",
        "a24",
        "chain_chunks",
        "sha256",
    }
    if set(residue) != required:
        raise ValueError("invalid residue fields")
    payload = {key: value for key, value in residue.items() if key != "sha256"}
    if (
        residue["sha256"] != digest(payload)
        or residue["version"] != RESIDUE_VERSION
    ):
        raise ValueError("corrupt residue")
    for key in ("id", "n", "b1", "b2", "sigma"):
        if (
            type(residue[key]) is not type(assignment[key])
            or residue[key] != assignment[key]
        ):
            raise ValueError("residue assignment changed")
    if residue["metadata"] != metadata:
        raise ValueError("residue schedule changed")
    n = assignment["n"]
    if not isinstance(residue["value"], list) or len(residue["value"]) != 2:
        raise ValueError("invalid residue point")
    for value in (*residue["value"], residue["a24"]):
        integer(value, "canonical residue integer", 0, n - 1)
    integer(residue["chain_chunks"], "chain chunks")
    if math.gcd(residue["value"][1], n) != 1:
        raise ValueError("boundary has unresolved nonunit/saturation")
    return json.loads(canonical(residue))


class Search:
    """Parent-owned curve/stage traversal and a bounded residue queue."""

    def __init__(self, manifest, settings, *, fixed_work=False):
        if type(fixed_work) is not bool:
            raise ValueError("fixed_work must be Boolean")
        self.manifest, self.settings = manifest, settings
        self.fixed_work = fixed_work
        self.attempts = [0] * len(manifest.jobs)
        self.continuation_attempts = [0] * len(manifest.jobs)
        self.failed_attempts = [0] * len(manifest.jobs)
        self.failed_continuation_attempts = [0] * len(manifest.jobs)
        self.boundaries = {}
        self.completed, self.active = {}, {}
        self.factor, self.retired, self.runs = None, False, 0

    def stage(self, position):
        if self.settings.stage_policy == "whole":
            return "curve"
        return "stage_two" if position in self.boundaries else "stage_one"

    def attempt(self, position):
        counts = (
            self.continuation_attempts
            if position in self.boundaries
            else self.attempts
        )
        return counts[position]

    def failures(self, position):
        counts = (
            self.failed_continuation_attempts
            if position in self.boundaries
            else self.failed_attempts
        )
        return counts[position]

    def record_failure(self, position):
        counts = (
            self.failed_continuation_attempts
            if position in self.boundaries
            else self.failed_attempts
        )
        counts[position] += 1

    def next(self, slot=None):
        remaining = [
            i for i in range(len(self.attempts)) if i not in self.completed
        ]
        if not remaining or self.retired:
            return None
        tier = self.manifest.jobs[remaining[0]][0]
        candidates = [
            i
            for i in remaining
            if self.manifest.jobs[i][0] == tier
            and i not in self.active
            and (
                slot is None
                or self.settings.assignment_policy != "static"
                or i % self.settings.workers == slot
            )
        ]
        if self.settings.stage_policy == "staged":
            # In-flight stage-one jobs reserve future queue positions too.
            occupied = len(self.boundaries) + sum(
                i not in self.boundaries for i in self.active
            )
            first = [i for i in candidates if i not in self.boundaries]
            second = [i for i in candidates if i in self.boundaries]
            candidates = (
                first + second
                if occupied < self.settings.residue_capacity
                else second
            )
        for position in candidates:
            if self.failures(position) >= self.settings.max_attempts:
                raise LeaseRefusedError("attempt_limit")
            return position
        return None

    def dispatchable_count(self):
        remaining = [
            i for i in range(len(self.attempts)) if i not in self.completed
        ]
        if not remaining or self.retired:
            return 0
        tier = self.manifest.jobs[remaining[0]][0]
        pending = [
            i
            for i in remaining
            if i not in self.active and self.manifest.jobs[i][0] == tier
        ]
        if self.settings.stage_policy == "whole":
            return len(pending)
        occupied = len(self.boundaries) + sum(
            i not in self.boundaries for i in self.active
        )
        continuations = sum(i in self.boundaries for i in pending)
        return continuations + min(
            len(pending) - continuations,
            max(0, self.settings.residue_capacity - occupied),
        )

    def dispatch(self, position, token):
        integer(position, "dispatch position", 0, len(self.attempts) - 1)
        integer(token, "lease token", 1)
        if position in self.completed or position in self.active:
            raise ValueError("duplicate dispatch")
        if self.failures(position) >= self.settings.max_attempts:
            raise LeaseRefusedError("attempt_limit")
        if self.stage(position) == "stage_one" and (
            len(self.boundaries)
            + sum(i not in self.boundaries for i in self.active)
            >= self.settings.residue_capacity
        ):
            raise LeaseRefusedError("residue_capacity")
        counts = (
            self.continuation_attempts
            if position in self.boundaries
            else self.attempts
        )
        counts[position] += 1
        self.active[position] = (token, counts[position])

    def _validate_result(self, position, token, message, stage):
        allowed = {
            "id",
            "token",
            "attempt",
            "work",
            "factor",
            "reason",
            "failure",
            "phase",
            "sigma",
            "prepared_schedule",
            "metadata_cache",
            "chain_chunks",
            "preparations",
            "kind",
            "cpu",
            "rss",
            "packet",
            "stage",
            "residue",
            "boundary_sha256",
            "cache_hits",
            "cache_evictions",
            "decode_hits",
            "decode_misses",
            "decode_refusals",
            "stage_two_hits",
            "stage_two_misses",
            "stage_two_refusals",
            "stage_two_build_bytes",
        }
        if not message.keys() <= allowed:
            raise ValueError("unbounded result record")
        integer(message["attempt"], "result attempt", 1)
        integer(message["token"], "result token", 1)
        for key in (
            "work",
            "rss",
            "chain_chunks",
            "preparations",
            "packet",
            "cache_hits",
            "cache_evictions",
            "decode_hits",
            "decode_misses",
            "decode_refusals",
            "stage_two_hits",
            "stage_two_misses",
            "stage_two_refusals",
            "stage_two_build_bytes",
        ):
            if key in message:
                integer(message[key], key)
        if "cpu" in message:
            finite(message["cpu"], "worker CPU")
        for key, bound in (("phase", 64), ("failure", 256)):
            if message.get(key) is not None and (
                not isinstance(message[key], str) or len(message[key]) > bound
            ):
                raise ValueError("unbounded result label")
        assignment = self.manifest.assignment(position)
        if (
            message.get("id") != assignment["id"]
            or message.get("token") != token
        ):
            raise ValueError("result belongs to another assignment")
        if message.get("stage", "curve") != stage:
            raise ValueError("result stage changed")
        if (
            message.get("sigma") is not None
            and message["sigma"] != assignment["sigma"]
        ):
            raise ValueError("curve sigma differs from manifest")
        factor = message.get("factor")
        if factor is not None and (
            type(factor) is not int
            or not 1 < factor < self.manifest.n
            or self.manifest.n % factor
        ):
            raise ValueError("worker returned an invalid proper divisor")
        if message.get("reason") not in (
            "complete",
            "boundary",
            "cancelled",
            "wall_limit",
            "cpu_limit",
            "error",
            "allocation_failure",
        ):
            raise ValueError("invalid worker outcome")
        if message["reason"] != "complete" and factor is not None:
            raise ValueError("unfinished result cannot publish a factor")
        if (
            message["reason"] != "boundary"
            and message.get("residue") is not None
        ):
            raise ValueError("unexpected residue")
        if stage == "stage_two":
            checksum = message.get("boundary_sha256")
            if not isinstance(checksum, str) or len(checksum) != 64:
                raise ValueError("missing continuation provenance")
        return factor

    def receive(self, position, token, message, *, publish=True):
        stage = self.stage(position)
        factor = self._validate_result(position, token, message, stage)
        if self.active.get(position) != (token, message["attempt"]):
            raise ValueError("duplicate or stale result")
        if (
            stage == "stage_two"
            and message["boundary_sha256"]
            != self.boundaries[position]["sha256"]
        ):
            raise ValueError("continuation boundary changed")
        boundary = None
        if message["reason"] == "boundary":
            if (
                stage != "stage_one"
                or message.get("phase") != "stage_two_setup"
            ):
                raise ValueError("handoff outside clean stage boundary")
            boundary = checked_residue(
                message.get("residue"),
                self.manifest.assignment(position),
                self.manifest.metadata,
            )
        del self.active[position]
        if message["reason"] in ("error", "allocation_failure"):
            self.record_failure(position)
        if not publish:
            # Draining settles costs and identities, but a refused batch
            # cannot publish a factor, boundary or completed search cursor.
            return
        if boundary is not None:
            self.boundaries[position] = boundary
        elif message["reason"] == "complete":
            self.boundaries.pop(position, None)
            self.completed[position] = dict(message)
            if factor is not None and publish and not self.retired:
                if self.factor is None:
                    self.factor = factor
                if not self.fixed_work:
                    self.retired = True

    def check_leases(self, receipts, saved):
        """Validate issued identities/attempts and observed saved charges."""
        if not isinstance(receipts, list):
            raise ValueError("invalid assignment history")
        identities = {
            digest([self.manifest.identity, tier, index]): position
            for position, (tier, index) in enumerate(self.manifest.jobs)
        }
        attempts, reported = {}, {}
        work, cpu = 0, 0.0
        for receipt in receipts:
            if not isinstance(receipt, dict):
                raise ValueError("invalid assignment receipt")
            integer(receipt["token"], "assignment token", 1)
            used = integer(receipt["work"], "reported work")
            consumed = finite(receipt["cpu"], "reported CPU")
            if receipt["status"] not in ("reported", "lost"):
                raise ValueError("checkpoint requires settled assignments")
            work += used
            cpu += consumed
            identity = receipt["id"]
            if identity not in identities:
                if identity not in (
                    "startup",
                    "shutdown",
                    "spawn_runtime_helper",
                ):
                    raise ValueError("unknown assignment receipt")
                continue
            position = identities[identity]
            stage = receipt["stage"]
            if stage not in ("curve", "stage_one", "stage_two"):
                raise ValueError("invalid receipt stage")
            attempt = integer(
                receipt["attempt"],
                "receipt attempt",
                1,
            )
            key = (position, stage, attempt)
            if key in attempts:
                raise ValueError("duplicate assignment attempt")
            attempts[key] = receipt
            if receipt["status"] == "reported":
                reported[key] = receipt
        for position, message in self.completed.items():
            key = (position, message.get("stage", "curve"), message["attempt"])
            receipt = reported.get(key)
            if receipt is None or receipt["token"] != message["token"]:
                raise ValueError("completed result lacks a matching receipt")
            if receipt["reason"] != message["reason"]:
                raise ValueError("completed receipt outcome changed")
        for position, residue in self.boundaries.items():
            receipt = reported.get(
                (position, "stage_one", self.attempts[position])
            )
            if (
                receipt is None
                or receipt["reason"] != "boundary"
                or receipt.get("boundary_sha256") != residue["sha256"]
            ):
                raise ValueError("boundary lacks a matching receipt")
        if saved["work"] < work or saved["cpu"] + 1e-9 < cpu:
            raise ValueError("resume would reset observed charges")

    @encode_checkpoint
    def snapshot(self):
        if self.active:
            raise ValueError("checkpoint requires drained leases")
        payload = dict(
            version=VERSION,
            runs=self.runs,
            manifest=self.manifest.identity,
            settings=asdict(self.settings),
            fixed_work=self.fixed_work,
            attempts=self.attempts,
            continuation_attempts=self.continuation_attempts,
            failed_attempts=self.failed_attempts,
            failed_continuation_attempts=self.failed_continuation_attempts,
            # Normalize before hashing: JSON decodes position keys as strings.
            boundaries={
                str(key): value for key, value in self.boundaries.items()
            },
            completed={
                str(key): value for key, value in self.completed.items()
            },
            factor=self.factor,
            retired=self.retired,
        )
        return dict(payload=payload, sha256=digest(payload))

    @classmethod
    @decode_checkpoint
    def restore(cls, checkpoint, manifest, settings, *, fixed_work=False):
        data = checkpoint["payload"]
        if checkpoint["sha256"] != digest(data) or data["version"] != VERSION:
            raise ValueError("corrupt/unsupported checkpoint")
        if data["manifest"] != manifest.identity or data["settings"] != asdict(
            settings
        ):
            raise ValueError(
                "resume identity/settings changed (including workers)"
            )
        result = cls(manifest, settings, fixed_work=fixed_work)
        if data["fixed_work"] != fixed_work:
            raise ValueError("resume traversal changed")
        result.runs = integer(data["runs"], "run count", 0)
        for key in ("attempts", "continuation_attempts"):
            values = data[key]
            if len(values) != len(result.attempts):
                raise ValueError("resume traversal changed")
            setattr(
                result,
                key,
                [integer(a, "dispatch ordinal", 0) for a in values],
            )
        for key in ("failed_attempts", "failed_continuation_attempts"):
            values = data[key]
            if not isinstance(values, list) or len(values) != len(
                result.attempts
            ):
                raise ValueError("retry counters changed shape")
            checked = [
                integer(value, "failed attempts", 0, settings.max_attempts)
                for value in values
            ]
            ordinals = (
                result.attempts
                if key == "failed_attempts"
                else result.continuation_attempts
            )
            if any(
                failures > issued
                for failures, issued in zip(checked, ordinals)
            ):
                raise ValueError("failures exceed issued attempts")
            setattr(result, key, checked)
        if len(data["boundaries"]) > settings.residue_capacity:
            raise ValueError("oversized boundary queue")
        for key, boundary in data["boundaries"].items():
            if (
                not isinstance(key, str)
                or not key.isdecimal()
                or str(int(key)) != key
            ):
                raise ValueError("noncanonical checkpoint position")
            position = integer(
                int(key), "boundary position", 0, len(result.attempts) - 1
            )
            if (
                settings.stage_policy != "staged"
                or result.attempts[position] == 0
                or position in result.boundaries
            ):
                raise ValueError("unissued/duplicate boundary")
            result.boundaries[position] = checked_residue(
                boundary, manifest.assignment(position), manifest.metadata
            )
        for key, message in data["completed"].items():
            if (
                not isinstance(key, str)
                or not key.isdecimal()
                or str(int(key)) != key
            ):
                raise ValueError("noncanonical checkpoint position")
            position = integer(
                int(key), "completed position", 0, len(result.attempts) - 1
            )
            stage = message.get("stage", "curve")
            allowed = (
                ("curve",)
                if settings.stage_policy == "whole"
                else ("stage_one", "stage_two")
            )
            if (
                position in result.completed
                or position in result.boundaries
                or message["reason"] != "complete"
                or stage not in allowed
            ):
                raise ValueError("duplicate/nonfinal checkpoint result")
            count = (
                result.continuation_attempts[position]
                if stage == "stage_two"
                else result.attempts[position]
            )
            if (
                not 1 <= message["attempt"] <= count
                or result.attempts[position] == 0
            ):
                raise ValueError("unissued checkpoint result")
            result._validate_result(position, message["token"], message, stage)
            result.completed[position] = dict(message)
        for position, count in enumerate(result.continuation_attempts):
            if count and (
                not result.attempts[position]
                or settings.stage_policy != "staged"
                or (
                    position not in result.boundaries
                    and position not in result.completed
                )
            ):
                raise ValueError(
                    "continuation has no boundary or final result"
                )
        factors = [
            m["factor"]
            for m in result.completed.values()
            if m["factor"] is not None
        ]
        if data["factor"] is not None and data["factor"] not in factors:
            raise ValueError("checkpoint factor has no acknowledged result")
        if type(data["retired"]) is not bool or (
            data["retired"] and not data["factor"]
        ):
            raise ValueError("invalid retirement")
        unfinished = [
            i for i in range(len(result.attempts)) if i not in result.completed
        ]
        if unfinished:
            tier = manifest.jobs[unfinished[0]][0]
            if any(
                count and manifest.jobs[i][0] > tier
                for i, count in enumerate(result.attempts)
            ):
                raise ValueError("checkpoint violates tier barrier")
        for position, message in result.completed.items():
            if (
                message.get("stage") == "stage_one"
                and result.continuation_attempts[position]
            ):
                raise ValueError("final stage one cannot follow continuation")
        result.factor, result.retired = data["factor"], data["retired"]
        return result


@removed_arguments(("limit", "frame_bytes"))
def frame(message):
    raw = canonical(message).encode("ascii")
    if len(raw) >= 2**32:
        raise ValueError("unsupported_encoding: frame length exceeds u32")
    return struct.pack("!I", len(raw)) + raw


class Decoder:
    """Incremental u32 frames; the byte buffer grows with received data."""

    @removed_arguments(("limit",))
    def __init__(self):
        self.buffer = bytearray()

    @property
    def limit(self):
        # This is the wire word domain, not a resource grant.
        return 2**32 - 1

    @property
    def complete(self):
        if len(self.buffer) < 4:
            return False
        size = struct.unpack("!I", self.buffer[:4])[0]
        return size == 0 or len(self.buffer) >= size + 4

    def feed(self, data):
        self.buffer.extend(data)
        if len(self.buffer) < 4:
            return None
        size = struct.unpack("!I", self.buffer[:4])[0]
        if size == 0:
            raise ValueError("invalid frame size")
        if len(self.buffer) < size + 4:
            return None
        stop = size + 4
        result = json.loads(bytes(self.buffer[4:stop]))
        del self.buffer[:stop]
        if type(result) is not dict:
            raise ValueError("message must be an object")
        return result
