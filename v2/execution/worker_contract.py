"""Pure B12 identities, framing and leases; no factorer or process imports."""

import hashlib
import json
import math
import random
import struct
from dataclasses import asdict, dataclass

SEARCH_VERSION = "b12-f4-v1"
VERSION = "b12-f4-v7"
RESIDUE_VERSION = "ecm-stage-boundary-v1"
RESULT_BYTES = 4096


def integer(value, name, low=0, high=2**63 - 1):
    """Reject booleans and unbounded protocol integers."""
    if type(value) is not int or not low <= value <= high:
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


@dataclass(frozen=True)
class WorkerSettings:
    """Finite pool ownership; acceleration choices remain explicitly opt-in."""

    workers: int = 1
    assignment_work: int = 32_000_000
    assignment_cpu: float = 30.0
    startup_work: int = 1_000_000
    startup_seconds: float = 10.0
    drain_seconds: float = 1.0
    kill_seconds: float = 1.0
    worker_memory: int = 64 * 2**20
    parent_memory: int = 32 * 2**20
    memory_bytes: int = 512 * 2**20
    frame_bytes: int = 16_384
    checkpoint_bytes: int = 2**20
    max_attempts: int = 3
    max_runs: int = 32
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
    residue_bytes: int = 8192
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
            ("max_runs", 1, 32),
            ("packet_curves", 1, 4),
            ("poll_every", 1, 64),
            ("metadata_slots", 1, 2),
            ("residue_capacity", 1, 64),
            ("residue_bytes", 1024, 16384),
        ):
            integer(getattr(self, name), name, low, high)
        for name in ("assignment_work", "startup_work"):
            integer(getattr(self, name), name, 1)
        for name in (
            "assignment_cpu",
            "startup_seconds",
            "drain_seconds",
            "kill_seconds",
        ):
            finite(getattr(self, name), name, 0.001)
        finite(self.parent_poll_seconds, "parent poll", 0.001)
        if self.parent_poll_seconds > 0.05:
            raise ValueError("parent poll exceeds cancellation bound")
        integer(self.frame_bytes, "frame_bytes", 1024, 65_536)
        integer(self.checkpoint_bytes, "checkpoint_bytes", 4096, 64 * 2**20)
        integer(self.worker_memory, "worker_memory", 2**20, 2**30)
        integer(self.parent_memory, "parent_memory", 2**20, 4 * 2**30)
        integer(self.memory_bytes, "memory_bytes", 2**20, 24 * 2**30)
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
        if (
            self.stage_policy == "staged"
            and self.residue_bytes + RESULT_BYTES > self.frame_bytes
        ):
            raise ValueError("stage result cannot fit the frame allowance")
        if self.owned_bytes > self.memory_bytes:
            raise MemoryError("simultaneous pool ownership refused")

    @property
    def shutdown_cpu(self):
        # A final interpreter teardown report is unavailable: retain in full.
        return self.drain_seconds + 2 * self.kill_seconds + 1

    @property
    def owned_bytes(self):
        contexts = (self.workers + 1) * self.metadata_slots
        residues = (
            self.residue_capacity if self.stage_policy == "staged" else 0
        )
        return (
            self.parent_memory
            + self.manifest_cache_bytes
            + contexts
            * (
                self.worker_memory
                + self.decoded_cache_bytes
                + self.packed_cache_bytes
                + self.stage_two_cache_bytes
                + 4096 * self.bind_context
            )
            + self.worker_memory  # Possible runtime-managed spawn helper.
            + self.workers
            * (12 * self.frame_bytes + RESULT_BYTES * self.packet_curves)
            + 4 * residues * self.residue_bytes
            + 3 * self.checkpoint_bytes
        )

    def checkpoint_bound(self, curves):
        stages = 2 if self.stage_policy == "staged" else 1
        residues = self.residue_capacity if stages == 2 else 0
        return (
            65536
            + curves * (RESULT_BYTES + 256 + stages * self.max_attempts * 1024)
            + self.max_runs * (self.workers * 2048 + 1024)
            + residues * self.residue_bytes
        )


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
        if not isinstance(lineage, str) or len(lineage) > 2048:
            raise ValueError("bounded cofactor lineage required")
        # Canonical copying prevents callers mutating the frozen schedule.
        self._options = canonical(options)
        if len(self._options) > 12_000:
            raise ValueError("arithmetic configuration exceeds transport cap")
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
        # Execution choice must not silently select a different curve stream.
        # Metadata and native snapshot identity still bind the actual executor.
        search_options = dict(options)
        search_options.pop("ecm_kernel", None)
        self.identity = digest(
            [SEARCH_VERSION, n, seed, lineage, search_options]
        )
        self.execution = None
        if options.get("ecm_kernel", "python") != "python":
            from ..ecm.native_packet import source_identity

            self.execution = digest([source_identity(), options["ecm_kernel"]])
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
    """A lease could not be admitted without exceeding the finite envelope."""


class Ledger:
    """Reserve before dispatch; refund only checked reports."""

    def __init__(self, work, cpu):
        self.work_limit = integer(work, "work allowance")
        self.cpu_limit = finite(cpu, "CPU allowance")
        self.work = 0
        self.cpu = 0.0
        self.parent_cpu = 0.0
        self.active = {}
        self.next_token = 0
        self.lost = 0
        self.receipts = []

    def charge(self, work=0, *, parent_cpu=None):
        integer(work, "controller work")
        if parent_cpu is not None:
            finite(parent_cpu, "parent CPU")
            if parent_cpu < self.parent_cpu:
                raise ValueError("parent CPU counter decreased")
            self.parent_cpu = parent_cpu
        if self.work + work > self.work_limit:
            raise LeaseRefusedError("work_limit")
        if self.cpu + self.parent_cpu > self.cpu_limit:
            raise LeaseRefusedError("cpu_limit")
        self.work += work

    def reserve(self, identity, work, cpu):
        integer(work, "lease work")
        finite(cpu, "lease CPU")
        self.charge(work)
        if self.cpu + self.parent_cpu + cpu > self.cpu_limit:
            self.work -= work
            raise LeaseRefusedError("cpu_limit")
        self.next_token += 1
        token = self.next_token
        self.active[token] = (identity, work, cpu)
        self.receipts.append(
            dict(
                token=token,
                id=identity,
                work_grant=work,
                cpu_grant=cpu,
                status="reserved",
            )
        )
        self.cpu += cpu
        return token

    def settle(self, token, identity, work, cpu):
        """Validate fully before refunding; a duplicate cannot refund twice."""
        integer(work, "reported work")
        finite(cpu, "reported CPU")
        expected, grant_work, grant_cpu = self.active[token]
        if identity != expected or work > grant_work or cpu > grant_cpu:
            raise ValueError("worker identity or lease overrun")
        del self.active[token]
        self.receipts[token - 1].update(status="reported", work=work, cpu=cpu)
        self.work -= grant_work - work
        self.cpu = max(0.0, self.cpu - (grant_cpu - cpu))

    def lose(self, token):
        self.active.pop(token)
        self.receipts[token - 1]["status"] = "lost"
        self.lost += 1


def checked_residue(residue, assignment, metadata, limit):
    """Validate a canonical boundary before any stage-two tables exist."""
    if not isinstance(residue, dict) or len(canonical(residue)) > limit:
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
        if (
            settings.checkpoint_bound(len(manifest.jobs))
            > settings.checkpoint_bytes
        ):
            raise MemoryError("checkpoint ownership refused before dispatch")
        self.fixed_work = fixed_work
        self.attempts = [0] * len(manifest.jobs)
        self.continuation_attempts = [0] * len(manifest.jobs)
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
            if self.attempt(position) >= self.settings.max_attempts:
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
        if self.attempt(position) >= self.settings.max_attempts:
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
        limit = RESULT_BYTES + (
            self.settings.residue_bytes
            if message.get("reason") == "boundary"
            else 0
        )
        if not message.keys() <= allowed or len(canonical(message)) > limit:
            raise ValueError("unbounded result record")
        integer(
            message["attempt"], "result attempt", 1, self.settings.max_attempts
        )
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
            "work_limit",
            "error",
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
                self.settings.residue_bytes,
            )
        del self.active[position]
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
        """Reject unpaid traversal or charges below retained lease evidence."""
        maximum = len(
            self.attempts
        ) * self.settings.max_attempts * 2 + self.settings.max_runs * (
            2 * self.settings.workers + 1
        )
        if not isinstance(receipts, list) or len(receipts) > maximum:
            raise ValueError("lease history exceeds ownership")
        identities = {
            digest([self.manifest.identity, tier, index]): position
            for position, (tier, index) in enumerate(self.manifest.jobs)
        }
        issued, acknowledged, phase_attempts = {}, {}, {}
        charged_work, charged_cpu = 0, []
        for receipt in receipts:
            if not isinstance(receipt, dict) or len(canonical(receipt)) > 1024:
                raise ValueError("invalid bounded lease receipt")
            integer(receipt["token"], "lease token", 1)
            work = integer(receipt["work_grant"], "work grant")
            cpu = finite(receipt["cpu_grant"], "CPU grant")
            identity = receipt["id"]
            if receipt["status"] == "reported":
                used, consumed = (
                    integer(receipt["work"], "used work"),
                    finite(receipt["cpu"], "used CPU"),
                )
                if (
                    used > work
                    or consumed > cpu
                    or "overrun_work" in receipt
                    or "overrun_cpu" in receipt
                ):
                    raise ValueError("invalid lease refund")
            elif receipt["status"] == "lost":
                used = max(
                    work,
                    integer(receipt.get("overrun_work", 0), "overrun work"),
                )
                consumed = max(
                    cpu, finite(receipt.get("overrun_cpu", 0), "overrun CPU")
                )
            else:
                raise ValueError("checkpoint has unsettled lease")
            charged_work += used
            charged_cpu.append(consumed)
            if identity not in identities:
                if identity not in (
                    "startup",
                    "shutdown",
                    "spawn_runtime_helper",
                ):
                    raise ValueError("foreign lease identity")
                continue
            position = identities[identity]
            stage = receipt.get("stage")
            allowed = (
                ("curve",)
                if self.settings.stage_policy == "whole"
                else ("stage_one", "stage_two")
            )
            if stage not in allowed:
                raise ValueError("missing or unsupported lease stage")
            attempt = integer(
                receipt.get("attempt"),
                "lease attempt",
                1,
                self.settings.max_attempts,
            )
            key = (position, stage, attempt)
            if key in issued:
                raise ValueError("duplicate assignment lease")
            issued[key] = receipt
            phase_attempts.setdefault((position, stage), set()).add(attempt)
            if receipt["status"] == "reported":
                acknowledged[key] = receipt
        for position in range(len(self.attempts)):
            first = (
                "curve"
                if self.settings.stage_policy == "whole"
                else "stage_one"
            )
            for stage, count in (
                (first, self.attempts[position]),
                ("stage_two", self.continuation_attempts[position]),
            ):
                actual = phase_attempts.get((position, stage), set())
                if actual != set(range(1, count + 1)):
                    raise ValueError("traversal and lease attempts disagree")
            if self.continuation_attempts[position]:
                first_ack = acknowledged.get(
                    (position, "stage_one", self.attempts[position])
                )
                if first_ack is None or first_ack.get("reason") != "boundary":
                    raise ValueError(
                        "continuation has no paid stage-one boundary"
                    )
                checksum = first_ack.get("boundary_sha256")
                for attempt in range(
                    1, self.continuation_attempts[position] + 1
                ):
                    row = issued[(position, "stage_two", attempt)]
                    if row.get("boundary_sha256") != checksum:
                        raise ValueError(
                            "continuation lease changed its boundary"
                        )
            message = self.completed.get(position)
            if message is not None:
                row = acknowledged.get(
                    (
                        position,
                        message.get("stage", "curve"),
                        message["attempt"],
                    )
                )
                if (
                    row is None
                    or row.get("reason") != "complete"
                    or row["token"] != message["token"]
                ):
                    raise ValueError("result has no paid acknowledgement")
            if position in self.boundaries:
                row = acknowledged.get(
                    (position, "stage_one", self.attempts[position])
                )
                if (
                    row is None
                    or row.get("reason") != "boundary"
                    or row.get("boundary_sha256")
                    != self.boundaries[position]["sha256"]
                ):
                    raise ValueError("boundary has no paid acknowledgement")
        try:
            total_cpu = math.fsum(charged_cpu)
        except OverflowError as error:
            raise ValueError("unbounded lease CPU total") from error
        if saved["work"] < charged_work or saved["cpu"] + 1e-9 < total_cpu:
            raise ValueError("resume would reset retained lease charges")

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
        if self.manifest.execution is not None:
            payload["execution"] = self.manifest.execution
        if len(canonical(payload)) > self.settings.checkpoint_bytes:
            raise MemoryError("checkpoint encoding exceeds cap")
        return dict(payload=payload, sha256=digest(payload))

    @classmethod
    def restore(cls, checkpoint, manifest, settings, *, fixed_work=False):
        if len(canonical(checkpoint)) > settings.checkpoint_bytes:
            raise ValueError("oversized checkpoint")
        data = checkpoint["payload"]
        if checkpoint["sha256"] != digest(data) or data["version"] != VERSION:
            raise ValueError("corrupt/unsupported checkpoint")
        if (
            data["manifest"] != manifest.identity
            or data["settings"] != asdict(settings)
            or data.get("execution") != manifest.execution
        ):
            raise ValueError(
                "resume identity/settings changed (including workers)"
            )
        result = cls(manifest, settings, fixed_work=fixed_work)
        if data["fixed_work"] != fixed_work:
            raise ValueError("resume traversal changed")
        result.runs = integer(data["runs"], "run count", 0, settings.max_runs)
        for key in ("attempts", "continuation_attempts"):
            values = data[key]
            if len(values) != len(result.attempts):
                raise ValueError("resume traversal changed")
            setattr(
                result,
                key,
                [
                    integer(a, "attempt", 0, settings.max_attempts)
                    for a in values
                ],
            )
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
                boundary,
                manifest.assignment(position),
                manifest.metadata,
                settings.residue_bytes,
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


def frame(message, limit):
    raw = canonical(message).encode("ascii")
    if len(raw) > limit:
        raise ValueError("oversized message")
    return struct.pack("!I", len(raw)) + raw


class Decoder:
    """Bounded incremental frames without blocking on partial headers."""

    def __init__(self, limit):
        self.limit, self.buffer = limit, bytearray()

    @property
    def complete(self):
        if len(self.buffer) < 4:
            return False
        size = struct.unpack("!I", self.buffer[:4])[0]
        return size == 0 or size > self.limit or len(self.buffer) >= size + 4

    def feed(self, data):
        self.buffer.extend(data)
        if len(self.buffer) > self.limit + 4:
            raise ValueError("transport capacity exceeded")
        if len(self.buffer) < 4:
            return None
        size = struct.unpack("!I", self.buffer[:4])[0]
        if size > self.limit or size == 0:
            raise ValueError("invalid frame size")
        if len(self.buffer) < size + 4:
            return None
        stop = size + 4
        result = json.loads(bytes(self.buffer[4:stop]))
        del self.buffer[: size + 4]
        if type(result) is not dict:
            raise ValueError("message must be an object")
        return result
