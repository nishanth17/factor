"""Pure B12 identities, framing and leases; no factorer or process imports."""

import hashlib
import json
import math
import random
import struct
from dataclasses import asdict, dataclass

VERSION = "b12-f4-v1"


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
    """Owned reservations, distinct from observed RSS or calibrated costs."""

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

    def __post_init__(self):
        if self.workers not in (1, 2, 4) or type(self.workers) is not int:
            raise ValueError("worker count must be 1, 2 or 4")
        integer(self.max_attempts, "max_attempts", 1, 3)
        integer(self.max_runs, "max_runs", 1, 32)
        for name in ("assignment_work", "startup_work"):
            integer(getattr(self, name), name, 1)
        for name in (
            "assignment_cpu",
            "startup_seconds",
            "drain_seconds",
            "kill_seconds",
        ):
            finite(getattr(self, name), name, 0.001)
        integer(self.frame_bytes, "frame_bytes", 1024, 65_536)
        integer(self.checkpoint_bytes, "checkpoint_bytes", 4096, 8 * 2**20)
        integer(self.worker_memory, "worker_memory", 2**20, 2**30)
        integer(self.parent_memory, "parent_memory", 2**20, 4 * 2**30)
        integer(self.memory_bytes, "memory_bytes", 2**20, 8 * 2**30)
        if self.assignment_policy not in ("dynamic", "static"):
            raise ValueError("unknown curve assignment policy")
        if self.owned_bytes > self.memory_bytes:
            raise MemoryError(
                "simultaneous worker/transport/checkpoint ownership refused"
            )

    @property
    def shutdown_cpu(self):
        # Retain this grant in full: interpreter teardown has no final report.
        return self.drain_seconds + 2 * self.kill_seconds + 1

    @property
    def owned_bytes(self):
        return (
            self.parent_memory
            + (self.workers + 1) * self.worker_memory
            + self.workers * 8 * self.frame_bytes
            + 3 * self.checkpoint_bytes
        )


@dataclass(frozen=True)
class EngagementRule:
    """Missing calibration stays serial; estimates are observable metadata."""

    schedule: str = ""
    cold_seconds: float | None = None
    warm_seconds: float | None = None
    seconds_per_curve: tuple = ()

    def __post_init__(self):
        # Invalid observations disable engagement; they must not enter hashes
        # as NaN/Infinity or silently become an optimistic cost estimate.
        for name in ("cold_seconds", "warm_seconds"):
            try:
                finite(getattr(self, name), name, 0.001)
            except (ValueError, TypeError):
                object.__setattr__(self, name, None)
        try:
            costs = tuple(self.seconds_per_curve)
            if len(costs) > 32:
                raise ValueError("too many cost rows")
            for cost in costs:
                finite(cost, "curve cost", 0.000001)
        except (ValueError, TypeError):
            costs = ()
        object.__setattr__(self, "seconds_per_curve", costs)
        if not isinstance(self.schedule, str) or len(self.schedule) > 64:
            object.__setattr__(self, "schedule", "")

    def engages(self, manifest, *, ready, remaining=None):
        threshold = self.warm_seconds if ready else self.cold_seconds
        try:
            finite(threshold, "threshold", 0.001)
            if self.schedule != manifest.schedule:
                return False
            if len(self.seconds_per_curve) != len(manifest.tiers):
                return False
            counts = (
                tuple(t[2] for t in manifest.tiers)
                if remaining is None
                else tuple(remaining)
            )
            if len(counts) != len(manifest.tiers):
                return False
            estimate = sum(
                finite(cost, "curve estimate", 0.000001)
                * integer(count, "remaining curves", 0, tier[2])
                for cost, count, tier in zip(
                    self.seconds_per_curve, counts, manifest.tiers
                )
            )
            finite(estimate, "total estimate")
        except (ValueError, TypeError, OverflowError):
            return False
        return estimate >= threshold


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
        self.schedule = digest([VERSION, options, n.bit_length()])
        self.identity = digest([VERSION, n, seed, lineage, options])
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
        tier, index = self.jobs[position]
        b1, b2, _ = self.tiers[tier]
        identity = digest([self.identity, tier, index])
        seed = int(digest([VERSION, "curve", identity]), 16)
        return dict(
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


class Search:
    """Parent-owned assignment and acknowledgement state."""

    def __init__(self, manifest, settings, *, fixed_work=False):
        if type(fixed_work) is not bool:
            raise ValueError("fixed_work must be Boolean")
        self.manifest, self.settings = manifest, settings
        if 65536 + len(manifest.jobs) * 4096 > settings.checkpoint_bytes:
            raise MemoryError("completed-result checkpoint ownership refused")
        self.fixed_work = fixed_work
        self.attempts = [0] * len(manifest.jobs)
        self.completed = {}
        self.active = {}
        self.factor = None
        self.retired = False
        self.runs = 0

    def next(self, slot):
        """Tier barrier prevents early dispatch at a higher numerical bound."""
        remaining = [
            i for i in range(len(self.attempts)) if i not in self.completed
        ]
        if not remaining or self.retired:
            return None
        tier = self.manifest.jobs[remaining[0]][0]
        for position in remaining:
            if self.manifest.jobs[position][0] != tier:
                break
            if position in self.active:
                continue
            if (
                self.settings.assignment_policy == "static"
                and position % self.settings.workers != slot
            ):
                continue
            if self.attempts[position] >= self.settings.max_attempts:
                raise LeaseRefusedError("attempt_limit")
            return position
        return None

    def dispatch(self, position, token):
        integer(position, "dispatch position", 0, len(self.attempts) - 1)
        integer(token, "lease token", 1)
        if position in self.completed or position in self.active:
            raise ValueError("duplicate dispatch")
        self.attempts[position] += 1
        self.active[position] = (token, self.attempts[position])

    def receive(self, position, token, message, *, publish=True):
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
            "chain_chunks",
            "preparations",
            "kind",
            "cpu",
            "rss",
        }
        if not message.keys() <= allowed or len(canonical(message)) > 2048:
            raise ValueError("unbounded result record")
        integer(
            message["attempt"], "result attempt", 1, self.settings.max_attempts
        )
        integer(message["token"], "result token", 1)
        for key in ("work", "rss", "chain_chunks", "preparations"):
            if key in message:
                integer(message[key], key)
        if "cpu" in message:
            finite(message["cpu"], "worker CPU")
        if "phase" in message and (
            not isinstance(message["phase"], str) or len(message["phase"]) > 64
        ):
            raise ValueError("unbounded phase label")
        if message.get("failure") is not None and (
            not isinstance(message["failure"], str)
            or len(message["failure"]) > 256
        ):
            raise ValueError("unbounded failure description")
        assignment = self.manifest.assignment(position)
        if self.active.get(position) != (token, message.get("attempt")):
            raise ValueError("duplicate or stale result")
        if (
            message.get("id") != assignment["id"]
            or message.get("token") != token
        ):
            raise ValueError("result belongs to another assignment")
        factor = message.get("factor")
        if factor is not None and (
            type(factor) is not int
            or not 1 < factor < self.manifest.n
            or self.manifest.n % factor
        ):
            raise ValueError("worker returned an invalid proper divisor")
        if message.get("reason") not in (
            "complete",
            "cancelled",
            "wall_limit",
            "cpu_limit",
            "work_limit",
            "error",
        ):
            raise ValueError("invalid worker outcome")
        if message["reason"] != "complete" and factor is not None:
            raise ValueError("unfinished result cannot publish a factor")
        del self.active[position]
        if message["reason"] == "complete":
            self.completed[position] = dict(message)
            if factor is not None and publish and not self.retired:
                if self.factor is None:
                    self.factor = factor
                if not self.fixed_work:
                    self.retired = True

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
            completed=self.completed,
            factor=self.factor,
            retired=self.retired,
        )
        raw = canonical(payload)
        if len(raw) > self.settings.checkpoint_bytes:
            raise MemoryError("checkpoint encoding exceeds cap")
        return dict(payload=payload, sha256=digest(payload))

    @classmethod
    def restore(cls, checkpoint, manifest, settings, *, fixed_work=False):
        if len(canonical(checkpoint)) > settings.checkpoint_bytes:
            raise ValueError("oversized checkpoint")
        data = checkpoint["payload"]
        if checkpoint["sha256"] != digest(data) or data["version"] != VERSION:
            raise ValueError("corrupt checkpoint")
        if data["manifest"] != manifest.identity or data["settings"] != asdict(
            settings
        ):
            raise ValueError(
                "resume identity/settings changed (including workers)"
            )
        result = cls(manifest, settings, fixed_work=fixed_work)
        if data["fixed_work"] != fixed_work or len(data["attempts"]) != len(
            result.attempts
        ):
            raise ValueError("resume traversal changed")
        result.runs = integer(data["runs"], "run count", 0, settings.max_runs)
        result.attempts = [
            integer(a, "attempt", 0, settings.max_attempts)
            for a in data["attempts"]
        ]
        for key, message in data["completed"].items():
            position = integer(
                int(key), "completed position", 0, len(result.attempts) - 1
            )
            if position in result.completed or message["reason"] != "complete":
                raise ValueError("duplicate/nonfinal checkpoint result")
            if not 1 <= message["attempt"] <= result.attempts[position]:
                raise ValueError("unissued checkpoint result")
            result.active[position] = (message["token"], message["attempt"])
            result.receive(position, message["token"], message, publish=False)
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
        result = json.loads(bytes(self.buffer[4:]))
        del self.buffer[: size + 4]
        if type(result) is not dict:
            raise ValueError("message must be an object")
        return result
