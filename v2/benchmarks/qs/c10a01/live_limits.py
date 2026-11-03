"""Durable cooperative grant extensions for a running benchmark worker."""

import fcntl
import json
import os
import time
from contextlib import contextmanager
from pathlib import Path

from ..c10a04.accounting import atomic
from . import plan_v4

STATE_BYTES = 256 * 1024
MAX_REVISIONS = 128
POLL_SECONDS = 0.2
EXTENSION_OVERHEAD = 10
CONTROL_HEADROOM = 16 * 1024


def allowances(calls, revision):
    return plan_v4.totals(
        list(calls), 12, 60, 30 + EXTENSION_OVERHEAD * revision
    )


def _read(output):
    path = Path(output) / "live.json"
    if path.stat().st_size > STATE_BYTES:
        raise RuntimeError("recorder_exhaustion: live state")
    state = json.loads(path.read_text())
    if state["schema"] != 1:
        raise ValueError("unsupported live grant schema")
    _validate(state)
    return state


def _validate(state):
    if not 0 <= state["revision"] <= MAX_REVISIONS:
        raise ValueError("invalid live revision capacity")
    replay = {
        name: (call["initial_seconds"], 0)
        for name, call in state["calls"].items()
    }
    if state["revision"] != len(state["journal"]):
        raise ValueError("live grant journal is incomplete")
    for revision, event in enumerate(state["journal"], 1):
        if event["revision"] != revision:
            raise ValueError("nonconsecutive live grant revision")
        for change in event["changes"]:
            previous, _ = replay[change["call"]]
            if (
                change["previous_seconds"] != previous
                or change["seconds"] <= previous
            ):
                raise ValueError("live grant journal changed")
            replay[change["call"]] = (change["seconds"], revision)
    for name, call in state["calls"].items():
        if type(call["seconds"]) is not int or call["seconds"] <= 0:
            raise ValueError("invalid finite live seconds")
        if (
            type(call["work"]) is not int
            or type(call["revision"]) is not int
            or type(call["acknowledged"]) is not int
        ):
            raise ValueError("exact integer grants and revisions required")
        if (call["seconds"], call["revision"]) != replay[name]:
            raise ValueError("live grant does not match its journal")
        if call["work"] != plan_v4.call_work(call["digits"], call["seconds"]):
            raise ValueError("live work allowance changed")
        if call["reservation_seconds"] != call["seconds"] + 5:
            raise ValueError("live time reservation changed")
        if not -1 <= call["acknowledged"] <= call["revision"]:
            raise ValueError("invalid live acknowledgement")
    if state["limits"] != allowances(
        state["calls"].values(), state["revision"]
    ):
        raise ValueError("live cumulative reservations changed")


def read(output):
    """Read one coherent atomic state; no mutation or acknowledgement."""
    return _read(output)


def _write(output, state):
    _validate(state)
    if len(json.dumps(state).encode()) > STATE_BYTES:
        raise RuntimeError("recorder_exhaustion: live state")
    atomic(Path(output) / "live.json", state)


@contextmanager
def transaction(output, identity):
    """Serialize worker, supervisor and explicit controller mutations."""
    output = Path(output)
    with (output / "live.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        state = _read(output)
        if state["identity"] != identity:
            raise ValueError("live source/release identity changed")
        yield state
        _write(output, state)


def initialize(output, identity, calls):
    """Reserve the entire finite stage before any worker starts."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "live.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if (output / "live.json").exists():
            raise RuntimeError("live stage already started; no implicit retry")
        rows = {}
        for call in calls:
            row = dict(call)
            row.update(
                status="pending",
                revision=0,
                acknowledged=-1,
                started=None,
                pid=None,
                used_at_ack=0,
                receipts={},
                initial_seconds=call["seconds"],
            )
            rows[row["id"]] = row
        if len(rows) != len(calls):
            raise ValueError("duplicate live call identity")
        state = dict(
            schema=1,
            identity=identity,
            revision=0,
            journal=[],
            calls=rows,
            limits=allowances(rows.values(), 0),
            status="running",
        )
        _write(output, state)
        return state


def extend(output, identity, names, seconds, expected_revision, reason):
    """Reserve a larger finite grant; a worker acknowledgement is separate."""
    if (
        type(seconds) is not int
        or not 0 < seconds < 2**63
        or type(expected_revision) is not int
    ):
        raise ValueError("finite integer seconds/revision required")
    if not isinstance(reason, str) or not 1 <= len(reason) <= 512:
        raise ValueError("a bounded explicit extension reason is required")
    if not names or len(names) != len(set(names)):
        raise ValueError("select distinct calls")
    with transaction(output, identity) as state:
        if state["status"] != "running":
            raise ValueError("stage is closed")
        if expected_revision != state["revision"]:
            raise ValueError("stale live grant revision")
        if state["revision"] >= MAX_REVISIONS:
            raise RuntimeError("live extension journal exhausted")
        now = time.monotonic()
        changed = []
        for name in names:
            call = state["calls"][name]
            if call["status"] not in ("pending", "warming", "measured"):
                raise ValueError("call no longer accepts extensions: " + name)
            if seconds <= call["seconds"]:
                raise ValueError("extensions must strictly increase limits")
            if call.get("extensions_closed"):
                raise ValueError("worker has closed extension admission")
            if call["started"] is not None:
                if now >= call["started"] + call["seconds"] - 0.02:
                    raise ValueError("old deadline passed; cannot resurrect")
            changed.append(
                dict(
                    call=name,
                    previous_seconds=call["seconds"],
                    seconds=seconds,
                    previous_work=call["work"],
                    work=plan_v4.call_work(call["digits"], seconds),
                )
            )
        state["revision"] += 1
        for change in changed:
            call = state["calls"][change["call"]]
            call.update(
                seconds=seconds,
                reservation_seconds=seconds + 5,
                work=change["work"],
                revision=state["revision"],
            )
        state["limits"] = allowances(
            state["calls"].values(), state["revision"]
        )
        state["journal"].append(
            dict(
                revision=state["revision"],
                requested_monotonic=now,
                reason=reason,
                changes=changed,
                stage_limits=state["limits"]["total"],
            )
        )
        # Leave room for acknowledgements, receipt references and failure
        # metadata; a full controller journal must not break the live worker.
        if len(json.dumps(state).encode()) > STATE_BYTES - CONTROL_HEADROOM:
            raise RuntimeError("live extension state capacity exhausted")
        return dict(revision=state["revision"], changes=changed)


class LiveBudget:
    """Apply durable grants before Budget checks in the worker thread."""

    def __init__(self, output, identity, name, budget):
        self.output = Path(output)
        self.identity = identity
        self.name = name
        self.budget = budget
        self.previous_cancel = budget.cancelled
        self.next_poll = 0.0
        self.signature = None
        self.closed = False
        with transaction(output, identity) as state:
            call = state["calls"][name]
            if call["status"] not in ("pending", "warming"):
                raise ValueError("measured call already started")
            call.update(
                status="measured",
                started=budget._wall_start - budget.prior_wall,
                pid=os.getpid(),
            )
            self._apply(call)
        budget.cancelled = self.poll

    def _apply(self, call):
        # The existing ledger, RNG, roots and clock origins stay in place.
        self.budget.seconds = call["seconds"]
        self.budget.cpu_seconds = call["seconds"]
        self.budget.work_limit = call["work"]
        call["acknowledged"] = call["revision"]
        call["used_at_ack"] = self.budget.used
        call["wall_at_ack"] = self.budget.wall_used
        call["cpu_at_ack"] = self.budget.cpu_used

    def poll(self):
        now = time.monotonic()
        near_time = (
            self.budget.wall_used >= self.budget.seconds - 0.02
            or self.budget.cpu_used >= self.budget.cpu_seconds - 0.02
        )
        near_limit = (
            near_time or self.budget.used >= self.budget.work_limit * 99 // 100
        )
        if near_time and not self.closed:
            # Serialize final admission with a controller's pending fsync.
            with transaction(self.output, self.identity) as state:
                call = state["calls"][self.name]
                if call["revision"] != call["acknowledged"]:
                    self._apply(call)
                if (
                    self.budget.wall_used >= self.budget.seconds - 0.02
                    or self.budget.cpu_used >= self.budget.cpu_seconds - 0.02
                ):
                    call["extensions_closed"] = True
                    self.closed = True
        if not self.closed and (now >= self.next_poll or near_limit):
            self.next_poll = now + POLL_SECONDS
            path = self.output / "live.json"
            stat = path.stat()
            signature = (stat.st_mtime_ns, stat.st_size, stat.st_ino)
            if signature != self.signature or near_limit:
                state = _read(self.output)
                if state["identity"] != self.identity:
                    raise ValueError("live release identity changed")
                call = state["calls"][self.name]
                if call["revision"] != call["acknowledged"]:
                    with transaction(self.output, self.identity) as state:
                        call = state["calls"][self.name]
                        if call["status"] != "measured":
                            raise RuntimeError("live call closed during work")
                        self._apply(call)
                self.signature = signature
        return bool(self.previous_cancel and self.previous_cancel())


def mark_phase(output, identity, name, phase):
    with transaction(output, identity) as state:
        call = state["calls"][name]
        if phase == "warming" and call["status"] != "pending":
            raise ValueError("warmup cannot be replayed")
        call["status"] = phase


def close_measurement(
    output, identity, name, receipt_name, receipt_sha, result
):
    """Close extension admission without resetting or refunding any charge."""
    with transaction(output, identity) as state:
        call = state["calls"][name]
        if call["status"] != "measured":
            raise ValueError("measurement is not live")
        if result["work"] > call["work"]:
            raise ValueError("reported work exceeds the live reservation")
        call["status"] = "measured_complete"
        call["receipts"]["measured"] = dict(
            name=receipt_name,
            sha256=receipt_sha,
            wall=max(call["seconds"] + 5, result["seconds"]),
            cpu=max(call["seconds"] + 5, result["cpu_seconds"]),
            work=result["work"],
        )
        call["reason"] = result["reason"]
        call["applied_revision"] = call["acknowledged"]


def watchdog_deadline(state, name, fallback_started):
    """A reserved pending extension also protects the parent's watchdog."""
    call = state["calls"][name]
    if call["started"] is None:
        return fallback_started + call["seconds"] + 80
    return call["started"] + call["seconds"] + 80


def check_capacity(state):
    """Include settled charges and still-admitted calls before launching."""
    measured = dict(wall=0, cpu=0, work=0)
    warm = dict(wall=0, cpu=0, work=0)
    for call in state["calls"].values():
        value = call["receipts"].get("measured")
        if value is None:
            value = dict(
                wall=call["seconds"] + 5,
                cpu=call["seconds"] + 5,
                work=call["work"],
            )
        for field in measured:
            measured[field] += value[field]
        for phase in ("before", "after"):
            value = call["receipts"].get(phase)
            if value is None and call["status"] != "complete":
                value = dict(wall=35, cpu=35, work=290_000_000_000)
            if value is not None:
                for field in warm:
                    warm[field] += value[field]
    for name, used in (("measured", measured), ("warmup", warm)):
        for field, value in used.items():
            if value > state["limits"]["phases"][name][field]:
                raise RuntimeError("live phase capacity exceeded: " + name)


def interrupt(output, identity, name, detail):
    with transaction(output, identity) as state:
        call = state["calls"][name]
        call["status"] = "interrupted"
        call["interruption"] = detail[:512]
        call["reserved_interrupted"] = dict(
            wall=call["seconds"] + 5,
            cpu=call["seconds"] + 5,
            work=call["work"],
        )
        state["status"] = "interrupted"
