"""Durable reservations and idempotent immutable result receipts."""

import fcntl
import hashlib
import json
import math
import os
import resource
import time
from contextlib import contextmanager
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic(path, value):
    """Publish one fsynced JSON object and persist the directory rename."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)
    descriptor = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def receipt(path, value):
    """Accept identical republication and refuse conflicting contents."""
    if path.exists():
        if json.loads(path.read_text()) != value:
            raise ValueError("receipt identity or content changed")
        return
    atomic(path, value)


@contextmanager
def output_lock(output):
    """Serialize all ledger mutations, including generation and selection."""
    output.mkdir(parents=True, exist_ok=True)
    with (output / "study.lock").open("a+") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


class Ledger:
    """Reserve work and time before execution; retain interrupted charges.

    A receipt is published before settling its reservation. Recovery can settle
    that narrow crash window from the bound receipt, without redoing the call.
    Missing measured receipts become censored cells. Diagnostic retries keep
    every earlier reservation charged.
    """

    def __init__(self, output, identity, limits):
        self.output = Path(output)
        self.path = self.output / "ledger.json"
        self.limits = limits
        self.data = (
            json.loads(self.path.read_text())
            if self.path.exists()
            else dict(identity=identity, limits=limits, calls=[])
        )
        if self.data["identity"] != identity or self.data["limits"] != limits:
            raise ValueError("ledger belongs to another freeze")
        self.recover()

    def used(self, phase=None, field="wall"):
        return sum(
            row[field]
            for row in self.data["calls"]
            if phase is None or row["phase"] == phase
        )

    def room(self, phase, seconds, work):
        return all(
            self.used(None, field) + amount <= self.limits["total"][field]
            and self.used(phase, field) + amount <= self.limits[phase][field]
            for field, amount in (
                ("wall", seconds),
                ("cpu", seconds),
                ("work", work),
            )
        )

    def write(self):
        atomic(self.path, self.data)

    def recover(self):
        changed = False
        for row in self.data["calls"]:
            path = self.output / row["receipt"]
            identity = row.get("identity", {})
            if row["status"] == "complete":
                if not path.exists() or row["sha256"] != digest(path):
                    raise ValueError("settled receipt changed or disappeared")
                continue
            if row["status"] == "failed" and "sample" in identity:
                raise RuntimeError(
                    "failed measurement needs a correctness audit"
                )
            if row["status"] not in ("reserved", "interrupted"):
                continue
            if path.exists():
                result = json.loads(path.read_text())
                if (
                    result["reservation"] != row["index"]
                    or result["identity"] != identity
                ):
                    raise ValueError("receipt reservation mismatch")
                self.settle(row, result)
                row["sha256"] = digest(path)
            elif "sample" in identity:
                # Lost measured calls remain censored. Retrying only
                # successful cells would bias both coverage and cost selection.
                result = dict(
                    id=identity["fixture"],
                    digits=identity["digits"],
                    sample=identity["sample"],
                    seed=identity["seed"],
                    arm=identity["arm"],
                    complete=False,
                    divisor=None,
                    factors=[],
                    certainty=[],
                    remaining=[identity["n"]],
                    reason="interrupted",
                    seconds=row["wall"],
                    cpu_seconds=row["cpu"],
                    work=row["work"],
                    stats={},
                    cap_seconds=identity["cap"],
                    capped_seconds=row["wall"],
                    configuration=None,
                    base={},
                    measurement_process=None,
                )
                wrapper = dict(
                    identity=identity,
                    reservation=row["index"],
                    result=result,
                    charge_wall=row["wall"],
                    charge_cpu=row["cpu"],
                    charge_work=row["work"],
                )
                receipt(path, wrapper)
                self.settle(row, wrapper)
                row["sha256"] = digest(path)
            else:
                row["status"] = "interrupted"
            changed = True
        if changed:
            self.write()

    @staticmethod
    def settle(row, result):
        for field in ("wall", "cpu", "work"):
            value = result["charge_" + field]
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0
            ):
                raise ValueError("invalid receipt resource charge")
        row.update(
            status="complete",
            wall=result["charge_wall"],
            cpu=result["charge_cpu"],
            work=result["charge_work"],
        )

    def call(self, phase, name, seconds, work, identity, function):
        path = self.output / (name + ".json")
        if path.exists():
            result = json.loads(path.read_text())
            if result["identity"] != identity:
                raise ValueError("saved call identity changed")
            index = result["reservation"]
            row = self.data["calls"][index]
            if row["receipt"] != path.name or row["status"] != "complete":
                raise ValueError("unaccounted saved receipt")
            if row["sha256"] != digest(path):
                raise ValueError("saved receipt changed")
            return result["result"]
        disk = sum(
            p.stat().st_size for p in self.output.rglob("*") if p.is_file()
        )
        if disk + 2**20 > self.limits["disk_bytes"]:
            raise RuntimeError("study disk allowance exhausted")
        if not self.room(phase, seconds, work):
            raise RuntimeError(
                "aggregate allowance cannot cover next reservation"
            )

        row = dict(
            index=len(self.data["calls"]),
            phase=phase,
            receipt=path.name,
            identity=identity,
            status="reserved",
            wall=seconds,
            cpu=seconds,
            work=work,
        )
        self.data["calls"].append(row)
        self.write()
        started, cpu = time.monotonic(), time.process_time()
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        child_before = children.ru_utime + children.ru_stime
        try:
            result = function()
            children = resource.getrusage(resource.RUSAGE_CHILDREN)
            spent_cpu = (
                time.process_time()
                - cpu
                + children.ru_utime
                + children.ru_stime
                - child_before
            )
            wrapper = dict(
                identity=identity,
                reservation=row["index"],
                result=result,
                actual_wall=time.monotonic() - started,
                actual_cpu=spent_cpu,
                charge_wall=max(seconds, time.monotonic() - started),
                charge_cpu=max(seconds, spent_cpu),
                charge_work=result.get("work", 0),
            )
            receipt(path, wrapper)
            self.settle(row, wrapper)
            row["sha256"] = digest(path)
            self.write()
            return result
        except BaseException as error:
            children = resource.getrusage(resource.RUSAGE_CHILDREN)
            spent_cpu = (
                time.process_time()
                - cpu
                + children.ru_utime
                + children.ru_stime
                - child_before
            )
            # No outcome is accepted on a failed call. Keep its whole reserved
            # work/time plus any cooperative overshoot; never refund on retry.
            row.update(
                status="interrupted"
                if isinstance(error, (KeyboardInterrupt, SystemExit))
                else "failed",
                error=type(error).__name__ + ": " + str(error),
                wall=max(seconds, time.monotonic() - started),
                cpu=max(seconds, spent_cpu),
            )
            self.write()
            raise
