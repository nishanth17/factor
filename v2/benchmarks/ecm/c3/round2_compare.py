"""Whole matched groups for the frozen C3 round-two comparison."""

import argparse
import contextlib
import json
import os
import signal
import subprocess
import time
import uuid
from unittest.mock import patch

from ..p52.p52_realistic import check_quiet
from . import c3_study as first
from . import round2_pilot as pilot
from . import round2_report as reporter
from . import round2_train as training
from .round2_corpus import load_training

CONTROL = first.INPUTS / "controls/c3_round2_comparison_v2.json"
ARMS = ("control", "fixed32", "fitted")
PHASE_SECONDS = 3600
WARM_GRANT = 60
REPORT_GRANT = 180


def atomic_save(path, value):
    """Replace the ledger atomically so interruption retains its prepayment."""
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


@contextlib.contextmanager
def deadline(seconds):
    """Bound the whole grant, including validation and capture overhead."""

    def expired(_signum, _frame):
        raise TimeoutError("comparison grant exhausted")

    previous = signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def verify_control():
    """Require unchanged, committed sources and the committed fitted table."""
    training.verify_control()
    protocol = json.loads(CONTROL.read_text())
    for name, digest in protocol["sha256"].items():
        if first.digest(first.ROOT / name) != digest:
            raise ValueError("comparison freeze changed: " + name)
    for path in (CONTROL, *(first.ROOT / n for n in protocol["sha256"])):
        relative = str(path.relative_to(first.ROOT))
        committed = subprocess.check_output(
            ["git", "show", "HEAD:" + relative], cwd=first.ROOT
        )
        if committed != path.read_bytes():
            raise ValueError(
                "comparison inputs must be committed: " + relative
            )
    assessment = json.loads(
        (first.INPUTS / "controls/c3_round2_assessment_v1.json").read_text()
    )
    for name, digest in assessment["sha256"].items():
        if first.digest(first.ROOT / name) != digest:
            raise ValueError("assessment freeze changed: " + name)
    selection = json.loads(training.SELECTION.read_text())
    if selection["protocol_sha256"] != first.digest(training.CONTROL):
        raise ValueError("selection belongs to another training freeze")
    return selection["bands"]


class Allowance:
    """A lease shares cumulative expense and retains interrupted grants."""

    def __init__(self, path, seconds, *, phase_reserve=0):
        self.path = path
        identity = first.digest(CONTROL)
        self.ledger = (
            json.loads(path.read_text())
            if path.exists()
            else dict(
                protocol_sha256=identity,
                active_seconds=0,
                leases=[],
                warmups=[],
            )
        )
        if self.ledger["protocol_sha256"] != identity:
            raise ValueError("continued comparison identity changed")
        if self.ledger.get("failed"):
            raise ValueError("failed comparison cannot resume or promote")
        pending = self.ledger.pop("pending", None)
        if pending and pending["kind"] != "warmup":
            raise ValueError(
                "interrupted matched group/report is inconclusive"
            )
        self.spent = self.ledger["active_seconds"]
        self.started = time.monotonic()
        self.limit = min(seconds, PHASE_SECONDS - self.spent - phase_reserve)
        self.prepaid = 0
        self.reserved_at = 0
        self.pending = None
        self.lease = uuid.uuid4().hex

    def elapsed(self):
        return time.monotonic() - self.started

    def fits(self, grant):
        # Leave two seconds for receipt/lock cleanup after the next full grant.
        return self.elapsed() + grant + 2 <= self.limit

    def save(self):
        charged = self.elapsed()
        if self.pending:
            charged = max(charged, self.reserved_at + self.prepaid)
        self.ledger["active_seconds"] = self.spent + charged
        if self.pending:
            self.ledger["pending"] = self.pending
        else:
            self.ledger.pop("pending", None)
        atomic_save(self.path, self.ledger)

    def reserve(self, kind, grant, identifier):
        if not self.fits(grant):
            return False
        self.prepaid = grant
        self.reserved_at = self.elapsed()
        self.pending = dict(kind=kind, identifier=identifier, lease=self.lease)
        self.save()
        return True

    def complete(self):
        self.prepaid = 0
        self.pending = None
        self.save()

    def fail(self, error):
        # Keep the full interrupted grant as conservative expense. Correctness
        # or capture failure cannot be discarded and retried into acceptance.
        self.ledger["failed"] = type(error).__name__ + ": " + str(error)
        self.save()

    def finish(self):
        self.ledger["leases"].append(
            dict(id=self.lease, wall=self.elapsed(), allowance=self.limit)
        )
        self.save()


def assignment(fixtures):
    """Preserve the original assignment and balanced arm-order rotation."""
    for sample, seed in enumerate(training.SEEDS):
        for index, fixture in enumerate(fixtures):
            offset = (sample + index) % len(ARMS)
            order = ARMS[offset:] + ARMS[:offset]
            if sample % 2:
                order = tuple(reversed(order))
            yield fixture, sample, seed, order


def group_name(fixture, sample):
    return f"{sample:02d}-{fixture['id']}"


def read_group(path, fixture, sample, seed):
    """A complete group has three immutable rows from the same lease."""
    if not path.exists():
        return None
    if not (path / "complete.json").exists():
        raise ValueError(
            "interrupted matched group is inconclusive: " + path.name
        )
    manifest = json.loads((path / "complete.json").read_text())
    begin = json.loads((path / "begin.json").read_text())
    expected = dict(
        id=fixture["id"],
        sample=sample,
        seed=seed,
        comparison_protocol_sha256=first.digest(CONTROL),
    )
    if any(begin.get(k) != v for k, v in expected.items()):
        raise ValueError("group assignment changed")
    if manifest["begin_sha256"] != first.digest(path / "begin.json"):
        raise ValueError("group admission record changed")
    if set(manifest["rows"]) != {a + ".json" for a in ARMS}:
        raise ValueError("group must contain exactly three arms")
    rows = []
    for arm in ARMS:
        name = arm + ".json"
        if first.digest(path / name) != manifest["rows"][name]:
            raise ValueError("matched row changed")
        row = json.loads((path / name).read_text())
        if (
            any(row.get(k) != v for k, v in expected.items())
            or row["arm"] != arm
            or row["lease"] != begin["lease"]
            or row["instrumented"]
            or row["protocol_sha256"] != first.digest(training.CONTROL)
        ):
            raise ValueError("incompatible matched row")
        rows.append(row)
    return rows


def capture_group(output, fixture, sample, seed, order, allowance, run):
    """Reserve all arms before the first one, never straddle process leases."""
    path = output / group_name(fixture, sample)
    if read_group(path, fixture, sample, seed) is not None:
        return True
    cap = 5 if first.band_of(fixture["n"]) == 30 else 30
    grant = len(ARMS) * (cap + 1)
    if not allowance.reserve("group", grant, path.name):
        return False
    path.mkdir()
    identity = dict(
        id=fixture["id"],
        sample=sample,
        seed=seed,
        lease=allowance.lease,
        comparison_protocol_sha256=first.digest(CONTROL),
    )
    first.save(path / "begin.json", dict(identity, order=order, grant=grant))
    with deadline(grant):
        for arm in order:
            check_quiet({os.getpid()})
            row = run(fixture, seed, arm)
            row.update(identity)
            row["protocol_sha256"] = first.digest(training.CONTROL)
            first.save(path / (arm + ".json"), row)
            print(
                sample,
                fixture["id"],
                arm,
                row["reason"],
                round(row["cpu"], 4),
                flush=True,
            )
        first.save(
            path / "complete.json",
            dict(
                begin_sha256=first.digest(path / "begin.json"),
                rows={
                    a + ".json": first.digest(path / (a + ".json"))
                    for a in ARMS
                },
            ),
        )
    allowance.complete()
    return True


def warmup(fixtures, allowance, run, control):
    """Warm both bands and force both historical/current relation paths."""
    for band in (30, 40):
        fixture = next(
            f for f in fixtures if f["id"] == f"r2_{band}_balanced_0"
        )
        for arm in (*ARMS, "old_no_ecm", "no_ecm"):
            if not allowance.reserve("warmup", WARM_GRANT, f"{band}:{arm}"):
                return False
            wall, cpu, count = time.monotonic(), time.process_time(), 0
            with deadline(WARM_GRANT):
                while (
                    min(time.monotonic() - wall, time.process_time() - cpu) < 3
                ):
                    if time.monotonic() - wall > 25:
                        raise TimeoutError("validated warmup did not finish")
                    check_quiet({os.getpid()})
                    if arm == "old_no_ecm":
                        with patch.dict(
                            pilot.POLICIES,
                            {"control": dict(tiers=(), mode=None)},
                        ):
                            row = training.run_one(
                                fixture, 17, "control", control
                            )
                    else:
                        row = run(fixture, 17, arm)
                    if not row["complete"]:
                        raise ValueError("warmup output incomplete")
                    count += 1
            allowance.ledger["warmups"].append(
                dict(
                    band=band,
                    arm=arm,
                    lease=allowance.lease,
                    wall=time.monotonic() - wall,
                    cpu=time.process_time() - cpu,
                    validated=count,
                )
            )
            allowance.complete()
    return True


def measure(output, lease_seconds):
    first.runtime()
    fitted = verify_control()
    if not 60 <= lease_seconds <= 1200:
        raise ValueError("measurement lease must be 60..1200 seconds")
    output.mkdir(parents=True, exist_ok=True)
    fixtures = load_training()
    with first.machine_window():
        allowance = Allowance(
            output / "ledger.json",
            lease_seconds,
            phase_reserve=REPORT_GRANT + 3,
        )
        try:
            control = first.baseline()

            def run(fixture, seed, arm):
                return training.run_one(
                    fixture, seed, arm, control, fitted=fitted
                )

            if not warmup(fixtures, allowance, run, control):
                return
            for fixture, sample, seed, order in assignment(fixtures):
                if not capture_group(
                    output, fixture, sample, seed, order, allowance, run
                ):
                    break
        except BaseException as error:
            allowance.fail(error)
            raise
        finally:
            allowance.finish()
    print(json.dumps(allowance.ledger), flush=True)


def report(output, destination):
    """Use the unchanged assessor under a charged, finite analysis grant."""
    first.runtime()
    verify_control()
    fixtures = load_training()
    with first.machine_window():
        allowance = Allowance(output / "ledger.json", REPORT_GRANT + 3)
        try:
            if not allowance.reserve("report", REPORT_GRANT, str(destination)):
                raise ValueError(
                    "comparison cap cannot fund frozen assessment"
                )
            with deadline(REPORT_GRANT):
                rows = []
                for fixture, sample, seed, _order in assignment(fixtures):
                    group = read_group(
                        output / group_name(fixture, sample),
                        fixture,
                        sample,
                        seed,
                    )
                    if group is None:
                        raise ValueError("incomplete comparison cannot select")
                    rows.extend(group)
                expected = {
                    group_name(f, s)
                    for f, s, _seed, _order in assignment(fixtures)
                }
                actual = {p.name for p in output.iterdir() if p.is_dir()}
                if actual != expected:
                    raise ValueError("unexpected group capture directory")
                groups = reporter.paired_groups(rows, fixtures)
                result = reporter.assess(groups)
                result["comparison_protocol_sha256"] = first.digest(CONTROL)
                result["complete_groups"] = len(expected)
                first.save(destination, result)
            allowance.complete()
        except BaseException as error:
            allowance.fail(error)
            raise
        finally:
            allowance.finish()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("measure", "report"))
    parser.add_argument("output", type=first.Path)
    parser.add_argument("--lease-seconds", type=int, default=600)
    parser.add_argument("--report", type=first.Path)
    args = parser.parse_args()
    if args.mode == "measure":
        measure(args.output, args.lease_seconds)
    else:
        if args.report is None:
            parser.error("report mode requires --report")
        report(args.output, args.report)


if __name__ == "__main__":
    main()
