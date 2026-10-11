"""Bounded, matched transfer screening of upstream ECM schedules."""

import argparse
import json
import os
import time
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import asdict
from statistics import mean
from unittest.mock import patch

from ....execution.ecm_presets import with_ecm_preset
from ..p52.p52_realistic import check_quiet
from . import c3_study as first
from . import round2_compare as matched
from . import round2_pilot as pilot
from . import round2_report as reporter
from . import round2_train as training
from .round2_corpus import load_training

CONTROL = first.INPUTS / "controls/c3_source_presets_v1.json"
PRESETS = (
    "yamaquasi_auto160",
    "yamaquasi_ecm64",
    "alpertron20",
    "gmp_ecm20",
)
ARMS = ("fixed32", "fitted", *PRESETS)
SEEDS = (
    17,
    43,
    89,
    113,
    151,
    181,
    223,
    269,
    307,
    347,
    389,
    431,
    479,
    523,
    569,
    617,
    661,
    709,
    751,
    797,
    839,
    887,
    929,
    977,
    1021,
    1069,
    1117,
)
CASE_IDS = tuple(
    f"r2_{band}_{kind}_0"
    for band, kinds in (
        (30, ("balanced", "uneven_12", "uneven_14")),
        (40, ("balanced", "uneven_12", "uneven_16")),
    )
    for kind in kinds
)


@contextmanager
def layout():
    """Reuse the frozen whole-group ledger in this isolated study process.

    Only benchmark identity and arm enumeration change inside this context;
    production sources, fitted tiers and the earlier comparison stay intact.
    Restoring the globals also makes imported helpers safe for ordinary tests.
    """
    with patch.multiple(matched, CONTROL=CONTROL, ARMS=ARMS):
        yield


def verify_control():
    """Require committed numerical presets and all measured dependencies."""
    matched.verify_control()
    protocol = json.loads(CONTROL.read_text())
    for name, digest in protocol["sha256"].items():
        if first.digest(first.ROOT / name) != digest:
            raise ValueError("source-preset freeze changed: " + name)
    for path in (CONTROL, *(first.ROOT / n for n in protocol["sha256"])):
        relative = str(path.relative_to(first.ROOT))
        committed = matched.subprocess.check_output(
            ["git", "show", "HEAD:" + relative], cwd=first.ROOT
        )
        if committed != path.read_bytes():
            raise ValueError("uncommitted source-preset input: " + relative)
    fixtures = {f["id"]: f for f in load_training()}
    return [fixtures[name] for name in CASE_IDS]


def run_one(fixture, seed, arm):
    """Apply the public helper to the identical protected C1 configuration."""
    original = pilot.configuration
    observed = []

    def configured(n, selected, engine):
        if selected in PRESETS:
            base = original(n, "protected32", engine)
            config = with_ecm_preset(base, selected)
        else:
            config = original(n, selected, engine)
        observed.append(config)
        return config

    with patch.object(pilot, "configuration", configured):
        if arm in PRESETS:
            row = pilot.run_one(fixture, seed, arm, None)
        else:
            fitted = json.loads(training.SELECTION.read_text())["bands"]
            row = training.run_one(fixture, seed, arm, None, fitted=fitted)
    if len(observed) != 1:
        raise ValueError("expected one root configuration")
    row["siqs_config"] = asdict(observed[0].siqs)
    return row


def assignment(fixtures, samples):
    for sample, seed in enumerate(SEEDS[:samples]):
        for index, fixture in enumerate(fixtures):
            offset = (sample + index) % len(ARMS)
            order = ARMS[offset:] + ARMS[:offset]
            if sample % 2:
                order = tuple(reversed(order))
            yield fixture, sample, seed, order


def target_samples(output):
    """Only the predeclared uncertainty rule can authorize extra samples."""
    ledger_path = output / "ledger.json"
    ledger = (
        json.loads(ledger_path.read_text()) if ledger_path.exists() else {}
    )
    if ledger and ledger["protocol_sha256"] != first.digest(CONTROL):
        raise ValueError("source-preset ledger identity changed")
    if (
        ledger.get("failed")
        or (ledger.get("pending") or {}).get("kind") == "report"
    ):
        raise ValueError("failed/interrupted assessment cannot resume")
    saved = ledger.get("assessments", {})
    paths = {p.name for p in output.glob("assessment-*.json")}
    if paths != {f"assessment-{int(count):02d}.json" for count in saved}:
        raise ValueError("unrecorded or missing preset assessment")
    expected = 9
    for count in (9, 18, 27):
        if str(count) not in saved:
            if any(int(c) >= count for c in saved):
                raise ValueError("preset assessments skipped a sample block")
            return expected
        if expected != count:
            raise ValueError("preset screening continued after its stop")
        path = output / f"assessment-{count:02d}.json"
        if first.digest(path) != saved[str(count)]:
            raise ValueError("preset assessment changed")
        assessment = json.loads(path.read_text())
        if (
            assessment["protocol_sha256"] != first.digest(CONTROL)
            or assessment["samples"] != count
        ):
            raise ValueError("assessment assignment changed")
        expected = assessment["next_samples"]
        if expected not in (None, count + 9) or (
            count == 27 and expected is not None
        ):
            raise ValueError("undeclared preset extension")
    return expected


def warmup(fixtures, allowance):
    """Warm each transferred tier schedule and force SIQS in both bands."""
    for band in (30, 40):
        fixture = next(
            f for f in fixtures if f["id"] == f"r2_{band}_balanced_0"
        )
        for arm in (*ARMS, "no_ecm"):
            if not allowance.reserve("warmup", 60, f"{band}:{arm}"):
                return False
            wall, cpu, count = time.monotonic(), time.process_time(), 0
            with matched.deadline(60):
                while (
                    min(time.monotonic() - wall, time.process_time() - cpu) < 3
                ):
                    if time.monotonic() - wall > 25:
                        raise TimeoutError(
                            "validated preset warmup did not finish"
                        )
                    check_quiet({os.getpid()})
                    row = run_one(fixture, 17, arm)
                    if not row["complete"]:
                        raise ValueError("preset warmup incomplete")
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
    fixtures = verify_control()
    if not 60 <= lease_seconds <= 1200:
        raise ValueError("measurement lease must be 60..1200 seconds")
    output.mkdir(parents=True, exist_ok=True)
    samples = target_samples(output)
    if samples is None:
        raise ValueError("source-preset screening has already stopped")
    with first.machine_window(), layout():
        allowance = matched.Allowance(
            output / "ledger.json", lease_seconds, phase_reserve=183
        )
        try:
            if not warmup(fixtures, allowance):
                return
            for fixture, sample, seed, order in assignment(fixtures, samples):
                if not matched.capture_group(
                    output, fixture, sample, seed, order, allowance, run_one
                ):
                    break
        except BaseException as error:
            allowance.fail(error)
            raise
        finally:
            allowance.finish()
    print(json.dumps(allowance.ledger), flush=True)


def assess(groups, samples, repetitions=10_000):
    """Screen complete costs; six inputs cannot establish a broad default."""
    results = {}
    for candidate in PRESETS:
        comparisons = {}
        for control in ("fixed32", "fitted"):
            paired = {
                key: {
                    subject: {
                        seed: {
                            "control": rows[control],
                            "fitted": rows[candidate],
                        }
                        for seed, rows in seeds.items()
                    }
                    for subject, seeds in subjects.items()
                }
                for key, subjects in groups.items()
            }
            comparisons[control] = reporter.improvement_interval(
                paired, "control", repetitions=repetitions
            )
        completion = all(
            mean(
                rows[candidate]["complete"]
                for seeds in subjects.values()
                for rows in seeds.values()
            )
            >= mean(
                rows[control]["complete"]
                for seeds in subjects.values()
                for rows in seeds.values()
            )
            - 0.05
            for subjects in groups.values()
            for control in ("fixed32", "fitted")
        )
        wall = {
            control: reporter.weighted_cost(groups, candidate, "wall")
            / reporter.weighted_cost(groups, control, "wall")
            for control in ("fixed32", "fitted")
        }
        positive_point = all(
            c["relative_improvement"] > 0 for c in comparisons.values()
        )
        stable = all(
            c["interval95"][0] > 0
            and c["interval95"][1] - c["interval95"][0] <= 0.20
            for c in comparisons.values()
        )
        results[candidate] = dict(
            comparisons=comparisons,
            wall_ratios=wall,
            completion_gate=completion,
            stable=stable,
            promising=positive_point
            and completion
            and max(wall.values()) <= 1.05,
        )
    extend = samples < 27 and any(
        r["promising"] and not r["stable"] for r in results.values()
    )
    eligible = [
        name for name, r in results.items() if r["promising"] and r["stable"]
    ]
    chosen = (
        min(eligible, key=lambda a: reporter.weighted_cost(groups, a, "cpu"))
        if eligible and not extend
        else None
    )
    return dict(
        protocol_sha256=first.digest(CONTROL),
        samples=samples,
        next_samples=(samples + 9 if extend else None),
        selected=chosen,
        candidates=results,
        status="Revealed six-input screen; no fresh or default acceptance",
    )


def report(output):
    first.runtime()
    fixtures = verify_control()
    samples = target_samples(output)
    if samples is None:
        raise ValueError("source-preset screening has already stopped")
    with first.machine_window(), layout():
        allowance = matched.Allowance(output / "ledger.json", 183)
        try:
            destination = output / f"assessment-{samples:02d}.json"
            if not allowance.reserve("report", 180, str(destination)):
                raise ValueError("source-preset cap cannot fund assessment")
            with matched.deadline(180):
                groups = defaultdict(lambda: defaultdict(dict))
                expected = set()
                for fixture, sample, seed, _order in assignment(
                    fixtures, samples
                ):
                    name = matched.group_name(fixture, sample)
                    expected.add(name)
                    rows = matched.read_group(
                        output / name, fixture, sample, seed
                    )
                    if rows is None:
                        raise ValueError(
                            "incomplete preset assignment cannot select"
                        )
                    # The relation bundle, total memory and reserve are matched
                    # across every arm; changed curve workspace is permitted.
                    reference = rows[0]
                    if any(
                        r["siqs_config"] != reference["siqs_config"]
                        or r["allocation"] != reference["allocation"]
                        or r["memory_cap"] != reference["memory_cap"]
                        for r in rows
                    ):
                        raise ValueError(
                            "preset changed protected SIQS bounds"
                        )
                    key = (first.band_of(fixture["n"]), fixture["kind"])
                    groups[key][fixture["id"]][sample] = {
                        r["arm"]: r for r in rows
                    }
                if {
                    p.name for p in output.iterdir() if p.is_dir()
                } != expected:
                    raise ValueError("unexpected preset groups")
                first.save(destination, assess(groups, samples))
                allowance.ledger.setdefault("assessments", {})[
                    str(samples)
                ] = first.digest(destination)
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
    parser.add_argument("--lease-seconds", type=int, default=1200)
    args = parser.parse_args()
    if args.mode == "measure":
        measure(args.output, args.lease_seconds)
    else:
        report(args.output)


if __name__ == "__main__":
    main()
