"""Frozen calibration, offline sequence fitting and warmed complete calls."""

import argparse
import contextlib
import json
import os
import subprocess
import time
from collections import Counter, defaultdict
from unittest.mock import patch

from .... import portfolio
from ..p52.p52_realistic import check_quiet
from . import c3_study as first
from . import round2_pilot as pilot
from .round2_corpus import load_training
from .round2_sequence import LADDERS, fit_bands

CONTROL = first.INPUTS / "controls/c3_round2_training_v1.json"
SELECTION = first.INPUTS / "controls/c3_round2_fitted_v1.json"
SEEDS = (17, 43, 89, 113, 151, 181, 223, 269, 307)
CALIBRATION_ARMS = ("no_ecm", "wide", "compact", "deep")


@contextlib.contextmanager
def root_trace(n, enabled):
    """Time whole curves at job/event boundaries, outside acceptance calls."""
    trace = dict(
        censored=False,
        curves=[],
        relation_cpu=None,
        relation_start=None,
        first_curve_start=None,
        hit_end=None,
    )
    if not enabled:
        yield trace
        return
    started = time.process_time()
    active = None
    new_job, event, relation = (
        portfolio.new_job,
        portfolio._event,
        portfolio.SIQSJob,
    )

    def begin(kind, modulus, seed, b1, b2):
        nonlocal active
        job = new_job(kind, modulus, seed, b1, b2)
        if kind == "ecm" and modulus == n:
            active = time.process_time()
            if trace["first_curve_start"] is None:
                trace["first_curve_start"] = active - started
        return job

    def record(state, config, **values):
        nonlocal active
        if values.get("n") == n:
            if values["stage"] == "ecm" and active is not None:
                ended = time.process_time()
                if values["outcome"] not in ("factor", "exhausted"):
                    trace["censored"] = True
                else:
                    trace["curves"].append(
                        dict(
                            cpu=ended - active,
                            outcome=values["outcome"],
                            b1=values["b1"],
                            b2=values["b2"],
                            seed=values["seed"],
                        )
                    )
                    if values["outcome"] == "factor":
                        trace["hit_end"] = ended
                active = None
            elif values["stage"] == "siqs":
                trace["relation_cpu"] = (
                    time.process_time() - started - trace["relation_start"]
                )
        return event(state, config, **values)

    def start_relation(modulus, *args, **kwargs):
        if modulus == n:
            trace["relation_start"] = time.process_time() - started
        return relation(modulus, *args, **kwargs)

    with (
        patch.object(portfolio, "new_job", begin),
        patch.object(portfolio, "_event", record),
        patch.object(portfolio, "SIQSJob", start_relation),
    ):
        yield trace
    trace["censored"] |= active is not None
    trace["recursive_cpu"] = (
        time.process_time() - trace["hit_end"] if trace["hit_end"] else 0
    )


def run_one(fixture, seed, arm, control, *, capture=False, fitted=None):
    """Only observable size and a compiled tier table configure dispatch."""
    if arm == "control":
        return pilot.run_one(fixture, seed, arm, control)
    if arm == "fixed32":
        tiers = ((2000, 147396, 32),)
    elif arm == "fitted":
        tiers = tuple(
            tuple(t) for t in fitted[str(first.band_of(fixture["n"]))]["tiers"]
        )
    elif arm == "no_ecm":
        tiers = ()
    else:
        tiers = LADDERS[arm]
    settings = {arm: dict(tiers=tiers, mode="campaign")}
    with patch.dict(pilot.POLICIES, settings):
        with root_trace(fixture["n"], capture) as trace:
            row = pilot.run_one(fixture, seed, arm, control)
    row["instrumented"] = capture
    if capture:
        row["trace"] = trace
    return row


def verify_control():
    protocol = json.loads(CONTROL.read_text())
    for name, digest in protocol["sha256"].items():
        if first.digest(first.ROOT / name) != digest:
            raise ValueError(
                "frozen round-two training source changed: " + name
            )
    return protocol


def row_path(output, fixture, seed, sample, arm):
    return output / f"{fixture['id']}-{seed}-{sample}-{arm}.json"


def measure(mode, output, *, lease_seconds):
    """Resume only completed captures and spend a finite prepaid phase cap."""
    first.runtime()
    protocol = verify_control()
    output.mkdir(parents=True, exist_ok=True)
    receipt_path = output / "ledger.json"
    ledger = (
        json.loads(receipt_path.read_text())
        if receipt_path.exists()
        else dict(
            mode=mode,
            protocol_sha256=first.digest(CONTROL),
            active_seconds=0,
            leases=[],
            warmups=[],
        )
    )
    if ledger["mode"] != mode or ledger["protocol_sha256"] != first.digest(
        CONTROL
    ):
        raise ValueError("incompatible continued phase")
    if not 60 <= lease_seconds <= 1200:
        raise ValueError("lease must be 60..1200 active seconds")
    fitted = None
    if mode == "calibrate":
        arms = CALIBRATION_ARMS
    else:
        selection = json.loads(SELECTION.read_text())
        if selection["protocol_sha256"] != first.digest(CONTROL):
            raise ValueError("selection belongs to another freeze")
        relative = str(SELECTION.relative_to(first.ROOT))
        committed = subprocess.check_output(
            ["git", "show", "HEAD:" + relative], cwd=first.ROOT
        )
        if committed != SELECTION.read_bytes():
            raise ValueError(
                "fitted candidate must be committed before comparison"
            )
        fitted = selection["bands"]
        arms = ("control", "fixed32", "fitted")
    fixtures = load_training()
    control = first.baseline()
    spent, started = ledger["active_seconds"], time.monotonic()
    allowance = min(lease_seconds, protocol["phase_seconds"][mode] - spent)
    stopped = None

    def available(reservation):
        return time.monotonic() - started + reservation <= allowance

    def save_ledger(prepaid=0):
        # A killed process leaves its full active-call grant charged. Completed
        # calls replace that conservative reservation with observed expense.
        ledger["active_seconds"] = spent + time.monotonic() - started + prepaid
        receipt_path.write_text(json.dumps(ledger, indent=2) + "\n")

    with first.machine_window():
        try:
            if not available(30 * len(arms)):
                raise RuntimeError("validated warmup does not fit phase")
            warm_fixture = next(
                f for f in fixtures if f["id"] == "r2_30_balanced_0"
            )
            for arm in arms:
                save_ledger(prepaid=30)
                wall, cpu, count = time.monotonic(), time.process_time(), 0
                while (
                    min(time.monotonic() - wall, time.process_time() - cpu) < 3
                ):
                    if time.monotonic() - wall > 25:
                        raise RuntimeError("warmup cap exhausted")
                    run_one(
                        warm_fixture,
                        17,
                        arm,
                        control,
                        capture=mode == "calibrate",
                        fitted=fitted,
                    )
                    count += 1
                ledger["warmups"].append(
                    dict(
                        arm=arm,
                        wall=time.monotonic() - wall,
                        cpu=time.process_time() - cpu,
                        validated=count,
                    )
                )
            for sample in range(9):
                seed = SEEDS[sample % len(SEEDS)]
                for fixture_index, fixture in enumerate(fixtures):
                    offset = (sample + fixture_index) % len(arms)
                    order = arms[offset:] + arms[:offset]
                    if sample % 2:
                        order = tuple(reversed(order))
                    for arm in order:
                        path = row_path(output, fixture, seed, sample, arm)
                        if path.exists():
                            row = json.loads(path.read_text())
                            if row["arm"] != arm or row["id"] != fixture["id"]:
                                raise ValueError("incompatible row capture")
                            continue
                        cap = 5 if first.band_of(fixture["n"]) == 30 else 30
                        if not available(cap + 1):
                            stopped = "phase_allowance"
                            break
                        check_quiet({os.getpid()})
                        save_ledger(prepaid=cap + 1)
                        row = run_one(
                            fixture,
                            seed,
                            arm,
                            control,
                            capture=mode == "calibrate",
                            fitted=fitted,
                        )
                        row["sample"] = sample
                        row["protocol_sha256"] = first.digest(CONTROL)
                        first.save(path, row)
                        save_ledger()
                        print(
                            mode,
                            sample,
                            fixture["id"],
                            arm,
                            row["reason"],
                            round(row["cpu"], 4),
                            flush=True,
                        )
                    if stopped:
                        break
                if stopped:
                    break
        finally:
            ledger["leases"].append(
                dict(
                    wall=time.monotonic() - started,
                    allowance=allowance,
                    stopped=stopped,
                )
            )
            save_ledger()
    return ledger


def fit(output):
    """Fit once only after the complete, verified nine-sample calibration."""
    protocol = verify_control()
    rows = [
        json.loads(p.read_text())
        for p in sorted(output.glob("*.json"))
        if p.name != "ledger.json"
    ]
    fixtures = load_training()
    expected = {
        (f["id"], sample, arm)
        for f in fixtures
        for sample in range(9)
        for arm in CALIBRATION_ARMS
    }
    actual = {(r["id"], r["sample"], r["arm"]) for r in rows}
    if actual != expected or len(rows) != len(expected):
        raise ValueError("incomplete calibration; cannot fit or select")
    if any(
        not r["complete"] or r["protocol_sha256"] != first.digest(CONTROL)
        for r in rows
    ):
        raise ValueError("incomplete or incompatible calibration outcome")
    if any(r["seed"] != SEEDS[r["sample"]] for r in rows):
        raise ValueError("calibration assignment changed")
    by_key = {(r["id"], r["sample"], r["arm"]): r for r in rows}
    strata = Counter((first.band_of(f["n"]), f["kind"]) for f in fixtures)
    records = defaultdict(list)
    for fixture in fixtures:
        band = first.band_of(fixture["n"])
        for sample in range(9):
            relation = by_key[(fixture["id"], sample, "no_ecm")]["trace"]
            for ladder in LADDERS:
                row = by_key[(fixture["id"], sample, ladder)]
                trace = row["trace"]
                if trace["censored"]:
                    raise ValueError(
                        "censored root curve cannot fit a sequence"
                    )
                curves = trace["curves"]
                if not curves:
                    continue
                if relation["relation_cpu"] is None:
                    raise ValueError("ECM entrant lacks matched SIQS cost")
                hit = next(
                    (
                        i + 1
                        for i, c in enumerate(curves)
                        if c["outcome"] == "factor"
                    ),
                    None,
                )
                records[ladder].append(
                    dict(
                        subject=fixture["id"],
                        seed=row["seed"],
                        band=band,
                        weight=1 / (9 * strata[(band, fixture["kind"])]),
                        curve_costs=[c["cpu"] for c in curves],
                        hit=hit,
                        qs_cost=relation["relation_cpu"],
                        recursive_cost=trace["recursive_cpu"],
                        setup=max(
                            0,
                            trace["first_curve_start"]
                            - relation["relation_start"],
                        ),
                    )
                )
    first.save(
        SELECTION,
        dict(
            protocol_sha256=first.digest(CONTROL),
            bands=fit_bands(records),
            calibration_sha256={
                p.name: first.digest(p) for p in sorted(output.glob("*.json"))
            },
            status="Fitted training candidate; not accepted policy",
            source_commit=protocol["source_commit"],
        ),
    )
    print(SELECTION)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("calibrate", "fit", "compare"))
    parser.add_argument("output", type=first.Path)
    parser.add_argument("--lease-seconds", type=int, default=600)
    args = parser.parse_args()
    if args.mode == "fit":
        fit(args.output)
    else:
        measure(args.mode, args.output, lease_seconds=args.lease_seconds)


if __name__ == "__main__":
    main()
