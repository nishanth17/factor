"""Machine-checkable v4 experiment design; this module never runs factoring."""

import argparse
import hashlib
import json
import sys
from dataclasses import asdict, replace
from pathlib import Path

from .configuration import (
    ATOMS,
    COUNTS,
    FACTORS,
    OWNED,
    PARTIALS,
    ROWS,
    RSS,
    TIERS,
)
from .configuration_v4 import configuration
from .control_v4 import BASELINE, PARENT_OWNED, START, configuration_values
from .inputs import CORPUS
from .work_model import RATE, charges

ROOT = Path(__file__).resolve().parents[4]
PLAN = Path(__file__).with_name("experiment_v4.json")
SOURCES = Path(__file__).with_name("experiment_sources_v4.json")
SETUP = 1_100_000_000_000
MARGIN = 5
POSITIONS = 131073
BLOCKS = 33
CAPS = {60: 30, 62: 45, 65: 60, 70: 90}
PILOT_CAPS = {60: 1800, 62: 2400, 65: 3600, 70: 7200}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def call_work(digits, seconds):
    """Finite loose guard, including all postprocessing and cache rebuilds."""
    f, r = COUNTS[digits], FACTORS[digits]
    bits = (3322 * digits + 999) // 1000
    cache = 100000 * (7 * (18 * r + 2) + f * (r + 7) * (18 * r + 1 + 324))
    pairs = 100000 * f * (18 * r + 3)
    window = f * (r + 1) * (bits + 8)
    return (
        4 * RATE * seconds
        + SETUP
        + cache
        + pairs
        + window
        + charges(digits)["postprocessing"]
    )


def diagnostic_work(seconds, setups):
    """Root/setup-only calls have no filtering or dependency solver."""
    # Additional bounded family traversal is below 100 billion charged units
    # for 9 samples * 8 families * 64 polynomials at the largest configuration.
    return 4 * RATE * seconds + setups * SETUP + 100_000_000_000


def config_values(digits, arm, *, multiplier=0, **options):
    if arm == "control":
        if multiplier or options:
            raise ValueError(
                "mainline control cannot accept challenger options"
            )
        return configuration_values(digits)
    if arm not in ("full", "nominee"):
        raise ValueError("unknown experiment arm")
    cfg = configuration(digits, "coupled", multiplier=multiplier, **options)
    cfg = replace(cfg)
    return asdict(cfg)


def cell(
    identity,
    phase,
    digits,
    balance,
    arm,
    *,
    seconds=None,
    polynomials=512,
    multiplier=0,
    kind="prefix",
    **options,
):
    seconds = CAPS[digits] if seconds is None else seconds
    fixture = f"c10a01_{digits}_{balance}"
    work = (
        diagnostic_work(seconds, 9 if kind == "setup" else 1)
        if kind in ("roots", "setup")
        else call_work(digits, seconds)
    )
    return dict(
        id=identity,
        phase=phase,
        fixture=fixture,
        digits=digits,
        arm=arm,
        kind=kind,
        multiplier=multiplier,
        seconds=seconds,
        reservation_seconds=seconds + MARGIN,
        work=work,
        tier=TIERS[digits],
        polynomials=polynomials if kind == "prefix" else 0,
        positions=polynomials * POSITIONS if kind == "prefix" else 0,
        blocks=polynomials * BLOCKS if kind == "prefix" else 0,
        configuration=config_values(
            digits, arm, multiplier=multiplier, **options
        ),
        samples=9 if kind in ("roots", "setup") else 1,
        families_per_sample=8 if kind == "roots" else 0,
        gray_steps_per_family=64 if kind == "roots" else 0,
    )


def totals(calls, warm_slots, analysis_seconds, orchestration_seconds):
    phases = {}
    for row in calls:
        entry = phases.setdefault(row["phase"], dict(wall=0, cpu=0, work=0))
        entry["wall"] += row["reservation_seconds"]
        entry["cpu"] += row["reservation_seconds"]
        entry["work"] += row["work"]
    # At most 32 fresh semantic warm/sentinel invocations per context. This
    # deliberately includes full postprocessing bounds; actual work is settled.
    warm_work = 32 * max(call_work(row["digits"], 20) for row in calls)
    phases["warmup"] = dict(
        wall=25 * warm_slots, cpu=25 * warm_slots, work=warm_slots * warm_work
    )
    phases["analysis"] = dict(
        wall=analysis_seconds, cpu=analysis_seconds, work=1_000_000_000_000
    )
    phases["orchestration"] = dict(
        wall=orchestration_seconds, cpu=orchestration_seconds, work=0
    )
    total = {
        field: sum(p[field] for p in phases.values())
        for field in ("wall", "cpu", "work")
    }
    return dict(
        phases=phases,
        total=total,
        warm_slots=warm_slots,
        warm_slot_work=warm_work,
        warm_slot_seconds=25,
    )


def build_plan():
    """Produce exact cells and sums without creating a job or a Budget."""
    calls = []
    for band, digits in enumerate((60, 62, 65, 70)):
        first = ("control", "full") if band % 2 == 0 else ("full", "control")
        for index, balance in enumerate(("balanced", "35_65")):
            for arm in first if index == 0 else first[::-1]:
                calls.append(
                    cell(
                        f"primary-{digits}-{balance}-{arm}",
                        "primary",
                        digits,
                        balance,
                        arm,
                    )
                )
    for digits, k in ((62, 13), (70, 1)):
        for policy in ("standard", "cached", "batch"):
            calls.append(
                cell(
                    f"roots-{digits}-{policy}",
                    "roots",
                    digits,
                    "balanced",
                    "full",
                    seconds=45,
                    multiplier=k,
                    kind="roots",
                    root_update_policy=policy,
                )
            )
        for policy in ("windowed_pair", "paired_flyer"):
            calls.append(
                cell(
                    f"population-{digits}-{policy}",
                    "population",
                    digits,
                    "balanced",
                    "full",
                    seconds=45,
                    multiplier=k,
                    assignment_policy=policy,
                )
            )
        for policy in ("parity2_norm", "parity2_wide"):
            calls.append(
                cell(
                    f"ranking-{digits}-{policy}",
                    "ranking",
                    digits,
                    "balanced",
                    "full",
                    seconds=45,
                    kind="setup",
                    multiplier_policy=policy,
                )
            )
    for form in ("standard", "parity2"):
        row = cell(
            f"forced-60-residue5-{form}",
            "forced",
            60,
            "balanced",
            "full",
            seconds=45,
            multiplier=5,
            polynomials=768,
            multiplier_policy="standard",
        )
        row["configuration"]["polynomial_policy"] = form
        calls.append(row)
    calls.append(
        cell(
            "forced-62-even",
            "forced",
            62,
            "35_65",
            "full",
            seconds=45,
            multiplier=2,
            polynomials=640,
        )
    )
    calls.append(
        cell(
            "forced-70-eligible",
            "forced",
            70,
            "balanced",
            "full",
            seconds=45,
            multiplier=31,
        )
    )
    for digits, order in ((60, (False, True)), (70, (True, False))):
        for observer in order:
            row = cell(
                f"observer-{digits}-{int(observer)}",
                "observer",
                digits,
                "balanced",
                "full",
                seconds=15,
                polynomials=32,
            )
            row["observer"] = observer
            calls.append(row)
    stages = {}
    for name, sizes, warms, analysis, orchestration in (
        ("A", (60, 62), 24, 120, 60),
        ("B", (65,), 12, 60, 30),
        ("C", (70,), 12, 60, 30),
    ):
        rows = []
        for digits in sizes:
            order = (
                ("control", "nominee")
                if digits in (60, 65)
                else ("nominee", "control")
            )
            for arm in order:
                row = cell(
                    f"pilot-{digits}-{arm}",
                    "measured",
                    digits,
                    "balanced",
                    arm,
                    seconds=PILOT_CAPS[digits],
                    kind="complete",
                )
                row["configuration_status"] = (
                    "frozen_control"
                    if arm == "control"
                    else "template_requires_nominee_freeze"
                )
                rows.append(row)
        stages[name] = dict(
            calls=rows, limits=totals(rows, warms, analysis, orchestration)
        )
    resources = []
    for d in (60, 62, 65, 70):
        resources.append(
            dict(
                digits=d,
                base_count=COUNTS[d],
                factor_count=FACTORS[d],
                rows=ROWS[d],
                atoms=ATOMS[d],
                partials=PARTIALS[d],
                nonzero_ceiling=ROWS[d] * (COUNTS[d] + 1),
                tier=TIERS[d],
                owned=OWNED[TIERS[d]],
                rss=RSS[TIERS[d]],
                parent_owned=PARENT_OWNED,
                parent_rss=256 * 2**20,
            )
        )
    return dict(
        identity="c10a01-v4-experiment-design-1",
        execution_status="design_only",
        implementation_commit="4b9a6163898a0b485436665d31c2a274e51ef643",
        control_commit=START,
        corpus_sha256=digest(CORPUS),
        control_snapshot_sha256=digest(BASELINE),
        source_manifest=SOURCES.name,
        seed=7,
        calibration=dict(
            calls=calls,
            limits=totals(calls, 40, 180, 260),
            expected_wall_seconds=[720, 1200],
        ),
        pilots=stages,
        resources=resources,
        accounting=dict(
            rate=RATE,
            headroom=4,
            setup=SETUP,
            original_calibration_ceiling=3600,
            original_calibration_charged=3124,
            original_calibration_work=138632524390,
            new_calibration_ceiling=3400,
            later_pilots=31600,
            combined_ceiling=38600,
            prior_charged_plus_new_ceiling=38124,
        ),
        recording=dict(
            total_study_bytes=512 * 2**20,
            new_output_bytes=64 * 2**20,
            raw_bytes=40 * 2**20,
            archive_bytes=24 * 2**20,
            receipt_bytes=512 * 1024,
            worker_stdout_bytes=64 * 1024,
            observer_records=8192,
            family_records=64,
            polynomial_records=4096,
        ),
        execution_prerequisites=[
            "Implement/test the bounded adapters in protocol_v4.md.",
            "Pin adapter sources and verify committed-only loading.",
            "Freeze one nominee and rate/CPU/RSS readiness before pilots.",
        ],
    )


def verify_plan(plan):
    """Reject altered inputs, configurations, omitted stages or budget sums."""
    if plan != build_plan():
        raise ValueError("design differs from its exact cells or accounting")
    calibration = plan["calibration"]
    assert len(calibration["calls"]) == 38
    assert calibration["limits"]["total"]["wall"] == 3400
    assert [
        plan["pilots"][s]["limits"]["total"]["wall"] for s in ("A", "B", "C")
    ] == [9200, 7600, 14800]
    assert (
        3600
        + 3400
        + sum(
            plan["pilots"][s]["limits"]["total"]["wall"]
            for s in ("A", "B", "C")
        )
        == 38600
    )
    for stage in [calibration, *plan["pilots"].values()]:
        assert (
            stage["limits"]["total"]["cpu"] == stage["limits"]["total"]["wall"]
        )
        assert len({c["id"] for c in stage["calls"]}) == len(stage["calls"])
        for row in stage["calls"]:
            assert row["configuration"]["memory_bytes"] > 0
            if row["kind"] == "prefix":
                assert row["positions"] == row["polynomials"] * POSITIONS
                assert row["blocks"] == row["polynomials"] * BLOCKS
            if row.get("observer"):
                assert row["blocks"] < plan["recording"]["observer_records"]
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check",))
    parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 is required")
    plan = json.loads(PLAN.read_text())
    verify_plan(plan)
    sources = json.loads(SOURCES.read_text())
    if sources["identity"] != plan["identity"]:
        raise RuntimeError("unexpected source manifest identity")
    runtime = sources["runtime"]
    if (
        sys.version != runtime["version"]
        or digest(Path(sys.executable).resolve()) != runtime["binary_sha256"]
    ):
        raise RuntimeError("runtime drift: prepare a new prospective release")
    for name, expected in sources["files"].items():
        if digest(ROOT / name) != expected:
            raise RuntimeError("design/source drift: " + name)
    print(
        json.dumps(
            dict(
                identity=plan["identity"],
                numerical_plan_verified=True,
                execution_ready=False,
                performance_calls=0,
                calibration_calls=38,
                calibration_seconds=3400,
                pilots={
                    s: plan["pilots"][s]["limits"]["total"]
                    for s in ("A", "B", "C")
                },
                combined_ceiling_seconds=38600,
                prerequisites=plan["execution_prerequisites"],
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
