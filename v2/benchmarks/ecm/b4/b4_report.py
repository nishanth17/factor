"""Summarize frozen B4 captures without changing selection or gate rules."""

import argparse
import json
import random
import statistics
from pathlib import Path

from ...support.paths import (
    source_path,
)
from ..p52.p52_a3 import interval
from . import b4_common as common
from . import b4_study as study


def median(capture):
    return statistics.median(
        sample["seconds"] for sample in capture["samples"]
    )


def completion_interval(control, candidate):
    differences = {}
    before = {
        (row["fixture"], row["seed"]): row
        for row in control["samples"][0]["rows"]
    }
    for row in candidate["samples"][0]["rows"]:
        key = row["fixture"]
        old = before[(key, row["seed"])]
        differences.setdefault(key, []).append(
            int(row["complete"]) - int(old["complete"])
        )
    means = [100 * statistics.mean(values) for values in differences.values()]
    generator = random.Random(2026100942)
    bootstrap = sorted(
        statistics.mean(generator.choices(means, k=len(means)))
        for _ in range(3000)
    )
    return statistics.mean(means), [bootstrap[75], bootstrap[2924]]


def full_rows(captures, protocol):
    rows = []
    for candidate in captures:
        if candidate["scope"] != "full" or candidate["arm"] == "baseline":
            continue
        control = next(
            c
            for c in captures
            if c["scope"] == "full"
            and c["arm"] == "baseline"
            and c["backend"] == candidate["backend"]
        )
        reduction = 100 * (1 - median(candidate) / median(control))
        timing = interval(
            [s["seconds"] for s in control["samples"]],
            [s["seconds"] for s in candidate["samples"]],
        )
        completion, uncertainty = completion_interval(control, candidate)
        regression = max(
            study.completions(control)[case]
            - study.completions(candidate)[case]
            for case in study.completions(control)
        )
        stable = study.stable(control, protocol) and study.stable(
            candidate, protocol
        )
        passed = (
            stable
            and regression <= 0.05
            and (
                (reduction >= 10 and timing[0] > 0)
                or (completion >= 10 and uncertainty[0] > 0)
            )
        )
        classes = []
        for case in sorted(study.completions(control)):
            before = [
                sum(r["seconds"] for r in sample["rows"] if r["case"] == case)
                for sample in control["samples"]
            ]
            after = [
                sum(r["seconds"] for r in sample["rows"] if r["case"] == case)
                for sample in candidate["samples"]
            ]
            classes.append(
                dict(
                    case=case,
                    control_seconds=statistics.median(before),
                    candidate_seconds=statistics.median(after),
                    reduction_percent=100
                    * (
                        1
                        - statistics.median(after) / statistics.median(before)
                    ),
                )
            )
        rows.append(
            dict(
                backend=candidate["backend"],
                arm=candidate["arm"],
                control_seconds=median(control),
                candidate_seconds=median(candidate),
                reduction_percent=reduction,
                conditional_timing_interval=timing,
                completion_change_points=completion,
                fixture_cluster_completion_interval=uncertainty,
                max_completion_regression_points=100 * regression,
                stable=stable,
                gate_passed=passed,
                classes=classes,
            )
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--confirmation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    common.require_runtime()
    study.verify_freeze()
    protocol, _ = common.inputs()
    training = json.loads(source_path(args.training).read_text())
    held_out = json.loads(source_path(args.confirmation).read_text())
    if training["freeze"] != held_out["freeze"]:
        raise ValueError("confirmation sources differ")
    for capture in held_out["captures"]:
        if capture["arm"] not in (
            "baseline",
            training["selected"][capture["backend"]],
        ):
            raise ValueError("confirmation violates training selection")
    report = dict(
        training_sha256=common.digest(args.training),
        confirmation_sha256=common.digest(args.confirmation),
        training=full_rows(training["captures"], protocol),
        confirmation=full_rows(held_out["captures"], protocol),
    )
    with args.output.open("x") as output:
        json.dump(report, output, indent=2)
        output.write("\n")
    print(json.dumps(report["confirmation"], indent=2))


if __name__ == "__main__":
    main()
