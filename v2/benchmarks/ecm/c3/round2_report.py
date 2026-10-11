"""Matched complete-call assessment; never rank instrumented calibration."""

import argparse
import json
import random
from collections import defaultdict
from statistics import mean

from . import c3_study as first
from . import round2_train as training
from .round2_corpus import load_training


def paired_groups(rows, fixtures):
    """Require the entire frozen assignment before computing a decision."""
    expected = {
        (f["id"], sample, arm)
        for f in fixtures
        for sample in range(9)
        for arm in ("control", "fixed32", "fitted")
    }
    keys = [(r["id"], r["sample"], r["arm"]) for r in rows]
    if set(keys) != expected or len(keys) != len(expected):
        raise ValueError("incomplete complete-call comparison")
    identity = {f["id"]: (first.band_of(f["n"]), f["kind"]) for f in fixtures}
    groups = defaultdict(lambda: defaultdict(dict))
    for row in rows:
        if row["instrumented"] or row["seed"] != training.SEEDS[row["sample"]]:
            raise ValueError("instrumented or unmatched observation")
        if row["protocol_sha256"] != first.digest(training.CONTROL):
            raise ValueError("comparison source identity changed")
        groups[identity[row["id"]]][row["id"]].setdefault(row["sample"], {})[
            row["arm"]
        ] = row
    return groups


def capped_cost(row, metric):
    """Early unresolved returns pay the full grant rather than look fast."""
    return (
        row[metric]
        if row["complete"]
        else max(row[metric], row["cap_seconds"])
    )


def weighted_cost(groups, arm, metric, generator=None):
    """Weight bands, kinds and subjects equally; resample paired seeds."""
    band_values = defaultdict(list)
    for (band, _), subjects in sorted(groups.items()):
        identifiers = sorted(subjects)
        chosen = (
            generator.choices(identifiers, k=len(identifiers))
            if generator
            else identifiers
        )
        values = []
        for identifier in chosen:
            samples = subjects[identifier]
            indices = sorted(samples)
            if generator:
                indices = generator.choices(indices, k=len(indices))
            values.append(
                mean(capped_cost(samples[i][arm], metric) for i in indices)
            )
        band_values[band].append(mean(values))
    return mean(mean(values) for values in band_values.values())


def improvement_interval(groups, control, *, repetitions=10_000):
    """Cluster by input and keep all treatment/control seed pairs together."""
    baseline = weighted_cost(groups, control, "cpu")
    candidate = weighted_cost(groups, "fitted", "cpu")
    generator = random.Random(2026101109)
    samples = []
    for _ in range(repetitions):
        # Replaying the same generator state pairs the input and seed draws.
        state = generator.getstate()
        old = weighted_cost(groups, control, "cpu", generator)
        generator.setstate(state)
        new = weighted_cost(groups, "fitted", "cpu", generator)
        samples.append((old - new) / old)
    samples.sort()
    return dict(
        baseline_cpu=baseline,
        fitted_cpu=candidate,
        relative_improvement=(baseline - candidate) / baseline,
        interval95=[
            samples[int(0.025 * repetitions)],
            samples[min(repetitions - 1, int(0.975 * repetitions))],
        ],
    )


def assess(groups, *, repetitions=10_000):
    """Apply the predeclared economics and class-completion gates."""
    comparisons = {
        arm: improvement_interval(groups, arm, repetitions=repetitions)
        for arm in ("control", "fixed32")
    }
    classes = {}
    for (band, kind), subjects in sorted(groups.items()):
        rows = [
            arms for samples in subjects.values() for arms in samples.values()
        ]
        classes[f"{band}:{kind}"] = {
            arm: dict(
                completion=mean(r[arm]["complete"] for r in rows),
                cpu=mean(capped_cost(r[arm], "cpu") for r in rows),
                wall=mean(capped_cost(r[arm], "wall") for r in rows),
            )
            for arm in ("control", "fixed32", "fitted")
        }
    completion_ok = all(
        values["fitted"]["completion"] >= values[arm]["completion"] - 0.05
        for values in classes.values()
        for arm in ("control", "fixed32")
    )
    wall_ratios = {
        arm: weighted_cost(groups, "fitted", "wall")
        / weighted_cost(groups, arm, "wall")
        for arm in ("control", "fixed32")
    }
    positive = all(c["interval95"][0] > 0 for c in comparisons.values())
    return dict(
        comparisons=comparisons,
        classes=classes,
        wall_ratios=wall_ratios,
        completion_gate=completion_ok,
        verdict=(
            "eligible for new fresh confirmation"
            if positive
            and completion_ok
            and all(r <= 1.05 for r in wall_ratios.values())
            else "retain fixed control; model has not earned promotion"
        ),
        status="Revealed training decision only; no fresh acceptance",
    )


def main():
    first.runtime()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("captures", type=first.Path)
    parser.add_argument("output", type=first.Path)
    args = parser.parse_args()
    training.verify_control()
    rows = [
        json.loads(p.read_text())
        for p in sorted(args.captures.glob("*.json"))
        if p.name != "ledger.json"
    ]
    groups = paired_groups(rows, load_training())
    first.save(args.output, assess(groups))


if __name__ == "__main__":
    main()
