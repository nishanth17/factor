"""Apply the frozen C3 paired input-subject acceptance rule."""

import argparse
import json
import random
import statistics
from collections import defaultdict
from pathlib import Path


def percentile(values, probability):
    values = sorted(values)
    position = probability * (len(values) - 1)
    left = int(position)
    weight = position - left
    return (
        values[left] * (1 - weight)
        + values[min(left + 1, len(values) - 1)] * weight
    )


def summarize(rows):
    """Keep timing replicates nested inside input/seed subjects."""
    cells = defaultdict(list)
    for row in rows:
        cells[row["id"], row["seed"], row["arm"]].append(row)
    medians, unstable = {}, []
    for key, samples in cells.items():
        wall = [r["wall"] for r in samples]
        quartiles = statistics.quantiles(wall, n=4)
        spread = (quartiles[2] - quartiles[0]) / statistics.median(wall)
        if spread > 0.15:
            unstable.append(
                dict(
                    id=key[0],
                    seed=key[1],
                    arm=key[2],
                    relative_iqr=spread,
                    samples=len(samples),
                )
            )
        medians[key] = dict(
            id=key[0],
            seed=key[1],
            arm=key[2],
            kind=samples[0]["kind"],
            band=samples[0]["band"],
            complete=sum(r["complete"] for r in samples) / len(samples),
            proper_factor=sum(r["proper_factor"] for r in samples)
            / len(samples),
            wall=statistics.median(wall),
            capped_cost=statistics.median(
                [
                    r["wall"] if r["complete"] else r["cap_seconds"]
                    for r in samples
                ]
            ),
            cpu=statistics.median(r["cpu"] for r in samples),
            work=statistics.median(r["work"] for r in samples),
            curves=statistics.median(r["curves"] for r in samples),
            fallback_started=sum(r["fallback_started"] for r in samples)
            / len(samples),
            owned_cap=max(r["memory_cap"] for r in samples),
            workspace_reserve=max(r["workspace_reserve"] for r in samples),
            fallback_owned_peak=max(r["fallback_owned_peak"] for r in samples),
            rss_high_water=max(r["rss"] for r in samples),
            checkpoint_bytes=max(r["checkpoint_bytes"] for r in samples),
            samples=len(samples),
            relative_iqr=spread,
        )
    return medians, unstable


def compare(medians, subjects, challenger):
    """Bootstrap independent inputs, retaining every matched seed per input."""
    pairs = {}
    for identity in subjects:
        seeds = sorted(
            {seed for name, seed, arm in medians if name == identity}
        )
        pairs[identity] = [
            (
                medians[identity, seed, "control"],
                medians[identity, seed, challenger],
            )
            for seed in seeds
        ]

    def saving(names):
        control = sum(
            a["capped_cost"] for name in names for a, b in pairs[name]
        )
        candidate = sum(
            b["capped_cost"] for name in names for a, b in pairs[name]
        )
        return 1 - candidate / control

    generator = random.Random(193003)
    distribution = [
        saving(generator.choices(subjects, k=len(subjects)))
        for _ in range(4000)
    ]
    cells = [pair for values in pairs.values() for pair in values]
    return dict(
        inputs=len(subjects),
        input_seeds=len(cells),
        saving=saving(subjects),
        interval95=[
            percentile(distribution, 0.025),
            percentile(distribution, 0.975),
        ],
        control_completion=statistics.mean(a["complete"] for a, b in cells),
        candidate_completion=statistics.mean(b["complete"] for a, b in cells),
        control_yield=statistics.mean(a["proper_factor"] for a, b in cells),
        candidate_yield=statistics.mean(b["proper_factor"] for a, b in cells),
    )


def analyze(capture, selected):
    data = json.loads(capture.read_text())
    if data.get("stopped") or data["mode"] != "confirm":
        raise ValueError("acceptance requires completed confirmation")
    medians, unstable = summarize(data["rows"])
    decisions = {}
    for band in (30, 40):
        arm = selected["bands"][str(band)]
        names = sorted(
            {r["id"] for r in medians.values() if r["band"] == band}
        )
        if arm == "control":
            decisions[str(band)] = dict(
                decision="retain",
                selected=arm,
                reason="training retained the control",
            )
            continue
        overall = compare(medians, names, arm)
        classes = {}
        for kind in sorted(
            {r["kind"] for r in medians.values() if r["band"] == band}
        ):
            subjects = sorted(
                {
                    r["id"]
                    for r in medians.values()
                    if r["band"] == band and r["kind"] == kind
                }
            )
            classes[kind] = compare(medians, subjects, arm)
        stable = not any(
            medians[key]["band"] == band
            for key in medians
            if medians[key]["relative_iqr"] > 0.15
        )
        completion_ok = all(
            value["candidate_completion"] >= value["control_completion"] - 0.05
            for value in classes.values()
        )
        promote = stable and completion_ok and overall["interval95"][0] > 0
        decisions[str(band)] = dict(
            decision="adopt-scoped" if promote else "retain",
            selected=arm,
            overall=overall,
            classes=classes,
            stable=stable,
            completion_gate=completion_ok,
        )
    return dict(
        decisions=decisions,
        unstable=unstable,
        rows=list(medians.values()),
        wall=data["wall"],
        cpu=data["cpu"],
        protocol_sha256=data["protocol_sha256"],
        limitations=(
            "Conditional input-subject bootstrap; small biased corpus. "
            "Single-input classes are regression controls, "
            "not population estimates."
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture", type=Path)
    parser.add_argument("selected", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = analyze(args.capture, json.loads(args.selected.read_text()))
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result["decisions"], indent=2))


if __name__ == "__main__":
    main()
