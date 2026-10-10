"""Summarize frozen C6 captures without filtering unsuccessful attempts."""

import argparse
import json
import random
import statistics
from pathlib import Path

from ...support.paths import (
    source_path,
)
from ..p41 import p41_campaign as control


def ratio_interval(candidate, baseline):
    """Repeat-timing uncertainty for a fixed cohort."""
    generator = random.Random(61045)
    samples = []
    for _ in range(2000):
        after = [generator.choice(candidate) for _ in candidate]
        before = [generator.choice(baseline) for _ in baseline]
        samples.append(statistics.median(after) / statistics.median(before))
    samples.sort()
    return samples[50], samples[1949]


def class_timings(capture, baseline):
    """Exploratory class breakdown; never a post-hoc promotion cohort."""
    classes = sorted(
        {r["case"].split("_")[0] for r in capture["samples"][0]["rows"]}
    )
    result = {}
    for name in classes:
        groups = []
        for arm in (capture, baseline):
            groups.append(
                [
                    sum(
                        r["seconds"]
                        for r in sample["rows"]
                        if r["case"].startswith(name + "_")
                    )
                    for sample in arm["samples"]
                ]
            )
        after, before = groups
        result[name] = dict(
            median_seconds=statistics.median(after),
            ratio=statistics.median(after) / statistics.median(before),
            ratio_interval=ratio_interval(after, before),
        )
    return result


def summarize(path):
    report = json.loads(source_path(path).read_text())
    indexed = {(r["backend"], r["arm"]): r for r in report["captures"]}
    corpus = {c["id"]: c for c in control.load_corpus()["fixtures"]}
    rows, pairs = [], 0
    for (backend, arm), capture in indexed.items():
        baseline = indexed[backend, "ladder"]
        times = [s["seconds"] for s in capture["samples"]]
        before = [s["seconds"] for s in baseline["samples"]]
        distinct = {}
        for sample in capture["samples"]:
            for result in sample["rows"]:
                control.validate_result(result, corpus[result["case"]])
                key = result["case"], result["seed"]
                outcome = {
                    field: result[field]
                    for field in (
                        "factor",
                        "cofactor",
                        "unresolved",
                        "stats",
                        "extra",
                        "point",
                        "timed_out",
                    )
                }
                if key in distinct and outcome != distinct[key]:
                    raise AssertionError("nondeterministic repeated outcome")
                distinct[key] = outcome
        if backend == "gmp":
            native = indexed["int", arm]["samples"][0]["rows"]
            for result in native:
                key = result["case"], result["seed"]
                for field, value in distinct[key].items():
                    if result[field] != value:
                        raise AssertionError(
                            (arm, key, field, "backend mismatch")
                        )
                pairs += 1
        rows.append(
            dict(
                backend=backend,
                arm=arm,
                scope=capture["scope"],
                median_seconds=statistics.median(times),
                ratio=statistics.median(times) / statistics.median(before),
                ratio_interval=ratio_interval(times, before),
                samples=len(times),
                relative_iqr=capture["relative_iqr"],
                stable=capture["relative_iqr"] <= 0.15
                and baseline["relative_iqr"] <= 0.15,
                successes=sum(
                    o["factor"] is not None for o in distinct.values()
                ),
                attempts=len(distinct),
                timeouts=sum(o["timed_out"] for o in distinct.values()),
                class_timings=class_timings(capture, baseline),
                by_case=[
                    dict(case=case, seed=seed, **outcome)
                    for (case, seed), outcome in distinct.items()
                ],
            )
        )
    return dict(source=str(path), backend_pairs=pairs, rows=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summaries = [summarize(path) for path in args.inputs]
    args.output.write_text(json.dumps(summaries, indent=2) + "\n")
    for report in summaries:
        print(report["source"], "backend pairs", report["backend_pairs"])
        for row in report["rows"]:
            print(
                row["backend"],
                row["arm"],
                row["median_seconds"],
                row["ratio"],
                row["ratio_interval"],
                row["stable"],
                row["successes"],
                row["attempts"],
            )


if __name__ == "__main__":
    main()
