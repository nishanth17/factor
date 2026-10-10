"""Summarize matched captures with timing and paired-input uncertainty."""

import argparse
import itertools
import json
import random
import statistics
from collections import defaultdict
from pathlib import Path

from ...support.paths import (
    source_path,
)


def interval(values):
    """Central 95% interval; small discrete populations stay coarse."""
    values = sorted(values)
    return [
        values[int(0.025 * len(values))],
        values[min(len(values) - 1, int(0.975 * len(values)))],
    ]


def summarize(before, after):
    """Never infer a successful-factor time ratio from censored attempts."""
    old = before["attempts"][-1]["samples"]
    new = after["attempts"][-1]["samples"]

    def inputs(samples):
        rows = defaultdict(list)
        for sample in samples:
            for row in sample["rows"]:
                rows[row["id"], row["seed"]].append(row)
        return rows

    old_rows, new_rows = inputs(old), inputs(new)
    if old_rows.keys() != new_rows.keys():
        raise ValueError("captures do not contain the same inputs and seeds")
    identities = sorted({key[0] for key in old_rows})
    by_input = {}
    complete = True

    for identity in identities:
        pairs = []

        for source in (old_rows, new_rows):
            groups = [v for k, v in source.items() if k[0] == identity]
            complete &= all(r["completed"] for g in groups for r in g)
            pairs.append(
                dict(
                    seconds=sum(
                        statistics.median(r["seconds"] for r in g)
                        for g in groups
                    ),
                    completion=statistics.mean(
                        statistics.mean(int(r["completed"]) for r in g)
                        for g in groups
                    ),
                )
            )

        by_input[identity] = pairs

    rng = random.Random(361042026)
    cohorts = (
        itertools.product(identities, repeat=len(identities))
        if len(identities) <= 6
        else (rng.choices(identities, k=len(identities)) for _ in range(10000))
    )
    ratios, deltas = [], []

    for cohort in cohorts:
        deltas.append(
            statistics.mean(
                by_input[i][1]["completion"] - by_input[i][0]["completion"]
                for i in cohort
            )
        )
        if complete:
            ratios.append(
                sum(by_input[i][1]["seconds"] for i in cohort)
                / sum(by_input[i][0]["seconds"] for i in cohort)
            )

    timing_ratios = []
    if complete:
        for _ in range(10000):
            timing_ratios.append(
                statistics.median(
                    s["seconds"] for s in rng.choices(new, k=len(new))
                )
                / statistics.median(
                    s["seconds"] for s in rng.choices(old, k=len(old))
                )
            )
    return dict(
        independent_inputs=len(identities),
        seeds=sorted({k[1] for k in old_rows}),
        before_median_seconds=before["median_seconds"],
        after_median_seconds=after["median_seconds"],
        stable=before["stable"] and after["stable"],
        all_timed_attempts_complete=complete,
        complete_median_ratio=(
            after["median_seconds"] / before["median_seconds"]
            if complete
            else None
        ),
        timing_ratio_interval=interval(timing_ratios) if complete else None,
        paired_input_ratio_interval=interval(ratios) if complete else None,
        completion_delta_interval=interval(deltas),
        by_input=by_input,
        note="Timing resampling and input-cluster resampling are separate. "
        "Tiny constructed cohorts do not establish broad scaling.",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new summary path")
    old = json.loads(source_path(args.before).read_text())
    new = json.loads(source_path(args.after).read_text())
    if old["corpus_sha256"] != new["corpus_sha256"]:
        raise ValueError("corpus hashes differ")
    results = {
        key: summarize(old["results"][key], new["results"][key])
        for key in old["results"].keys() & new["results"].keys()
    }
    args.output.write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
