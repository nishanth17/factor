"""Summarize paired C6 captures; reruns are not new factor trials."""

import argparse
import hashlib
import json
import statistics
from pathlib import Path

from . import c6_fast_study, p41_campaign


def summarize(capture, settings, shape=None, seed=None, metric="seconds"):
    pairs = []
    outcomes = {}
    for sample in capture["samples"]:
        times = {}
        for arm in ("ladder", "candidate"):
            rows = [
                row
                for row in sample[arm]["rows"]
                if (shape is None or row["case"].startswith(shape))
                and (seed is None or row["seed"] == seed)
            ]
            times[arm] = sum(row[metric] for row in rows)
            current = dict(
                attempts=len(rows),
                splits=sum(row["factor"] is not None for row in rows),
                timeouts=sum(row["timed_out"] for row in rows),
                unresolved=[
                    dict(
                        case=row["case"], seed=row["seed"], n=row["unresolved"]
                    )
                    for row in rows
                    if row["factor"] is None
                ],
                curves=sum(row["stats"]["curves"] for row in rows),
                stage_two_calls=sum(
                    row["stats"]["stage_two_calls"] for row in rows
                ),
                block_replays=sum(
                    row["extra"].get("block_replays", 0) for row in rows
                ),
            )
            if arm in outcomes and outcomes[arm] != current:
                raise ValueError("nonrepeatable validated outcomes")
            outcomes[arm] = current
        pairs.append(times)
    ratios = [p["candidate"] / p["ladder"] for p in pairs]
    middle = len(ratios) // 2
    return dict(
        median_ratio=statistics.median(ratios),
        ratio_95_interval=c6_fast_study.interval(ratios, settings),
        seconds={
            arm: statistics.median(p[arm] for p in pairs)
            for arm in ("ladder", "candidate")
        },
        outcomes=outcomes,
        chronological_ratio_medians=[
            statistics.median(ratios[:middle]),
            statistics.median(ratios[middle:]),
        ],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("captures", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {}
    corpora = {
        confirmation: {
            c["id"]: c for c in c6_fast_study.corpus(confirmation)["fixtures"]
        }
        for confirmation in (False, True)
    }
    for path in args.captures:
        data = json.loads(path.read_text())
        rows = []
        for capture in data["captures"]:
            for sample in capture["samples"]:
                for arm in ("ladder", "candidate"):
                    for result in sample[arm]["rows"]:
                        p41_campaign.validate_result(
                            result,
                            corpora[capture["confirmation"]][result["case"]],
                        )
            row = {
                key: capture[key]
                for key in ("arm", "backend", "scope", "confirmation")
            }
            row["samples"] = len(capture["samples"])
            row["preparation_seconds"] = capture["preparation_seconds"]
            if capture["scope"].endswith("reuse"):
                # Preparation includes both programs. Charging all of it to
                # the candidate is a conservative first-cohort estimate;
                # warmed execution plus this cost is not a cold process run.
                charged = [
                    (
                        sample["candidate"]["seconds"]
                        + capture["preparation_seconds"]
                    )
                    / sample["ladder"]["seconds"]
                    for sample in capture["samples"]
                ]
                row["first_cohort_preparation_charged"] = dict(
                    median_ratio=statistics.median(charged),
                    ratio_95_interval=c6_fast_study.interval(
                        charged, data["protocol"]
                    ),
                )
            row["stable"] = capture["summary"]["stable"]
            row["all"] = summarize(capture, data["protocol"])
            row["cpu"] = summarize(
                capture, data["protocol"], metric="cpu_seconds"
            )
            first_rows = capture["samples"][0]["candidate"]["rows"]
            row["seeds"] = {
                str(seed): summarize(capture, data["protocol"], seed=seed)
                for seed in sorted({r["seed"] for r in first_rows})
            }
            row["input_classes"] = {
                case: summarize(capture, data["protocol"], shape=case)
                for case in sorted({r["case"] for r in first_rows})
            }
            if capture["scope"].startswith("campaign"):
                row["classes"] = {
                    shape: summarize(capture, data["protocol"], shape)
                    for shape in ("balanced", "small10")
                }
            rows.append(row)
        report[path.name] = dict(
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(), rows=rows
        )
    with args.output.open("x") as output:
        json.dump(report, output, indent=2)
        output.write("\n")


if __name__ == "__main__":
    main()
