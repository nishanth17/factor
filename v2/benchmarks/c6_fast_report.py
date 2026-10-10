"""Summarize paired C6 captures; reruns are not new factor trials."""

import argparse
import hashlib
import json
import statistics
from pathlib import Path

from . import c6_fast_study, p41_campaign


def summarize(capture, settings, shape=None):
    pairs = []
    outcomes = {}
    for sample in capture["samples"]:
        times = {}
        for arm in ("ladder", "candidate"):
            rows = [
                row
                for row in sample[arm]["rows"]
                if shape is None or row["case"].startswith(shape)
            ]
            times[arm] = sum(row["seconds"] for row in rows)
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
    return dict(
        median_ratio=statistics.median(ratios),
        ratio_95_interval=c6_fast_study.interval(ratios, settings),
        seconds={
            arm: statistics.median(p[arm] for p in pairs)
            for arm in ("ladder", "candidate")
        },
        outcomes=outcomes,
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
            row["stable"] = capture["summary"]["stable"]
            row["all"] = summarize(capture, data["protocol"])
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
