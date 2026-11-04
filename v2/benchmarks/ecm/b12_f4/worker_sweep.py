"""Summarize fixed-curve worker screens without implying portfolio adoption."""

import argparse
import json
import statistics
from pathlib import Path


def summarize(directory):
    """Keep censoring and cold process lifetime beside warm call times."""
    directory = Path(directory)
    rows = []
    for path in sorted(directory.glob("*.receipt.json")):
        receipt = json.loads(path.read_text())
        cell = receipt["cell"]
        result = directory / (cell["id"] + ".capture") / "result.json"
        data = json.loads(result.read_text()) if result.exists() else None
        complete = (
            receipt["exit"] == 0
            and receipt["reason"] is None
            and receipt["owned_group_exited"]
            and data is not None
            and data["complete"]
        )
        row = dict(
            id=cell["id"],
            input=cell["input"],
            arm=cell["arm"],
            workers=cell["workers"],
            mode=cell["mode"],
            observer=cell["observer"],
            complete=complete,
            stop=receipt["reason"],
            cli_wall=receipt["wall"],
            cpu_observed=receipt["cpu"],
            cpu_charged=receipt["cpu_charge"],
            signature=receipt["signature"],
        )
        if complete:
            calls = data["calls"]
            row.update(
                samples=len(calls),
                median=statistics.median(calls),
                mean=statistics.mean(calls),
                minimum=min(calls),
                maximum=max(calls),
                cv=statistics.pstdev(calls) / statistics.mean(calls),
                warmup=data["warmup_seconds"],
                startup=data["stats"][0]["startup_seconds"],
                preparations=data["lifetime"]["preparations"],
                owned_bytes=data["stats"][0]["owned_bytes"],
                parent_cpu_per_call=statistics.mean(
                    item["parent_cpu"] for item in data["stats"]
                ),
                worker_cpu_per_call=statistics.mean(
                    item["reported_child_cpu"] for item in data["stats"]
                ),
                worker_hwm_sum=max(
                    item["worker_peak_rss_sum"] for item in data["stats"]
                ),
                sampled_rss_peak=max(
                    (
                        item["total"]
                        for item in data["rss"] or []
                        if item["complete"]
                    ),
                    default=None,
                ),
            )
        rows.append(row)
    return rows


def select(rows):
    """Select exploratory nominees separately for each size/work/pool state."""
    groups = {}
    for row in rows:
        if row["observer"] != "none" or not row["complete"]:
            continue
        groups.setdefault((row["input"], row["mode"]), []).append(row)
    nominees = []
    for (case, mode), arms in sorted(groups.items()):
        metric = "cli_wall" if mode == "cold" else "median"
        best = min(arms, key=lambda row: row[metric])
        nominees.append(
            dict(
                input=case,
                mode=mode,
                arm=best["arm"],
                workers=best["workers"],
                seconds=best[metric],
                metric=metric,
                scope="exposed fixed-work input; exploratory selection",
            )
        )
    return nominees


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = summarize(args.directory)
    with args.output.open("x") as stream:
        json.dump(dict(rows=rows, nominees=select(rows)), stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
