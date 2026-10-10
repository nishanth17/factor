"""Separate cold-start records and prespecified paired cost summaries."""

import json
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path

from .b1_calibration import save


def cold_start(fixture, arm, selected, corpus, output, seconds):
    """Parent holds the performance lease; child includes fresh interpreter."""
    path = output / f"cold-{fixture['digits']}-{arm}.json"
    command = [
        sys.executable,
        "-B",
        "-m",
        "v2.benchmarks.c1_implementation",
        "cold-worker",
        "--output",
        str(path),
        "--selected",
        str(selected),
        "--corpus",
        str(corpus),
        "--fixture",
        fixture["id"],
        "--arm",
        arm,
    ]
    started = time.perf_counter()
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=seconds + 30,
        check=True,
    )
    if result.stdout or result.stderr:
        raise RuntimeError("unexpected cold worker output")
    save(
        path.with_name(path.stem + "-startup.json"),
        dict(
            process_seconds=time.perf_counter() - started,
            worker=str(path.name),
            classification="cold diagnostic; not a warmed sample",
        ),
    )


def interval(values):
    values = sorted(values)
    return [values[int((len(values) - 1) * p)] for p in (0.025, 0.975)]


def summarize(directory):
    """Resample inputs, then paired sample indices; retain capped failures."""
    groups = {}
    for path in sorted(Path(directory).glob("c1_fresh_*.json")):
        if "warmup" in path.name or "stability" in path.name:
            continue
        row = json.loads(path.read_text())
        if "sample" not in row:
            continue
        groups.setdefault(row["digits"], {}).setdefault(
            row["id"], {}
        ).setdefault(row["arm"], []).append(row)
    report = {}
    for digits, cases in groups.items():
        band = {}
        for arm in sorted(
            {a for case in cases.values() for a in case if a != "slp"}
        ):
            paired, per_input = [], {}
            for identity, arms in cases.items():
                if arm not in arms or "slp" not in arms:
                    continue
                reference = {r["sample"]: r for r in arms["slp"]}
                candidate = {r["sample"]: r for r in arms[arm]}
                indices = sorted(reference.keys() & candidate.keys())
                pairs = [(reference[i], candidate[i]) for i in indices]
                paired.append(pairs)
                stability_path = Path(directory) / f"{identity}-stability.json"
                stability = (
                    json.loads(stability_path.read_text())
                    if stability_path.exists()
                    else None
                )
                per_input[identity] = dict(
                    paired_count=len(pairs),
                    stability=stability,
                    arms={
                        name: dict(
                            complete=sum(r["complete"] for r in rows),
                            count=len(rows),
                            median_seconds=statistics.median(
                                r["seconds"] for r in rows
                            ),
                            capped_mean_seconds=statistics.mean(
                                r["capped_seconds"] for r in rows
                            ),
                            reasons=sorted({r["reason"] for r in rows}),
                            max_workspace=max(
                                r["stats"].get("workspace_bytes", 0)
                                for r in rows
                            ),
                        )
                        for name, rows in (
                            ("slp", [p[0] for p in pairs]),
                            (arm, [p[1] for p in pairs]),
                        )
                    },
                )
            generator = random.Random(380454 + digits)
            effects, completions = [], []
            for _ in range(10000):
                old_cost = new_cost = old_complete = new_complete = count = 0
                for _ in paired:
                    pairs = generator.choice(paired)
                    for _ in pairs:
                        old, new = generator.choice(pairs)
                        old_cost += old["capped_seconds"]
                        new_cost += new["capped_seconds"]
                        old_complete += old["complete"]
                        new_complete += new["complete"]
                        count += 1
                effects.append(1 - new_cost / old_cost)
                completions.append((new_complete - old_complete) / count)
            flat = [p for pairs in paired for p in pairs]
            band[arm] = dict(
                inputs=per_input,
                capped_cost_benefit=1
                - sum(p[1]["capped_seconds"] for p in flat)
                / sum(p[0]["capped_seconds"] for p in flat),
                benefit_95_interval=interval(effects),
                completion_gain_95_interval=interval(completions),
                limitation=(
                    "Two independently generated inputs per band; paired "
                    "repeated runs do not enlarge that input population."
                ),
            )
        report[str(digits)] = band
    return report
