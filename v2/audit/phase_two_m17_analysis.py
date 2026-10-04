"""Summarize paired confirmation without counting seeds as new inputs."""

import argparse
import gzip
import hashlib
import json
import random
import statistics
from collections import defaultdict
from pathlib import Path


def read_capture(path):
    """Read a JSON capture, retaining the hash of its original bytes."""
    data = path.read_bytes()
    if path.suffix == ".gz":
        data = gzip.decompress(data)
    return json.loads(data), hashlib.sha256(data).hexdigest()


def interval(values):
    """Return the empirical central 95% bootstrap interval."""
    ordered = sorted(values)
    return [ordered[len(ordered) // 40], ordered[-len(ordered) // 40 - 1]]


def compare_cell(rows, resamples, generator):
    """Pair identical input/seed/repetition assignments within one mode.

    Completion includes all outcomes. Successful timing ratios use only
    pairs where both engines succeeded, and are conditional diagnostics.
    Bootstrap completion by input, preserving within-input dependence.
    """
    assignments = defaultdict(dict)
    for row in rows:
        key = row["id"], row["seed"], row["repetition"]
        if row["engine"] in assignments[key]:
            raise AssertionError("duplicate assignment")
        assignments[key][row["engine"]] = row
    clusters = defaultdict(list)
    successful_ratios = []
    for (identity, _, _), pair in assignments.items():
        if set(pair) != {"m12", "bounded"}:
            raise AssertionError("unmatched confirmation assignment")
        old, new = pair["m12"], pair["bounded"]
        clusters[identity].append(int(new["success"]) - int(old["success"]))
        if old["success"] and new["success"]:
            successful_ratios.append(
                old["total_seconds"] / new["total_seconds"]
            )
    differences = [statistics.mean(values) for values in clusters.values()]
    bootstrap = [
        100 * statistics.mean(generator.choices(differences, k=len(clusters)))
        for _ in range(resamples)
    ]
    result = {
        "independent_inputs": len(clusters),
        "paired_assignments": len(assignments),
        "completion_delta_percentage_points": 100
        * statistics.mean(differences),
        "completion_delta_95pct_input_bootstrap": interval(bootstrap),
        "both_successful_pairs": len(successful_ratios),
        "conditional_successful_pair_median_time_ratio": statistics.median(
            successful_ratios
        )
        if successful_ratios
        else None,
    }
    for engine in ("m12", "bounded"):
        selected = [row for row in rows if row["engine"] == engine]
        medians = [
            statistics.median(
                row["total_seconds"]
                for row in selected
                if row["repetition"] == repetition
            )
            for repetition in sorted({row["repetition"] for row in selected})
        ]
        result[engine] = {
            "completion_fraction": statistics.mean(
                row["success"] for row in selected
            ),
            "median_seconds_all_outcomes": statistics.median(
                row["total_seconds"] for row in selected
            ),
            "repetition_median_seconds": medians,
            "peak_rss_bytes": max(row["peak_rss_bytes"] for row in selected),
            "cpu_seconds": sum(row["cpu_seconds"] for row in selected),
            "stop_reasons": {
                reason: sum(row["reason"] == reason for row in selected)
                for reason in sorted({row["reason"] for row in selected})
            },
        }
    return result


def analyze(path, resamples=2000, seed=20261003):
    """Validate accounting and summarize warm comparisons and cold costs."""
    capture, digest = read_capture(path)
    samples = capture["samples"]
    failures = sum(not row["reconstructs"] for row in samples)
    memory_failures = sum(row["memory_exceeded"] for row in samples)
    if failures or memory_failures:
        raise AssertionError("confirmation failed correctness/memory gates")
    cells = defaultdict(list)
    for row in samples:
        if row["temperature"] == "warm":
            cells[row["band"], row["mode"]].append(row)
    generator = random.Random(seed)
    comparisons = []
    for (band, mode), rows in sorted(cells.items()):
        comparisons.append(
            {
                "band": band,
                "mode": mode,
                **compare_cell(rows, resamples, generator),
            }
        )
    return {
        "milestone": "M17",
        "capture": str(path),
        "original_capture_sha256": digest,
        "corpus": capture["corpus"],
        "corpus_sha256": capture["corpus_sha256"],
        "sample_count": len(samples),
        "reconstruction_failures": failures,
        "memory_failures": memory_failures,
        "bootstrap_resamples": resamples,
        "bootstrap_seed": seed,
        "comparisons": comparisons,
        "cold_summaries": [
            row for row in capture["summaries"] if row["temperature"] == "cold"
        ],
        "warmups": capture["warmups"],
        "warmup_scope": capture.get("warmup_scope", "representatives"),
        "limits": {
            key: capture[key]
            for key in (
                "wall_cap_seconds",
                "cpu_cap_seconds",
                "rss_cap_bytes",
                "work_cap_bounded_only",
            )
        },
        "limitations": [
            "Three fresh bands; no claim of broad portfolio leadership",
            "Median costs include censored/unfinished outcomes",
            "Successful timing ratios are conditional on shared success",
            "Cold observations cover only one input per band and five seeds",
            "RSS includes interpreter/JIT and is measured, not OS enforced",
        ],
    }


def main():
    """Write a fresh, compact decision artifact from preserved raw samples."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("select a fresh --output filename")
    result = analyze(args.capture)
    with args.output.open("x") as output:
        output.write(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
