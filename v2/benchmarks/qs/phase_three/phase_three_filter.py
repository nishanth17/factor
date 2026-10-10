"""Matched filtering control and checkpoint continuation after M31 repairs."""

import argparse
import hashlib
import json
import statistics
import time
from pathlib import Path
from unittest.mock import patch

from ....qs.linear_algebra import filter_matrix
from ...suites.phase_one import environment
from ...support.paths import (
    source_path,
)
from .phase_three_audit import SNAPSHOT, load_arm
from .phase_three_large import CORPUS, run_one
from .phase_three_siqs import _config


def stability(rows):
    """Recheck drift and spread after any extension."""
    times = [row["seconds"] for row in rows]
    quartiles = statistics.quantiles(times, n=4)
    drift = abs(
        statistics.median(times[:3]) / statistics.median(times[-3:]) - 1
    )
    spread = (quartiles[2] - quartiles[0]) / statistics.median(times)
    return dict(
        completed=sum(row["complete"] for row in rows),
        samples=len(rows),
        median_seconds=statistics.median(times),
        median_drift=drift,
        relative_iqr=spread,
        stable=drift <= 0.15 and spread <= 0.2,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    corpus = json.loads(source_path(CORPUS).read_text())
    frozen = json.loads(source_path(args.frozen).read_text())
    assert (
        hashlib.sha256(source_path(CORPUS).read_bytes()).hexdigest()
        == frozen["corpus_sha256"]
    )
    source = environment()
    previous = load_arm(
        "_m31_filter_reference"
    ).qs.linear_algebra.filter_matrix
    # The owned M26 function is identical to M31's prior filter at the
    # small control's row count. Its 4096-row cap is disclosed separately.
    baseline = json.loads(source_path(SNAPSHOT).read_text())
    old_source = baseline["source"]["v2/qs/linear_algebra.py"]
    fixture = next(
        row
        for row in corpus["fixtures"]
        if row.get("performance_representative") and row["digits"] == 30
    )
    config = _config(frozen["configs"]["30"])
    rows, warmups, repeats = [], [], {}
    journal_path = args.output.with_suffix(".jsonl")
    checkpoint_dir = args.output.with_suffix(".checkpoints")
    with journal_path.open("x") as journal:

        def record(row):
            journal.write(json.dumps(row) + "\n")
            journal.flush()
            print(
                row.get("control"),
                row.get("cap_seconds"),
                row.get("complete"),
                row.get("reason"),
                round(row.get("seconds", 0), 3),
                flush=True,
            )

        for label, function in (
            ("prior_filter", previous),
            ("incidence_filter", filter_matrix),
        ):
            with patch("v2.qs.pipeline.filter_matrix", function):
                began, calls = time.perf_counter(), 0
                while time.perf_counter() - began < 3:
                    run_one(fixture, 7, config, "siqs")
                    calls += 1
                warmups.append(
                    dict(
                        control=label,
                        seconds=time.perf_counter() - began,
                        calls=calls,
                    )
                )
                samples = []

                for index in range(9):
                    row = run_one(
                        fixture,
                        7,
                        config,
                        "siqs",
                        checkpoint_dir=checkpoint_dir,
                    )
                    row.update(control=label, sample=index)
                    record(row)
                    samples.append(row)
                    rows.append(row)

                if not stability(samples)["stable"]:
                    began = time.perf_counter()
                    while time.perf_counter() - began < 5:
                        run_one(fixture, 7, config, "siqs")
                    samples = []

                    for index in range(15):
                        row = run_one(
                            fixture,
                            7,
                            config,
                            "siqs",
                            checkpoint_dir=checkpoint_dir,
                        )
                        row.update(control=label, sample=index, extension=True)
                        record(row)
                        samples.append(row)
                        rows.append(row)

                repeats[label] = stability(samples)

        fixture = next(
            row
            for row in corpus["fixtures"]
            if row.get("performance_representative") and row["digits"] == 50
        )
        config = _config(frozen["configs"]["50"])
        encoded = source_path(args.checkpoint).read_bytes()
        checkpoint = json.loads(encoded)
        parent_sha = hashlib.sha256(encoded).hexdigest()

        for ceiling in (2400, 3600):
            row = run_one(
                fixture,
                7,
                config,
                "siqs",
                seconds=ceiling,
                checkpoint=checkpoint,
                checkpoint_dir=checkpoint_dir,
            )
            row.update(
                control="incidence_checkpoint_continuation",
                parent_checkpoint=parent_sha,
            )
            record(row)
            rows.append(row)
            if row["complete"] or row["reason"] not in (
                "wall_limit",
                "cpu_limit",
                "work_limit",
            ):
                break

            metadata = row["checkpoint"]
            encoded = (
                source_path(args.output.parent / metadata["path"])
            ).read_bytes()
            parent_sha = hashlib.sha256(encoded).hexdigest()
            assert parent_sha == metadata["sha256"]
            checkpoint = json.loads(encoded)

        result = dict(
            corpus_sha256=frozen["corpus_sha256"],
            frozen_sha256=hashlib.sha256(
                source_path(args.frozen).read_bytes()
            ).hexdigest(),
            parent_checkpoint_sha256=hashlib.sha256(
                source_path(args.checkpoint).read_bytes()
            ).hexdigest(),
            previous_filter_sha256=hashlib.sha256(
                old_source.encode()
            ).hexdigest(),
            previous_row_cap=4096,
            environment=source,
            warmups=warmups,
            rows=rows,
            repeat_stability=repeats,
            scope=(
                "Inspected fixed representative. Filtering-only warmed "
                "control; cumulative continuation changes code and is "
                "not fresh held-out promotion evidence."
            ),
        )
        assert environment()["source_sha256"] == source["source_sha256"]
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2)

    print("Filtering acceptance complete", flush=True)


if __name__ == "__main__":
    main()
