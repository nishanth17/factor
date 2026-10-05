"""Validate cooperative stop latency with warmed worker pools."""

import argparse
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

from ..budget import Budget
from ..qs.parallel import CollectionPool, ParallelConfig, ParallelSIQSJob
from .performance_audit import CORPUS, fingerprint, measure, verify_corpus


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new capture path")
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 is required")
    corpus = json.loads(CORPUS.read_text())
    verify_corpus(corpus)
    fixture = next(
        f
        for f in corpus["fixtures"]
        if f["split"] == "training" and f["band"] == "30d"
    )
    root = Path(__file__).resolve().parents[2]
    before = fingerprint(root)
    config = ParallelConfig(
        base_bound=400,
        half_width=8192,
        factor_count=3,
        family_count=64,
        pool_size=16,
        batch_width=256,
        max_batch_atoms=256,
        assignment_work=50_000_000,
    )
    results = {}

    for mode, workers in (
        ("serial", 1),
        ("thread", 1),
        ("thread", 2),
        ("thread", 4),
        ("process", 1),
        ("process", 2),
        ("process", 4),
    ):
        with CollectionPool(mode, workers) as pool:
            for resource in ("wall", "cpu"):

                def run():
                    budget = Budget(
                        work_limit=10**10,
                        seconds=0.02 if resource == "wall" else 5,
                        cpu_seconds=0.02 if resource == "cpu" else 5,
                    )
                    start = time.perf_counter()
                    job = ParallelSIQSJob(
                        fixture["n"], seed=7, config=config, budget=budget
                    )

                    result = job.run(pool=pool, fixed_work=True)
                    elapsed = time.perf_counter() - start
                    if (result.divisor or 1) * result.cofactor != fixture["n"]:
                        raise AssertionError("limit result reconstruction")
                    if result.reason != resource + "_limit":
                        raise AssertionError("unexpected limit result")
                    if budget.used > budget.work_limit or pool.lock.locked():
                        raise AssertionError("undrained pool or work overdraw")
                    return [
                        dict(
                            seconds=elapsed,
                            cpu_seconds=budget.cpu_used,
                            work_used=budget.used,
                            reason=result.reason,
                            stats=result.stats,
                        )
                    ]

                results[f"{mode}_{workers}/{resource}"] = measure(
                    run, SimpleNamespace(warmup_seconds=3, repetitions=9)
                )

    if fingerprint(root) != before:
        raise AssertionError("source changed during measurement")
    args.output.write_text(
        json.dumps(
            dict(
                source_sha256=before,
                fixture=fixture["id"],
                seconds=0.02,
                scope="Warmed cooperative limits, excluding pool startup.",
                results=results,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
