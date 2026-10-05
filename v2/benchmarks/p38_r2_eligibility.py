"""Count whole-interval prime/power eligibility for a future CRT scheduler."""

import argparse
import importlib
import json
import sys
from dataclasses import asdict
from pathlib import Path

from .p38_r2 import HERE, ROOT, WORK, digest, modules, selected_config
from .performance_audit import fingerprint, verify_corpus


def probe(mods, fixture, seed):
    siqs = mods["qs.siqs"]
    selected = selected_config(siqs, fixture["band"], "current")
    budget = mods["budget"].Budget(work_limit=WORK, seconds=30, cpu_seconds=30)
    job = siqs.SIQSJob(fixture["n"], seed=seed, config=selected, budget=budget)
    job._setup()
    polynomial, roots = job._next_step()
    worker = siqs.SieveCollector(
        polynomial,
        job.base,
        config=selected.collector,
        precomputed_roots=roots,
        budget=budget,
    )
    lo, hi = -selected.half_width, selected.half_width + 1
    worker._set_plan_interval(lo, hi)
    maximum = worker._plan_window[2]
    powers = importlib.import_module(siqs.__package__ + ".power_sieve")
    count = dict(
        primes=0,
        power_moduli=0,
        power_hits=0,
        fallback_marks=0,
        all_marks=0,
        all_hits=0,
    )
    for root in roots:
        count["primes"] += root.prime > hi - lo
        for modulus, residues, weight in powers.prime_power_roots(
            polynomial,
            root,
            maximum,
            1,
            budget,
            lo=lo,
            hi=hi,
        ):
            hits = sum(
                len(range(lo + (r - lo) % modulus, hi, modulus))
                for r in residues
            )
            count["all_marks"] += 1
            count["all_hits"] += hits
            if weight > 1:
                count["fallback_marks"] += 1
            elif modulus > max(root.prime, hi - lo):
                count["power_moduli"] += 1
                count["power_hits"] += hits
    return dict(
        id=fixture["id"],
        seed=seed,
        polynomial=asdict(polynomial),
        base_cardinality=len(job.base.entries),
        interval=hi - lo,
        block=selected.collector.block_width,
        counts=count,
        work=budget.used,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    if args.output.exists():
        parser.error("choose a new capture path")
    corpus_path = HERE / "p38_r2_training_corpus.json"
    corpus = json.loads(corpus_path.read_text())
    verify_corpus(corpus)
    before = fingerprint(ROOT)
    mods = modules(importlib.import_module("v2"))
    results = [
        probe(mods, fixture, seed)
        for fixture in corpus["fixtures"]
        for seed in corpus["seeds"]
    ]
    if fingerprint(ROOT) != before:
        raise RuntimeError("source changed during eligibility capture")
    args.output.write_text(
        json.dumps(
            dict(
                source_sha256=before,
                corpus_sha256=digest(corpus_path),
                driver_sha256=digest(Path(__file__)),
                results=results,
                scope="First polynomial per seeded workload. Counts are "
                "bounded eligibility diagnostics, without family-wide "
                "performance evidence.",
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
