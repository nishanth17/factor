"""Matched finite candidate throughput against the frozen M12 batch loops."""

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from ...common import utils
from ...execution import stage_jobs
from ...execution.budget import Budget
from ...execution.schedules import SieveContext
from ...portfolio import PortfolioConfig
from ..support.snapshot_loader import load_stage_jobs
from .phase_one import _measure_case, environment


def run(args):
    """Consume identical candidate assignments and verify final state/work."""
    measured_environment = environment()
    baseline = load_stage_jobs()
    config = PortfolioConfig(
        max_input_bits=512,
        pm1_b1=50,
        pm1_b2=10000,
        ecm_tiers=((50, 10000, 1),),
    )
    n = 1000000000039 * 1000000000061
    seed = 20261003

    def candidate(action, kind):
        """Finish a fixed candidate, consuming its schedule and output."""
        job = stage_jobs.new_job(kind, n, seed, 50, 10000)
        context = SieveContext(config.max_hi, segment_size=config.segment_size)
        budget = Budget(seconds=None, cpu_seconds=None)
        while not job["done"]:
            action(job, budget, context, config)
        if job["factor"] is not None and not utils.valid_divisor(
            job["factor"], n
        ):
            raise AssertionError("invalid candidate factor")
        encoded = json.dumps(job, sort_keys=True, separators=(",", ":"))
        return {
            "factor": job["factor"],
            "work": budget.used,
            "state_sha256": hashlib.sha256(encoded.encode()).hexdigest(),
        }

    rows = []

    for kind in ("rho", "pm1", "ecm"):
        expected = candidate(baseline.advance_job, kind)
        row = _measure_case(
            f"fixed_{kind}_candidate",
            {
                "v2_native": lambda: candidate(stage_jobs.advance_job, kind),
                "m12": lambda: candidate(baseline.advance_job, kind),
            },
            lambda answer: answer == expected,
            args.repetitions,
            args.warmup_seconds,
        )
        row["expected"] = expected
        rows.append(row)
        print(kind, "done", flush=True)

    if environment()["source_sha256"] != measured_environment["source_sha256"]:
        raise RuntimeError("source changed during measurements")
    return {
        "milestone": "phase_two_m13_batch_controls",
        "environment": measured_environment,
        "config": asdict(config),
        "n": n,
        "seed": seed,
        "repetitions": args.repetitions,
        "warmup_seconds": args.warmup_seconds,
        "benchmarks": rows,
        "limitations": [
            "Fixed candidate throughput; exhaustion is an explicit outcome",
            "No exhausted candidate is a successful factorization baseline",
            "These utility ratios do not establish portfolio speedups",
        ],
    }


def main():
    """Save new controls without overwriting earlier measurements."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmup-seconds", type=float, default=10)
    parser.add_argument("--repetitions", type=int, default=15)
    args = parser.parse_args()
    if args.output.exists() or args.warmup_seconds < 3 or args.repetitions < 1:
        parser.error("use fresh output, >=3s warmup, and positive samples")
    args.output.write_text(json.dumps(run(args), indent=2) + "\n")


if __name__ == "__main__":
    main()
