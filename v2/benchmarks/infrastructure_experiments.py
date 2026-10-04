"""Training-only phase-two cutoff, chunk, schedule, and sieve experiments."""

import argparse
import hashlib
import json
import resource
import sys
import time
from array import array
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

from v2 import ecm, prime_sieve, utils
from v2.budget import Budget
from v2.factor import factorize
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.schedules import ScheduleCache, SieveContext
from v2.stage_jobs import advance_job, new_job

from .phase_one import _measure_case, _reference_primes, environment
from .sieve_candidates import integer_bitset, presieved, wheel_thirty


def run(args):
    """Validate every warmed and timed output before recording samples."""
    measured_environment = environment()
    config = PortfolioConfig(
        max_input_bits=512,
        pm1_b1=10,
        pm1_b2=200,
        ecm_tiers=((50, 1000, 3),),
    )
    rows = []

    def measure(name, candidates, validate):
        """Retain raw samples without assigning ratios to incomplete work."""
        row = _measure_case(
            name, candidates, validate, args.repetitions, args.warmup_seconds
        )
        rows.append(row)
        print(name, "done", flush=True)

    def job(kind, chunk, seed):
        """Drive exactly one candidate, retaining work and explicit failure."""
        candidate = new_job(
            kind, 1000000000039 * 1000000000061, seed, 50, 1000
        )
        context = SieveContext(config.max_hi, segment_size=config.segment_size)
        budget = Budget(seconds=2, cpu_seconds=2)
        while not candidate["done"]:
            advance_job(
                candidate, budget, context, replace(config, chunk_size=chunk)
            )
        return candidate["factor"]

    for kind in ("pm1", "ecm"):
        measure(
            f"{kind}_chunk_candidates",
            {
                "v2_native": lambda kind=kind: job(kind, 16, 0),
                **{
                    f"chunk_{size}": lambda size=size, kind=kind: job(
                        kind, size, 0
                    )
                    for size in (1, 4, 8, 32, 64)
                },
            },
            lambda value: value is None
            or utils.valid_divisor(value, 1000000000039 * 1000000000061),
        )
    corpus_bytes = (
        Path(__file__).with_name("phase_two_corpus.json").read_bytes()
    )
    corpus = json.loads(corpus_bytes)
    fixtures = [
        fixture
        for fixture in corpus["fixtures"]
        if fixture["split"] == "training"
        and fixture["band"] in ("powers", "random_small", "close_small")
    ]
    seeds = corpus["seeds"]

    def portfolio_runs(candidate_config):
        """Consume result and checkpoint output under identical budgets."""
        return [
            factorize_bounded(
                fixture["n"],
                seed=seed,
                config=candidate_config,
                budget=Budget(
                    work_limit=2_000_000, seconds=0.2, cpu_seconds=0.2
                ),
            )
            for fixture in fixtures
            for seed in seeds
        ]

    def validate_runs(runs):
        """Validate successful factors independently; retain exhausted runs."""
        for run, fixture in zip(
            runs, (fixture for fixture in fixtures for _ in seeds)
        ):
            if run.result.reconstruct() != fixture["n"]:
                return False
            if run.result.complete and {
                factor.value: factor.exponent for factor in run.result.factors
            } != dict(fixture["factors"]):
                return False
        return len(runs) == len(fixtures) * len(seeds)

    candidates = {"v2_native": lambda: portfolio_runs(config)}
    candidates.update(
        {
            f"trial_{bound}": lambda bound=bound: portfolio_runs(
                replace(config, trial_bound=bound)
            )
            for bound in (1000, 5000, 100000)
        }
    )
    measure("complete_portfolio_trial_cutoffs", candidates, validate_runs)
    outcomes = {
        label: [
            {
                "complete": run.result.complete,
                "reason": run.reason,
                "work": run.work_used,
            }
            for run in function()
        ]
        for label, function in candidates.items()
    }
    rho_candidates = {"v2_native": lambda: portfolio_runs(config)}
    rho_candidates.update(
        {
            f"batch_{size}": lambda size=size: portfolio_runs(
                replace(config, rho_batch=size)
            )
            for size in (32, 128, 256)
        }
    )
    measure("complete_portfolio_rho_batches", rho_candidates, validate_runs)
    measure(
        "complete_portfolio_context_controls",
        {
            "v2_native": lambda: portfolio_runs(config),
            "rolling": lambda: portfolio_runs(replace(config, rolling=True)),
            "cached": lambda: portfolio_runs(
                replace(config, schedule_cache_bytes=65536)
            ),
        },
        validate_runs,
    )
    # ECM tier grids use common inputs/seeds and remain experimental.
    for b1, b2 in ((2000, 147396), (11000, 500000), (50000, 2000000)):
        n = 1009 * 1000000000039

        def curves(bound=b1, upper=b2):
            """Measure fixed candidates and validate every returned divisor."""
            return [
                ecm.factorize_ecm(
                    n, seed=seed, b1=bound, b2=upper, max_curves=1
                )
                for seed in seeds
            ]

        measure(
            f"ecm_tier_{b1}_{b2}",
            {"v2_native": curves},
            lambda values: len(values) == 5
            and all(
                value is None or utils.valid_divisor(value, n)
                for value in values
            ),
        )
    for hi in (1000, 10000, 200000):
        expected = _reference_primes(hi)
        plain = SieveContext(hi, segment_size=4096)
        rolling = SieveContext(hi, segment_size=4096, rolling=True)
        measure(
            f"prime_consumption_{hi}",
            {
                "v2_native": lambda hi=hi: prime_sieve.segmented_sieve(2, hi),
                "context": lambda hi=hi, plain=plain: list(
                    plain.primes(2, hi)
                ),
                "rolling": lambda hi=hi, rolling=rolling: list(
                    rolling.primes(2, hi)
                ),
                "wheel_30": lambda hi=hi: wheel_thirty(2, hi),
                "pre_sieve": lambda hi=hi: presieved(2, hi),
                "integer_bitset": lambda hi=hi: integer_bitset(2, hi),
                "packed_output": lambda hi=hi: array(
                    "Q", prime_sieve.segmented_sieve(2, hi)
                ),
            },
            lambda values, expected=expected: list(values) == expected,
        )
    expected = _reference_primes(10000)
    context = SieveContext(10000, segment_size=4096)
    cache = ScheduleCache(context, cache_bytes=100000)
    cache_path = Path(__file__).with_name(".phase_two_schedule_control.bin")
    if cache_path.exists():
        raise FileExistsError("temporary schedule path already exists")
    disk_started = time.perf_counter()
    packed = array("Q", context.primes(2, 10000)).tobytes()
    cache_path.write_bytes(packed)
    disk_creation_seconds = time.perf_counter() - disk_started
    disk_hash = hashlib.sha256(packed).hexdigest()

    def disk_schedule():
        """Include disk reads, checksums, decoding, and consumption."""
        data = cache_path.read_bytes()
        if hashlib.sha256(data).hexdigest() != disk_hash:
            raise ValueError("corrupt packed schedule")
        values = array("Q")
        values.frombytes(data)
        return list(values)

    try:
        measure(
            "schedule_reuse_10000",
            {
                "v2_native": lambda: list(context.primes(2, 10000)),
                "bounded_ram": lambda: list(cache.values(2, 10000)),
                "cold_ram": lambda: list(
                    ScheduleCache(
                        SieveContext(10000), cache_bytes=100000
                    ).values(2, 10000)
                ),
                "packed_disk": disk_schedule,
            },
            lambda values: values == expected,
        )
    finally:
        cache_path.unlink()
    # Full factoring control: useful primes are consumed by actual ECM stages.
    factor_cases = [1009 * 1013, 1009 * 1000003, 1019 * 10007]

    def factorizations(sieve):
        """Include setup/output costs and require complete answers."""
        with patch.object(prime_sieve, "segmented_sieve", sieve):
            return [
                factorize(
                    n, level=1, seed=7, ecm_b1=50, ecm_b2=10000, ecm_curves=16
                )
                for n in factor_cases
            ]

    measure(
        "complete_ecm_sieve_controls",
        {
            "v2_native": lambda: factorizations(prime_sieve.segmented_sieve),
            "pre_sieve": lambda: factorizations(presieved),
            "wheel_30": lambda: factorizations(wheel_thirty),
            "integer_bitset": lambda: factorizations(integer_bitset),
        },
        lambda answers: len(answers) == len(factor_cases)
        and all(
            answer.complete and answer.reconstruct() == n
            for answer, n in zip(answers, factor_cases)
        ),
    )
    if environment()["source_sha256"] != measured_environment["source_sha256"]:
        raise RuntimeError("source changed during measurements")
    return {
        "milestone": "phase_two_infrastructure_training",
        "environment": measured_environment,
        "config": asdict(config),
        "corpus_sha256": hashlib.sha256(corpus_bytes).hexdigest(),
        "fixture_ids": [fixture["id"] for fixture in fixtures],
        "seeds": seeds,
        "repetitions": args.repetitions,
        "warmup_seconds": args.warmup_seconds,
        "peak_rss_platform_units": resource.getrusage(
            resource.RUSAGE_SELF
        ).ru_maxrss,
        "rss_units": "bytes" if sys.platform == "darwin" else "KiB",
        "cutoff_outcomes": outcomes,
        "packed_disk_sha256": disk_hash,
        "packed_disk_byteorder": sys.byteorder,
        "disk_creation_seconds": disk_creation_seconds,
        "benchmarks": rows,
        "limitations": [
            "Training only; incomplete outputs receive no speed ratios",
            "Compare process RSS separately in the isolated runner",
            "ECM tier probes do not select production tiers",
            "Disk creation is outside reuse timing; no disk default",
            "No native C crossover or threading threshold is imported",
        ],
    }


def main():
    """Capture new experiment evidence without replacing earlier results."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--repetitions", type=int, default=9)
    args = parser.parse_args()
    if args.output.exists() or args.repetitions < 1 or args.warmup_seconds < 3:
        parser.error("use a fresh output, positive samples, and >=3s warmup")
    args.output.write_text(json.dumps(run(args), indent=2) + "\n")


if __name__ == "__main__":
    main()
