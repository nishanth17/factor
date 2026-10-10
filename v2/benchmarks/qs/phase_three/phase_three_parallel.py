"""P3.6 fixed-family throughput and complete-factor parallelism experiments."""

import argparse
import hashlib
import json
import platform
import random
import resource
import statistics
import sys
import time
from dataclasses import asdict, replace
from math import prod
from pathlib import Path

from ....common import utils
from ....execution.budget import Budget, BudgetExhaustedError
from ....qs.parallel import (
    CollectionPool,
    ParallelConfig,
    ParallelSIQSJob,
    _rss,
)
from ....qs.sieve_collector import SieveConfig
from ....qs.siqs import SIQSConfig, SIQSJob
from ...support.paths import (
    BENCHMARK_ROOT,
    PACKAGE_ROOT,
    source_path,
)
from .phase_three_sss import sources

CORPUS = BENCHMARK_ROOT / "inputs/corpora/phase_three_p36_corpus.json"
SEEDS = (7, 29)
WORK = 2_000_000_000
SECONDS = 5
RSS = 1024 * 2**20
ARMS = (
    ("serial", 1),
    ("thread", 2),
    ("thread", 4),
    ("process", 1),
    ("process", 2),
    ("process", 4),
)


def configuration(band, bound=None):
    return replace(
        ParallelConfig(),
        base_bound=bound or (200 if band == "small" else 1000),
        half_width=256 if band == "small" else 512,
        factor_count=1 if band == "small" else 3,
        family_count=4,
        pool_size=16,
        max_batch_atoms=512,
        assignment_work=50_000_000,
    )


def reference_config(config):
    return SIQSConfig(
        base_bound=config.base_bound,
        half_width=config.half_width,
        factor_count=config.factor_count,
        family_count=config.family_count,
        pool_size=config.pool_size,
        max_stalled=64,
        max_trivial=4096,
        memory_bytes=config.parent_memory_bytes,
        collector=config.collector,
    )


def process_cpu():
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


def run_one(
    fixture, seed, config, *, pool=None, native=False, fixed=False, cold=False
):
    """Include setup, transfer, verification, extraction and classification."""
    started, cpu_started = time.perf_counter(), process_cpu()
    budget = Budget(work_limit=WORK, seconds=SECONDS, cpu_seconds=SECONDS)
    own_pool = None
    if cold:
        own_pool = CollectionPool(pool.mode, pool.workers)
        pool = own_pool
    try:
        if native:
            job = SIQSJob(
                fixture["n"],
                seed=seed,
                config=reference_config(config),
                budget=budget,
            )

            result = job.run()
        else:
            job = ParallelSIQSJob(
                fixture["n"], seed=seed, config=config, budget=budget
            )

            result = job.run(pool=pool, fixed_work=fixed)

        remaining, factors, labels = [fixture["n"]], [], []
        reason = result.reason
        if result.divisor is not None:
            if not utils.valid_divisor(result.divisor, fixture["n"]):
                raise AssertionError("invalid proper split")
            children = sorted((result.divisor, result.cofactor))

            try:
                for child in children:
                    budget.consume(child.bit_length() * 32)
                    label = utils.classify_prime(
                        child, rng=random.Random(seed)
                    )

                    if label is utils.Primality.COMPOSITE:
                        raise AssertionError("composite terminal factor")
                    labels.append(label.value)

                factors, remaining = children, []
            except BudgetExhaustedError:
                labels = []
                reason = "classification_" + budget.reason

        if prod(factors) * prod(remaining) != fixture["n"]:
            raise AssertionError("unfinished result does not reconstruct")
        if factors and factors != fixture["factors"]:
            raise AssertionError("split disagrees with independent corpus")
    finally:
        if own_pool is not None:
            own_pool.close()

    peak = _rss() + (
        sum(pool.state["rss"]) if pool and pool.mode == "process" else 0
    )
    if peak > RSS:
        raise MemoryError(
            "aggregate process high-water bound exceeds RSS gate"
        )
    return dict(
        id=fixture["id"],
        seed=seed,
        completed=not remaining,
        fixed_work_complete=not fixed or result.stats["schedule_complete"],
        factors=factors,
        remaining=remaining,
        certainty=labels,
        reason=reason,
        seconds=time.perf_counter() - started,
        cpu_seconds=budget.cpu_used,
        observed_cpu_seconds=process_cpu() - cpu_started
        if cold
        else budget.cpu_used,
        work=budget.used,
        rss_high_water_sum=peak,
        stats=result.stats,
    )


def measure(
    fixtures,
    config,
    mode,
    workers,
    args,
    *,
    fixed=False,
    native=False,
    cold=False,
):
    """Validate every warmup/sample, extending noisy measurements."""
    attempts = []
    with CollectionPool(mode, workers) as pool:

        def cohort():
            return [
                run_one(
                    f,
                    seed,
                    config,
                    pool=pool,
                    native=native,
                    fixed=fixed,
                    cold=cold,
                )
                for f in fixtures
                for seed in SEEDS
            ]

        for attempt in range(3):
            warm_start, calls = time.perf_counter(), 0
            while time.perf_counter() - warm_start < max(
                3, args.warmup_seconds, 5 if attempt else 3
            ):
                cohort()
                calls += 1

            warm_seconds = time.perf_counter() - warm_start
            samples = []

            for _ in range(max(9 if attempt == 0 else 15, args.repetitions)):
                started = time.perf_counter()
                rows = cohort()
                samples.append(
                    dict(seconds=time.perf_counter() - started, rows=rows)
                )

            times = [s["seconds"] for s in samples]
            median = statistics.median(times)
            q1, _, q3 = statistics.quantiles(times, n=4)
            drift = abs(
                statistics.median(times[:3]) / statistics.median(times[-3:])
                - 1
            )
            relative_iqr = (q3 - q1) / median
            stable = drift <= 0.15 and relative_iqr <= 0.2
            attempts.append(
                dict(
                    warmup_seconds=warm_seconds,
                    warmup_calls=calls,
                    samples=samples,
                    median_seconds=median,
                    relative_iqr=relative_iqr,
                    drift=drift,
                    stable=stable,
                )
            )
            if stable:
                break

    return dict(
        attempts=attempts,
        median_seconds=median,
        stable=stable,
        completion=[sum(r["completed"] for r in s["rows"]) for s in samples],
        work_completion=[
            sum(r["fixed_work_complete"] for r in s["rows"]) for s in samples
        ],
        rss_observation=(
            "Sum of parent/worker process-lifetime high-water marks; "
            "conservative bound, not simultaneous sampled peak"
        ),
        scope="cold job/pool startup and shutdown"
        if cold
        else "warmed reused executor",
    )


def ratio_interval(reference, challenger):
    """Bootstrap sample ratios; small-cohort uncertainty remains."""
    a = [s["seconds"] for s in reference["attempts"][-1]["samples"]]
    b = [s["seconds"] for s in challenger["attempts"][-1]["samples"]]
    size = min(len(a), len(b))
    generator = random.Random(3604)
    ratios = []

    for _ in range(2000):
        indices = [generator.randrange(size) for _ in range(size)]
        ratios.append(
            statistics.median(b[i] for i in indices)
            / statistics.median(a[i] for i in indices)
        )

    ratios.sort()
    return [ratios[50], ratios[1949]]


def main():
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("use PyPy implementing Python 3.11")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase", choices=("train", "warm", "cold", "all"), default="all"
    )
    parser.add_argument(
        "--band", choices=("small", "medium", "all"), default="all"
    )
    parser.add_argument("--training", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    args = parser.parse_args()
    corpus = json.loads(source_path(CORPUS).read_text())

    for fixture in corpus["fixtures"]:
        if prod(fixture["factors"]) != fixture["n"] or any(
            utils.classify_prime(p) != utils.Primality.PROVEN
            for p in fixture["factors"]
        ):
            raise AssertionError("invalid independent corpus")

    before = sources()
    before["qs/parallel.py"] = hashlib.sha256(
        source_path(PACKAGE_ROOT.joinpath("qs/parallel.py")).read_bytes()
    ).hexdigest()
    before["benchmarks/phase_three_parallel.py"] = hashlib.sha256(
        source_path(Path(__file__)).read_bytes()
    ).hexdigest()
    output = dict(
        schema=1,
        runtime=platform.python_implementation()
        + " "
        + platform.python_version(),
        source_sha256=before,
        corpus_sha256=hashlib.sha256(
            source_path(CORPUS).read_bytes()
        ).hexdigest(),
        seeds=SEEDS,
        limits=dict(
            work=WORK,
            wall=SECONDS,
            aggregate_cpu=SECONDS,
            owned_memory=512 * 2**20,
            observed_rss=RSS,
        ),
        gil="PyPy 3.11 GIL; no native backend release assumed",
        results={},
        training={},
    )
    if args.training:
        loaded = json.loads(source_path(args.training).read_text())
        if loaded["corpus_sha256"] != output["corpus_sha256"]:
            raise ValueError("training corpus mismatch")
        output["training"] = loaded["training"]

    bands = ("small", "medium") if args.band == "all" else (args.band,)

    def save():
        args.output.write_text(json.dumps(output, indent=2) + "\n")

    for band in bands:
        if args.phase in ("train", "all"):
            fixtures = [
                f
                for f in corpus["fixtures"]
                if f["band"] == band and f["split"] == "training"
            ]
            candidates = {}
            for bound in (200, 400) if band == "small" else (400, 1000):
                print(f"training {band} bound={bound}", flush=True)
                candidates[str(bound)] = measure(
                    fixtures, configuration(band, bound), "serial", 1, args
                )

            selected = min(
                candidates,
                key=lambda b: (
                    -min(candidates[b]["completion"]),
                    candidates[b]["median_seconds"],
                ),
            )
            output["training"][band] = dict(
                config=asdict(configuration(band, int(selected))),
                candidates=candidates,
                policy=(
                    "Training completion then complete-time median; "
                    "frozen across execution modes"
                ),
            )
            save()

        if args.phase == "train":
            continue
        if band not in output["training"]:
            raise ValueError("provide --training before held-out evaluation")
        options = dict(output["training"][band]["config"])
        options["collector"] = SieveConfig(**options["collector"])
        config = ParallelConfig(**options)
        fixtures = [
            f
            for f in corpus["fixtures"]
            if f["band"] == band and f["split"] == "held_out"
        ]
        output["results"][band] = {}
        if args.phase in ("warm", "all"):
            for fixed in (False, True):
                group = "fixed_work" if fixed else "first_factor"
                comparisons = {}
                if not fixed:
                    print(f"{band} native serial first factor", flush=True)
                    comparisons["native_serial"] = measure(
                        fixtures, config, "serial", 1, args, native=True
                    )

                for mode, workers in ARMS:
                    print(f"{band} {group} {mode}/{workers}", flush=True)
                    comparisons[f"{mode}_{workers}"] = measure(
                        fixtures, config, mode, workers, args, fixed=fixed
                    )
                    output["results"][band][group] = comparisons
                    save()

                baseline = comparisons[
                    "serial_1" if fixed else "native_serial"
                ]
                for name, capture in comparisons.items():
                    capture["ratio_to_baseline"] = (
                        capture["median_seconds"] / baseline["median_seconds"]
                    )
                    capture["ratio_ci95"] = ratio_interval(baseline, capture)

                save()

        if args.phase in ("cold", "all"):
            comparisons = {}

            for mode, workers in ARMS:
                print(f"{band} cold {mode}/{workers}", flush=True)
                comparisons[f"{mode}_{workers}"] = measure(
                    fixtures, config, mode, workers, args, cold=True
                )
                output["results"][band]["cold_first_factor"] = comparisons
                save()

    after = sources()
    after["qs/parallel.py"] = hashlib.sha256(
        source_path(PACKAGE_ROOT.joinpath("qs/parallel.py")).read_bytes()
    ).hexdigest()
    after["benchmarks/phase_three_parallel.py"] = hashlib.sha256(
        source_path(Path(__file__)).read_bytes()
    ).hexdigest()
    output["source_changed"] = before != after
    if output["source_changed"]:
        raise RuntimeError("experiment sources changed during execution")
    output["promotion"] = (
        "Experimental evaluator; throughput alone cannot promote."
    )
    save()


if __name__ == "__main__":
    main()
