"""P3.5 trained SSS/SSSf challengers versus matched single-prime SIQS."""

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
from math import isfinite, prod
from pathlib import Path

from ....common import utils
from ....execution.budget import Budget, BudgetExhaustedError
from ....qs import SieveConfig, SIQSConfig, SIQSJob
from ....qs.polynomial import a_target
from ....qs.sss import SSSConfig, SSSJob
from ...suites.build_phase_two_corpus import verify_certificates
from ...support.paths import (
    BENCHMARK_ROOT,
    PACKAGE_ROOT,
    source_path,
)

CORPUS = BENCHMARK_ROOT / "inputs/corpora/phase_three_p34_corpus.json"
MEMORY = 64 * 1024 * 1024
RSS_LIMIT = 512 * 1024 * 1024
WORK = 2_000_000_000


def sources():
    """Identify every loaded arithmetic and relation-engine source byte."""
    root = PACKAGE_ROOT
    paths = list((root / "qs").glob("*.py"))
    paths += [
        root / name for name in ("utils.py", "budget.py", "prime_sieve.py")
    ]
    paths.append(Path(__file__))
    return {
        str(p.relative_to(root)): hashlib.sha256(
            source_path(p).read_bytes()
        ).hexdigest()
        for p in sorted(paths)
    }


def rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def run_one(fixture, seed, config, seconds):
    """Time setup through factor classification; hide the known factors."""
    budget = Budget(work_limit=WORK, seconds=seconds, cpu_seconds=seconds)
    started, cpu = time.perf_counter(), time.process_time()
    job_type = SSSJob if isinstance(config, SSSConfig) else SIQSJob
    job = job_type(fixture["n"], seed=seed, config=config, budget=budget)

    result = job.run()
    n, split = fixture["n"], result.divisor
    factors, remaining, labels = [], [n], []
    reason = result.reason
    if split is not None:
        if not utils.valid_divisor(split, n):
            raise AssertionError("invalid proper divisor")
        children = sorted((split, n // split))

        try:
            for child in children:
                budget.consume(child.bit_length() * 32)
                label = utils.classify_prime(child, rng=random.Random(seed))

                if label is utils.Primality.COMPOSITE:
                    raise AssertionError("semiprime child is composite")
                labels.append(label.value)

            factors, remaining = children, []
        except BudgetExhaustedError:
            labels = []
            reason = "classification_" + budget.reason

    if prod(factors) * prod(remaining) != n:
        raise AssertionError("result does not reconstruct")
    if factors and factors != fixture["factors"]:
        raise AssertionError(
            "split differs from the independently certified factors"
        )
    return dict(
        id=fixture["id"],
        seed=seed,
        completed=not remaining,
        factors=factors,
        remaining=remaining,
        certainty=labels,
        reason=reason,
        seconds=time.perf_counter() - started,
        cpu_seconds=time.process_time() - cpu,
        work_used=budget.used,
        stats=result.stats,
    )


def measure(fixtures, seeds, config, seconds, args):
    """Validate warmups and all samples; extend noisy timing measurements."""

    def cohort():
        return [
            run_one(fixture, seed, config, seconds)
            for fixture in fixtures
            for seed in seeds
        ]

    attempts = []

    for attempt in range(3):
        started = time.perf_counter()
        warmup = max(args.warmup_seconds, 3 if attempt == 0 else 5)
        warm_calls = 0
        while time.perf_counter() - started < warmup:
            cohort()
            warm_calls += 1
        warm_elapsed = time.perf_counter() - started
        samples = []
        repetitions = max(args.repetitions, 9 if attempt == 0 else 15)

        for _ in range(repetitions):
            wall, cpu = time.perf_counter(), time.process_time()
            rows = cohort()
            peak = rss_bytes()
            if peak > RSS_LIMIT:
                raise MemoryError(
                    "comparison process exceeds the declared RSS limit"
                )
            samples.append(
                dict(
                    seconds=time.perf_counter() - wall,
                    cpu_seconds=time.process_time() - cpu,
                    peak_rss_bytes=peak,
                    rows=rows,
                )
            )

        times = [sample["seconds"] for sample in samples]
        median = statistics.median(times)
        q1, _, q3 = statistics.quantiles(times, n=4)
        drift = abs(
            statistics.median(times[:3]) / statistics.median(times[-3:]) - 1
        )
        iqr = (q3 - q1) / median
        stable = drift <= 0.15 and iqr <= 0.2
        attempts.append(
            dict(
                warmup_seconds=warm_elapsed,
                warmup_calls=warm_calls,
                samples=samples,
                median_seconds=median,
                drift=drift,
                relative_iqr=iqr,
                stable=stable,
            )
        )
        if stable:
            break

    return dict(
        attempts=attempts,
        stable=stable,
        median_seconds=median,
        completion_per_sample=[
            sum(row["completed"] for row in s["rows"]) for s in samples
        ],
    )


def collector_config():
    return SieveConfig(
        division="roots",
        residual_bound=10000,
        max_atoms=8192,
        max_relations=4096,
        max_partials=1024,
        memory_bytes=MEMORY,
    )


def train(corpus, args):
    fixtures = [
        f
        for f in corpus["fixtures"]
        if (f["split"] == "training" and f["band"] == "small")
    ]
    candidates = []

    for bound in (400, 1000):
        for size in (3, 6):
            config = SSSConfig(
                base_bound=bound,
                selection_size=size,
                search_rounds=128,
                collector=collector_config(),
            )
            result = measure(fixtures, (7,), config, 5, args)
            candidates.append(dict(config=asdict(config), measurement=result))
            print(
                "training", bound, size, result["median_seconds"], flush=True
            )

    eligible = [
        c
        for c in candidates
        if min(c["measurement"]["completion_per_sample"]) == len(fixtures)
    ]
    if not eligible:
        eligible = candidates
    best = min(
        eligible,
        key=lambda c: (
            -min(c["measurement"]["completion_per_sample"]),
            c["measurement"]["median_seconds"],
        ),
    )
    return dict(
        candidates=candidates, frozen=best["config"], held_out_used=False
    )


def deserialize_config(data):
    data = dict(data)
    data["collector"] = SieveConfig(**data["collector"])
    return SSSConfig(**data)


def band_bound(band):
    """Predeclare bounds from observable digit bands, without known factors."""
    return {
        "balanced_30d": 5000,
        "balanced_40d": 8000,
        "balanced_50d": 22000,
        "balanced_60d": 82000,
        "balanced_70d": 82000,
        "balanced_80d": 82000,
    }[band]


def train_siqs(corpus, band, args):
    """Freeze a control with training probes; elapsed times are diagnostics."""
    if band == "small":
        fixtures = [
            f
            for f in corpus["fixtures"]
            if f["split"] == "training" and f["band"] == band
        ]
        candidates = []

        for bound in (200, 400):
            for count in (1, 3):
                config = SIQSConfig(
                    base_bound=bound,
                    factor_count=count,
                    half_width=256,
                    max_half_width=512,
                    family_count=64,
                    pool_size=64,
                    memory_bytes=MEMORY,
                    collector=replace(
                        collector_config(),
                        score_policy="powers",
                        division="bucket",
                    ),
                )
                measurement = measure(fixtures, (7,), config, 5, args)
                candidates.append(
                    dict(config=asdict(config), measurement=measurement)
                )
                print(
                    "SIQS small training",
                    bound,
                    count,
                    measurement["median_seconds"],
                    flush=True,
                )

        best = min(
            candidates,
            key=lambda c: (
                -min(c["measurement"]["completion_per_sample"]),
                c["measurement"]["median_seconds"],
            ),
        )
        return dict(
            frozen=best["config"],
            candidates=candidates,
            band=band,
            source_sha256=sources(),
            policy="repeated small training; no held-out tuning",
        )

    fixture = next(
        f
        for f in corpus["fixtures"]
        if (f["split"] == "training" and f["band"] == band)
    )
    bound = band_bound(band)
    target = a_target(fixture["n"], 512)
    minimum = next((i for i in range(1, 9) if bound**i >= target), 8)
    counts = sorted({max(1, minimum - 1), minimum, min(8, minimum + 1)})
    rows = []

    for count in counts:
        for width in (512, 4096):
            config = SIQSConfig(
                base_bound=bound,
                factor_count=count,
                half_width=width,
                max_half_width=4096,
                family_count=64,
                pool_size=64,
                memory_bytes=MEMORY,
                collector=replace(
                    collector_config(),
                    score_policy="powers",
                    division="bucket",
                ),
            )
            row = run_one(fixture, 7, config, 5)
            rows.append(dict(config=asdict(config), row=row))
            print(
                "SIQS training",
                band,
                count,
                width,
                row["completed"],
                row["reason"],
                row["seconds"],
                flush=True,
            )

    best = min(
        rows,
        key=lambda r: (
            not r["row"]["completed"],
            r["row"]["seconds"]
            if r["row"]["completed"]
            else (-int(r["row"]["stats"].get("relations", 0))),
        ),
    )
    return dict(
        frozen=best["config"],
        rows=rows,
        band=band,
        source_sha256=sources(),
        policy="training probes; no repeated timing or held-out tuning",
    )


def confidence(control, challenger):
    """Bootstrap matched cohort times; report completion separately."""
    left = [s["seconds"] for s in control["attempts"][-1]["samples"]]
    right = [s["seconds"] for s in challenger["attempts"][-1]["samples"]]
    generator = random.Random(350035)
    ratios = sorted(
        statistics.median(generator.choices(right, k=len(right)))
        / statistics.median(generator.choices(left, k=len(left)))
        for _ in range(2000)
    )
    return dict(
        median_ratio=statistics.median(right) / statistics.median(left),
        ratio_95_interval=[ratios[50], ratios[1949]],
        includes_exhausted_outcomes=True,
    )


def main():
    if platform.python_implementation() != "PyPy" or sys.version_info[:2] != (
        3,
        11,
    ):
        raise RuntimeError("benchmark requires PyPy implementing Python 3.11")

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--training", type=Path)
    parser.add_argument("--siqs-training", type=Path)
    parser.add_argument("--arms", nargs="+", choices=("siqs", "sss", "sssf"))
    parser.add_argument("--train-only", action="store_true")
    parser.add_argument("--train-siqs-only", action="store_true")
    parser.add_argument("--diagnostic-only", action="store_true")
    parser.add_argument("--band", default="small")
    parser.add_argument("--large-seconds", type=float, default=5)
    parser.add_argument("--limit-inputs", type=int, default=0)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    args = parser.parse_args()
    if (
        args.repetitions < 9
        or args.warmup_seconds < 3
        or not isfinite(args.warmup_seconds)
        or args.large_seconds <= 0
        or not isfinite(args.large_seconds)
        or args.limit_inputs < 0
    ):
        parser.error("need >=9 repetitions, >=3s warmup and positive budgets")

    corpus = json.loads(source_path(CORPUS).read_text())
    verify_certificates(corpus["certificates"])
    for fixture in corpus["fixtures"]:
        if prod(fixture["factors"]) != fixture["n"]:
            raise ValueError("corpus factors do not reconstruct")
    manifest = sources()
    output = dict(
        schema=1,
        runtime=sys.version,
        implementation="PyPy",
        source_sha256=manifest,
        corpus_sha256=hashlib.sha256(
            source_path(CORPUS).read_bytes()
        ).hexdigest(),
        resources=dict(
            cores=1,
            work=WORK,
            owned_memory=MEMORY,
            observed_process_rss_limit=RSS_LIMIT,
        ),
        timing_scope="setup through reconstructed/classified results",
    )
    if args.train_siqs_only:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(
                train_siqs(corpus, args.band, args),
                indent=2,
            )
            + "\n"
        )
        return

    if args.training:
        training = json.loads(source_path(args.training).read_text())[
            "training"
        ]
    else:
        training = train(corpus, args)
    output["training"] = training
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    if args.train_only:
        return
    base = deserialize_config(training["frozen"])
    fixtures = [
        f
        for f in corpus["fixtures"]
        if (f["split"] == "held_out" and f["band"] == args.band)
    ]
    if args.limit_inputs:
        fixtures = fixtures[: args.limit_inputs]
    if not fixtures:
        raise ValueError("no inputs for the requested band")
    seconds = 5 if args.band == "small" else args.large_seconds
    if args.band != "small":
        base = replace(
            base, base_bound=band_bound(args.band), search_rounds=4096
        )
    configs = {
        "siqs": SIQSConfig(
            base_bound=base.base_bound,
            half_width=512,
            factor_count=3,
            family_count=64,
            pool_size=32,
            memory_bytes=MEMORY,
            collector=collector_config(),
        ),
        "sss": base,
        "sssf": replace(
            base,
            mode="sssf",
            selection_size=7,
            filter_divisor=2,
            filter_bound=10 ** max(1, len(str(fixtures[0]["n"])) // 2 - 1),
        ),
    }
    if args.siqs_training:
        control_training = json.loads(
            source_path(args.siqs_training).read_text()
        )
        if control_training.get("band", args.band) != args.band:
            raise ValueError("SIQS settings were trained for another band")
        data = dict(control_training["frozen"])
        data["collector"] = SieveConfig(**data["collector"])
        configs["siqs"] = SIQSConfig(**data)
        for name, digest in control_training["source_sha256"].items():
            if name.startswith("qs/") and manifest.get(name) != digest:
                raise ValueError("SIQS training and evaluation sources differ")
        output["siqs_training"] = control_training

    if args.arms:
        configs = {name: configs[name] for name in args.arms}
    output.update(
        band=args.band,
        ids=[f["id"] for f in fixtures],
        seeds=corpus["seeds"],
        seconds_per_input=seconds,
        configs={name: asdict(config) for name, config in configs.items()},
        comparisons={},
        measurement=(
            "single diagnostic cohort; no timing comparison"
            if args.diagnostic_only
            else "validated warmup and repeated cohorts"
        ),
    )
    for name, config in configs.items():
        if args.diagnostic_only:
            rows = [
                run_one(fixture, seed, config, seconds)
                for fixture in fixtures
                for seed in corpus["seeds"]
            ]
            peak = rss_bytes()
            if peak > RSS_LIMIT:
                raise MemoryError("diagnostic exceeds the declared RSS limit")
            result = dict(rows=rows, peak_rss_bytes=peak)
            output["comparisons"][name] = result
            args.output.write_text(json.dumps(output, indent=2) + "\n")
            print(
                args.band,
                name,
                "diagnostic complete",
                sum(row["completed"] for row in rows),
                "/",
                len(rows),
                flush=True,
            )
            continue

        result = measure(fixtures, corpus["seeds"], config, seconds, args)
        output["comparisons"][name] = result
        args.output.write_text(json.dumps(output, indent=2) + "\n")
        print(
            args.band,
            name,
            result["median_seconds"],
            result["completion_per_sample"],
            "stable",
            result["stable"],
            flush=True,
        )

    if "siqs" in output["comparisons"] and not args.diagnostic_only:
        output["relative_to_siqs"] = {
            name: confidence(
                output["comparisons"]["siqs"], output["comparisons"][name]
            )
            for name in ("sss", "sssf")
            if name in output["comparisons"]
        }

    output["sources_changed_during_run"] = sources() != manifest
    output["promotion"] = (
        "retain challenger; evaluate all declared classes before promotion"
    )
    args.output.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
