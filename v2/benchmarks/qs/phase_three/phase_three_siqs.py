"""Full P3.4 tuning, frozen held-out coverage, cap accounting and decisions."""

import argparse
import cProfile
import hashlib
import json
import math
import random
import resource
import statistics
import subprocess
import sys
import time
from dataclasses import asdict, replace
from math import prod
from pathlib import Path

from ....common import utils
from ....execution.budget import Budget, BudgetExhaustedError
from ....portfolio import PortfolioConfig, factorize_bounded
from ....qs import SieveConfig, SIQSConfig, SIQSJob
from ...suites.build_phase_two_corpus import verify_certificates
from ...suites.phase_one import environment
from ...support.paths import (
    BENCHMARK_ROOT,
    source_path,
)
from .phase_three_reference import _rss_bytes

CORPUS = BENCHMARK_ROOT / "inputs/corpora/phase_three_p34_corpus.json"
MEMORY = 64 * 1024 * 1024
SMALL_SECONDS = 10
LARGE_SECONDS = 0.2
WORK = 200_000_000


def base_config():
    return SIQSConfig(
        base_bound=200,
        half_width=256,
        family_count=8,
        memory_bytes=MEMORY,
        collector=SieveConfig(
            division="bucket",
            residual_bound=1,
            max_atoms=4096,
            max_relations=4096,
            max_partials=4096,
        ),
    )


def _validate(row, fixture):
    if prod(row["factors"]) * prod(row["remaining"]) != fixture["n"]:
        raise AssertionError("factoring result does not reconstruct")
    if any(value not in fixture["factors"] for value in row["factors"]):
        raise AssertionError("returned terminal factor differs from proof")
    if row["completed"] and sorted(row["factors"]) != fixture["factors"]:
        raise AssertionError("incorrect complete split")
    if row["split"] is not None and row["split"] not in fixture["factors"]:
        raise AssertionError("invalid proper divisor")


def run_one(fixture, seed, config, *, ecm=False, seconds=None):
    n = fixture["n"]
    seconds = (
        seconds
        if seconds is not None
        else (SMALL_SECONDS if fixture["band"] == "small" else LARGE_SECONDS)
    )
    budget = Budget(work_limit=WORK, seconds=seconds, cpu_seconds=seconds)
    started, cpu = time.perf_counter(), time.process_time()
    if ecm:
        portfolio = PortfolioConfig(
            trial_bound=2,
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=((2000, 147396, 2),),
            memory_bytes=MEMORY,
        )

        result = factorize_bounded(
            n, seed=seed, config=portfolio, budget=budget
        )
        factors = [
            f.value for f in result.result.factors for _ in range(f.exponent)
        ]
        remaining = list(result.result.remaining)
        labels = [f.certainty.value for f in result.result.factors]
        reason = result.reason
        stats = dict(
            events=result.events,
            stage_seconds=result.checkpoint["payload"]["state"][
                "stage_seconds"
            ],
        )
        values = factors + remaining
        split = next((v for v in values if utils.valid_divisor(v, n)), None)
    else:
        job = SIQSJob(n, seed=seed, config=config, budget=budget)

        result = job.run()
        split = result.divisor
        factors, remaining, labels = [], [n], []
        reason = result.reason
        stats = result.stats
        if split is not None:
            children = sorted((split, n // split))

            try:
                for child in children:
                    budget.consume(child.bit_length() * 32)
                    label = utils.classify_prime(
                        child, rng=random.Random(seed)
                    )

                    if label is utils.Primality.COMPOSITE:
                        raise AssertionError("semiprime child was composite")
                    labels.append(label.value)

                factors, remaining = children, []
            except BudgetExhaustedError:
                labels = []
                reason = "classification_" + budget.reason

    row = dict(
        id=fixture["id"],
        band=fixture["band"],
        seed=seed,
        reason=reason,
        completed=not remaining,
        split=split,
        factors=factors,
        remaining=remaining,
        certainty=labels,
        seconds=time.perf_counter() - started,
        cpu_seconds=time.process_time() - cpu,
        work_used=budget.used,
        stats=stats,
    )
    _validate(row, fixture)
    return row


def cohort(fixtures, seed, config, *, ecm=False, seconds=None):
    return [
        run_one(fixture, seed, config, ecm=ecm, seconds=seconds)
        for fixture in fixtures
    ]


def _measure(fixtures, seed, config, args, *, ecm=False, seconds=None):
    attempts = []

    for attempt in range(2):
        began = time.perf_counter()
        warmup = (
            args.warmup_seconds
            if attempt == 0
            else max(5, args.warmup_seconds)
        )
        warm_calls = 0
        while time.perf_counter() - began < warmup:
            cohort(fixtures, seed, config, ecm=ecm, seconds=seconds)
            warm_calls += 1
        warm_elapsed = time.perf_counter() - began
        samples = []

        for _ in range(
            args.repetitions if attempt == 0 else max(15, args.repetitions)
        ):
            began, cpu = time.perf_counter(), time.process_time()
            rows = cohort(fixtures, seed, config, ecm=ecm, seconds=seconds)
            samples.append(
                dict(
                    seconds=time.perf_counter() - began,
                    cpu_seconds=time.process_time() - cpu,
                    peak_rss_bytes=_rss_bytes(),
                    rows=rows,
                )
            )

        times = [s["seconds"] for s in samples]
        median = statistics.median(times)
        quartiles = statistics.quantiles(times, n=4)
        drift = abs(
            statistics.median(times[:3]) / statistics.median(times[-3:]) - 1
        )
        iqr = (quartiles[2] - quartiles[0]) / median
        stable = drift <= 0.15 and iqr <= 0.2
        attempts.append(
            dict(
                warmup_seconds=warm_elapsed,
                warmup_calls=warm_calls,
                samples=samples,
                median_seconds=median,
                stable=stable,
                drift=drift,
                relative_iqr=iqr,
            )
        )
        if stable:
            break

    return dict(
        attempts=attempts,
        median_seconds=median,
        stable=stable,
        timing_scope="All attempts; unfinished inputs are censored outcomes",
    )


def training(corpus, args):
    small = [
        f
        for f in corpus["fixtures"]
        if f["split"] == "training" and f["band"] == "small"
    ]
    large = [
        f
        for f in corpus["fixtures"]
        if f["split"] == "training" and f["band"] != "small"
    ]
    base = base_config()
    choices = [("baseline", base)]

    for field, values in (
        ("base_bound", (100, 500)),
        ("factor_count", (1, 2, 4)),
        ("half_width", (128, 512)),
        ("pool_size", (8, 32)),
        ("multiplier", (0,)),
        ("growth_steps", (1,)),
        ("diverse", (False,)),
    ):
        choices.extend(
            (f"{field}_{v}", replace(base, **{field: v})) for v in values
        )

    for field, values in (
        ("block_width", (128, 512)),
        ("residual_bound", (500, 10000)),
        ("threshold_extra", (1,)),
        ("score_policy", ("candidate",)),
    ):
        choices.extend(
            (
                f"{field}_{v}",
                replace(base, collector=replace(base.collector, **{field: v})),
            )
            for v in values
        )

    records = []

    for name, config in choices:
        measurement = _measure(small, 7, config, args)
        repetitions = [
            sample["rows"] for sample in measurement["attempts"][-1]["samples"]
        ]
        completed = sum(r["completed"] for rows in repetitions for r in rows)
        seconds = statistics.median(
            sum(r["seconds"] for r in rows) for rows in repetitions
        )
        records.append(
            dict(
                name=name,
                config=asdict(config),
                completed=completed,
                completion_rate=statistics.mean(
                    r["completed"] for rows in repetitions for r in rows
                ),
                median_seconds=seconds,
                repetitions=repetitions,
                measurement=measurement,
            )
        )

    selected = max(
        records, key=lambda r: (r["completion_rate"], -r["median_seconds"])
    )
    chosen = next(c for name, c in choices if name == selected["name"])
    large_choices = [
        ("base_1000_k3", replace(base, base_bound=1000, half_width=512)),
        (
            "base_2000_k4",
            replace(base, base_bound=2000, half_width=512, factor_count=4),
        ),
        (
            "base_2000_k6",
            replace(base, base_bound=2000, half_width=512, factor_count=6),
        ),
        (
            "base_2000_k8",
            replace(base, base_bound=2000, half_width=512, factor_count=8),
        ),
        (
            "scored",
            replace(
                base,
                base_bound=2000,
                half_width=512,
                factor_count=6,
                multiplier=0,
            ),
        ),
        (
            "recovery",
            replace(
                base,
                base_bound=2000,
                half_width=512,
                factor_count=6,
                growth_steps=1,
            ),
        ),
    ]
    large_records = []

    for name, config in large_choices:
        measurement = _measure(large, 7, config, args, seconds=0.08)
        rows = [
            row
            for sample in measurement["attempts"][-1]["samples"]
            for row in sample["rows"]
        ]
        large_records.append(
            dict(
                name=name,
                config=asdict(config),
                rows=rows,
                measurement=measurement,
                completed=sum(r["completed"] for r in rows),
                completion_rate=statistics.mean(r["completed"] for r in rows),
                mean_relations=statistics.mean(
                    r["stats"].get("relations", 0) for r in rows
                ),
                relations=sum(r["stats"].get("relations", 0) for r in rows),
                seconds=sum(r["seconds"] for r in rows),
            )
        )

    selected_large = max(
        large_records,
        key=lambda r: (
            r["completion_rate"],
            r["mean_relations"],
            -r["measurement"]["median_seconds"],
        ),
    )
    chosen_large = next(
        c for name, c in large_choices if name == selected_large["name"]
    )
    return (
        chosen,
        chosen_large,
        dict(
            small=records,
            large=large_records,
            selected_small=selected["name"],
            selected_large=selected_large["name"],
            policy="Small complete/time; large complete/yield/effort.",
        ),
    )


def _config(values):
    values = dict(values)
    values["collector"] = SieveConfig(**values["collector"])
    return SIQSConfig(**values)


def bootstrap(records, before, after):
    values = {}

    for label in (before, after):
        by_input = {}

        for record in records:
            if record["arm"] == label:
                for sample in record["measurement"]["attempts"][-1]["samples"]:
                    for row in sample["rows"]:
                        by_input.setdefault(
                            (row["id"], row["seed"]), []
                        ).append(row)

        values[label] = by_input

    keys = sorted(set(values[before]) & set(values[after]))
    generator = random.Random(3034)
    ratios, differences = [], []

    for _ in range(2000):
        picked = [generator.choice(keys) for _ in keys]
        totals = {}
        completions = {}

        for label in (before, after):
            totals[label] = sum(
                statistics.median(r["seconds"] for r in values[label][k])
                for k in picked
            )
            completions[label] = sum(
                statistics.mean(r["completed"] for r in values[label][k])
                for k in picked
            ) / len(picked)

        ratios.append(totals[after] / totals[before])
        differences.append(completions[after] - completions[before])

    ratios.sort()
    differences.sort()
    return dict(
        median_ratio_interval=[ratios[50], ratios[1949]],
        completion_difference_interval=[differences[50], differences[1949]],
        unit="Paired input/seed bootstrap on this fixed corpus.",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--cold", choices=("qs", "mpqs", "siqs", "ecm"))
    args = parser.parse_args()
    corpus = json.loads(source_path(CORPUS).read_text())
    frozen_path = args.output.with_suffix(".frozen.json")
    if args.cold:
        frozen = json.loads(source_path(frozen_path).read_text())
        held = [f for f in corpus["fixtures"] if f["split"] == "held_out"]
        representatives = [
            next(f for f in held if f["band"] == band)
            for band in (
                "small",
                "balanced_30d",
                "balanced_40d",
                "balanced_50d",
                "balanced_60d",
                "balanced_70d",
                "balanced_80d",
            )
        ]
        started = time.perf_counter()
        rows = []

        for f in representatives:
            config = _config(
                frozen["small"] if f["band"] == "small" else frozen["large"]
            )
            config = replace(
                config, mode=args.cold if args.cold != "ecm" else "siqs"
            )
            rows.append(run_one(f, 7, config, ecm=args.cold == "ecm"))

        usage = resource.getrusage(resource.RUSAGE_SELF)
        print(
            json.dumps(
                dict(
                    rows=rows,
                    work_seconds=time.perf_counter() - started,
                    process_cpu_seconds=usage.ru_utime + usage.ru_stime,
                    peak_rss_bytes=_rss_bytes(),
                )
            )
        )
        return

    if (
        args.output.exists()
        or frozen_path.exists()
        or not math.isfinite(args.warmup_seconds)
        or args.warmup_seconds < 3
        or args.repetitions < 9
    ):
        parser.error(
            "new output, three-second warmup and nine samples required"
        )

    verify_certificates(corpus["certificates"])
    measured = environment()
    small_config, large_config, tuning = training(corpus, args)
    frozen = dict(
        small=asdict(small_config),
        large=asdict(large_config),
        corpus_sha256=hashlib.sha256(
            source_path(CORPUS).read_bytes()
        ).hexdigest(),
        source_sha256=measured["source_sha256"],
        budgets=dict(
            work=WORK,
            small_seconds=SMALL_SECONDS,
            large_seconds=LARGE_SECONDS,
            training_large_seconds=0.08,
            memory_bytes=MEMORY,
        ),
        training=tuning,
    )
    frozen_path.write_text(json.dumps(frozen, indent=2) + "\n")
    print(
        "frozen",
        tuning["selected_small"],
        tuning["selected_large"],
        flush=True,
    )
    held = [f for f in corpus["fixtures"] if f["split"] == "held_out"]
    records = []

    for cohort_name in ("small", "large"):
        fixtures = [
            f
            for f in held
            if (f["band"] == "small") == (cohort_name == "small")
        ]
        selected = small_config if cohort_name == "small" else large_config
        arms = [
            ("qs", replace(selected, mode="qs")),
            ("mpqs", replace(selected, mode="mpqs")),
            ("siqs", selected),
            ("ecm", selected),
        ]
        if cohort_name == "small":
            arms.extend(
                [
                    (
                        "fresh_store",
                        replace(base_config(), shared_relations=False),
                    ),
                    ("shared_store", base_config()),
                    ("simple", replace(selected, diverse=False)),
                    ("recovery", replace(selected, growth_steps=1)),
                    ("h1", replace(selected, multiplier=1)),
                    ("scored", replace(selected, multiplier=0)),
                ]
            )
        else:
            # Matched heuristic/recovery probes retain capped outcomes.
            for label, options in (
                ("h1", dict(multiplier=1)),
                ("scored", dict(multiplier=0)),
                ("simple", dict(diverse=False)),
                ("recovery", dict(growth_steps=1)),
            ):
                rows = [
                    run_one(f, 7, replace(selected, **options))
                    for f in fixtures
                ]
                frozen.setdefault("large_controls", []).append(
                    dict(arm=label, rows=rows)
                )

        for seed in corpus["seeds"]:
            for arm, config in arms:
                measurement = _measure(
                    fixtures, seed, config, args, ecm=arm == "ecm"
                )
                records.append(
                    dict(
                        cohort=cohort_name,
                        arm=arm,
                        seed=seed,
                        config=asdict(config),
                        measurement=measurement,
                    )
                )
                completed = sum(
                    row["completed"]
                    for sample in measurement["attempts"][-1]["samples"]
                    for row in sample["rows"]
                )
                print(
                    cohort_name,
                    arm,
                    seed,
                    "ms",
                    round(measurement["median_seconds"] * 1000, 3),
                    "completed",
                    completed,
                    "stable",
                    measurement["stable"],
                    flush=True,
                )

    cold = []

    for arm in ("qs", "mpqs", "siqs", "ecm"):
        for _ in range(9):
            started = time.perf_counter()

            child = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "v2.benchmarks.qs.phase_three.phase_three_siqs",
                    "--output",
                    str(args.output),
                    "--cold",
                    arm,
                ],
                capture_output=True,
                text=True,
                check=True,
                timeout=30,
            )
            value = json.loads(child.stdout)
            value.update(
                arm=arm, lifecycle_seconds=time.perf_counter() - started
            )
            cold.append(value)

    checkpoint_samples = []
    fixture = next(f for f in held if f["band"] == "small")
    warm_started = time.perf_counter()
    calls = 0

    def roundtrip():
        job = SIQSJob(
            fixture["n"],
            config=small_config,
            budget=Budget(work_limit=WORK, seconds=10, cpu_seconds=10),
        )
        job.run(max_blocks=1)
        snapshot = job.checkpoint()
        before = snapshot["resources"]["work_used"]

        resumed = SIQSJob.from_checkpoint(
            snapshot,
            budget=Budget(work_limit=WORK, seconds=10, cpu_seconds=10),
        )

        result = resumed.run()

        if result.divisor is None or result.divisor not in fixture["factors"]:
            raise AssertionError("checkpoint roundtrip failed factoring")
        return dict(
            bytes=len(snapshot["blob"]) + 1024,
            rebuild_work=resumed.budget.used - before,
            stats=result.stats,
        )

    while time.perf_counter() - warm_started < 3:
        roundtrip()
        calls += 1
    checkpoint_warmup = time.perf_counter() - warm_started
    for _ in range(9):
        started = time.perf_counter()
        value = roundtrip()
        value["seconds"] = time.perf_counter() - started
        checkpoint_samples.append(value)

    profiles = []

    for label, config in (
        ("qs", replace(small_config, mode="qs")),
        ("siqs", small_config),
    ):
        profile = cProfile.Profile()
        profile.runcall(run_one, fixture, 7, config)
        profiles.append(
            dict(
                arm=label,
                stage_functions=[
                    dict(
                        function=getattr(
                            entry.code, "co_name", str(entry.code)
                        ),
                        calls=entry.callcount,
                        total_seconds=entry.totaltime,
                        own_seconds=entry.inlinetime,
                    )
                    for entry in profile.getstats()
                    if getattr(entry.code, "co_name", "")
                    in (
                        "build_factor_base",
                        "_make_step",
                        "_sieve",
                        "_divide",
                        "_admit",
                        "classify_prime",
                        "verify_atomic",
                        "combine_relations",
                        "filter_matrix",
                        "extract_dependency",
                        "step",
                        "modular_square_roots",
                    )
                ],
                entries=[
                    dict(
                        function=str(entry.code),
                        calls=entry.callcount,
                        total_seconds=entry.totaltime,
                    )
                    for entry in sorted(
                        profile.getstats(),
                        key=lambda e: e.totaltime,
                        reverse=True,
                    )[:20]
                ],
            )
        )

    assert measured["source_sha256"] == environment()["source_sha256"]
    small_records = [r for r in records if r["cohort"] == "small"]
    intervals = {
        label: bootstrap(small_records, before, after)
        for label, before, after in (
            ("shared_vs_fresh", "fresh_store", "shared_store"),
            ("siqs_vs_qs", "qs", "siqs"),
            ("scored_vs_h1", "h1", "scored"),
            ("recovery_vs_fixed", "siqs", "recovery"),
            ("diverse_vs_simple", "simple", "siqs"),
        )
    }
    output = dict(
        milestone="M30 / complete bounded P3.4",
        environment=measured,
        command=sys.orig_argv,
        corpus_sha256=frozen["corpus_sha256"],
        frozen=frozen,
        records=records,
        confidence=intervals,
        cold_samples=cold,
        checkpoint_samples=checkpoint_samples,
        checkpoint_warmup=dict(seconds=checkpoint_warmup, calls=calls),
        profiles=profiles,
        limitations=[
            "Large results are capped; failures have no success-time ratio.",
            "Runtime probable labels stay distinct from independent proofs.",
            "Warm RSS is shared; cold RSS belongs to individual children.",
            "Finite tested settings are not a universal parameter optimum.",
        ],
    )
    args.output.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
