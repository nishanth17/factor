"""Long-budget P3.4 training, certified large coverage, and varied portfolios.

Diagnostic training screens are not timing promotion evidence. Matched
input/seed outcomes report completion under the declared
caps; failed attempts cannot establish a successful-factorization speed ratio.
A separate repeated 30-digit control measures the optional power sieve.
"""

import argparse
import cProfile
import hashlib
import json
import resource
import statistics
import subprocess
import sys
import time
from collections import Counter
from dataclasses import asdict, replace
from math import prod
from pathlib import Path

from ....execution.budget import Budget
from ....portfolio import PortfolioConfig, factorize_bounded
from ....qs import SieveConfig, SIQSConfig
from ....qs.polynomial import a_target
from ...suites.build_phase_two_corpus import verify_certificates
from ...suites.phase_one import environment
from ...support.paths import (
    BENCHMARK_ROOT,
    source_path,
)
from .phase_three_reference import _rss_bytes
from .phase_three_siqs import _config

CORPUS = BENCHMARK_ROOT / "inputs/corpora/phase_three_p34_large_v3_corpus.json"
MEMORY = 1024 * 1024 * 1024
WORK = 10**13
TIERS = ((2000, 147396, 32), (11000, 1000000, 32), (50000, 5000000, 64))


def allowance(digits):
    return 30 if digits <= 40 else 60 if digits <= 60 else 300


def configuration(n, bound, width, residual, multiplier=1):
    target = a_target(n, width)
    count = next((k for k in range(1, 9) if (bound // 2) ** k >= target), 8)
    return SIQSConfig(
        base_bound=bound,
        half_width=width,
        max_half_width=width,
        factor_count=count,
        family_count=64,
        pool_size=32,
        max_stalled=8192,
        max_trivial=4096,
        row_excess=32,
        batch_width=4096,
        memory_bytes=768 * 1024 * 1024,
        checkpoint_bytes=16 * 1024 * 1024,
        multiplier=multiplier,
        collector=SieveConfig(
            block_width=4096,
            division="bucket",
            score_policy="powers",
            residual_bound=residual,
            max_atoms=65536,
            max_relations=32768,
            max_partials=32768,
        ),
    )


def validate(row, fixture):
    if prod(row["factors"]) * prod(row["remaining"]) != fixture["n"]:
        raise AssertionError("every partial result must reconstruct")
    if any(p not in fixture["factors"] for p in row["factors"]):
        raise AssertionError("unproved terminal factor")
    if row["complete"] and sorted(row["factors"]) != fixture["factors"]:
        raise AssertionError("wrong complete prime factorization")


def run_one(
    fixture,
    seed,
    config,
    arm,
    seconds=None,
    *,
    varied=False,
    checkpoint=None,
    checkpoint_dir=None,
):
    seconds = allowance(fixture["digits"]) if seconds is None else seconds
    if varied:
        portfolio = PortfolioConfig(
            memory_bytes=MEMORY,
            rho_attempts=4,
            rho_evaluations=100000,
            pm1_attempts=2,
            pm1_b1=11000,
            pm1_b2=1000000,
            ecm_tiers=((2000, 147396, 2), (11000, 1000000, 4)),
            fermat_steps=100000,
            siqs=config if arm == "fallback" else None,
        )
    else:
        portfolio = PortfolioConfig(
            memory_bytes=MEMORY,
            trial_bound=2,
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=TIERS if arm == "ecm" else (),
            fermat_steps=0,
            siqs=None if arm == "ecm" else replace(config, mode=arm),
        )

    print("starting", fixture["id"], arm, seed, "cap", seconds, flush=True)
    budget = Budget(work_limit=WORK, seconds=seconds, cpu_seconds=seconds)
    began, cpu = time.perf_counter(), time.process_time()

    result = factorize_bounded(
        fixture["n"],
        seed=seed,
        config=portfolio,
        budget=budget,
        checkpoint=checkpoint,
    )
    factors = [
        f.value for f in result.result.factors for _ in range(f.exponent)
    ]
    state = result.checkpoint["payload"]["state"]
    stats = {}
    for event in result.events:
        if event.get("stage") == "siqs":
            stats = event.get("stats", {})
    current = state.get("current")
    if current and current.get("siqs_checkpoint"):
        stats = json.loads(current["siqs_checkpoint"]["blob"])["stats"]
    row = dict(
        id=fixture["id"],
        kind=fixture["kind"],
        digits=fixture["digits"],
        arm=arm,
        seed=seed,
        cap_seconds=seconds,
        seconds=time.perf_counter() - began,
        cpu_seconds=time.process_time() - cpu,
        work_used=budget.used,
        complete=result.result.complete,
        reason=result.reason,
        factors=factors,
        remaining=list(result.result.remaining),
        certainty=[f.certainty.value for f in result.result.factors],
        stats=stats,
        stage_seconds=state["stage_seconds"],
        events=result.events,
        peak_rss_bytes=_rss_bytes(),
    )
    validate(row, fixture)
    row["total_wall_used"] = budget.wall_used
    row["total_cpu_used"] = budget.cpu_used
    row["resumed"] = checkpoint is not None
    if checkpoint_dir is not None and not row["complete"]:
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
        began_write = time.perf_counter()
        encoded = json.dumps(result.checkpoint, separators=(",", ":")).encode()
        digest = hashlib.sha256(encoded).hexdigest()
        path = checkpoint_dir / (
            fixture["id"]
            + "."
            + arm
            + "."
            + str(seed)
            + "."
            + digest[:16]
            + ".json"
        )
        with path.open("xb") as output:
            output.write(encoded)
        row["checkpoint"] = dict(
            path=str(path.relative_to(checkpoint_dir.parent)),
            sha256=digest,
            bytes=len(encoded),
            write_seconds=time.perf_counter() - began_write,
        )
        row["seconds"] += row["checkpoint"]["write_seconds"]
    return row


def summary(rows):
    return dict(
        attempts=len(rows),
        completed=sum(r["complete"] for r in rows),
        reasons=dict(Counter(r["reason"] for r in rows)),
        median_attempt_seconds=statistics.median(r["seconds"] for r in rows),
        successful_median_seconds=statistics.median(
            r["seconds"] for r in rows if r["complete"]
        )
        if any(r["complete"] for r in rows)
        else None,
        relations=sum(r["stats"].get("relations", 0) for r in rows),
        partials=sum(r["stats"].get("partials", 0) for r in rows),
        timing_scope=(
            "All attempts including failures; success times "
            "are descriptive and may have different completed "
            "subsets."
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--phase",
        choices=(
            "train",
            "balanced",
            "varied",
            "controls",
            "cold",
            "continuation",
            "all",
        ),
        default="all",
    )
    parser.add_argument("--cold-child", choices=("qs", "mpqs", "siqs", "ecm"))
    parser.add_argument(
        "--bands", nargs="+", type=int, choices=(30, 40, 50, 60, 70, 80)
    )
    parser.add_argument(
        "--training-capture",
        type=Path,
        help="Freeze inherited choices after fresh calibration",
    )
    args = parser.parse_args()
    checkpoint_dir = args.output.with_suffix(".checkpoints")
    corpus = json.loads(source_path(CORPUS).read_text())
    frozen_path = args.output.with_suffix(".frozen.json")
    if args.cold_child:
        frozen = json.loads(source_path(frozen_path).read_text())
        fixture = next(
            f
            for f in corpus["fixtures"]
            if f["kind"] == "balanced"
            and f["digits"] == 30
            and f["split"] == "held_out"
        )
        row = run_one(
            fixture, 7, _config(frozen["configs"]["30"]), args.cold_child
        )
        usage = resource.getrusage(resource.RUSAGE_SELF)
        print(
            json.dumps(
                dict(
                    row=row,
                    process_cpu_seconds=usage.ru_utime + usage.ru_stime,
                    peak_rss_bytes=_rss_bytes(),
                )
            )
        )
        return

    verify_certificates(corpus["certificates"])
    corpus_sha = hashlib.sha256(source_path(CORPUS).read_bytes()).hexdigest()
    source = environment()
    phases = (
        ("train", "balanced", "varied", "continuation", "controls", "cold")
        if args.phase == "all"
        else (args.phase,)
    )

    for phase in phases:
        journal_path = args.output.with_name(
            args.output.stem + "." + phase + ".jsonl"
        )
        if journal_path.exists():
            raise ValueError(
                "preserve each prior phase; choose a fresh output"
            )
        with journal_path.open("x") as journal:

            def record(row):
                journal.write(json.dumps(row) + "\n")
                journal.flush()
                print(
                    phase,
                    row.get("id"),
                    row.get("arm"),
                    row.get("seed"),
                    row.get("complete"),
                    row.get("reason"),
                    round(row.get("seconds", 0), 3),
                    "rows",
                    row.get("stats", {}).get("relations", 0),
                    "partials",
                    row.get("stats", {}).get("partials", 0),
                    flush=True,
                )

            if phase == "train":
                if frozen_path.exists():
                    raise ValueError("freeze is immutable")
                records, configs = {}, {}
                inherited = None
                if args.training_capture is not None:
                    encoded = source_path(args.training_capture).read_bytes()
                    inherited = json.loads(encoded)
                    assert inherited["corpus_sha256"] == corpus_sha
                    inherited_sha = hashlib.sha256(encoded).hexdigest()

                for digits in (30, 40, 50, 60, 70, 80):
                    fixtures = [
                        f
                        for f in corpus["fixtures"]
                        if f["kind"] == "balanced"
                        and f["digits"] == digits
                        and f["split"] == "training"
                    ]
                    choices = []
                    if inherited is not None:
                        config = replace(
                            _config(inherited["configs"][str(digits)]),
                            memory_bytes=768 * 1024 * 1024,
                            checkpoint_bytes=16 * 1024 * 1024,
                        )
                        rows = []

                        for fixture in fixtures:
                            row = run_one(
                                fixture,
                                7,
                                config,
                                "siqs",
                                seconds=15 if digits <= 40 else 30,
                                checkpoint_dir=checkpoint_dir,
                            )
                            row["training_candidate"] = "inherited_calibration"
                            record(row)
                            rows.append(row)

                        configs[str(digits)] = asdict(config)
                        records[str(digits)] = dict(
                            rows=rows, summary=summary(rows)
                        )
                        continue

                    grid = (
                        (
                            ("small", 10000, 8192, 10**7, 1),
                            ("small_lp", 10000, 8192, 10**10, 1),
                            ("moderate", 50000, 32768, 10**8, 1),
                        )
                        if digits <= 40
                        else (
                            ("moderate", 50000, 32768, 10**8, 1),
                            ("large_lp", 100000, 65536, 10**10, 1),
                            ("scored", 100000, 65536, 10**10, 0),
                        )
                    )

                    for label, bound, width, residual, multiplier in grid:
                        config = configuration(
                            fixtures[0]["n"],
                            bound,
                            width,
                            residual,
                            multiplier,
                        )
                        rows = []

                        for fixture in fixtures:
                            # A diagnostic screen includes generous time and
                            # all stages, but supplies no promotion statistic.
                            row = run_one(
                                fixture,
                                7,
                                config,
                                "siqs",
                                seconds=15 if digits <= 40 else 30,
                                checkpoint_dir=checkpoint_dir,
                            )
                            row["training_candidate"] = label
                            record(row)
                            rows.append(row)

                        choices.append(
                            dict(
                                name=label,
                                config=asdict(config),
                                rows=rows,
                                summary=summary(rows),
                            )
                        )

                    chosen = max(
                        choices,
                        key=lambda c: (
                            c["summary"]["completed"],
                            c["summary"]["relations"]
                            / max(0.001, sum(r["seconds"] for r in c["rows"])),
                            c["summary"]["partials"]
                            / max(0.001, sum(r["seconds"] for r in c["rows"])),
                        ),
                    )
                    configs[str(digits)] = chosen["config"]
                    records[str(digits)] = choices
                    print("selected", digits, chosen["name"], flush=True)

                frozen = dict(
                    corpus_sha256=corpus_sha,
                    environment=source,
                    configs=configs,
                    training=records,
                    budgets=dict(
                        work=WORK,
                        owned_memory=MEMORY,
                        wall_cpu_seconds_by_digits={
                            str(d): allowance(d)
                            for d in (30, 40, 50, 60, 70, 80)
                        },
                    ),
                    ecm_tiers=TIERS,
                    inherited_training_sha256=(
                        inherited_sha if inherited is not None else None
                    ),
                    calibration_policy=(
                        "Inherited algorithm choices with enlarged "
                        "memory/checkpoint allowances are frozen "
                        "after two fresh training checks per band; original "
                        "capture hash is retained. No default promotion."
                    ),
                    policy=(
                        "Single-pass diagnostic selection; completion then "
                        "checked-row yield/elapsed. Optional parameters "
                        "only, no default promotion."
                    ),
                    profile_policy=(
                        "Nine matched case/seed outcomes per balanced 30-60 "
                        "band; one fixed representative at 70/80 digits. "
                        "Independent case coverage and capped exploration, "
                        "not nine timing repeats of one call."
                    ),
                )
                with frozen_path.open("x") as stream:
                    json.dump(frozen, stream, indent=2)
            else:
                frozen = json.loads(source_path(frozen_path).read_text())
                assert frozen["corpus_sha256"] == corpus_sha
                assert (
                    frozen["environment"]["source_sha256"]
                    == source["source_sha256"]
                )
                rows = []
                if phase == "balanced":
                    for digits in args.bands or (30, 40, 50, 60, 70, 80):
                        fixtures = [
                            f
                            for f in corpus["fixtures"]
                            if f["kind"] == "balanced"
                            and f["digits"] == digits
                            and f["split"] == "held_out"
                        ]
                        config = _config(frozen["configs"][str(digits)])
                        warm_fixture = next(
                            f
                            for f in corpus["fixtures"]
                            if f["kind"] == "balanced"
                            and f["digits"] == digits
                            and f["split"] == "training"
                        )
                        arms = (
                            ("siqs", "ecm")
                            if digits >= 70
                            else ("qs", "mpqs", "siqs", "ecm")
                        )

                        for arm in arms:
                            warm = time.perf_counter()
                            calls = 0
                            while time.perf_counter() - warm < 3:
                                run_one(
                                    warm_fixture, 7, config, arm, seconds=3
                                )
                                calls += 1

                            record(
                                dict(
                                    type="warmup",
                                    digits=digits,
                                    arm=arm,
                                    seconds=time.perf_counter() - warm,
                                    calls=calls,
                                )
                            )
                            pairs = [
                                (f, s) for s in (7, 29) for f in fixtures
                            ] + [(fixtures[0], 47)]
                            # Large-band exploration uses a fixed case at
                            # seed 7, with no timing promotion claim. Nine
                            # repeats are in the separate timing control.
                            if digits >= 70:
                                pairs = [(f, 7) for f in fixtures[:1]]
                            for fixture, seed in pairs:
                                row = run_one(
                                    fixture,
                                    seed,
                                    config,
                                    arm,
                                    checkpoint_dir=checkpoint_dir,
                                )
                                record(row)
                                rows.append(row)
                elif phase == "varied":
                    fixtures = [
                        f
                        for f in corpus["fixtures"]
                        if f["kind"] != "balanced" and f["split"] == "held_out"
                    ]

                    for arm in ("ecm", "fallback"):
                        warm = next(
                            f
                            for f in corpus["fixtures"]
                            if f["kind"] == "balanced"
                            and f["digits"] == 30
                            and f["split"] == "training"
                        )
                        began = time.perf_counter()

                        while time.perf_counter() - began < 3:
                            run_one(
                                warm,
                                7,
                                _config(frozen["configs"]["30"]),
                                arm,
                                seconds=3,
                                varied=True,
                            )

                        for fixture in fixtures:
                            digits = min(
                                80,
                                max(30, 10 * ((fixture["digits"] + 9) // 10)),
                            )
                            config = _config(frozen["configs"][str(digits)])

                            for seed in (7, 29):
                                row = run_one(
                                    fixture,
                                    seed,
                                    config,
                                    arm,
                                    varied=True,
                                    checkpoint_dir=checkpoint_dir,
                                )
                                record(row)
                                rows.append(row)
                elif phase == "continuation":
                    prior_rows = [
                        json.loads(line)
                        for line in source_path(
                            args.output.with_name(
                                args.output.stem + ".balanced.jsonl"
                            )
                        )
                        .read_text()
                        .splitlines()
                    ]
                    fixture = next(
                        f
                        for f in corpus["fixtures"]
                        if f.get("performance_representative")
                        and f["digits"] == 50
                    )
                    prior = next(
                        r
                        for r in prior_rows
                        if r.get("id") == fixture["id"]
                        and r.get("arm") == "siqs"
                        and r.get("seed") == 7
                    )
                    config = _config(frozen["configs"]["50"])

                    for total_seconds in (600, 1800):
                        if prior["complete"] or prior["reason"] not in (
                            "wall_limit",
                            "cpu_limit",
                            "work_limit",
                        ):
                            break

                        metadata = prior["checkpoint"]
                        path = args.output.parent / metadata["path"]
                        encoded = source_path(path).read_bytes()
                        assert (
                            hashlib.sha256(encoded).hexdigest()
                            == metadata["sha256"]
                        )
                        row = run_one(
                            fixture,
                            7,
                            config,
                            "siqs",
                            seconds=total_seconds,
                            checkpoint=json.loads(encoded),
                            checkpoint_dir=checkpoint_dir,
                        )
                        row["continuation_policy"] = (
                            "Fixed 50-digit representative; "
                            "increased cumulative limits"
                        )
                        row["parent_checkpoint"] = metadata["sha256"]
                        record(row)
                        rows.append(row)
                        prior = row
                elif phase == "controls":
                    fixture = next(
                        f
                        for f in corpus["fixtures"]
                        if f.get("performance_representative")
                        and f["digits"] == 30
                    )
                    base = _config(frozen["configs"]["30"])

                    for label, options in (
                        ("powers", {}),
                        ("upper_bound", {"score_policy": "candidate"}),
                    ):
                        config = replace(
                            base, collector=replace(base.collector, **options)
                        )
                        warm = time.perf_counter()
                        while time.perf_counter() - warm < 3:
                            run_one(fixture, 7, config, "siqs")
                        samples = []

                        for index in range(9):
                            row = run_one(fixture, 7, config, "siqs")
                            row["control"] = label
                            row["sample"] = index
                            record(row)
                            samples.append(row)
                            rows.append(row)

                        times = [r["seconds"] for r in samples]
                        quartiles = statistics.quantiles(times, n=4)
                        stable = (
                            abs(
                                statistics.median(times[:3])
                                / statistics.median(times[-3:])
                                - 1
                            )
                            <= 0.15
                            and (quartiles[2] - quartiles[0])
                            / statistics.median(times)
                            <= 0.2
                        )
                        if not stable:
                            warm = time.perf_counter()
                            while time.perf_counter() - warm < 5:
                                run_one(fixture, 7, config, "siqs")
                            for index in range(15):
                                row = run_one(fixture, 7, config, "siqs")
                                row["control"] = label
                                row["extension"] = True
                                row["sample"] = index
                                record(row)
                                rows.append(row)

                    for digits in (50, 60, 80):
                        fixture = next(
                            f
                            for f in corpus["fixtures"]
                            if f.get("performance_representative")
                            and f["digits"] == digits
                        )
                        base = _config(frozen["configs"][str(digits)])

                        for label, changes in (
                            ("simple", dict(diverse=False)),
                            (
                                "width_recovery",
                                dict(
                                    growth_steps=1,
                                    max_half_width=min(
                                        499999, 2 * base.half_width
                                    ),
                                ),
                            ),
                        ):
                            row = run_one(
                                fixture, 7, replace(base, **changes), "siqs"
                            )
                            row["control"] = label
                            record(row)
                            rows.append(row)

                    profiler = cProfile.Profile()
                    profiler.enable()
                    run_one(fixture, 7, base, "siqs", seconds=10)
                    profiler.disable()
                    profiles = [
                        dict(
                            function=str(s.code),
                            calls=s.callcount,
                            own_seconds=s.inlinetime,
                            total_seconds=s.totaltime,
                        )
                        for s in profiler.getstats()
                    ]
                    with args.output.with_name(
                        args.output.stem + ".profiles.json"
                    ).open("x") as stream:
                        json.dump(profiles, stream, indent=2)
                elif phase == "cold":
                    for arm in ("qs", "mpqs", "siqs", "ecm"):
                        for index in range(9):
                            began = time.perf_counter()

                            child = subprocess.run(
                                [
                                    sys.executable,
                                    "-m",
                                    __spec__.name,
                                    "--output",
                                    str(args.output),
                                    "--cold-child",
                                    arm,
                                ],
                                capture_output=True,
                                text=True,
                                check=True,
                            )
                            value = json.loads(
                                child.stdout.strip().splitlines()[-1]
                            )
                            row = value["row"]
                            row.update(
                                cold=True,
                                sample=index,
                                lifecycle_seconds=time.perf_counter() - began,
                                process_cpu_seconds=value[
                                    "process_cpu_seconds"
                                ],
                                peak_rss_bytes=value["peak_rss_bytes"],
                            )
                            record(row)
                            rows.append(row)

                result = dict(
                    phase=phase,
                    corpus_sha256=corpus_sha,
                    environment=source,
                    rows=rows,
                    summary={
                        arm: summary([r for r in rows if r["arm"] == arm])
                        for arm in sorted({r["arm"] for r in rows})
                    },
                )
                with args.output.with_name(
                    args.output.stem + "." + phase + ".json"
                ).open("x") as stream:
                    json.dump(result, stream, indent=2)

        print("phase_complete", phase, flush=True)


if __name__ == "__main__":
    main()
