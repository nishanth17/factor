"""Train, freeze, then confirm finite SIQS capacity across sub-100-digit bands.

Training screens select controls, not speedups. Confirmation measures full
calls under shared allowances, with validated warmup and noise extension.
Censored attempts never imply successful-factorization time ratios.
"""

import argparse
import cProfile
import hashlib
import json
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from collections import Counter
from dataclasses import asdict, replace
from math import isqrt, prod
from pathlib import Path

from .. import utils
from ..budget import Budget
from ..portfolio import PortfolioConfig, factorize_bounded
from ..qs import SieveConfig, SIQSConfig
from .build_p38_r1_corpus import BANDS
from .build_phase_two_corpus import verify_certificates
from .phase_three_reference import _rss_bytes

ROOT = Path(__file__).resolve().parents[2]
WORK = 10**13
MEMORY = 512 * 2**20
PRESETS = ((1000, 8192), (3000, 32768), (10000, 65536))
ARMS = ("siqs", "mpqs", "ecm", "ecm_siqs")


def source_hashes():
    """Pin loaded arithmetic and this driver; exclude local result captures."""
    paths = sorted((ROOT / "v2").glob("*.py"))
    paths += sorted((ROOT / "v2/qs").glob("*.py"))
    paths += [
        Path(__file__),
        Path(__file__).with_name("build_p38_r1_corpus.py"),
    ]
    return {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in paths
    }


def configuration(digits, preset, policy):
    """Select only from public size/configuration, never known factor sizes."""
    bound, width = PRESETS[preset]
    # Upper bands use a declared larger base, still within the same owned cap.
    if digits >= 60:
        bound *= 10
    target = isqrt(2 * 10**digits) // width
    count = next(k for k in range(1, 33) if (bound // 3) ** k >= target)
    rows = min(8192, max(2048, bound // 8))
    return SIQSConfig(
        base_bound=bound,
        half_width=width,
        max_half_width=width,
        factor_count=count,
        family_count=100000,
        pool_size=max(32, count * 4),
        assignment_policy=policy,
        polynomials_per_family=64,
        max_stalled=100000,
        max_trivial=100000,
        row_excess=32,
        batch_width=4096,
        memory_bytes=384 * 2**20,
        checkpoint_bytes=16 * 2**20,
        collector=SieveConfig(
            block_width=4096,
            score_policy="powers",
            division="bucket",
            residual_bound=min(10**12, bound**2),
            max_relations=rows,
            max_partials=rows,
            max_atoms=min(65536, rows * 4),
        ),
    )


def decode_config(values):
    values = dict(values)
    values["collector"] = SieveConfig(**values["collector"])
    return SIQSConfig(**values)


def validate(row, fixture):
    """Validate every proper split, terminal factor, label and cofactor."""
    n = fixture["n"]
    if prod(row["factors"]) * prod(row["remaining"]) != n:
        raise AssertionError("result does not reconstruct the complete input")
    if any(p not in fixture["factors"] for p in row["factors"]):
        raise AssertionError(
            "terminal factor disagrees with independent proof"
        )
    if row["complete"] and sorted(row["factors"]) != fixture["factors"]:
        raise AssertionError("complete result has incorrect multiplicity")
    if row["first_factor"] is not None and not utils.valid_divisor(
        row["first_factor"], n
    ):
        raise AssertionError("invalid proper split")
    if any(
        label not in ("proven_prime", "probable_prime")
        for label in row["certainty"]
    ):
        raise AssertionError("invalid certainty label")


def run_one(fixture, seed, config, arm, seconds):
    """Measure setup through classification under one shared budget."""
    if arm not in ARMS:
        raise ValueError("unknown R1 arm")
    if arm == "mpqs":
        config = replace(
            config,
            mode="mpqs",
            assignment_policy="reference",
            external_coefficients=True,
            factor_count=1,
            polynomials_per_family=1,
        )
    portfolio = PortfolioConfig(
        trial_bound=1000,
        rho_attempts=0,
        pm1_attempts=1,
        pm1_b1=1000,
        pm1_b2=10000,
        fermat_steps=1000,
        ecm_tiers=((2000, 147396, 8), (11000, 1000000, 8))
        if arm in ("ecm", "ecm_siqs")
        else (),
        siqs=None if arm == "ecm" else config,
        memory_bytes=MEMORY,
    )
    allowance = Budget(work_limit=WORK, seconds=seconds, cpu_seconds=seconds)
    start, cpu = time.perf_counter(), time.process_time()
    result = factorize_bounded(
        fixture["n"], seed=seed, config=portfolio, budget=allowance
    )
    elapsed, cpu_elapsed = (
        time.perf_counter() - start,
        time.process_time() - cpu,
    )
    factors = [
        f.value for f in result.result.factors for _ in range(f.exponent)
    ]
    remaining = list(result.result.remaining)
    state = result.checkpoint["payload"]["state"]
    stats = {}
    for event in result.events:
        if event.get("stage") == "siqs":
            stats = event.get("stats", {})
    current = state.get("current")
    if current and current.get("siqs_checkpoint"):
        stats = json.loads(current["siqs_checkpoint"]["blob"])["stats"]
    first = next(
        (
            p
            for p in factors + remaining
            if utils.valid_divisor(p, fixture["n"])
        ),
        None,
    )
    row = dict(
        id=fixture["id"],
        kind=fixture["kind"],
        digits=fixture["digits"],
        smaller_factor_digits=fixture["smaller_factor_digits"],
        seed=seed,
        arm=arm,
        seconds=elapsed,
        cpu_seconds=cpu_elapsed,
        work_used=allowance.used,
        cap_seconds=seconds,
        reason=result.reason,
        complete=result.result.complete,
        first_factor=first,
        factors=factors,
        remaining=remaining,
        certainty=[f.certainty.value for f in result.result.factors],
        stats=stats,
        stage_seconds=state["stage_seconds"],
        peak_process_rss_bytes=_rss_bytes(),
        events=result.events,
    )
    validate(row, fixture)
    return row


def summarize(rows):
    """Keep successful and censored populations explicit."""
    successes = [r["seconds"] for r in rows if r["complete"]]
    return dict(
        attempts=len(rows),
        completed=len(successes),
        split_cases=sum(r["first_factor"] is not None for r in rows),
        median_attempt_seconds=statistics.median(r["seconds"] for r in rows),
        successful_median_seconds=statistics.median(successes)
        if successes
        else None,
        reasons=dict(Counter(r["reason"] for r in rows)),
        verified_rows=sum(r["stats"].get("relations", 0) for r in rows),
        owned_peak_bytes=max(
            (r["stats"].get("workspace_bytes", 0) for r in rows), default=0
        ),
        note="Unfinished times are censored, not time-to-factor ratios. "
        "Process RSS is a lifetime peak, separate from owned reservations.",
    )


def measure(fixtures, config, arm, seconds, seeds, repetitions=9, warmup=3):
    """Warm validated calls; extend noisy cohorts to fifteen samples."""
    attempts = []
    for attempt in range(3):
        started, warm_calls = time.perf_counter(), 0
        target_warm = warmup if attempt == 0 else max(5, warmup)
        while time.perf_counter() - started < target_warm:
            for fixture in fixtures:
                for seed in seeds:
                    run_one(fixture, seed, config, arm, seconds)
                    warm_calls += 1
        warm_seconds = time.perf_counter() - started
        samples = []
        for _ in range(repetitions if attempt == 0 else max(15, repetitions)):
            samples.append(
                [
                    run_one(f, seed, config, arm, seconds)
                    for f in fixtures
                    for seed in seeds
                ]
            )
        totals = [sum(r["seconds"] for r in sample) for sample in samples]
        median = statistics.median(totals)
        mad = statistics.median(abs(t - median) for t in totals)
        stable = (
            mad <= median * 0.1 and max(totals) - min(totals) <= median * 0.35
        )
        attempts.append(
            dict(
                warmup_seconds=warm_seconds,
                warmup_calls=warm_calls,
                samples=samples,
                median_seconds=median,
                mad_seconds=mad,
                stable=stable,
                summary=summarize([r for sample in samples for r in sample]),
            )
        )
        if stable:
            break
    return dict(
        attempts=attempts,
        accepted=attempts[-1],
        timing_scope="Warmed calls include setup and serialization. "
        "Coarse stage clocks are present in every arm; cProfile is separate.",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase",
        choices=("train", "confirm", "screen", "profile", "cold", "once"),
        required=True,
    )
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--frozen", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--bands", type=int, nargs="+", default=BANDS, choices=BANDS
    )
    parser.add_argument("--kinds", default="balanced")
    parser.add_argument("--arms", default=",".join(ARMS))
    parser.add_argument(
        "--seconds",
        type=float,
        default=10,
        help="training/screen cap; confirmation uses frozen caps",
    )
    parser.add_argument("--confirmation-small-seconds", type=float, default=10)
    parser.add_argument("--confirmation-large-seconds", type=float, default=1)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    if (
        min(
            args.seconds,
            args.confirmation_small_seconds,
            args.confirmation_large_seconds,
        )
        <= 0
        or args.repetitions < 9
        or args.warmup_seconds < 3
    ):
        parser.error(
            "positive time, three seconds warmup and nine samples are required"
        )
    if args.output.exists():
        parser.error("preserve existing output bytes")
    corpus = json.loads(args.corpus.read_text())
    verify_certificates(corpus["certificates"])
    for fixture in corpus["fixtures"]:
        if prod(fixture["factors"]) != fixture["n"]:
            raise ValueError("corrupt corpus factorization")
        if any(
            str(p) not in corpus["certificates"] for p in fixture["factors"]
        ):
            raise ValueError("missing independent prime proof")
    hashes = source_hashes()
    data = dict(
        schema=1,
        phase=args.phase,
        source_sha256=hashes,
        corpus_sha256=hashlib.sha256(args.corpus.read_bytes()).hexdigest(),
        runtime=sys.version,
        machine=platform.machine(),
        work_limit=WORK,
        memory_bytes=MEMORY,
        seconds=args.seconds,
        results={},
    )
    fixtures = [
        f
        for f in corpus["fixtures"]
        if f["digits"] in args.bands
        and (args.kinds == "all" or f["kind"] in args.kinds.split(","))
    ]
    if not fixtures:
        parser.error("no matching fixtures")
    arms = args.arms.split(",")
    if any(arm not in ARMS for arm in arms):
        parser.error("unknown arm")

    def save():
        if hashes != source_hashes():
            raise RuntimeError("source changed during capture")
        args.output.write_text(json.dumps(data, indent=2) + "\n")

    if args.phase == "train":
        if (
            corpus["split"] != "training"
            or args.frozen is None
            or args.frozen.exists()
        ):
            parser.error(
                "training requires a training corpus and a fresh --frozen path"
            )
        selected = {}
        for digits in args.bands:
            group = [
                f
                for f in fixtures
                if f["digits"] == digits and f["kind"] == "balanced"
            ]
            choices = []
            for preset in range(len(PRESETS)):
                for policy in ("nearest", "flyer"):
                    config = configuration(digits, preset, policy)
                    rows = [
                        run_one(f, seed, config, "siqs", args.seconds)
                        for f in group
                        for seed in corpus["seeds"]
                    ]
                    summary = summarize(rows)
                    choices.append(
                        (
                            summary["completed"],
                            -summary["median_attempt_seconds"]
                            if summary["completed"]
                            else summary["verified_rows"],
                            summary["verified_rows"]
                            if summary["completed"]
                            else -summary["median_attempt_seconds"],
                            preset,
                            policy,
                            config,
                        )
                    )
                    data["results"][f"{digits}/{preset}/{policy}"] = dict(
                        rows=rows, summary=summary
                    )
                    print(digits, preset, policy, summary, flush=True)
                    save()
            winner = max(choices, key=lambda row: row[:3])
            selected[str(digits)] = asdict(winner[-1])
        frozen = dict(
            schema=1,
            source_sha256=hashes,
            configs=selected,
            training_corpus_sha256=data["corpus_sha256"],
            training_capture_sha256=hashlib.sha256(
                args.output.read_bytes()
            ).hexdigest(),
            training_seconds=args.seconds,
            seeds=corpus["seeds"],
            confirmation_seconds={
                str(d): args.confirmation_small_seconds
                if d == 30
                else args.confirmation_large_seconds
                for d in args.bands
            },
            work_limit=WORK,
            memory_bytes=MEMORY,
            decision="Rank completion then time; censored ties use rows. "
            "Exploratory censored bands establish no performance crossover.",
        )
        with args.frozen.open("x") as stream:
            json.dump(frozen, stream, indent=2)
            stream.write("\n")
        return
    if args.frozen is None:
        parser.error("this phase requires frozen configurations")
    frozen_bytes = args.frozen.read_bytes()
    frozen = json.loads(frozen_bytes)
    if frozen["source_sha256"] != hashes:
        parser.error("frozen source differs; retrain before new confirmation")
    if (
        args.phase != "screen"
        and corpus.get("frozen_sha256")
        != hashlib.sha256(frozen_bytes).hexdigest()
    ):
        parser.error("confirmation corpus was not generated from this freeze")
    if corpus["seeds"] != frozen["seeds"]:
        parser.error("confirmation seeds differ from the frozen control")
    data["frozen_sha256"] = hashlib.sha256(frozen_bytes).hexdigest()
    data["confirmation_seconds"] = frozen["confirmation_seconds"]
    for digits in args.bands:
        config = decode_config(frozen["configs"][str(digits)])
        seconds = (
            args.seconds
            if args.phase == "screen"
            else frozen["confirmation_seconds"][str(digits)]
        )
        group = [f for f in fixtures if f["digits"] == digits]
        if not group:
            continue
        for arm in arms:
            key = f"{digits}/{arm}"
            if args.phase == "confirm":
                result = measure(
                    group,
                    config,
                    arm,
                    seconds,
                    corpus["seeds"],
                    args.repetitions,
                    args.warmup_seconds,
                )
            elif args.phase in ("screen", "once"):
                rows = [
                    run_one(f, seed, config, arm, seconds)
                    for f in group
                    for seed in corpus["seeds"]
                ]
                result = dict(
                    rows=rows, summary=summarize(rows), timing_evidence=False
                )
            elif args.phase == "profile":
                profiler = cProfile.Profile()
                row = profiler.runcall(
                    run_one,
                    group[0],
                    corpus["seeds"][0],
                    config,
                    arm,
                    seconds,
                )
                profiler.dump_stats(str(args.output) + f".{digits}.{arm}.prof")
                result = dict(row=row, timing_evidence=False)
            else:
                times, outcomes = [], []
                for _ in range(args.repetitions):
                    with tempfile.TemporaryDirectory() as directory:
                        output = Path(directory) / "child.json"
                        command = [
                            sys.executable,
                            "-m",
                            "v2.benchmarks.p38_r1_capacity",
                            "--phase",
                            "once",
                            "--corpus",
                            str(args.corpus.resolve()),
                            "--frozen",
                            str(args.frozen.resolve()),
                            "--output",
                            str(output),
                            "--bands",
                            str(digits),
                            "--kinds",
                            args.kinds,
                            "--arms",
                            arm,
                            "--seconds",
                            str(seconds),
                        ]
                        started = time.perf_counter()
                        subprocess.run(
                            command, check=True, cwd=ROOT, capture_output=True
                        )
                        times.append(time.perf_counter() - started)
                        outcomes.append(
                            json.loads(output.read_text())["results"][key]
                        )
                result = dict(
                    seconds=times,
                    median_seconds=statistics.median(times),
                    outcomes=outcomes,
                    scope="Cold process including import and output.",
                )
            data["results"][key] = result
            print(key, "finished", flush=True)
            save()


if __name__ == "__main__":
    main()
