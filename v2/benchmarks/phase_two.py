"""Capped corpus runner: complete/factor-one, isolated cold and warmed runs."""

import argparse
import hashlib
import json
import math
import os
import resource
import signal
import statistics
import subprocess
import sys
import time
from contextlib import ExitStack
from dataclasses import asdict
from pathlib import Path

from v2 import ecm, pollard_pm1, pollard_rho, portfolio, utils
from v2.budget import Budget
from v2.factor import factorize, factorize_bf
from v2.portfolio import PortfolioConfig, factorize_bounded

from .build_phase_two_corpus import verify_certificates
from .phase_one import environment
from .snapshot_loader import load_stage_jobs

CORPUS = (
    Path(__file__).parent / "inputs/corpora/phase_two_complete_corpus.json"
)


class SampleTimeoutError(Exception):
    """The external wall/CPU timer stopped either comparison arm."""


def _timeout(signum, frame):
    """Interrupt Python/native work using the Unix benchmark watchdog."""
    raise SampleTimeoutError(
        "wall_limit" if signum == signal.SIGALRM else ("cpu_limit")
    )


def _rss_bytes():
    """Convert the platform's process high-water RSS to bytes."""
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def _baseline_one(n, seed, config):
    """Run actual factor-one work, exposing composite children explicitly."""
    if utils.is_prime(n, rng=utils.resolve_rng(seed)):
        return [], [n]
    factors, remainder = factorize_bf(n, bound=config.trial_bound)
    if factors:
        return factors, [remainder] if remainder > 1 else []
    divisor = pollard_rho.factorize_rho(
        n,
        seed=seed,
        max_attempts=config.rho_attempts,
        max_evaluations=config.rho_evaluations,
        batch_size=config.rho_batch,
    )
    if divisor is None:
        divisor = pollard_pm1.factorize_pm1(
            n,
            b1=config.pm1_b1,
            b2=config.pm1_b2,
            max_attempts=config.pm1_attempts,
        )

    for b1, b2, curves in config.ecm_tiers:
        if divisor is not None:
            break
        divisor = ecm.factorize_ecm(
            n, seed=seed, b1=b1, b2=b2, max_curves=curves
        )

    return [], [divisor, n // divisor] if divisor else [n]


def _sample(request):
    """Time one validated answer under identical external time/CPU caps."""
    config = PortfolioConfig(**request["config"])
    n, seed = request["n"], request["seed"]
    started = time.perf_counter()
    cpu_started = time.process_time()
    watchdog_active = True

    def timed_out(signum, frame):
        """Latch the first expiry before raising, including during cleanup."""
        nonlocal watchdog_active
        if watchdog_active:
            watchdog_active = False
            _timeout(signum, frame)

    signal.signal(signal.SIGALRM, timed_out)
    signal.signal(signal.SIGPROF, timed_out)
    signal.setitimer(signal.ITIMER_REAL, request["seconds"])
    signal.setitimer(signal.ITIMER_PROF, request["cpu_seconds"])
    factors, remaining, certainty = [], [n], []
    complete = False
    reason = "exhausted"
    stage_seconds = {}
    events = []
    work_used = None

    try:
        if request["engine"] in ("bounded", "m12"):
            run = factorize_bounded(
                n,
                seed=seed,
                config=config,
                stop_after_split=request["mode"] == "factor_one",
                budget=Budget(
                    work_limit=request["work_limit"],
                    seconds=request["seconds"],
                    cpu_seconds=request["cpu_seconds"],
                ),
            )
            answer = run.result
            reason = run.reason
            work_used = run.work_used
            events = list(run.events)
            stage_seconds = run.checkpoint["payload"]["state"]["stage_seconds"]
            factors = [(item.value, item.exponent) for item in answer.factors]
            certainty = [item.certainty.value for item in answer.factors]
            remaining, complete = list(answer.remaining), answer.complete
        elif request["mode"] == "factor_one":
            factors, remaining = _baseline_one(n, seed, config)
            certainty = [utils.Primality.PROVEN.value for _ in factors]
            reason = (
                "factor_found"
                if any(
                    utils.valid_divisor(value, n)
                    for value in remaining + [p for p, _ in factors]
                )
                else "exhausted"
            )
        else:
            b1, b2, curves = (
                config.ecm_tiers[0] if config.ecm_tiers else (2, 2, 0)
            )

            answer = factorize(
                n,
                seed=seed,
                trial_bound=config.trial_bound,
                rho_attempts=config.rho_attempts,
                rho_evaluations=config.rho_evaluations,
                ecm_curves=curves,
                ecm_b1=b1,
                ecm_b2=b2,
            )
            factors = [(item.value, item.exponent) for item in answer.factors]
            certainty = [item.certainty.value for item in answer.factors]
            remaining, complete = list(answer.remaining), answer.complete
            reason = "complete" if complete else "exhausted"
    except SampleTimeoutError as error:
        reason = str(error)
    finally:
        # Both timers can expire together. A second pending signal must not
        # raise while disarming the other timer or reconstructing the result.
        watchdog_active = False
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.setitimer(signal.ITIMER_PROF, 0)

    elapsed = time.perf_counter() - started
    cpu = time.process_time() - cpu_started
    reconstructed = 1
    for value, exponent in factors:
        reconstructed *= value**exponent
    for value in remaining:
        reconstructed *= value
    if reconstructed != n:
        raise AssertionError("sample lost a cofactor")
    rss = _rss_bytes()
    memory_exceeded = rss > request["rss_cap_bytes"]
    split = any(
        utils.valid_divisor(value, n)
        for value in (remaining + [p for p, _ in factors])
    )
    return {
        "elapsed_seconds": elapsed,
        "cpu_seconds": cpu,
        "peak_rss_bytes": rss,
        "memory_exceeded": memory_exceeded,
        "reason": reason,
        "complete": complete,
        "success": (complete if request["mode"] == "complete" else split)
        and not memory_exceeded,
        "factors": factors,
        "remaining": remaining,
        "certainty": certainty,
        "reconstructs": True,
        "stage_seconds": stage_seconds,
        "events": events,
        "work_used": work_used,
    }


def _worker():
    """Process requests sequentially so warmed measurements reuse the JIT."""
    baseline_ready = False

    for line in sys.stdin:
        request = json.loads(line)
        candidates = request.get("requests", [request])
        if not baseline_ready and any(
            item.get("engine") == "m12" for item in candidates
        ):
            # A dedicated control worker loads once outside operation timing.
            # Cold lifecycle samples still include verification/compilation.
            portfolio.advance_job = load_stage_jobs().advance_job
            baseline_ready = True

        if "warmup" in request:
            start = time.perf_counter()
            calls = 0

            while time.perf_counter() - start < request["warmup"]:
                for item, fixture in zip(
                    request["requests"], request["fixtures"]
                ):
                    _validate(_sample(item), fixture)
                    calls += 1

            result = {"seconds": time.perf_counter() - start, "calls": calls}
        else:
            result = _sample(request)

        print(json.dumps(result), flush=True)


class Worker:
    """Own a single-core subprocess, counting startup separately."""

    def __init__(self):
        self.process = subprocess.Popen(
            [
                sys.executable,
                "-u",
                "-m",
                "v2.benchmarks.phase_two",
                "--worker",
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        )

    def call(self, request):
        """Send one request; worker output is data, never instructions."""
        self.process.stdin.write(json.dumps(request) + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        if not line:
            error = self.process.stderr.read()
            raise RuntimeError(f"benchmark worker failed: {error}")
        return json.loads(line)

    def close(self):
        """Include shutdown in cold samples and close worker pipes."""
        self.process.stdin.close()
        self.process.wait(timeout=10)
        self.process.stdout.close()
        self.process.stderr.close()


def _validate(row, fixture):
    """Verify complete and partial answers against hidden certificates."""
    expected = dict(fixture["factors"])

    if not row["reconstructs"]:
        raise AssertionError("invalid reconstruction")
    for prime, exponent in row["factors"]:
        if prime not in expected or not 0 < exponent <= expected[prime]:
            raise AssertionError("invalid terminal factor")
    if row["complete"] and dict(row["factors"]) != expected:
        raise AssertionError("complete answer disagrees with hidden oracle")


def _quantile(values, percentile):
    """Nearest-rank percentile across all observations, including timeouts."""
    values = sorted(values)
    return values[max(0, math.ceil(percentile * len(values)) - 1)]


def summarize(rows):
    """Retain censored outcomes; never summarize successful times alone."""
    groups = {}
    for row in rows:
        key = (row["engine"], row["mode"], row["temperature"], row["band"])
        groups.setdefault(key, []).append(row)
    summaries = []

    for key, samples in groups.items():
        times = [sample["total_seconds"] for sample in samples]
        successes = sum(sample["success"] for sample in samples)
        summaries.append(
            {
                "engine": key[0],
                "mode": key[1],
                "temperature": key[2],
                "band": key[3],
                "observations": len(samples),
                "successes": successes,
                "completion_fraction": successes / len(samples),
                "median_seconds": statistics.median(times),
                "p90_seconds": _quantile(times, 0.90),
                "p95_seconds": _quantile(times, 0.95),
                "min_seconds": min(times),
                "max_seconds": max(times),
                "censored_timeouts": sum(
                    sample["reason"] in ("wall_limit", "cpu_limit")
                    for sample in samples
                ),
                "work_exhaustions": sum(
                    sample["reason"] == "work_limit" for sample in samples
                ),
                "cpu_seconds": sum(
                    sample["cpu_seconds"] for sample in samples
                ),
                "peak_rss_bytes": max(
                    sample["peak_rss_bytes"] for sample in samples
                ),
                "feasibility": "no completion under these caps"
                if not successes
                else "some completions; no promotion inference",
            }
        )

    return summaries


def run(args):
    """Evaluate fixed corpus IDs/seeds, with hidden factorization oracles."""
    measured_environment = environment()
    corpus_path = args.corpus
    corpus_bytes = corpus_path.read_bytes()
    corpus = json.loads(corpus_bytes)
    verify_certificates(corpus["certificates"])
    for fixture in corpus["fixtures"]:
        for prime, _ in fixture["factors"]:
            if str(prime) not in corpus["certificates"] and not all(
                prime % d for d in range(2, isqrt_integer(prime) + 1)
            ):
                raise AssertionError("uncertified oracle factor")

    fixtures = [
        fixture
        for fixture in corpus["fixtures"]
        if fixture["split"] == args.split
        and (not args.bands or fixture["band"] in args.bands.split(","))
    ]
    if not fixtures:
        raise ValueError("no matching fixtures")
    if args.config:
        config = PortfolioConfig(**json.loads(args.config.read_text()))
    else:
        config = PortfolioConfig(
            trial_bound=args.trial_bound,
            max_input_bits=512,
            ecm_tiers=((2000, 147396, args.curves),),
        )

    rows, warmups = [], []
    seeds = corpus["seeds"][: args.seeds]
    engines = args.engines.split(",")
    modes = args.modes.split(",")

    def request(engine, fixture, seed, mode):
        """Keep hidden oracle factors outside the timed algorithm request."""
        return {
            "n": fixture["n"],
            "seed": seed,
            "engine": engine,
            "mode": mode,
            "config": asdict(config),
            "work_limit": args.work_limit,
            "seconds": args.seconds,
            "cpu_seconds": args.cpu_seconds,
            "rss_cap_bytes": args.rss_mib * 1024 * 1024,
        }

    representatives = {}
    for fixture in fixtures:
        representatives.setdefault(fixture["band"], fixture)
    warmup_inputs = (
        fixtures
        if args.warmup_scope == "corpus"
        else list(representatives.values())
    )
    with ExitStack() as cleanup:
        workers = {}

        for engine in engines:
            worker = Worker()
            cleanup.callback(worker.close)
            workers[engine] = worker
            warmup_fixtures = [
                fixture for mode in modes for fixture in warmup_inputs
            ]
            warmup = worker.call(
                {
                    "warmup": args.warmup_seconds,
                    "requests": [
                        request(engine, fixture, seeds[0], mode)
                        for mode in modes
                        for fixture in warmup_inputs
                    ],
                    "fixtures": warmup_fixtures,
                }
            )
            warmups.append({"engine": engine, **warmup})

        # Startup and validated warmup precede the retained execution samples.
        # Only one worker executes at a time. Rotate comparison order while
        # preserving each interpreter's warmed process and JIT traces.
        for repetition in range(args.repetitions):
            offset = repetition % len(engines)

            for engine in engines[offset:] + engines[:offset]:
                worker = workers[engine]

                for mode in modes:
                    for fixture in fixtures:
                        for seed in seeds:
                            started = time.perf_counter()
                            row = worker.call(
                                request(engine, fixture, seed, mode)
                            )
                            total_seconds = time.perf_counter() - started
                            _validate(row, fixture)
                            row.update(
                                id=fixture["id"],
                                band=fixture["band"],
                                seed=seed,
                                engine=engine,
                                mode=mode,
                                temperature="warm",
                                repetition=repetition,
                                total_seconds=total_seconds,
                            )
                            rows.append(row)

                    print(engine, mode, repetition, "done", flush=True)

    # Cold startup includes imports, one operation, output, and shutdown.
    # Each band/mode has one representative under all declared seeds.
    for engine in engines:
        for mode in modes:
            for fixture in representatives.values():
                for seed in seeds:
                    rows.append(
                        _cold_sample(
                            request(engine, fixture, seed, mode), fixture
                        )
                    )

        print(engine, "cold samples done", flush=True)

    if environment()["source_sha256"] != measured_environment["source_sha256"]:
        raise RuntimeError("source changed during measurements")
    return {
        "milestone": "phase_two",
        "environment": measured_environment,
        "corpus_sha256": hashlib.sha256(corpus_bytes).hexdigest(),
        "corpus": str(corpus_path),
        "m12_snapshot_sha256": hashlib.sha256(
            CORPUS.parents[1]
            .joinpath("benchmarks/inputs/baselines/m12_source_snapshot.json")
            .read_bytes()
        ).hexdigest()
        if "m12" in engines
        else None,
        "split": args.split,
        "config": asdict(config),
        "seeds": seeds,
        "wall_cap_seconds": args.seconds,
        "cpu_cap_seconds": args.cpu_seconds,
        "rss_cap_bytes": args.rss_mib * 1024 * 1024,
        "work_cap_bounded_only": args.work_limit,
        "repetitions": args.repetitions,
        "warmups": warmups,
        "warmup_scope": args.warmup_scope,
        "competitors": json.loads(
            (
                CORPUS.parent.parent / "controls/phase_two_competitors.json"
            ).read_text()
        ),
        "limitations": [
            "RSS cap is a measured acceptance gate, not a hard OS limit",
            "Unix wall/CPU watchdog enforces matching time caps on both arms",
            "Legacy interruption preserves only the untouched original input",
            "Work units apply only to bounded; wall/CPU/RSS caps are matched",
            "Worker RSS is cumulative and includes warmup",
            "Cold samples include process lifecycle and IPC",
            "Cold runs cover one representative per band/mode and all seeds",
            "Comparison order rotates by repetition; active work is serial",
            "Baseline full factoring retains its rho cutoff and omits p-1",
            "Competitor revisions are pinned but not installed/executed here",
            "Failed outputs receive no speed ratios; defaults are unpromoted",
        ],
        "summaries": summarize(rows),
        "samples": rows,
    }


def _cold_sample(request, fixture):
    """Include parent CPU, worker startup, IPC, and shutdown in a cold run."""
    started = time.perf_counter()
    cpu_started = time.process_time()
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    child_cpu_started = children.ru_utime + children.ru_stime
    cold = Worker()
    try:
        row = cold.call(request)
        _validate(row, fixture)
    finally:
        cold.close()

    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    row["operation_cpu_seconds"] = row["cpu_seconds"]
    row["cpu_seconds"] = (
        time.process_time()
        - cpu_started
        + children.ru_utime
        + children.ru_stime
        - child_cpu_started
    )
    row.update(
        id=fixture["id"],
        band=fixture["band"],
        seed=request["seed"],
        engine=request["engine"],
        mode=request["mode"],
        temperature="cold",
        repetition=0,
        total_seconds=time.perf_counter() - started,
    )
    return row


def isqrt_integer(n):
    """Exact oracle root kept separate from Factor's helper dispatch."""
    from math import isqrt

    return isqrt(n)


def main():
    """Run an isolated benchmark or service its private worker protocol."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument(
        "--split", choices=("training", "held_out"), default="held_out"
    )
    parser.add_argument("--bands")
    parser.add_argument("--config", type=Path)
    parser.add_argument("--engines", default="phase_one,bounded")
    parser.add_argument("--modes", default="complete,factor_one")
    parser.add_argument("--seeds", type=int, default=5)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument(
        "--warmup-scope",
        choices=("representatives", "corpus"),
        default="representatives",
        help="warm one input per band or every selected input before timing",
    )
    parser.add_argument("--seconds", type=float, default=0.05)
    parser.add_argument("--cpu-seconds", type=float, default=0.05)
    parser.add_argument("--rss-mib", type=int, default=256)
    parser.add_argument("--work-limit", type=int, default=2_000_000)
    parser.add_argument("--trial-bound", type=int, default=25_000)
    parser.add_argument("--curves", type=int, default=2)
    args = parser.parse_args()
    if args.worker:
        _worker()
        return
    if args.output is None or args.output.exists():
        parser.error("select a fresh --output filename")
    if any(
        not math.isfinite(value) or value <= 0
        for value in (args.seconds, args.cpu_seconds, args.warmup_seconds)
    ):
        parser.error("time caps and warmup must be finite and positive")

    if not 1 <= args.seeds <= 5 or args.repetitions < 1 or args.rss_mib < 1:
        parser.error("invalid seeds, repetitions, or RSS cap")
    if not set(args.engines.split(",")) <= {"phase_one", "bounded", "m12"}:
        parser.error("unknown engine")
    if not set(args.modes.split(",")) <= {"complete", "factor_one"}:
        parser.error("unknown mode")
    args.output.write_text(json.dumps(run(args), indent=2) + "\n")


if __name__ == "__main__":
    main()
