"""Frozen A6 p-1 stages, continuations and portfolio measurements."""

import argparse
import contextlib
import fcntl
import hashlib
import json
import math
import os
import platform
import random
import statistics
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from v2 import portfolio, stage_jobs, utils
from v2.benchmarks.build_phase_two_corpus import verify_certificates
from v2.budget import Budget
from v2.pm1_bounded import PM1Config, factorize_pm1_bounded
from v2.schedules import ScheduleCache, SieveContext

INPUTS = Path(__file__).parent / "inputs/corpora"
PROTOCOL = INPUTS / "a6_pm1_protocol.json"


@contextlib.contextmanager
def performance_window():
    """Use the cross-chat lock and fail closed on competing heavy processes."""
    owner = Path("/private/tmp/factor-performance-owner.json")
    with open("/private/tmp/factor-performance.lock", "a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner.write_text(
            json.dumps(
                {"owner": "A6", "pid": os.getpid(), "started": time.time()}
            )
        )
        try:
            assert_quiet()
            yield
            assert_quiet()
        finally:
            owner.unlink(missing_ok=True)


def assert_quiet():
    processes = subprocess.check_output(
        ["ps", "-axo", "pid=,command="], text=True
    )
    for line in processes.splitlines():
        fields = line.strip().split(None, 1)
        if len(fields) != 2 or int(fields[0]) == os.getpid():
            continue
        command = fields[1]
        executable = Path(command.split()[0]).name.lower()
        interpreter = executable.startswith(("pypy", "python"))
        heavy = any(
            name in command
            for name in (
                "v2.benchmarks",
                "unittest",
                "test_a6",
                "test_c6",
                "test_b4",
            )
        )
        if interpreter and heavy:
            raise RuntimeError(
                "competing heavy process during timing: " + line
            )


def require_runtime():
    if platform.python_implementation() != "PyPy" or sys.version_info[:2] != (
        3,
        11,
    ):
        raise RuntimeError("A6 requires PyPy implementing Python 3.11")


def budget():
    return Budget(work_limit=2_000_000, seconds=30, cpu_seconds=30)


def stage(n, b1, b2, *, base=2, chunk=16, context=None, gaps="cached"):
    config = portfolio.PortfolioConfig(
        trial_bound=2,
        rho_attempts=0,
        pm1_b1=b1,
        pm1_b2=b2,
        ecm_tiers=(),
        chunk_size=chunk,
        gcd_batch=64,
        segment_size=256,
    )
    ledger = budget()
    if context is None:
        ledger.consume((math.isqrt(b2) + 1) // 2)
        context = SieveContext(b2 + 1, segment_size=256)
    job = stage_jobs.new_job("pm1", n, base - 2, b1, b2)
    while not job["done"]:
        if gaps == "direct" and job["phase"] in ("stage_two", "term_replay"):
            direct_stage_two(job, ledger, context, config)
        else:
            stage_jobs.advance_job(job, ledger, context, config)
    divisor = job["factor"]
    if divisor is not None and not utils.valid_divisor(divisor, n):
        raise AssertionError("invalid complete-stage output")
    return {"factor": divisor, "work": ledger.used, "phase": job["phase"]}


def direct_stage_two(job, ledger, context, config):
    """Independent a**q oracle with the same relation batches and replay."""
    if job["phase"] == "term_replay" or len(job["terms"]) == config.gcd_batch:
        stage_jobs._batch_check(job, ledger)
        return
    prime = stage_jobs.peek_prime(job["cursor"], context, ledger)
    if prime is None:
        stage_jobs._batch_check(job, ledger)
        if not job["terms"] and not job["done"]:
            stage_jobs._finish(job)
        return
    cursor = job["cursor"]
    count = min(
        config.gcd_batch - len(job["terms"]),
        len(cursor["values"]) - cursor["index"],
    )
    start = cursor["index"]
    stop = start + count
    primes = cursor["values"][start:stop]
    ledger.consume(sum(prime.bit_length() + 1 for prime in primes))
    for prime in primes:
        term = (pow(job["value"], prime, job["n"]) - 1) % job["n"]
        job["terms"].append(term)
        job["product"] = job["product"] * term % job["n"]
    cursor["index"] += count


def stage_fixtures(split="screen"):
    data = json.loads((INPUTS / "p43_size_corpus.json").read_text())
    return [
        f
        for f in data["fixtures"]
        if f["split"] == split and f["digits"] in (20, 50, 100)
    ]


def verify_inputs():
    """Check independent prime proofs and all selected product identities."""
    for name in (
        "p43_size_corpus.json",
        "phase_two_m15_independent_corpus.json",
    ):
        data = json.loads((INPUTS / name).read_text())
        verify_certificates(data["certificates"])
        for fixture in data["fixtures"]:
            factors = fixture["factors"]
            value = (
                math.prod(p**e for p, e in factors)
                if isinstance(factors[0], list)
                else math.prod(factors)
            )
            if value != fixture["n"]:
                raise AssertionError("invalid frozen product identity")


def profile():
    started = time.monotonic()
    while time.monotonic() - started < 3:
        for fixture in stage_fixtures():
            stage(fixture["n"], 2000, 20000)
    records = []
    for fixture in stage_fixtures():
        counters = {}

        def wrap(name, function):
            def measured(*args, **kwargs):
                start = time.process_time()
                result = function(*args, **kwargs)
                entry = counters.setdefault(name, {"calls": 0, "cpu": 0.0})
                entry["calls"] += 1
                entry["cpu"] += time.process_time() - start
                return result

            return measured

        with contextlib.ExitStack() as stack:
            for name in ("pow", "gcd"):
                stack.enter_context(
                    patch.object(
                        stage_jobs, name, wrap(name, getattr(stage_jobs, name))
                    )
                )
            stack.enter_context(
                patch.object(
                    utils,
                    "prime_power",
                    wrap("prime_power", utils.prime_power),
                )
            )
            stack.enter_context(
                patch.object(
                    SieveContext,
                    "prime_segment",
                    wrap("sieve", SieveContext.prime_segment),
                )
            )
            stack.enter_context(
                patch.object(
                    stage_jobs,
                    "_stage_one",
                    wrap("stage_one_inclusive", stage_jobs._stage_one),
                )
            )
            stack.enter_context(
                patch.object(
                    stage_jobs,
                    "_stage_two",
                    wrap("stage_two_inclusive", stage_jobs._stage_two),
                )
            )
            start = time.process_time()
            outcome = stage(fixture["n"], 2000, 20000)
            total = time.process_time() - start
        records.append(
            {
                "id": fixture["id"],
                "total_cpu": total,
                "instrumented": counters,
                "outcome": outcome,
            }
        )
    return records


def _measure_arms(functions, samples=9):
    """Alternate frozen arm order and extend the entire comparison together."""
    generator = random.Random(2026100906)
    loops, values, warmups = {}, {name: [] for name in functions}, {}
    for name, function in functions.items():
        start, count = time.monotonic(), 0
        while time.monotonic() - start < 3:
            function()
            count += 1
        warmups[name] = count
        loops[name] = 1
        while True:
            start = time.process_time()
            for _ in range(loops[name]):
                function()
            if time.process_time() - start >= 0.08:
                break
            loops[name] *= 2
    target = samples
    while len(next(iter(values.values()))) < target:
        order = list(functions)
        generator.shuffle(order)
        for name in order:
            assert_quiet()
            start = time.process_time()
            for _ in range(loops[name]):
                functions[name]()
            values[name].append((time.process_time() - start) / loops[name])
            assert_quiet()
        if len(next(iter(values.values()))) == samples:
            for series in values.values():
                q = statistics.quantiles(series, n=4)
                if (q[2] - q[0]) / statistics.median(series) > 0.15:
                    target = max(target, 27)
    return {
        name: {
            "cpu_samples": series,
            "median_cpu": statistics.median(series),
            "loops": loops[name],
            "warmup_iterations": warmups[name],
            "extended": target > samples,
        }
        for name, series in values.items()
    }


def stage_screen(samples=9, split="screen"):
    records = []
    for fixture in stage_fixtures(split):
        n = fixture["n"]
        for b1, b2 in ((2000, 20000), (11000, 100000)):
            functions = {}
            for chunk in (1, 16, 64):
                functions[str(chunk)] = lambda chunk=chunk: stage(
                    n, b1, b2, chunk=chunk
                )
            functions["direct_q"] = lambda: stage(n, b1, b2, gaps="direct")
            measurements = _measure_arms(functions, samples)
            for arm, function in functions.items():
                records.append(
                    {
                        "id": fixture["id"],
                        "bounds": [b1, b2],
                        "arm": arm,
                        "outcome": function(),
                        **measurements[arm],
                    }
                )
    return records


def reuse_screen(samples=9):
    records = []
    for fixture in stage_fixtures():

        def run(reuse):
            context = SieveContext(20001, segment_size=256)
            if reuse:
                context = ScheduleCache(context, cache_bytes=262144)
            outcomes = [
                stage(fixture["n"], 2000, 20000, base=base, context=context)
                for base in (2, 3, 4)
            ]
            return {
                "stages": outcomes,
                "retained_bytes": getattr(context, "used_bytes", 0),
                "hits": getattr(context, "hits", 0),
            }

        functions = {
            str(reuse): lambda reuse=reuse: run(reuse)
            for reuse in (False, True)
        }
        measurements = _measure_arms(functions, samples)
        for arm, function in functions.items():
            records.append(
                {
                    "id": fixture["id"],
                    "reuse": arm,
                    "outcomes": function(),
                    **measurements[arm],
                }
            )
    return records


def continuations(samples=9):
    records = []
    for fixture in stage_fixtures():
        n = fixture["n"]
        config = PM1Config(
            bounds=((500, 5000), (2000, 20000)), segment_size=256
        )
        final = replace(config, bounds=((2000, 20000),))
        first = replace(config, bounds=((500, 5000),))
        from v2.pm1_bounded import _advance, _initial

        state = _initial(n, 2, config)
        context = SieveContext(20001, segment_size=256)
        ledger = budget()
        while not state["done"] and state["rung"] == 0:
            _advance(state, ledger, context, config)
        boundary_actions = state["steps"]

        def run(arm):
            prior_work = 0
            if arm == "fresh_final":
                result = factorize_pm1_bounded(
                    n, config=final, budget=budget()
                )
            elif arm == "fresh_each":
                result = factorize_pm1_bounded(
                    n, config=first, budget=budget()
                )
                if result.divisor is None and result.reason == "exhausted":
                    prior_work = result.work_used
                    # Independent fresh stages share the campaign's total
                    # allowance, just as in-memory and resumed execution do.
                    remaining = Budget(
                        work_limit=2_000_000 - prior_work,
                        seconds=max(0, 30 - result.wall_seconds),
                        cpu_seconds=max(0, 30 - result.cpu_seconds),
                    )
                    result = factorize_pm1_bounded(
                        n, config=final, budget=remaining
                    )
            elif arm == "in_memory":
                result = factorize_pm1_bounded(
                    n, config=config, budget=budget()
                )
            else:
                paused = factorize_pm1_bounded(
                    n,
                    config=config,
                    budget=budget(),
                    max_actions=boundary_actions,
                )
                result = factorize_pm1_bounded(
                    n,
                    config=config,
                    budget=budget(),
                    checkpoint=paused.checkpoint,
                )
            if result.result.reconstruct() != n:
                raise AssertionError("continuation lost a cofactor")
            return {
                "factor": result.divisor,
                "reason": result.reason,
                "work": prior_work + result.work_used,
                "verify_work": result.verification_work,
                "bytes": len(json.dumps(result.checkpoint).encode()),
            }

        functions = {
            arm: lambda arm=arm: run(arm)
            for arm in ("fresh_final", "fresh_each", "in_memory", "checkpoint")
        }
        measurements = _measure_arms(functions, samples)
        for arm, function in functions.items():
            records.append(
                {
                    "id": fixture["id"],
                    "arm": arm,
                    "outcome": function(),
                    **measurements[arm],
                }
            )
    return records


def portfolio_screen(samples=9):
    data = json.loads(
        (INPUTS / "phase_two_m15_independent_corpus.json").read_text()
    )
    fixtures = [f for f in data["fixtures"] if len(str(f["n"])) == 20][:12]
    records = []

    def run(attempts):
        config = portfolio.PortfolioConfig(
            trial_bound=100,
            rho_attempts=1,
            rho_evaluations=512,
            pm1_attempts=attempts,
            pm1_b1=2000,
            pm1_b2=20000,
            ecm_tiers=((50, 1000, 2),),
            segment_size=256,
        )
        completed, splits = 0, 0
        for fixture in fixtures:
            outcome = portfolio.factorize_bounded(
                fixture["n"],
                seed=7,
                config=config,
                budget=Budget(work_limit=500000, seconds=5, cpu_seconds=5),
            )
            assert outcome.result.reconstruct() == fixture["n"]
            completed += outcome.result.complete
            splits += len(outcome.result.remaining) != 1 or bool(
                outcome.result.factors
            )
        return {
            "complete": completed,
            "split": splits,
            "attempts": len(fixtures),
        }

    functions = {
        str(attempts): lambda attempts=attempts: run(attempts)
        for attempts in (0, 1)
    }
    measurements = _measure_arms(functions, samples)
    for arm, function in functions.items():
        records.append(
            {
                "pm1_attempts": int(arm),
                "outcome": function(),
                **measurements[arm],
            }
        )
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mode",
        choices=(
            "profile",
            "stages",
            "reuse",
            "continuation",
            "portfolio",
            "cold",
        ),
    )
    parser.add_argument(
        "--split", default="screen", choices=("screen", "confirmation")
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require_runtime()
    functions = {
        "profile": profile,
        "stages": lambda: stage_screen(split=args.split),
        "reuse": reuse_screen,
        "continuation": continuations,
        "portfolio": portfolio_screen,
    }
    with performance_window():
        verify_inputs()
        if args.mode == "cold":
            samples = []
            for _ in range(9):
                start = time.monotonic()
                subprocess.run(
                    [
                        sys.executable,
                        "-c",
                        "from v2.pm1_bounded import "
                        "factorize_pm1_bounded, PM1Config; "
                        "r=factorize_pm1_bounded(618533,"
                        "config=PM1Config(bounds=((10,200),))); "
                        "assert r.result.reconstruct()==618533 and r.divisor",
                    ],
                    check=True,
                    capture_output=True,
                )
                samples.append(time.monotonic() - start)
            records = {"startup_wall_samples": samples}
        else:
            records = functions[args.mode]()
    output = {
        "runtime": sys.version,
        "platform": platform.platform(),
        "pid": os.getpid(),
        "mode": args.mode,
        "protocol_sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
        "runner_sha256": hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        "source_sha256": {
            str(path.relative_to(Path(__file__).parents[2])): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in sorted(Path(__file__).parents[1].glob("*.py"))
        },
        "input_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (
                PROTOCOL,
                INPUTS / "p43_size_corpus.json",
                INPUTS / "phase_two_m15_independent_corpus.json",
            )
        },
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
