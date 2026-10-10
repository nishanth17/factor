"""B3 equal-schedule follow-up with a common-prefix cooperative pause."""

import argparse
import fcntl
import json
import os
import resource
import signal
import subprocess
import sys
import time
from math import prod
from pathlib import Path

from . import b3_production as study
from .a6_pm1 import assert_quiet, require_runtime
from .build_b3_inputs import write_new
from .build_c6_fast_inputs import LimitedRandom
from .build_phase_two_corpus import certified_prime, verify_certificates

PROTOCOL = study.INPUTS / "controls/b3_coverage_protocol.json"
CORPUS = study.INPUTS / "corpora/b3_coverage_confirmation.json"
WORK_LIMIT = 8_000_000
SCOPES = ("portfolio", "resume_equal")


def primary_inputs():
    """Read immutable initial inputs; its trial source remains historical."""
    settings = json.loads(study.PROTOCOL.read_text())
    for path in (study.CORPUS, study.BASELINE):
        name = str(path.relative_to(study.ROOT))
        if study.digest(path) != settings["sha256"][name]:
            raise ValueError("primary B3 input changed")
    corpus = json.loads(study.CORPUS.read_text())
    verify_certificates(corpus["certificates"])
    return settings, corpus


def freeze():
    """Add untouched inputs without altering primary sources or captures."""
    if PROTOCOL.exists() or CORPUS.exists():
        raise ValueError("coverage follow-up is already frozen")
    original, original_corpus = primary_inputs()
    generator, certificates, fixtures = LimitedRandom(108841), {}, []

    def prime(digits):
        for _ in range(64):
            value = certified_prime(
                (10**digits - 1).bit_length(), generator, certificates
            )
            if len(str(value)) == digits:
                return value
        raise RuntimeError("coverage prime draw cap")

    for digits in (40, 50, 60, 70, 80):
        for shape in ("balanced", "small10"):
            target = digits // 2 if shape == "balanced" else 10
            for _ in range(64):
                factors = [prime(target), prime(digits - target)]
                if (
                    len(str(prod(factors))) == digits
                    and len(set(factors)) == 2
                ):
                    break
            else:
                raise RuntimeError("coverage composite draw cap")
            fixtures.append(
                dict(
                    id=f"{shape}_{digits}d",
                    n=prod(factors),
                    factors=factors,
                    factor_digits=[target, digits - target],
                )
            )
    for index, sizes in enumerate(((6, 6, 8), (6, 10, 20))):
        factors = [prime(size) for size in sizes]
        fixtures.append(
            dict(
                id=f"recursive_{index}",
                n=prod(factors),
                factors=factors,
                factor_digits=list(sizes),
            )
        )
    verify_certificates(certificates)
    old = {f["n"] for f in original_corpus["fixtures"]}
    for path in CORPUS.parent.glob("*.json"):
        data = json.loads(path.read_text())
        if isinstance(data, dict):
            old.update(
                f["n"]
                for f in data.get("fixtures", [])
                if isinstance(f, dict) and "n" in f
            )
    values = {f["n"] for f in fixtures}
    if old & values or len(values) != len(fixtures):
        raise ValueError("coverage confirmation inputs overlap")
    data = dict(
        generation_seed=108841, certificates=certificates, fixtures=fixtures
    )
    if len(json.dumps(data).encode()) > 16777216:
        raise MemoryError("coverage corpus output cap")
    write_new(CORPUS, data)
    write_new(
        PROTOCOL,
        dict(
            parent_sha256=study.digest(study.PROTOCOL),
            production_sha256={
                name: study.digest(study.ROOT / name)
                for name in original["sha256"]
            },
            routing=(
                "Explicit reuse, B1=2000, chunk16, tier >=8, "
                "40-80 digits; off default. Unsupported small recursive "
                "cofactors retain B4."
            ),
            source_sha256=study.digest(Path(__file__)),
            corpus_sha256=study.digest(CORPUS),
            work=WORK_LIMIT,
            seconds=20,
            cpu_seconds=20,
            memory_bytes=33554432,
            b1=2000,
            b2=147396,
            curves=8,
            training_seeds=original["training_seeds"],
            confirmation_seeds=[62677, 70596],
            scopes=list(SCOPES),
            sampling=original["sampling"],
            max_relative_iqr=0.15,
            bootstrap_seed=193001,
            reason=(
                "Primary equal-cap runs can end early under conservative "
                "chain reservations. Keep those captures, but require eight "
                "completed root curves on every unresolved matched-schedule "
                "case before a full-run speed claim. Eight million units "
                "allow both existing eight-curve campaigns and the at-most-"
                "three-factor recursive fixtures without changing time/RAM "
                "caps, sources, chains, bounds or routing. "
                "No candidate tuning."
            ),
            pause=(
                "Both arms cooperatively cancel after the first certified "
                "16-prime chunk on balanced inputs, then JSON roundtrip and "
                "resume under the same cumulative cap. Earlier splits fail "
                "this equal-prefix diagnostic rather than replacing inputs."
            ),
            confirmation=(
                "New independent twelve-input corpus and repeated seeds "
                "frozen before any follow-up timings. Primary data "
                "and decisions remain separate. Revised positive CI, "
                "CPU, chronological halves and stability gates apply."
            ),
            generation_limits=dict(
                random_draws=500000,
                decimal_attempts=64,
                wall_seconds=60,
                output_bytes=16777216,
            ),
        ),
    )


def inputs(split):
    original, original_corpus = primary_inputs()
    settings = json.loads(PROTOCOL.read_text())
    if (
        settings["parent_sha256"] != study.digest(study.PROTOCOL)
        or settings["source_sha256"] != study.digest(Path(__file__))
        or settings["corpus_sha256"] != study.digest(CORPUS)
    ):
        raise ValueError("coverage follow-up source/input changed")
    for name, expected in settings["production_sha256"].items():
        if study.digest(study.ROOT / name) != expected:
            raise ValueError("routed production source/input changed: " + name)
    corpus = json.loads(CORPUS.read_text())
    verify_certificates(corpus["certificates"])
    for fixture in corpus["fixtures"]:
        if prod(fixture["factors"]) != fixture["n"]:
            raise ValueError("coverage corpus reconstruction failed")
    fixtures = (
        corpus["fixtures"]
        if split == "confirmation"
        else [
            f for f in original_corpus["fixtures"] if f["split"] == "training"
        ]
    )
    return settings, fixtures


def allowance(engine, cancelled=None):
    return engine.Budget(
        work_limit=WORK_LIMIT, seconds=20, cpu_seconds=20, cancelled=cancelled
    )


def equal_prefix_pause(engine, fixture, config, seed):
    """Trigger the public cancellation callback at the same scalar prefix."""
    stopped = [False]
    original = engine.advance_job

    def advance(job, budget, context, selected):
        original(job, budget, context, selected)
        if (
            job["kind"] == "ecm"
            and not job["done"]
            and job["phase"] == "stage_one"
            and not job["powers"]
            and job["cursor"]["index"] == 16
        ):
            stopped[0] = True

    engine.advance_job = advance
    try:
        partial = engine.factorize_bounded(
            fixture["n"],
            seed=seed,
            config=config,
            budget=allowance(engine, lambda: stopped[0]),
        )
    finally:
        engine.advance_job = original
    study.validate(partial, fixture)
    if partial.reason != "cancelled":
        raise AssertionError("equal-prefix pause was not reached")
    current = partial.checkpoint["payload"]["state"]["current"]
    if (
        current["attempt"] != 0
        or current["job"]["powers"]
        or current["job"]["cursor"]["index"] != 16
    ):
        raise AssertionError("pause did not retain the first certified chunk")
    return json.loads(json.dumps(partial.checkpoint))


def cohort(fixtures, seeds, backend, candidate, scope):
    engine = study.portfolio if candidate else study.baseline()
    config = study.configuration(engine, backend, candidate)
    rows = []
    for fixture in fixtures:
        if scope == "resume_equal" and not fixture["id"].startswith(
            "balanced"
        ):
            continue
        for seed in seeds:
            started = time.perf_counter()
            if scope == "resume_equal":
                checkpoint = equal_prefix_pause(engine, fixture, config, seed)
                run = engine.factorize_bounded(
                    fixture["n"],
                    config=config,
                    checkpoint=checkpoint,
                    budget=allowance(engine),
                )
            else:
                run = engine.factorize_bounded(
                    fixture["n"],
                    seed=seed,
                    config=config,
                    budget=allowance(engine),
                )
            row = study.validate(run, fixture)
            root_curves = sum(
                event.get("stage") == "ecm" and event.get("n") == fixture["n"]
                for event in run.events
            )
            if not run.result.complete and (
                run.reason != "exhausted" or root_curves != 8
            ):
                raise AssertionError("unequal completed curve schedule")
            rows.append(
                dict(
                    row,
                    seed=seed,
                    root_curves=root_curves,
                    certified_pause_prefix=16
                    if scope == "resume_equal"
                    else None,
                    seconds=time.perf_counter() - started,
                )
            )
    return rows


def worker(args):
    settings, fixtures = inputs(args.split)
    seeds = settings[args.split + "_seeds"]
    signatures = [None, None]

    def measure(arm):
        start, cpu = time.perf_counter(), time.process_time()
        rows = cohort(fixtures, seeds, args.backend, bool(arm), args.scope)
        elapsed, used = time.perf_counter() - start, time.process_time() - cpu
        signature = json.dumps(
            [{k: v for k, v in row.items() if k != "seconds"} for row in rows],
            sort_keys=True,
        )
        if signatures[arm] is not None and signatures[arm] != signature:
            raise AssertionError("coverage outcome/work changed")
        signatures[arm] = signature
        return dict(wall=elapsed, cpu=used, rows=rows)

    pairs, warmups = [], []
    for seconds, count in settings["sampling"]:
        for arm in (0, 1):
            started, rounds = time.perf_counter(), 0
            while time.perf_counter() - started < seconds:
                measure(arm)
                rounds += 1
            warmups.append(
                dict(
                    arm=arm,
                    seconds=time.perf_counter() - started,
                    cohorts=rounds,
                )
            )
        while len(pairs) < count:
            pair = [None, None]
            for arm in (0, 1) if len(pairs) % 2 == 0 else (1, 0):
                pair[arm] = measure(arm)
            pairs.append(pair)
        summary = study.summarize(pairs, settings)
        if summary["stable"]:
            break
    return dict(
        scope=args.scope,
        split=args.split,
        backend=args.backend,
        pairs=pairs,
        warmups=warmups,
        summary=summary,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", action="store_true")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--scope", choices=SCOPES, default="portfolio")
    parser.add_argument(
        "--split", choices=("training", "confirmation"), default="training"
    )
    parser.add_argument(
        "--backend", choices=("python-int", "gmpy2-mpz"), default="python-int"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require_runtime()
    if args.freeze:

        def expired(signum, frame):
            raise TimeoutError("coverage generation wall limit")

        signal.signal(signal.SIGALRM, expired)
        signal.alarm(60)
        try:
            freeze()
        finally:
            signal.alarm(0)
        return
    if args.worker:
        print(json.dumps(worker(args)))
        return
    if args.output is None or args.output.exists():
        parser.error("choose a new output path")
    inputs(args.split)
    owner = Path("/private/tmp/factor-performance-owner.json")
    with open("/private/tmp/factor-performance.lock", "a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert_quiet()
        owner.write_text(
            json.dumps(dict(owner="B3 equal schedule", pid=os.getpid()))
        )
        try:
            command = [
                sys.executable,
                "-B",
                "-m",
                __spec__.name,
                "--worker",
                "--scope",
                args.scope,
                "--split",
                args.split,
                "--backend",
                args.backend,
            ]
            result = subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
                timeout=1200,
            )
            data = json.loads(result.stdout)
            inputs(args.split)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(data, indent=2) + "\n")
            print(json.dumps(data["summary"]))
        finally:
            owner.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
