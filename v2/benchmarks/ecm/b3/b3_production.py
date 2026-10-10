"""Frozen current-mainline versus bounded production chain integration."""

import argparse
import fcntl
import hashlib
import json
import os
import random
import resource
import statistics
import subprocess
import sys
import time
from math import prod
from pathlib import Path

from .... import portfolio
from ....common import arithmetic, utils
from ....ecm import core as ecm
from ....ecm.chains import MIN_MEMORY_BYTES, ChainPlans
from ....execution.schedules import SieveContext
from ...pm1.a6.a6_pm1 import assert_quiet, require_runtime
from ...suites.build_phase_two_corpus import verify_certificates
from ...support.paths import (
    BENCHMARK_ROOT,
    REPOSITORY_ROOT,
    matches_source_pin,
    source_path,
)
from ...support.prac_oracle import affine_multiply, matches, twist_point
from ..b4.b4_common import FrozenFinder

ROOT = REPOSITORY_ROOT
INPUTS = BENCHMARK_ROOT / "inputs"
PROTOCOL = INPUTS / "controls/b3_protocol_v3.json"
CORPUS = INPUTS / "corpora/b3_corpus.json"
BASELINE = INPUTS / "baselines/b3_mainline.json"
SCOPES = (
    "fresh1",
    "reuse2",
    "reuse4",
    "reuse8",
    "reuse16",
    "eviction",
    "portfolio",
    "resume",
)


def digest(path):
    return hashlib.sha256(source_path(path).read_bytes()).hexdigest()


def inputs():
    settings = json.loads(source_path(PROTOCOL).read_text())
    for name, expected in settings["sha256"].items():
        if not matches_source_pin(ROOT / name, expected):
            raise ValueError("frozen B3 source/input changed: " + name)
    corpus = json.loads(source_path(CORPUS).read_text())
    verify_certificates(corpus["certificates"])
    for fixture in corpus["fixtures"]:
        if prod(fixture["factors"]) != fixture["n"]:
            raise ValueError("corpus reconstruction failed")
    return settings, corpus


def baseline():
    name = "_b3_control"
    data = json.loads(source_path(BASELINE).read_text())
    for path, source in data["source"].items():
        if hashlib.sha256(source.encode()).hexdigest() != data["sha256"][path]:
            raise ValueError("corrupt B3 mainline snapshot")
    if name not in sys.modules:
        sys.meta_path.insert(0, FrozenFinder(name, data["source"]))
    return __import__(name + ".portfolio", fromlist=["portfolio"])


def allowance(engine=portfolio, work=2_000_000):
    return engine.Budget(work_limit=work, seconds=20, cpu_seconds=20)


def configuration(engine, backend, candidate, curves=8):
    options = dict(
        backend=backend,
        trial_bound=5,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=((2000, 147396, curves),),
        chunk_size=16,
        segment_size=128,
        max_input_bits=512,
        memory_bytes=32 * 1024**2,
        ecm_program_bytes=262144,
    )
    if candidate:
        options.update(
            ecm_chain_mode="reuse", ecm_chain_bytes=MIN_MEMORY_BYTES
        )
    return engine.PortfolioConfig(**options)


def validate(run, fixture):
    if run.result.reconstruct() != fixture["n"]:
        raise AssertionError("portfolio reconstruction failed")
    known = fixture["factors"]
    for factor in run.result.factors:
        if factor.value not in known:
            raise AssertionError("unknown terminal factor")
        if factor.certainty.value == utils.Primality.COMPOSITE.value:
            raise AssertionError("invalid certainty")
    return dict(
        fixture=fixture["id"],
        complete=run.result.complete,
        factors=[
            (f.value, f.exponent, f.certainty.value)
            for f in run.result.factors
        ],
        unresolved=list(run.result.remaining),
        reason=run.reason,
        work=run.checkpoint["payload"]["work_used"],
    )


def portfolio_cohort(fixtures, seeds, backend, candidate, resume):
    engine = portfolio if candidate else baseline()
    config = configuration(engine, backend, candidate)
    rows = []
    for fixture in fixtures:
        for seed in seeds:
            row_started = time.perf_counter()
            if resume:
                partial = engine.factorize_bounded(
                    fixture["n"],
                    seed=seed,
                    config=config,
                    budget=allowance(engine, 250000),
                )
                validate(partial, fixture)
                # JSON packing/unpacking and output checks belong to the
                # timed complete resumed call, including failed preparations.
                checkpoint = json.loads(json.dumps(partial.checkpoint))
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
            rows.append(
                dict(
                    validate(run, fixture),
                    seed=seed,
                    seconds=time.perf_counter() - row_started,
                )
            )
    return rows


def stage_cohort(fixtures, seeds, backend, candidate, scope, targets):
    engine = portfolio if candidate else baseline()
    jobs_module = __import__(
        engine.__package__ + ".stage_jobs", fromlist=["stage_jobs"]
    )
    programs_module = __import__(
        engine.__package__ + ".ecm_programs", fromlist=["ecm_programs"]
    )
    config = configuration(engine, backend, False)
    integer = arithmetic.get_backend(backend).integer
    curves = (
        3
        if scope == "eviction"
        else int(scope.replace("fresh", "").replace("reuse", ""))
    )
    rows = []
    for fixture in fixtures:
        if not fixture["id"].startswith("balanced"):
            continue
        for seed in seeds:
            context = SieveContext(147397, segment_size=128)
            store = programs_module.ECMPrograms(context, memory_bytes=262144)
            if candidate:
                store.chains = ChainPlans(
                    MIN_MEMORY_BYTES, backend, (1999, 2000)
                )
            for curve in range(curves):
                row_started = time.perf_counter()
                bound = 1999 if scope == "eviction" and curve == 1 else 2000
                ledger = allowance(engine)
                job = jobs_module.new_job(
                    "ecm", integer(fixture["n"]), seed + curve, bound, 147396
                )
                while not job["done"] and job["phase"] in (
                    "setup",
                    "stage_one",
                    "replay",
                ):
                    jobs_module.advance_job(job, ledger, store, config)
                factor = job["factor"]
                if factor is not None and not utils.valid_divisor(
                    factor, fixture["n"]
                ):
                    raise AssertionError("invalid stage factor")
                unresolved = (
                    fixture["n"] if factor is None else fixture["n"] // factor
                )
                if (factor or 1) * unresolved != fixture["n"]:
                    raise AssertionError("stage reconstruction failed")
                # Unit-Z output is nondegenerate modulo every hidden prime.
                point = job.get("value")
                if not job["done"]:
                    for prime, target in targets[(fixture["id"], seed, curve)]:
                        if not matches(point, target, prime):
                            raise AssertionError(
                                "independent affine stage mismatch"
                            )
                    if arithmetic.gcd(point[1], job["n"]) != 1:
                        raise AssertionError("uncertified stage continuation")
                rows.append(
                    dict(
                        fixture=fixture["id"],
                        seed=seed,
                        curve=curve,
                        factor=factor,
                        unresolved=unresolved,
                        complete=factor is not None,
                        work=ledger.used,
                        seconds=time.perf_counter() - row_started,
                        misses=getattr(
                            getattr(store, "chains", None), "misses", 0
                        ),
                        evictions=getattr(
                            getattr(store, "chains", None), "evictions", 0
                        ),
                    )
                )
    return rows


def summarize(pairs, settings):
    ratios = [
        candidate["wall"] / control["wall"] for control, candidate in pairs
    ]
    generator = random.Random(settings["bootstrap_seed"])
    medians = sorted(
        statistics.median(generator.choices(ratios, k=len(ratios)))
        for _ in range(4000)
    )
    low, high = medians[100], medians[3899]

    def spread(values):
        ordered = sorted(values)
        quartiles = statistics.quantiles(ordered, n=4, method="inclusive")
        return (quartiles[2] - quartiles[0]) / statistics.median(ordered)

    variability = max(
        spread(ratios),
        *(spread([p[i]["wall"] for p in pairs]) for i in (0, 1)),
    )
    classes = {}
    for fixture in sorted({r["fixture"] for r in pairs[0][0]["rows"]}):
        outcomes = []
        for arm in (0, 1):
            rows = [
                r for r in pairs[0][arm]["rows"] if r["fixture"] == fixture
            ]
            outcomes.append(sum(r["complete"] for r in rows) / len(rows))
        times = [
            statistics.median(
                sum(
                    r["seconds"]
                    for r in pair[arm]["rows"]
                    if r["fixture"] == fixture
                )
                for pair in pairs
            )
            for arm in (0, 1)
        ]
        classes[fixture] = dict(completion=outcomes, median_seconds=times)
    regression = max(
        (
            row["completion"][0] - row["completion"][1]
            for row in classes.values()
        ),
        default=0,
    )
    cpu_ratio = statistics.median(c["cpu"] / b["cpu"] for b, c in pairs)
    midpoint = len(ratios) // 2
    halves = [
        statistics.median(ratios[:midpoint]),
        statistics.median(ratios[midpoint:]),
    ]
    return dict(
        samples=len(pairs),
        control_seconds=statistics.median(p[0]["wall"] for p in pairs),
        candidate_seconds=statistics.median(p[1]["wall"] for p in pairs),
        saving=1 - statistics.median(ratios),
        saving_ci95=[1 - high, 1 - low],
        cpu_saving=1 - cpu_ratio,
        half_savings=[1 - r for r in halves],
        relative_iqr=variability,
        stable=variability <= settings["max_relative_iqr"],
        completion_by_fixture=classes,
        max_completion_regression=regression,
        accepted=high < 1
        and cpu_ratio < 1
        and max(halves) < 1
        and variability <= settings["max_relative_iqr"]
        and regression <= 0.05,
    )


def worker(args):
    settings, corpus = inputs()
    fixtures = [f for f in corpus["fixtures"] if f["split"] == args.split]
    seeds = settings[args.split + "_seeds"]
    targets = {}
    if args.scope not in ("portfolio", "resume"):
        from math import lcm

        curves = (
            3
            if args.scope == "eviction"
            else int(args.scope.replace("fresh", "").replace("reuse", ""))
        )
        for fixture in fixtures:
            if not fixture["id"].startswith("balanced"):
                continue
            for seed in seeds:
                for curve in range(curves):
                    bound = (
                        1999
                        if args.scope == "eviction" and curve == 1
                        else 2000
                    )
                    sigma = random.Random(seed + curve).randrange(6, 2**63)
                    setup = ecm.setup_curve(fixture["n"], sigma)
                    expected = []
                    if setup.point is not None:
                        scalar = lcm(*range(1, bound + 1))
                        for prime in fixture["factors"]:
                            point, curve_a, curve_b = twist_point(
                                setup.point, setup.a24, prime
                            )
                            expected.append(
                                (
                                    prime,
                                    affine_multiply(
                                        scalar, point, prime, curve_a, curve_b
                                    ),
                                )
                            )
                    targets[(fixture["id"], seed, curve)] = expected
    functions = []
    for candidate in (False, True):
        if args.scope in ("portfolio", "resume"):
            functions.append(
                lambda candidate=candidate: portfolio_cohort(
                    fixtures,
                    seeds,
                    args.backend,
                    candidate,
                    args.scope == "resume",
                )
            )
        else:
            functions.append(
                lambda candidate=candidate: stage_cohort(
                    fixtures,
                    seeds,
                    args.backend,
                    candidate,
                    args.scope,
                    targets,
                )
            )
    warmups, signatures = [], [None, None]

    def measure(arm):
        start, cpu = time.perf_counter(), time.process_time()
        rows = arithmetic.canonical(functions[arm]())
        elapsed, cpu_used = (
            time.perf_counter() - start,
            time.process_time() - cpu,
        )
        signature = json.dumps(
            [{k: v for k, v in row.items() if k != "seconds"} for row in rows],
            sort_keys=True,
        )
        if signatures[arm] is not None and signatures[arm] != signature:
            raise AssertionError("nondeterministic factoring outcome/work")
        signatures[arm] = signature
        return dict(wall=elapsed, cpu=cpu_used, rows=rows)

    if args.cold:
        return dict(cold=True, arm=args.arm, observation=measure(args.arm))
    pairs = []
    for seconds, count in settings["sampling"]:
        for arm in (0, 1):
            start, rounds = time.perf_counter(), 0
            while time.perf_counter() - start < seconds:
                measure(arm)
                rounds += 1
            warmups.append(
                dict(
                    arm=arm,
                    seconds=time.perf_counter() - start,
                    cohorts=rounds,
                )
            )
        while len(pairs) < count:
            order = (0, 1) if len(pairs) % 2 == 0 else (1, 0)
            pair = [None, None]
            for arm in order:
                pair[arm] = measure(arm)
            pairs.append(pair)
        summary = summarize(pairs, settings)
        if summary["stable"]:
            break
    return dict(
        scope=args.scope,
        split=args.split,
        backend=args.backend,
        warmups=warmups,
        pairs=pairs,
        summary=summary,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--scope", choices=SCOPES, default="fresh1")
    parser.add_argument(
        "--split", choices=("training", "confirmation"), default="training"
    )
    parser.add_argument(
        "--backend", choices=("python-int", "gmpy2-mpz"), default="python-int"
    )
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--cold", action="store_true")
    parser.add_argument("--arm", type=int, choices=(0, 1), default=0)
    args = parser.parse_args()
    require_runtime()
    if args.worker:
        print(json.dumps(worker(args)))
        return
    if args.output is None or args.output.exists():
        parser.error("choose a new output path")
    inputs()
    owner = Path("/private/tmp/factor-performance-owner.json")
    with open("/private/tmp/factor-performance.lock", "a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert_quiet()
        owner.write_text(
            json.dumps(
                dict(
                    owner="B3 timing",
                    pid=os.getpid(),
                    scope=args.scope,
                    split=args.split,
                    worktree=str(ROOT),
                )
            )
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
            if args.cold:
                captures = []
                for index in range(9):
                    for arm in (0, 1) if index % 2 == 0 else (1, 0):
                        started = time.perf_counter()
                        result = subprocess.run(
                            command + ["--cold", "--arm", str(arm)],
                            check=True,
                            capture_output=True,
                            text=True,
                            timeout=1200,
                        )
                        captures.append(
                            dict(
                                arm=arm,
                                seconds=time.perf_counter() - started,
                                result=json.loads(result.stdout),
                            )
                        )
                args.output.parent.mkdir(parents=True, exist_ok=True)
                args.output.write_text(
                    json.dumps(dict(cold=captures), indent=2) + "\n"
                )
                print(
                    json.dumps(
                        dict(
                            cold=True,
                            scope=args.scope,
                            backend=args.backend,
                            medians=[
                                statistics.median(
                                    c["seconds"]
                                    for c in captures
                                    if c["arm"] == arm
                                )
                                for arm in (0, 1)
                            ],
                        )
                    )
                )
                return
            result = subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
                timeout=1200,
            )
            data = json.loads(result.stdout)
            inputs()
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(data, indent=2) + "\n")
            print(
                json.dumps(
                    dict(
                        scope=args.scope,
                        split=args.split,
                        backend=args.backend,
                        **data["summary"],
                    )
                )
            )
        finally:
            owner.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
