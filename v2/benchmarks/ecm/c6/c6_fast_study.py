"""Paired C6 optimized-executor screen and untouched confirmation."""

import argparse
import cProfile
import hashlib
import json
import os
import platform
import random
import resource
import signal
import statistics
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

from .... import constants
from ....common import prime_sieve, utils
from ....ecm import core as ecm
from ...suites.build_phase_two_corpus import verify_certificates
from ...support.paths import (
    BENCHMARK_ROOT,
    source_path,
)
from ...support.prac_oracle import matches
from ..p41 import p41_campaign
from ..p52 import p52_realistic
from . import c6_chains, c6_fast, c6_study
from .build_c6_fast_inputs import HELDOUT, PROTOCOL
from .build_c6_inputs import require_runtime

SELECTION = BENCHMARK_ROOT / "inputs/controls/c6_fast_selection.json"


def protocol():
    data = json.loads(source_path(PROTOCOL).read_text())
    for name, digest in data["sha256"].items():
        actual = hashlib.sha256(
            (source_path(c6_study.ROOT / name)).read_bytes()
        ).hexdigest()
        if actual != digest:
            raise ValueError("frozen source changed: " + name)
    return data


def corpus(confirmation=False):
    if not confirmation:
        return p41_campaign.load_corpus()
    data = json.loads(source_path(HELDOUT).read_text())
    verify_certificates(data["certificates"])
    for case in data["fixtures"]:
        p, q = case["factors"]
        assert p != q and p * q == case["n"]
        assert len(str(case["n"])) == case["digits"]
        assert case["factor_digits"] == [len(str(p)), len(str(q))]
    return data


def construct(arm, backend):
    if arm == "ladder":
        return backend.ecm.stage_one_scalar(2000)
    if arm.startswith("strict/"):
        family = arm.split("/")[1]
        return c6_chains.build_program(
            2000, "compact" if family == "prac" else "lucas", backend
        )
    _, family, mode, batch = arm.split("/")
    return c6_fast.build_program(2000, family, mode, backend, int(batch))


def attempt(case, seed, arm, backend, settings, scope, program, oracle):
    """Charge construction, execution, recovery and output validation."""

    def expired(signum, frame):
        raise p41_campaign.CampaignTimeoutError

    handlers = [
        signal.signal(sig, expired) for sig in (signal.SIGALRM, signal.SIGPROF)
    ]
    started, cpu = time.perf_counter(), time.process_time()
    stats, extra = ecm.EcmStats(), {}
    factor, point, timed_out = None, None, False
    b1, b2 = settings["b1"], settings["b2"]
    for timer in (signal.ITIMER_REAL, signal.ITIMER_PROF):
        signal.setitimer(timer, settings["seconds_per_attempt"])
    try:
        n = backend.integer(case["n"])
        campaign = scope.startswith("campaign")
        if campaign and arm == "ladder":
            factor = backend.ecm.factorize_ecm(
                n,
                b1=b1,
                b2=b2,
                seed=seed,
                max_curves=settings["curves"],
                stats=stats,
                _known_composite=True,
            )
        else:
            if program is None:
                program = construct(arm, backend)
            rng = utils.resolve_rng(seed, None)
            primes = None
            for _ in range(settings["curves"] if campaign else 1):
                stats.curves += 1
                setup = backend.ecm.setup_curve(
                    n, rng.randint(6, constants.MAX_RANDOM_ECM)
                )
                if setup.factor is not None:
                    factor = setup.factor
                    break
                if setup.retry:
                    stats.setup_retries += 1
                    continue
                stats.stage_one_calls += 1
                if arm == "ladder":
                    point = backend.ecm.scalar_multiply(
                        program, *setup.point, n, setup.a24
                    )
                    divisor = backend.gcd(point[1], n)
                    if 1 < divisor < n:
                        factor = divisor
                    elif divisor == n:
                        point = None
                elif arm.startswith("strict/"):
                    extra.setdefault("prime_power_replays", 0)
                    extra.setdefault("prime_units_replayed", 0)
                    point, factor = c6_study.stage_one(
                        setup.point, n, setup.a24, program, backend, extra
                    )
                else:
                    point, factor = program(setup.point, n, setup.a24, extra)
                if factor is not None:
                    break
                if point is None:
                    stats.stage_one_saturations += 1
                    continue
                if not campaign:
                    break
                if primes is None:
                    primes = prime_sieve.segmented_sieve(b1 + 1, b2 + 1)
                stats.stage_two_calls += 1
                factor, saturated = backend.ecm.stage_two(
                    point, n, setup.a24, b1, primes
                )
                if saturated:
                    stats.stage_two_saturations += 1
                if factor is not None:
                    break
    except p41_campaign.CampaignTimeoutError:
        timed_out = True
    finally:
        for timer in (signal.ITIMER_REAL, signal.ITIMER_PROF):
            signal.setitimer(timer, 0)
        for sig, handler in zip((signal.SIGALRM, signal.SIGPROF), handlers):
            signal.signal(sig, handler)
    factor = int(factor) if factor is not None else None
    result = dict(
        case=case["id"],
        seed=seed,
        arm=arm,
        digits=case["digits"],
        factor_digits=case["factor_digits"],
        factor=factor,
        cofactor=case["n"] // factor if factor else None,
        unresolved=None if factor else case["n"],
        timed_out=timed_out,
        stats=asdict(stats),
        extra=extra,
        point=[int(v) for v in point] if point is not None else None,
    )
    p41_campaign.validate_result(result, case)
    if not timed_out and point is not None and factor is None and oracle:
        for modulus, expected in oracle:
            if not matches(point, expected, modulus):
                raise AssertionError("invalid independent stage action")
    result["seconds"] = time.perf_counter() - started
    result["cpu_seconds"] = time.process_time() - cpu
    return result


def interval(ratios, settings):
    generator = random.Random(settings["bootstrap_seed"])
    estimates = sorted(
        statistics.median(generator.choices(ratios, k=len(ratios)))
        for _ in range(settings["bootstrap_replicates"])
    )
    return estimates[len(estimates) // 40], estimates[
        -len(estimates) // 40 - 1
    ]


def summarize(capture, settings):
    ratios = [
        sample["candidate"]["seconds"] / sample["ladder"]["seconds"]
        for sample in capture["samples"]
    ]
    spreads = [
        c6_study.spread([s[arm]["seconds"] for s in capture["samples"]])
        for arm in ("candidate", "ladder")
    ]
    return dict(
        ratio=statistics.median(ratios),
        interval=interval(ratios, settings),
        relative_iqr=max(spreads),
        stable=max(spreads) <= settings["max_relative_iqr"],
        candidate_seconds=statistics.median(
            s["candidate"]["seconds"] for s in capture["samples"]
        ),
        ladder_seconds=statistics.median(
            s["ladder"]["seconds"] for s in capture["samples"]
        ),
    )


def worker(args):
    settings = protocol()
    data = corpus(args.confirmation)
    cases = [c for c in data["fixtures"] if c["id"] in settings["cases"]]
    if args.scope.startswith("stage"):
        cases = [c for c in cases if c["id"].startswith("balanced")]
    seeds = settings["confirmation_seeds" if args.confirmation else "seeds"]
    backend = (
        p41_campaign.gmp_backend()
        if args.backend == "gmp"
        else p41_campaign.PYTHON_BACKEND
    )
    controls = {}
    if args.scope.startswith("stage"):
        for case in cases:
            for seed in seeds:
                controls[case["id"], seed] = c6_study.affine_controls(
                    case["n"], tuple(case["factors"]), seed, settings["b1"]
                )
    prepared = time.perf_counter()
    arms = (args.arm,) if args.single else ("ladder", args.arm)
    programs = {
        arm: construct(arm, backend) if args.scope.endswith("reuse") else None
        for arm in arms
    }
    preparation = time.perf_counter() - prepared

    def sample(arm):
        rows = [
            attempt(
                case,
                seed,
                arm,
                backend,
                settings,
                args.scope,
                programs[arm],
                controls.get((case["id"], seed)),
            )
            for seed in seeds
            for case in cases
        ]
        return dict(seconds=sum(row["seconds"] for row in rows), rows=rows)

    elapsed = dict.fromkeys(arms, 0.0)
    warmups = dict.fromkeys(arms, 0)
    while min(elapsed.values()) < args.warmup:
        for arm in arms:
            result = sample(arm)
            elapsed[arm] += result["seconds"]
            warmups[arm] += 1
    if args.profile is not None:
        profiler = cProfile.Profile()
        profiler.runcall(sample, args.arm)
        profiler.dump_stats(str(args.profile))
        return dict(profile=str(args.profile), warmup_seconds=elapsed)
    samples = []
    for index in range(args.samples):
        result = {}
        for arm in arms[:: 1 if index % 2 == 0 else -1]:
            result[
                "ladder"
                if arm == "ladder" and not args.single
                else "candidate"
            ] = sample(arm)
        samples.append(result)
    result = dict(
        arm=args.arm,
        backend=args.backend,
        scope=args.scope,
        confirmation=args.confirmation,
        samples=samples,
        warmup_seconds=elapsed,
        validated_warmups=warmups,
        preparation_seconds=preparation,
        maxrss_bytes_macos=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    )
    if not args.single and args.samples >= 9:
        result["summary"] = summarize(result, settings)
    return result


def launch(args, arm, backend, warmup, count, *, single=False, profile=None):
    command = [
        sys.executable,
        "-B",
        "-m",
        __spec__.name,
        "--worker",
        "--arm",
        arm,
        "--backend",
        backend,
        "--scope",
        args.scope,
        "--warmup",
        str(warmup),
        "--samples",
        str(count),
    ]
    if args.confirmation:
        command.append("--confirmation")
    if single:
        command.append("--single")
    if profile is not None:
        command += ["--profile", str(profile)]
    started = time.perf_counter()
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    try:
        while True:
            try:
                stdout, stderr = process.communicate(timeout=1)
                break
            except subprocess.TimeoutExpired:
                p52_realistic.check_quiet({os.getpid(), process.pid})
                if time.perf_counter() - started > 1200:
                    raise RuntimeError("C6 worker exceeded finite limit")
        p52_realistic.check_quiet({os.getpid(), process.pid})
        if process.returncode:
            raise RuntimeError(stderr)
    except BaseException:
        process.terminate()
        process.communicate()
        raise
    return json.loads(stdout), time.perf_counter() - started


def select(path):
    data = json.loads(source_path(path).read_text())
    settings = protocol()
    result = dict(
        screen_sha256=hashlib.sha256(
            source_path(path).read_bytes()
        ).hexdigest(),
        protocol_sha256=hashlib.sha256(
            source_path(PROTOCOL).read_bytes()
        ).hexdigest(),
        candidates={},
    )
    for backend in settings["backends"]:
        choices = []
        for family in settings["families"]:
            rows = [
                r
                for r in data["captures"]
                if r["backend"] == backend
                and r["arm"].startswith("fast/" + family + "/")
                and r["summary"]["stable"]
            ]
            if not rows:
                raise RuntimeError("no stable candidate for " + family)
            best = min(r["summary"]["ratio"] for r in rows)
            tied = [r for r in rows if r["summary"]["ratio"] <= best * 1.01]
            winner = min(
                tied,
                key=lambda r: (
                    c6_fast.MODES.index(r["arm"].split("/")[2]),
                    int(r["arm"].split("/")[3]),
                ),
            )
            choices.append(winner["arm"])
        result["candidates"][backend] = choices
    with SELECTION.open("x") as output:
        json.dump(result, output, indent=2)
        output.write("\n")
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--arm", default="fast/lucas/tuple/16")
    parser.add_argument("--backend", choices=("int", "gmp"), default="int")
    parser.add_argument(
        "--scope",
        choices=(
            "stage_reuse",
            "stage_fresh",
            "campaign_reuse",
            "campaign_fresh",
        ),
        default="stage_reuse",
    )
    parser.add_argument("--warmup", type=float, default=3)
    parser.add_argument("--samples", type=int, default=9)
    parser.add_argument("--single", action="store_true")
    parser.add_argument("--confirmation", action="store_true")
    parser.add_argument("--screen", action="store_true")
    parser.add_argument("--cold", action="store_true")
    parser.add_argument("--profiles", action="store_true")
    parser.add_argument("--profile", type=Path)
    parser.add_argument("--select", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require_runtime()
    if args.select:
        select(args.select)
        return
    if args.worker:
        print(json.dumps(worker(args)))
        return
    settings = protocol()
    if args.screen and args.confirmation:
        raise ValueError("held-out data cannot be used for screening")
    if args.screen:
        candidates = [
            f"fast/{family}/{mode}/{batch}"
            for family in settings["families"]
            for mode in settings["modes"]
            for batch in settings["batches"]
        ] + ["strict/prac", "strict/lucas"]
        assignments = dict.fromkeys(settings["backends"], candidates)
        selection_hash = None
    else:
        selected = json.loads(source_path(SELECTION).read_text())
        if (
            selected["protocol_sha256"]
            != hashlib.sha256(source_path(PROTOCOL).read_bytes()).hexdigest()
        ):
            raise ValueError("selection belongs to a different protocol")
        assignments = selected["candidates"]
        selection_hash = hashlib.sha256(
            source_path(SELECTION).read_bytes()
        ).hexdigest()
    report = dict(
        protocol=settings,
        selection_sha256=selection_hash,
        runtime=sys.version,
        platform=platform.platform(),
        captures=[],
        attempts=[],
        cold=[],
        profiles=[],
    )
    with (
        c6_study.performance_window("optimized C6 " + args.scope),
        args.output.open("x") as output,
    ):
        started = time.perf_counter()

        def save():
            output.seek(0)
            json.dump(report, output, indent=2)
            output.truncate()
            output.flush()

        save()
        for backend, arms in assignments.items():
            ordered = list(arms)
            random.Random(6104604 if backend == "int" else 6104605).shuffle(
                ordered
            )
            if args.cold:
                ordered = ["ladder"] + ordered
            for arm in ordered:
                if (
                    time.perf_counter() - started
                    > settings["limits"]["phase_seconds"]
                ):
                    raise RuntimeError(
                        "phase cap reached; retain partial evidence"
                    )
                if args.profiles:
                    profile = args.output.parent / (
                        backend + "-" + arm.replace("/", "-") + ".pstats"
                    )
                    capture, _ = launch(
                        args, arm, backend, 3, 0, profile=profile
                    )
                    report["profiles"].append(capture)
                elif args.cold:
                    for _ in range(9):
                        capture, elapsed = launch(
                            args, arm, backend, 0, 1, single=True
                        )
                        report["cold"].append(
                            dict(capture=capture, seconds=elapsed)
                        )
                        save()
                else:
                    for warmup, count in settings["sampling"]:
                        capture, _ = launch(args, arm, backend, warmup, count)
                        report["attempts"].append(capture)
                        save()
                        if capture["summary"]["stable"]:
                            break
                    report["captures"].append(capture)
                    print(
                        backend,
                        arm,
                        json.dumps(capture["summary"]),
                        count,
                        flush=True,
                    )
                save()
        protocol()


if __name__ == "__main__":
    main()
