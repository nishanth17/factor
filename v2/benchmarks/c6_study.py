"""Frozen C6 stage/campaign comparisons; never changes production routing."""

import argparse
import fcntl
import hashlib
import json
import os
import platform
import signal
import statistics
import subprocess
import sys
import time
from contextlib import contextmanager
from dataclasses import asdict
from functools import lru_cache
from pathlib import Path

from .. import constants, ecm, prac, prime_sieve, utils
from . import c6_chains as chains
from . import p41_campaign as control
from . import p52_realistic
from .prac_oracle import affine_multiply, matches, twist_point

ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = Path(__file__).parent / "inputs/controls/c6_protocol.json"
ARMS = ("ladder", "checked", "compact", "lucas", "rolling")


@contextmanager
def performance_window(phase):
    """Share the machine lease with B4/A6, failing closed on contention."""
    owner = Path("/private/tmp/factor-performance-owner.json")
    with open("/private/tmp/factor-performance.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p52_realistic.check_quiet({os.getpid()})
        owner.write_text(
            json.dumps(
                dict(
                    owner="C6 Lucas chains",
                    pid=os.getpid(),
                    phase=phase,
                    worktree=str(ROOT),
                )
            )
        )
        try:
            yield
        finally:
            if (
                owner.exists()
                and json.loads(owner.read_text()).get("pid") == os.getpid()
            ):
                owner.unlink()


def freeze():
    """Pin controls and selection rules before observing C6 timings."""
    sources = [str(p.relative_to(ROOT)) for p in (ROOT / "v2").glob("*.py")]
    sources += [
        "v2/benchmarks/" + name
        for name in (
            "c6_chains.py",
            "c6_costs.py",
            "c6_study.py",
            "build_c6_inputs.py",
            "p41_campaign.py",
            "p41_gmp.py",
            "prac_oracle.py",
        )
    ]
    sources += [
        str(chains.DATA.relative_to(ROOT)),
        str(control.CORPUS.relative_to(ROOT)),
    ]
    sources += [
        str(p.relative_to(ROOT))
        for p in (ROOT / "v2/benchmarks/inputs/upstream/c6_gmp_ecm").iterdir()
        if p.is_file()
    ]
    corpus = control.load_corpus()
    data = dict(
        baseline_commit="bcf5f3d1e57304694b48ba6e7ef8b4ea2ffd0db0",
        sha256={
            name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in sources
        },
        arms=ARMS,
        backends=["int", "gmp"],
        b1=2000,
        b2=147396,
        curves=8,
        seconds_per_attempt=20,
        batch_size=128,
        seeds=corpus["seeds"][:2],
        confirmation_seeds=corpus["seeds"][2:4],
        cases=[
            f"{shape}_{digits}d"
            for digits in (40, 50, 60, 70, 80)
            for shape in ("balanced", "small10")
        ],
        sampling=[[3, 9], [5, 18], [8, 27]],
        max_relative_iqr=0.15,
        promotion="At least 10% full-stage and complete-campaign gain on each "
        "backend independently, stable timing and fresh confirmation; "
        "no proper-factor/validity regression. Otherwise retain ladder.",
        search_gate="Only if compact Lucas beats compact PRAC by >=5% in "
        "stable stage reuse and is within 10% of ladder, consider "
        "bounded offline continued-fraction search; otherwise stop.",
        limits=dict(
            records=512,
            steps=512,
            live_points=16,
            bytecode_bytes=1048576,
            bound=2000,
            generator_seconds=60,
            generator_data_bytes=536870912,
        ),
        scope="Finite inspected corpus; known-composite ECM, no portfolio "
        "or B3 checkpoint integration. Generation charged separately.",
    )
    PROTOCOL.write_text(json.dumps(data, indent=2) + "\n")


def protocol():
    data = json.loads(PROTOCOL.read_text())
    for name, digest in data["sha256"].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
            raise ValueError("frozen source changed: " + name)
    return data


def stage_one(point, n, a24, program, backend, extra):
    for prime, power, action, unit in program:
        original = point
        try:
            point = chains.apply(action, power, point, n, a24, backend)
            if backend.gcd(point[1], n) != n:
                continue
        except prac.NonunitPointError as result:
            if result.factor is not None:
                return None, result.factor
        extra["prime_power_replays"] += 1
        point, remaining = original, power
        while remaining > 1:
            extra["prime_units_replayed"] += 1
            try:
                point = chains.apply(unit, prime, point, n, a24, backend)
            except prac.NonunitPointError as result:
                return None, result.factor
            if backend.gcd(point[1], n) == n:
                return None, None
            remaining //= prime
    return point, None


def attempt(case, seed, arm, backend, settings, scope, program=None):
    """Charge conversion, setup, construction, dispatch and recovery."""

    def expired(signum, frame):
        raise control.CampaignTimeoutError

    previous = [
        signal.signal(s, expired) for s in (signal.SIGALRM, signal.SIGPROF)
    ]
    started, cpu = time.perf_counter(), time.process_time()
    stats, factor = ecm.EcmStats(), None
    extra = dict(prime_power_replays=0, prime_units_replayed=0)
    timed_out, point = False, None
    b1, b2 = settings["b1"], settings["b2"]
    for timer in (signal.ITIMER_REAL, signal.ITIMER_PROF):
        signal.setitimer(timer, settings["seconds_per_attempt"])
    try:
        n = backend.integer(case["n"])
        if arm == "ladder" and scope == "campaign":
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
            if arm != "ladder" and program is None:
                program = chains.build_program(b1, arm, backend)
            scalar = None
            if arm == "ladder":
                scalar = (
                    program
                    if program is not None
                    else backend.ecm.stage_one_scalar(b1)
                )
            generator, primes = utils.resolve_rng(seed, None), None
            curves = settings["curves"] if scope == "campaign" else 1
            for _ in range(curves):
                stats.curves += 1
                setup = backend.ecm.setup_curve(
                    n, generator.randint(6, constants.MAX_RANDOM_ECM)
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
                        scalar, *setup.point, n, setup.a24
                    )
                    divisor = backend.gcd(point[1], n)
                    if 1 < divisor < n:
                        factor = divisor
                    elif divisor == n:
                        point = None
                else:
                    point, factor = stage_one(
                        setup.point, n, setup.a24, program, backend, extra
                    )
                if factor is not None:
                    break
                if point is None:
                    stats.stage_one_saturations += 1
                    continue
                if scope != "campaign":
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
    except control.CampaignTimeoutError:
        timed_out = True
    finally:
        for timer in (signal.ITIMER_REAL, signal.ITIMER_PROF):
            signal.setitimer(timer, 0)
        for sig, handler in zip((signal.SIGALRM, signal.SIGPROF), previous):
            signal.signal(sig, handler)
    factor = int(factor) if factor is not None else None
    result = dict(
        case=case["id"],
        digits=case["digits"],
        factor_digits=case["factor_digits"],
        seed=seed,
        arm=arm,
        seconds=time.perf_counter() - started,
        cpu_seconds=time.process_time() - cpu,
        factor=factor,
        cofactor=case["n"] // factor if factor else None,
        unresolved=None if factor else case["n"],
        timed_out=timed_out,
        stats=asdict(stats),
        extra=extra,
        point=[int(v) for v in point] if point is not None else None,
    )
    control.validate_result(result, case)
    if (
        not timed_out
        and scope != "campaign"
        and point is not None
        and factor is None
    ):
        for modulus, expected in affine_controls(
            case["n"], tuple(case["factors"]), seed, b1
        ):
            if not matches(point, expected, modulus):
                raise AssertionError("invalid stage-one projective action")
    return result


@lru_cache(maxsize=128)
def affine_controls(n, factors, seed, bound):
    """Independent full-coordinate prime-field controls outside the timer."""
    generator = utils.resolve_rng(seed, None)
    setup = ecm.setup_curve(n, generator.randint(6, constants.MAX_RANDOM_ECM))
    scalar = ecm.stage_one_scalar(bound)
    controls = []
    for modulus in factors:
        point, curve_a, curve_b = twist_point(setup.point, setup.a24, modulus)
        expected = affine_multiply(scalar, point, modulus, curve_a, curve_b)
        controls.append((modulus, expected))
    return tuple(controls)


def spread(values):
    quartiles = statistics.quantiles(values, n=4)
    return (quartiles[2] - quartiles[0]) / statistics.median(values)


def worker(args):
    settings = protocol()
    corpus = control.load_corpus()
    cases = [c for c in corpus["fixtures"] if c["id"] in settings["cases"]]
    if args.scope != "campaign":
        cases = [c for c in cases if c["id"].startswith("balanced")]
    backend = (
        control.gmp_backend()
        if args.backend == "gmp"
        else control.PYTHON_BACKEND
    )
    preparation = time.perf_counter()
    program = None
    if args.scope == "stage_reuse":
        program = (
            backend.ecm.stage_one_scalar(2000)
            if args.arm == "ladder"
            else chains.build_program(2000, args.arm, backend)
        )
    preparation = time.perf_counter() - preparation

    def sample():
        rows = [
            attempt(c, seed, args.arm, backend, settings, args.scope, program)
            for seed in settings["seeds"]
            for c in cases
        ]
        return dict(seconds=sum(r["seconds"] for r in rows), rows=rows)

    # Every warmup output is validated. Seeds are frozen and all outcomes kept.
    start, runs = time.perf_counter(), 0
    while time.perf_counter() - start < args.warmup:
        sample()
        runs += 1
    warmed = time.perf_counter() - start
    samples = [sample() for _ in range(args.samples)]
    return dict(
        backend=args.backend,
        arm=args.arm,
        scope=args.scope,
        warmup_seconds=warmed,
        validated_warmups=runs,
        preparation_seconds=preparation,
        samples=samples,
        censored=any(r["timed_out"] for s in samples for r in s["rows"]),
        median=statistics.median(s["seconds"] for s in samples),
        relative_iqr=spread([s["seconds"] for s in samples])
        if len(samples) > 1
        else None,
    )


def launch(args, arm, backend, scope, warmup, samples):
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
        scope,
        "--warmup",
        str(warmup),
        "--samples",
        str(samples),
    ]
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
                if time.perf_counter() - started > 900:
                    raise RuntimeError("worker exceeded 900-second cap")
        p52_realistic.check_quiet({os.getpid(), process.pid})
        if process.returncode:
            raise RuntimeError(stderr)
    except BaseException:
        process.terminate()
        process.communicate()
        raise
    return json.loads(stdout), time.perf_counter() - started


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", action="store_true")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument(
        "--scope",
        choices=("stage_fresh", "stage_reuse", "campaign"),
        default="campaign",
    )
    parser.add_argument("--arm", choices=ARMS, default="ladder")
    parser.add_argument("--backend", choices=("int", "gmp"), default="int")
    parser.add_argument("--warmup", type=float, default=3)
    parser.add_argument("--samples", type=int, default=9)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cold", action="store_true")
    args = parser.parse_args()
    from .build_c6_inputs import require_runtime

    require_runtime()
    if args.freeze:
        freeze()
        return
    if args.worker:
        print(json.dumps(worker(args)))
        return
    settings = protocol()
    report = dict(
        protocol=settings,
        runtime=sys.version,
        executable=sys.executable,
        platform=platform.platform(),
        captures=[],
        attempts=[],
        cold=[],
    )
    with performance_window(args.scope), args.output.open("x") as output:

        def save():
            output.seek(0)
            json.dump(report, output, indent=2)
            output.truncate()
            output.flush()

        save()
        # Reverse candidate ordering by backend to reduce monotonic drift.
        for backend in settings["backends"]:
            for arm in ARMS[:: 1 if backend == "int" else -1]:
                for warmup, samples in settings["sampling"]:
                    capture, _ = launch(
                        args, arm, backend, args.scope, warmup, samples
                    )
                    report["attempts"].append(capture)
                    save()
                    if capture["relative_iqr"] <= settings["max_relative_iqr"]:
                        break
                report["captures"].append(capture)
                save()
                print(
                    backend,
                    args.scope,
                    arm,
                    round(capture["median"], 6),
                    capture["relative_iqr"],
                    len(capture["samples"]),
                    flush=True,
                )
                if args.cold:
                    for _ in range(9):
                        capture, elapsed = launch(
                            args, arm, backend, args.scope, 0, 1
                        )
                        report["cold"].append(
                            dict(capture=capture, seconds=elapsed)
                        )
                        save()
        protocol()


if __name__ == "__main__":
    main()
