"""Frozen B4 comparisons under a machine-wide exclusive performance window."""

import argparse
import fcntl
import json
import os
import platform
import random
import resource
import statistics
import subprocess
import sys
import time
from math import gcd
from pathlib import Path

from ...support.paths import (
    REPOSITORY_ROOT,
    matches_source_pins,
    source_path,
)
from ..p52 import p52_a3, p52_realistic
from . import b4_common as common
from . import b4_kernels as kernels

FREEZE = common.INPUTS / "controls/b4_freeze.json"
ROOT = REPOSITORY_ROOT


def source_hashes():
    paths = [
        Path(__file__),
        Path(kernels.__file__),
        Path(common.__file__),
        Path(p52_a3.__file__),
        Path(p52_realistic.__file__),
        ROOT / "v2/arithmetic.py",
        ROOT / "v2/utils.py",
        ROOT / "v2/benchmarks/build_phase_two_corpus.py",
        ROOT / "v2/benchmarks/prac_oracle.py",
        common.CONTROL,
        common.CORPUS,
        common.PROTOCOL,
    ]
    return {str(path.relative_to(ROOT)): common.digest(path) for path in paths}


def freeze(profile):
    report = json.loads(source_path(profile).read_text())
    if not report["rows"] or report["protocol_sha256"] != common.digest(
        common.PROTOCOL
    ):
        raise ValueError("baseline profile must precede candidate freeze")
    data = dict(
        source_hashes=source_hashes(),
        arms=list(kernels.ARMS),
        baseline_profile_sha256=common.digest(profile),
        selection_reason="Point arithmetic is a measured component; test "
        "five bounded syntax/fusion/reduction/normalization transfers. "
        "Leave reducers, backends, schedules and curve families unchanged.",
    )
    with FREEZE.open("x") as output:
        json.dump(data, output, indent=2)
        output.write("\n")


def verify_freeze():
    data = json.loads(source_path(FREEZE).read_text())
    if not matches_source_pins(data["source_hashes"]) or data["arms"] != list(
        kernels.ARMS
    ):
        raise ValueError("sources changed after B4 freeze")
    return data


def check_point(actual, expected, n):
    """Equality requires primitive pairs; degeneracy is a separate outcome."""
    x, z = map(int, actual)
    u, v = map(int, expected)
    ideals = (gcd(gcd(x, z), n), gcd(gcd(u, v), n))
    if ideals != (1, 1):
        if ideals[0] != ideals[1] or gcd(z, n) != gcd(v, n):
            raise AssertionError("exceptional ideals differ")
        return "degenerate"
    if (x * v - u * z) % n or gcd(z, n) != gcd(v, n):
        raise AssertionError("projective actions differ")
    return "primitive"


def kernel_oracles(protocol):
    ecm = kernels.ecm_module(kernels.engine("baseline"))
    cases = []
    for bits in protocol["kernel_bits"]:
        n = 2**bits - 59
        for sigma in protocol["kernel_sigmas"]:
            setup = ecm.setup_curve(n, sigma)
            if setup.point is None:
                cases.append((n, sigma, None))
            else:
                expected = ecm.scalar_multiply(
                    protocol["kernel_scalar"], *setup.point, n, setup.a24
                )
                cases.append((n, sigma, expected))
    return cases


def measure_kernels(engine, protocol, backend, oracles):
    ecm = kernels.ecm_module(engine)
    integer = (
        __import__(engine.__package__ + ".arithmetic", fromlist=["arithmetic"])
        .get_backend(backend)
        .integer
    )
    rows = []
    for n, sigma, expected in oracles:
        # Entry conversion and Suyama setup are paid each time, including its
        # checked inverse. A normalization inverse is additional, never reused.
        modulus = integer(n)
        setup = ecm.setup_curve(modulus, sigma)
        if setup.point is None:
            rows.append([n.bit_length(), sigma, "setup_exception"])
            continue
        point = ecm.scalar_multiply(
            protocol["kernel_scalar"], *setup.point, modulus, setup.a24
        )
        status = check_point(point, expected, n)
        rows.append([n.bit_length(), sigma, status])
    return rows


def measure_stages(engine, protocol, backend, fixtures):
    ecm = kernels.ecm_module(engine)
    package = __import__(
        engine.__package__, fromlist=["arithmetic", "prime_sieve", "utils"]
    )
    integer = package.arithmetic.get_backend(backend).integer
    b1, b2 = protocol["stage_bounds"]
    rows = []
    for fixture in fixtures:
        for seed in protocol["stage_seeds"]:
            start = time.perf_counter()
            n = integer(fixture["n"])
            sigma = random.Random(seed).randint(
                6, package.constants.MAX_RANDOM_ECM
            )
            setup = ecm.setup_curve(n, sigma)
            if setup.point is None:
                rows.append(
                    dict(
                        fixture=fixture["id"],
                        seed=seed,
                        stage_one=time.perf_counter() - start,
                        stage_two=0,
                        factor=int(setup.factor or 0),
                        status="setup_exception",
                    )
                )
                continue
            scalar = ecm.stage_one_scalar(b1, backend=backend)
            point = ecm.scalar_multiply(scalar, *setup.point, n, setup.a24)
            divisor = int(package.arithmetic.gcd(point[1], n))
            # Prime factors from the corpus validate divisors but never guide
            # execution. The readable formula validates scalar action outside
            # the timing window through the precomputed stage oracle below.
            first = time.perf_counter() - start
            start = time.perf_counter()
            result, saturated = None, divisor == n
            if divisor == 1:
                primes = package.prime_sieve.segmented_sieve(b1 + 1, b2 + 1)
                result, saturated = ecm.stage_two(
                    point, n, setup.a24, b1, primes
                )
            elif divisor < n:
                result = divisor
            second = time.perf_counter() - start
            if result is not None and not package.utils.valid_divisor(
                result, n
            ):
                raise AssertionError("invalid stage divisor")
            rows.append(
                dict(
                    fixture=fixture["id"],
                    seed=seed,
                    stage_one=first,
                    stage_two=second,
                    factor=int(result or 0),
                    status="saturated" if saturated else "valid",
                    point=[int(c) for c in point],
                )
            )
    return rows


def measure_full(engine, protocol, backend, fixtures):
    rows = []
    for fixture in fixtures:
        for seed in protocol["seeds"]:
            started = time.perf_counter()
            config = common.config(engine, protocol, fixture["case"], backend)
            run = engine.factorize_bounded(
                fixture["n"],
                seed=seed,
                config=config,
                budget=engine.Budget(
                    work_limit=protocol["work_limit"],
                    seconds=protocol["seconds"],
                    cpu_seconds=protocol["cpu_seconds"],
                ),
            )
            row = common.validate_run(run, fixture)
            # _pack has serialized and checked canonical checkpoints already;
            # include a public JSON roundtrip and complete output validation.
            json.loads(json.dumps(run.checkpoint))
            row.update(
                seed=seed,
                case=fixture["case"],
                seconds=time.perf_counter() - started,
            )
            if run.reason in ("wall_limit", "cpu_limit"):
                raise AssertionError("censored sample cannot be accepted")
            rows.append(row)
    return rows


def signature(rows):
    return [
        {
            k: v
            for k, v in row.items()
            if k not in ("seconds", "stage_one", "stage_two", "point")
        }
        if isinstance(row, dict)
        else row
        for row in rows
    ]


def worker(args):
    common.require_runtime()
    verify_freeze()
    protocol, corpus = common.inputs()
    engine = kernels.engine(args.arm)
    stage_fixtures = [f for f in corpus["fixtures"] if f["split"] == "stages"]
    fixtures = [f for f in corpus["fixtures"] if f["split"] == args.phase]
    oracles = kernel_oracles(protocol) if args.scope == "kernel" else None
    expected_stages = None
    if args.scope == "stages":
        expected_stages = measure_stages(
            kernels.engine("baseline"), protocol, "python-int", stage_fixtures
        )

    def measure():
        started, cpu = time.perf_counter(), time.process_time()
        if args.scope == "full":
            rows = measure_full(engine, protocol, args.backend, fixtures)
        elif args.scope == "kernel":
            rows = measure_kernels(engine, protocol, args.backend, oracles)
        else:
            rows = measure_stages(
                engine, protocol, args.backend, stage_fixtures
            )
            for row, expected in zip(rows, expected_stages):
                if (
                    row["factor"] != expected["factor"]
                    or row["status"] != expected["status"]
                ):
                    raise AssertionError("stage outcomes differ")
                if "point" in row:
                    fixture = next(
                        f for f in stage_fixtures if f["id"] == row["fixture"]
                    )
                    check_point(row["point"], expected["point"], fixture["n"])
        return dict(
            seconds=time.perf_counter() - started,
            cpu_seconds=time.process_time() - cpu,
            rows=rows,
        )

    start, warm_count, expected = time.perf_counter(), 0, None
    while time.perf_counter() - start < args.warmup:
        sample = measure()
        current = signature(sample["rows"])
        if expected is not None and current != expected:
            raise AssertionError("warmup outcomes changed")
        expected = current
        warm_count += 1
    warm_elapsed = time.perf_counter() - start
    samples = []
    for _ in range(args.samples):
        sample = measure()
        current = signature(sample["rows"])
        if expected is not None and current != expected:
            raise AssertionError("sample outcomes changed")
        expected = current
        samples.append(sample)
    return dict(
        arm=args.arm,
        backend=args.backend,
        phase=args.phase,
        scope=args.scope,
        warmup_seconds=warm_elapsed,
        warmup_cohorts=warm_count,
        samples=samples,
        rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    )


def launch(args, arm, backend, scope, warmup, samples):
    p52_realistic.check_quiet({os.getpid()})
    command = [
        sys.executable,
        "-B",
        "-m",
        "v2.benchmarks.ecm.b4.b4_study",
        "--worker",
        "--phase",
        args.phase,
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
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    try:
        while True:
            try:
                stdout, stderr = process.communicate(timeout=1)
                break
            except subprocess.TimeoutExpired:
                p52_realistic.check_quiet({os.getpid(), process.pid})
                if time.perf_counter() - started > 600:
                    raise TimeoutError("bounded worker exceeded 600 seconds")
        if process.returncode:
            raise RuntimeError(stderr)
        return json.loads(stdout), time.perf_counter() - started
    finally:
        if process.poll() is None:
            process.terminate()
        process.communicate()


def stable(capture, protocol):
    return (
        p52_realistic.spread([s["seconds"] for s in capture["samples"]])
        <= protocol["relative_iqr_limit"]
    )


def completions(capture):
    rows = capture["samples"][0]["rows"]
    return {
        case: statistics.mean(
            row["complete"] for row in rows if row["case"] == case
        )
        for case in {row["case"] for row in rows}
    }


def select(captures, protocol):
    selected = {}
    for backend in protocol["backends"]:
        candidates = [
            c
            for c in captures
            if c["backend"] == backend and c["scope"] == "full"
        ]
        control = next(c for c in candidates if c["arm"] == "baseline")
        baseline = completions(control)
        eligible = [
            c
            for c in candidates
            if c["arm"] != "baseline"
            and stable(c, protocol)
            and all(
                completions(c)[case] >= value - 0.05
                for case, value in baseline.items()
            )
        ]
        selected[backend] = (
            min(
                eligible,
                key=lambda c: statistics.median(
                    s["seconds"] for s in c["samples"]
                ),
            )["arm"]
            if eligible and stable(control, protocol)
            else None
        )
    return selected


def summaries(captures):
    rows = []
    for candidate in captures:
        if candidate["arm"] == "baseline":
            continue
        control = next(
            c
            for c in captures
            if c["arm"] == "baseline"
            and c["backend"] == candidate["backend"]
            and c["scope"] == candidate["scope"]
        )
        before = [s["seconds"] for s in control["samples"]]
        after = [s["seconds"] for s in candidate["samples"]]
        rows.append(
            dict(
                backend=candidate["backend"],
                scope=candidate["scope"],
                arm=candidate["arm"],
                reduction_percent=100
                * (1 - statistics.median(after) / statistics.median(before)),
                conditional_timing_interval=p52_a3.interval(before, after),
            )
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-profile", type=Path)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument(
        "--phase", choices=("training", "held_out"), default="training"
    )
    parser.add_argument(
        "--scope", choices=("kernel", "stages", "full"), default="full"
    )
    parser.add_argument(
        "--scopes",
        nargs="+",
        choices=("kernel", "stages", "full"),
        default=["kernel", "stages", "full"],
    )
    parser.add_argument("--arm", choices=kernels.ARMS, default="baseline")
    parser.add_argument(
        "--backend", choices=("python-int", "gmpy2-mpz"), default="python-int"
    )
    parser.add_argument("--warmup", type=float, default=3)
    parser.add_argument("--samples", type=int, default=9)
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--cold", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    common.require_runtime()
    if args.freeze_profile:
        freeze(args.freeze_profile)
        return
    if args.worker:
        print(json.dumps(worker(args)))
        return
    pins = verify_freeze()
    protocol, _ = common.inputs()
    selection = None
    if args.phase == "held_out":
        data = json.loads(source_path(args.selection).read_text())
        if data["freeze"] != pins or data["phase"] != "training":
            raise ValueError(
                "selection must come from matching frozen training"
            )
        selection = data["selected"]
    report = dict(
        phase=args.phase,
        freeze=pins,
        runtime=sys.version,
        executable=sys.executable,
        platform=platform.platform(),
        attempts=[],
        captures=[],
        cold=[],
    )
    p52_realistic.check_quiet({os.getpid()})
    with open(protocol["machine_lock"], "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner = Path("/private/tmp/factor-performance-owner.json")
        owner.write_text(
            json.dumps(
                dict(
                    owner="B4 arithmetic kernels",
                    pid=os.getpid(),
                    phase=args.phase,
                )
            )
        )
        try:
            with args.output.open("x") as output:

                def save():
                    output.seek(0)
                    json.dump(report, output, indent=2)
                    output.truncate()
                    output.flush()

                save()
                tasks = [
                    (backend, scope, arm)
                    for backend in protocol["backends"]
                    for scope in args.scopes
                    for arm in (
                        kernels.ARMS
                        if selection is None
                        else ("baseline", selection[backend])
                    )
                    if arm is not None
                ]
                random.Random(2026100942).shuffle(tasks)
                for backend, scope, arm in tasks:
                    for warmup, samples in protocol["sampling"]:
                        capture, _ = launch(
                            args, arm, backend, scope, warmup, samples
                        )
                        report["attempts"].append(capture)
                        save()
                        if stable(capture, protocol):
                            break
                    report["captures"].append(capture)
                    save()
                    print(
                        backend,
                        scope,
                        arm,
                        len(capture["samples"]),
                        "stable"
                        if stable(capture, protocol)
                        else "inconclusive",
                        flush=True,
                    )
                    if args.cold and scope == "full":
                        for _ in range(9):
                            cold, elapsed = launch(
                                args, arm, backend, scope, 0, 1
                            )
                            report["cold"].append(
                                dict(seconds=elapsed, capture=cold)
                            )
                            save()
                verify_freeze()
                report.update(summary=summaries(report["captures"]))
                if selection is None:
                    report["selected"] = select(report["captures"], protocol)
                save()
        finally:
            owner.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
