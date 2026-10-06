"""Reproduce A4 correctness and diagnostic PRAC costs, without promotion.

Run with PyPy 3.11. Stage measurements include setup, schedule construction,
cache misses, dispatch and checked recovery. These are standalone stage-one
campaigns, not the budgeted portfolio: B3 owns that integration and its gate.
"""

import argparse
import json
import platform
import statistics
import subprocess
import sys
import time
from math import gcd, prod
from pathlib import Path

from v2 import ecm, prac, prime_sieve, utils

from .phase_one import environment
from .prac_oracle import (
    affine_add,
    affine_multiply,
    historical_points,
    matches,
    twist_point,
)

SIGMAS = (6, 11, 23)
MODULI = (
    (2**31 - 1, 2**61 - 1),
    (2**61 - 1, 2**127 - 1),
    (2**61 - 1, 2**255 - 19),
)


def powers_for(bound):
    return tuple(
        utils.prime_power(p, bound) for p in prime_sieve.prime_sieve(bound + 1)
    )


def validate_historical():
    fixtures = list(historical_points())
    expected = [None] * len(fixtures)
    checks = exceptional = 0
    examples = []
    for scalar in range(1001):
        chain = prac.get_chain(scalar)
        for index, (modulus, curve_a, point) in enumerate(fixtures):
            if scalar:
                expected[index] = affine_add(
                    expected[index], point, modulus, curve_a
                )
            a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
            try:
                prac._execute(
                    chain,
                    (point[0], 1),
                    modulus,
                    a24,
                    ecm.point_add,
                    ecm.point_double,
                )
            except prac._ExceptionalDifferenceError:
                exceptional += 1
                if len(examples) < 16:
                    examples.append((scalar, modulus, curve_a, point))
            actual = ecm.multiply_prac(scalar, point[0], 1, modulus, a24)
            if not matches(actual, expected[index], modulus):
                raise AssertionError((scalar, modulus, point, actual))
            checks += 1
    return {
        "affine_comparisons": checks,
        "mismatches": 0,
        "exceptional_chains_recovered": exceptional,
        "degenerate_points_counted_as_equal": 0,
        "recovery_examples": examples,
    }


def stage_sample(mode, bound, primes):
    """Three fixed curves with the exact lcm action; include every setup."""
    n = prod(primes)
    powers = powers_for(bound)
    scalars = powers
    if mode == "ladder_chunks":
        chunks = []
        for start in range(0, len(powers), 16):
            end = start + 16
            chunks.append(prod(powers[start:end]))
        scalars = tuple(chunks)
    multiply = (
        ecm.multiply_prac if mode == "prac_powers" else ecm.scalar_multiply
    )
    outputs = []
    for sigma in SIGMAS:
        setup = ecm.setup_curve(n, sigma)
        if setup.factor or setup.retry:
            outputs.append((sigma, None, setup.factor, setup.retry))
            continue
        point, factor, retry = setup.point, None, False
        for scalar in scalars:
            try:
                point = multiply(scalar, *point, n, setup.a24)
            except prac.NonunitPointError as result:
                factor, retry = result.factor, result.factor is None
                break
            divisor = gcd(point[1], n)
            if divisor != 1:
                factor = divisor if divisor < n else None
                retry = divisor == n
                break
        outputs.append((sigma, point, factor, retry))
    return outputs


def stage_control(bound, primes):
    controls = {}
    for sigma in SIGMAS:
        setup = ecm.setup_curve(prod(primes), sigma)
        if setup.factor or setup.retry:
            raise AssertionError("diagnostic fixture must have a usable curve")
        results = []
        for prime in primes:
            point, curve_a, curve_b = twist_point(
                setup.point, setup.a24, prime
            )
            for scalar in powers_for(bound):
                point = affine_multiply(scalar, point, prime, curve_a, curve_b)
            results.append(point)
        controls[sigma] = results
    return controls


def validate_stage(outputs, primes, controls):
    for sigma, point, factor, retry in outputs:
        if factor is not None:
            if not utils.valid_divisor(factor, prod(primes)):
                raise AssertionError("invalid factor")
        elif retry:
            # These fixed, nonsaturating fixtures deliberately isolate the
            # full stage cost. Saturation/replay belongs to other controls.
            raise AssertionError("unexpected saturation in diagnostic fixture")
        elif not all(
            matches(point, expected, prime)
            for prime, expected in zip(primes, controls[sigma])
        ):
            raise AssertionError("invalid stage-one output")


def measure(call, validate, repetitions, warmup_seconds):
    start = time.perf_counter()
    warmups = 0
    while time.perf_counter() - start < warmup_seconds:
        validate(call())
        warmups += 1
    warmup_elapsed = time.perf_counter() - start
    samples = []
    # Retain nine samples, extending noisy groups rather than claiming a win
    # from a single timing. Validation is outside each execution timer.
    for _ in range(repetitions):
        started = time.perf_counter()
        output = call()
        samples.append(time.perf_counter() - started)
        validate(output)
    median = statistics.median(samples)
    if (max(samples) - min(samples)) / median > 0.25:
        for _ in range(repetitions):
            started = time.perf_counter()
            output = call()
            samples.append(time.perf_counter() - started)
            validate(output)
    return {
        "seconds": samples,
        "median_seconds": statistics.median(samples),
        "min_seconds": min(samples),
        "max_seconds": max(samples),
        "warmup_seconds": warmup_elapsed,
        "validated_warmups": warmups,
    }


def construction_sample(bound):
    prac.clear_cache()
    return tuple(prac.get_chain(scalar) for scalar in powers_for(bound))


def kernel_sample(operation, fixtures, n, a24):
    output = None
    for _ in range(32):
        for point, twice in fixtures:
            output = (
                ecm.point_double(*point, n, a24)
                if operation == "double"
                else ecm.point_add(*point, *twice, *point, n)
            )
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--warmup-seconds", type=float, default=3.0)
    args = parser.parse_args()
    if platform.python_implementation() != "PyPy" or sys.version_info[:2] != (
        3,
        11,
    ):
        parser.error("requires PyPy implementing Python 3.11")
    if args.repetitions < 9 or args.warmup_seconds < 3:
        parser.error("requires at least nine samples and three seconds warmup")
    report = {
        "environment": environment(),
        "validation": validate_historical(),
        "construction": [],
        "stages": [],
        "kernels": [],
        "cold": [],
        "limitations": [
            "A4 diagnostic; no production/default promotion",
            "No portfolio, stage-two or whole-factor speed claim",
            "Large B1 can churn the finite scalar cache",
            "Full stage timings include checked execution costs",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        args.output.write_text(json.dumps(report, indent=2) + "\n")

    save()
    examples = report["validation"]["recovery_examples"]

    def recovery_sample():
        return [
            ecm.multiply_prac(k, point[0], 1, n, (a + 2) * pow(4, -1, n) % n)
            for k, n, a, point in examples
        ]

    recovery_controls = [
        affine_multiply(k, point, n, a) for k, n, a, point in examples
    ]

    def validate_recovery(outputs):
        for point, expected, (_, n, _, _) in zip(
            outputs, recovery_controls, examples
        ):
            if not matches(point, expected, n):
                raise AssertionError("recovery mismatch")

    report["recovery"] = measure(
        recovery_sample,
        validate_recovery,
        args.repetitions,
        args.warmup_seconds,
    )
    save()
    for bound in (128, 1000, 2000):
        row = {
            "b1": bound,
            **measure(
                lambda: construction_sample(bound),
                lambda records: [
                    prac.verify_chain(record) for record in records
                ],
                args.repetitions,
                args.warmup_seconds,
            ),
        }
        records = construction_sample(bound)
        row.update(
            records=len(records),
            instructions=sum(len(r.instructions) for r in records),
            abstract_prac_cost=sum(r.cost() for r in records),
            abstract_ladder_cost=sum(
                prac._binary_chain(r.scalar).cost() for r in records
            ),
        )
        report["construction"].append(row)
        save()
        print("construction", bound, row["median_seconds"], flush=True)

    for primes in MODULI:
        n = prod(primes)
        setup = ecm.setup_curve(n, 6)
        fixtures = [
            (
                ecm.scalar_multiply(k, *setup.point, n, setup.a24),
                ecm.scalar_multiply(2 * k, *setup.point, n, setup.a24),
            )
            for k in range(1, 33)
        ]
        for operation in ("add", "double"):
            expected = ecm.scalar_multiply(
                96 if operation == "add" else 64, *setup.point, n, setup.a24
            )

            def validate_kernel(actual):
                if (
                    gcd(gcd(*actual), n) != 1
                    or (actual[0] * expected[1] - actual[1] * expected[0]) % n
                ):
                    raise AssertionError("kernel mismatch")

            row = {
                "bits": n.bit_length(),
                "operation": operation,
                "operations_per_sample": 1024,
                **measure(
                    lambda: kernel_sample(operation, fixtures, n, setup.a24),
                    validate_kernel,
                    args.repetitions,
                    args.warmup_seconds,
                ),
            }
            report["kernels"].append(row)
            save()
        for bound in (128, 1000, 2000):
            controls = stage_control(bound, primes)
            for mode in ("ladder_chunks", "ladder_powers", "prac_powers"):
                prac.clear_cache()
                started = time.perf_counter()
                first = stage_sample(mode, bound, primes)
                first_seconds = time.perf_counter() - started
                validate_stage(first, primes, controls)
                row = {
                    "bits": n.bit_length(),
                    "b1": bound,
                    "mode": mode,
                    "curves": len(SIGMAS),
                    "first_call_seconds": first_seconds,
                    **measure(
                        lambda: stage_sample(mode, bound, primes),
                        lambda outputs: validate_stage(
                            outputs, primes, controls
                        ),
                        args.repetitions,
                        args.warmup_seconds,
                    ),
                }
                row["cache"] = prac.cache_info()._asdict()
                row["outcomes"] = [
                    {
                        "sigma": sigma,
                        "factor": factor,
                        "retry": retry,
                        "full_stage_completed": factor is None and not retry,
                    }
                    for sigma, _, factor, retry in first
                ]
                report["stages"].append(row)
                save()
                print(
                    "stage",
                    n.bit_length(),
                    bound,
                    mode,
                    row["median_seconds"],
                    flush=True,
                )

    # Fresh-process samples include PyPy startup/import, first chain and exit.
    command = (
        "import json; from v2 import ecm; "
        "print(json.dumps(ecm.multiply_prac(1009,3,1,1009,2)))"
    )
    expected = affine_multiply(1009, (3, 293), 1009, 6)
    for _ in range(args.repetitions):
        started = time.perf_counter()
        child = subprocess.run(
            [sys.executable, "-c", command],
            check=True,
            text=True,
            capture_output=True,
        )
        elapsed = time.perf_counter() - started
        if not matches(json.loads(child.stdout), expected, 1009):
            raise AssertionError("cold-process mismatch")
        report["cold"].append(elapsed)
    save()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
