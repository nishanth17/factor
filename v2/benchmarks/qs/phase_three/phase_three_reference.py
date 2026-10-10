"""P3.1 costs and oracle comparisons; no SIQS factorization promotion."""

import argparse
import hashlib
import json
import math
import random
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

from ....execution.budget import Budget
from ....portfolio import PortfolioConfig, factorize_bounded
from ....qs import (
    build_factor_base,
    collect_block,
    combine_relations,
    modular_square_roots,
    mpqs_polynomial,
    polynomial_roots,
    qs_polynomial,
    verify_atomic,
)
from ....qs.polynomial import Polynomial
from ....tests.test_qs import reference_positions, reference_primes
from ...suites.phase_one import environment
from ...support.paths import (
    REPOSITORY_ROOT,
    source_path,
)

CORPUS_SEED = 20261003
ROOT = REPOSITORY_ROOT
BASELINE = (
    ROOT / "v2/benchmarks/inputs/controls/phase_two_m17_frozen_baseline.json"
)


def _rss_bytes():
    """Return process high-water RSS, including JIT and untimed warmup."""
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def _budget():
    """Declare finite per-call controls instead of inheriting M17's 50 ms."""
    return Budget(work_limit=20_000_000, seconds=5, cpu_seconds=5)


def _collect(n, multiplier, a, b, residual_bound):
    """Measure setup through verified collection under one shared allowance."""
    budget = _budget()
    base = build_factor_base(
        n, multiplier=multiplier, bound=100, budget=budget
    ).factor_base
    polynomial = Polynomial(n, multiplier, a, b)
    result = collect_block(
        polynomial,
        base,
        -128,
        129,
        residual_bound=residual_bound,
        budget=budget,
    )

    if result.reason != "complete" or result.divisor is not None:
        raise AssertionError("reference collection did not complete")
    signature = {
        atom.position: (atom.sign, atom.exponents, atom.residual)
        for atom in result.relations
    }
    return signature, result.workspace_bytes, budget.used


def _cases():
    """Freeze small diagnostic inputs and independent expected values."""
    generator = random.Random(CORPUS_SEED)
    root_inputs = [
        (generator.getrandbits(64), prime)
        for prime in reference_primes(150)
        for _ in range(4)
    ]

    def reference_roots():
        """Enumerate residues without using the modular-root algorithm."""
        return [
            tuple(x for x in range(prime) if x * x % prime == value % prime)
            for value, prime in root_inputs
        ]

    cases = [
        (
            "modular_roots_140",
            {
                "qs": lambda: [
                    modular_square_roots(v, p) for v, p in root_inputs
                ],
                "exhaustive_oracle": reference_roots,
            },
            reference_roots(),
            10,
        )
    ]
    setup_inputs = [(10403, h) for h in (1, 3, 9)] + [(1022117, 1)]

    def reference_bases():
        """Enumerate every prime/root directly for the same setup inputs."""
        answer = []

        for n, h in setup_inputs:
            entries = []

            for prime in reference_primes(100):
                roots = tuple(
                    x for x in range(prime) if x * x % prime == n * h % prime
                )
                if roots:
                    entries.append((prime, roots))

            answer.append(tuple(entries))

        return answer

    def native_bases():
        """Include primality/root validation and factor-base construction."""
        answer = []

        for n, h in setup_inputs:
            base = build_factor_base(
                n, multiplier=h, bound=100, budget=_budget()
            ).factor_base
            answer.append(
                tuple((e.prime, e.square_roots) for e in base.entries)
            )

        return answer

    cases.append(
        (
            "factor_base_setup_4",
            {
                "qs": native_bases,
                "exhaustive_oracle": reference_bases,
            },
            reference_bases(),
            10,
        )
    )
    base = build_factor_base(10403, bound=100).factor_base
    polynomials = [
        qs_polynomial(base),
        Polynomial(10403, 1, 7, 1),
        Polynomial(10403, 1, 49, 8),
    ]

    def reference_polynomial_roots():
        """Enumerate the normalized polynomial at every residue."""
        return [
            tuple(
                x
                for x in range(entry.prime)
                if polynomial.value(x) % entry.prime == 0
            )
            for polynomial in polynomials
            for entry in base.entries
        ]

    def native_polynomial_roots():
        """Use cached target roots and exceptional normalized branches."""
        answers = []

        for polynomial in polynomials:
            for entry in base.entries:
                answer = polynomial_roots(polynomial, base, entry)
                answers.append(
                    tuple(range(entry.prime))
                    if answer.all_positions
                    else answer.roots
                )

        return answers

    cases.append(
        (
            "normalized_polynomial_roots",
            {
                "qs": native_polynomial_roots,
                "exhaustive_oracle": reference_polynomial_roots,
            },
            reference_polynomial_roots(),
            10,
        )
    )
    widths = (1, 16, 128)
    expected_lifts = []

    for width in widths:
        target = max(1, math.isqrt(2 * base.n_prime) // width)
        prime = min(
            (
                entry.prime
                for entry in base.entries
                if entry.prime != 2 and base.n_prime % entry.prime
            ),
            key=lambda p: (abs(p * p - target), p),
        )
        a = prime * prime
        b = next(b for b in range(a) if b * b % a == base.n_prime % a)
        expected_lifts.append((a, b, (b * b - base.n_prime) // a))

    def native_lifts():
        """Include integer target, nearest-A selection, inverse and lift."""
        polynomials = [
            mpqs_polynomial(base, width, budget=_budget()) for width in widths
        ]
        return [
            (polynomial.a, polynomial.b, polynomial.c)
            for polynomial in polynomials
        ]

    cases.append(
        (
            "mpqs_integer_selection_lift",
            {"qs": native_lifts},
            expected_lifts,
            10,
        )
    )
    fixtures = []

    for name, n, h, a, b, residual in (
        ("qs_full", 10403, 1, 1, 102, 1),
        ("mpqs_full", 10403, 1, 49, 8, 1),
        ("qs_partial", 1022117, 3, 1, 1752, 500),
    ):
        local_base = build_factor_base(n, multiplier=h, bound=100).factor_base
        polynomial = Polynomial(n, h, a, b)
        expected = reference_positions(
            polynomial, local_base, -128, 129, residual
        )

        if not expected:
            raise AssertionError("diagnostic corpus must produce useful atoms")
        fixtures.append(
            {
                "name": name,
                "n": n,
                "h": h,
                "a": a,
                "b": b,
                "lo": -128,
                "hi": 129,
                "factor_base_bound": 100,
                "residual_bound": residual,
                "verified_yield": len(expected),
            }
        )
        cases.append(
            (
                name + "_setup_to_verified_collection",
                {
                    "qs": lambda n=n, h=h, a=a, b=b, r=residual: _collect(
                        n, h, a, b, r
                    )[0],
                    "trial_factor_oracle": lambda p=polynomial,
                    fb=local_base,
                    r=residual: reference_positions(p, fb, -128, 129, r),
                },
                expected,
                3,
            )
        )

    atoms = collect_block(
        qs_polynomial(base), base, -128, 129, residual_bound=500
    ).relations
    cases.append(
        (
            "atomic_verification",
            {
                "qs": lambda: [
                    verify_atomic(
                        atom, base, residual_bound=500, budget=_budget()
                    )
                    for atom in atoms
                ],
            },
            [True] * len(atoms),
            10,
        )
    )
    pairs, seen = [], {}

    for atom in atoms:
        if atom.residual > 1 and math.gcd(atom.residual, base.n) == 1:
            if atom.residual in seen:
                pairs.append((seen.pop(atom.residual), atom))
            else:
                seen[atom.residual] = atom

    if not pairs:
        raise AssertionError("combination corpus must contain matched atoms")

    def combinations():
        """Include both repeated atomic verification and combined checking."""
        return [
            combine_relations(pair, base, budget=_budget()).relation
            for pair in pairs
        ]

    expected_combinations = combinations()

    for relation, pair in zip(expected_combinations, pairs):
        actual = math.prod(atom.u**2 - base.n_prime for atom in pair)
        reconstructed = relation.sign * pair[0].residual ** 2
        for prime, exponent in relation.exponents:
            reconstructed *= prime**exponent
        if actual != reconstructed:
            raise AssertionError(
                "independent two-norm combination check failed"
            )

    cases.append(
        (
            "checked_partial_combination",
            {"qs": combinations},
            expected_combinations,
            3,
        )
    )

    control_inputs = [
        25013 * 25031,
        1000003 * 1000033,
        1009**3 * 1013**2,
        -(2**16 * 3**6 * 101),
        2147483647,
    ]
    config = PortfolioConfig()

    def complete_control():
        """Measure complete factoring with unchanged frozen production code."""
        answers = []
        for n in control_inputs:
            run = factorize_bounded(n, seed=7, config=config, budget=_budget())

            if not run.result.complete or run.result.reconstruct() != n:
                raise AssertionError("complete control failed reconstruction")
            answers.append(run.result.reconstruct())

        return answers

    cases.append(
        (
            "unchanged_complete_factorization_control",
            {
                "m17_runtime": complete_control,
            },
            control_inputs,
            3,
        )
    )
    corpus = {
        "seed": CORPUS_SEED,
        "root_inputs": root_inputs,
        "setup_inputs": setup_inputs,
        "collection": fixtures,
        "verification_atoms": len(atoms),
        "combination_pairs": len(pairs),
        "complete_inputs": control_inputs,
        "complete_seed": 7,
        "scope": "small diagnostic costs; no tuned or held-out promotion",
    }
    return cases, corpus


def _measure(function, expected, iterations, repetitions, warmup_seconds):
    """Validate warmups/samples and retain adaptive stabilization attempts."""
    attempts = []

    for attempt in range(2):
        started, calls = time.perf_counter(), 0
        duration = warmup_seconds if attempt == 0 else max(5, warmup_seconds)
        while time.perf_counter() - started < duration:
            if function() != expected:
                raise AssertionError("invalid warmup output")
            calls += 1
        warmup = time.perf_counter() - started
        samples = []
        count = repetitions if attempt == 0 else max(15, repetitions)

        for _ in range(count):
            started, cpu_started = time.perf_counter(), time.process_time()
            for _ in range(iterations):
                if function() != expected:
                    raise AssertionError("invalid measured output")
            samples.append(
                {
                    "seconds_per_call": (time.perf_counter() - started)
                    / iterations,
                    "cpu_seconds_per_call": (time.process_time() - cpu_started)
                    / iterations,
                    "peak_rss_bytes": _rss_bytes(),
                    "correct": True,
                }
            )

        times = [sample["seconds_per_call"] for sample in samples]
        median = statistics.median(times)
        quartiles = statistics.quantiles(times, n=4)
        drift = abs(
            statistics.median(times[:3]) / statistics.median(times[-3:]) - 1
        )
        relative_iqr = (quartiles[2] - quartiles[0]) / median
        stable = drift <= 0.15 and relative_iqr <= 0.20
        attempts.append(
            {
                "warmup_seconds": warmup,
                "warmup_calls": calls,
                "samples": samples,
                "median_seconds": median,
                "min_seconds": min(times),
                "max_seconds": max(times),
                "relative_iqr": relative_iqr,
                "early_late_median_drift": drift,
                "stable": stable,
            }
        )
        if stable:
            break

    return {
        "iterations_per_sample": iterations,
        "attempts": attempts,
        "median_seconds": attempts[-1]["median_seconds"],
        "stable": attempts[-1]["stable"],
    }


def _cold_worker():
    """Emit a consumed verified collection for isolated lifecycle timing."""
    signature, workspace, work = _collect(10403, 1, 1, 102, 1)
    return {
        "signature": signature,
        "workspace_bytes": workspace,
        "work_used": work,
        "peak_rss_bytes": _rss_bytes(),
    }


def main():
    """Capture finite warm/cold costs and exact source/runtime provenance."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--cold-worker", action="store_true")
    args = parser.parse_args()
    if args.cold_worker:
        print(json.dumps(_cold_worker()))
        return
    if args.output is None:
        parser.error("--output is required")
    if (
        not math.isfinite(args.warmup_seconds)
        or args.warmup_seconds < 3
        or args.repetitions < 9
    ):
        parser.error("use at least three seconds warmup and nine samples")

    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("supported runtime is PyPy Python 3.11")
    measured_environment = environment()
    frozen = json.loads(source_path(BASELINE).read_text())
    unchanged = all(
        hashlib.sha256((source_path(ROOT / path)).read_bytes()).hexdigest()
        == digest
        for path, digest in frozen["production_source_sha256"].items()
    )

    if not unchanged:
        raise AssertionError("M17 production control changed")
    cases, corpus = _cases()
    rows = []

    for name, candidates, expected, iterations in cases:
        results = {}

        for label, function in candidates.items():
            results[label] = _measure(
                function,
                expected,
                iterations,
                args.repetitions,
                args.warmup_seconds,
            )
            print(
                name,
                label,
                f"{results[label]['median_seconds'] * 1000:.3f} ms",
                f"stable={results[label]['stable']}",
                flush=True,
            )

        rows.append({"name": name, "candidates": results})

    cold = []
    expected_signature = json.loads(
        json.dumps(
            next(
                expected
                for name, _, expected, _ in cases
                if name == "qs_full_setup_to_verified_collection"
            )
        )
    )

    for _ in range(args.repetitions):
        started = time.perf_counter()
        cpu_started = resource.getrusage(resource.RUSAGE_CHILDREN)

        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "v2.benchmarks.qs.phase_three.phase_three_reference",
                "--cold-worker",
            ],
            capture_output=True,
            text=True,
            check=True,
            timeout=30,
        )
        response = json.loads(completed.stdout)
        elapsed = time.perf_counter() - started
        cpu_finished = resource.getrusage(resource.RUSAGE_CHILDREN)

        if response["signature"] != expected_signature:
            raise AssertionError("cold worker returned an invalid signature")
        cold.append(
            {
                "lifecycle_seconds": elapsed,
                "cpu_seconds": cpu_finished.ru_utime
                + cpu_finished.ru_stime
                - cpu_started.ru_utime
                - cpu_started.ru_stime,
                "correct": True,
                **response,
            }
        )

    if environment()["source_sha256"] != measured_environment["source_sha256"]:
        raise AssertionError("measured source changed during capture")
    encoded_corpus = json.dumps(corpus, sort_keys=True).encode()
    args.output.write_text(
        json.dumps(
            {
                "milestone": "M20 / P3.1",
                "environment": measured_environment,
                "command": sys.orig_argv,
                "corpus": corpus,
                "corpus_sha256": hashlib.sha256(encoded_corpus).hexdigest(),
                "production_m17_hashes_unchanged": unchanged,
                "pipeline_or_control_input_limits": {
                    "work": 20_000_000,
                    "wall_seconds": 5,
                    "cpu_seconds": 5,
                    "owned_bytes": 8_388_608,
                },
                "jit_policy": "default JIT; validated workload warmup",
                "timing_policy": (
                    "include output consumption and validation; "
                    "extra warmup/samples on drift"
                ),
                "rss_policy": (
                    "high-water includes JIT, oracles and warmup; "
                    "owned bytes do not imply an RSS cap"
                ),
                "promotion": (
                    "none; no previous QS implementation or full extraction"
                ),
                "results": rows,
                "cold_samples": cold,
                "cold_median_seconds": statistics.median(
                    sample["lifecycle_seconds"] for sample in cold
                ),
                "cold_scope": (
                    "isolated harness startup/imports, setup, "
                    "collection and output consumption"
                ),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
