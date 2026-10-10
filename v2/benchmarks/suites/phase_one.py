"""Measure native v2 and optionally emulated v1 on the same interpreter.

Every timing sample retains its output validation. Invalid legacy work is
recorded but never receives a speedup ratio. This small Phase 1 corpus is not
the broad, held-out portfolio benchmark planned for Phase 2.
"""

import argparse
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
from bisect import bisect_right
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

from ...common import prime_sieve, utils
from ...ecm import core as ecm
from ...factor import factorize
from ..support.paths import (
    REPOSITORY_ROOT,
    source_path,
)

CORPUS_SEED = 20261003


def environment():
    """Identify the runtime and measured sources, excluding dependencies."""
    root = REPOSITORY_ROOT
    paths = sorted(
        path
        for path in (root / "v2").rglob("*.py")
        if not {"audit", ".venv"}.intersection(
            path.relative_to(root / "v2").parts
        )
    )
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True
    ).strip()
    jit_defaults = None
    translation = None
    if hasattr(sys, "pypy_version_info"):
        import pypyjit

        jit_defaults = getattr(pypyjit, "defaults", None)
        translation = getattr(sys, "pypy_translation_info", None)
    return {
        "python": sys.version,
        "implementation": platform.python_implementation(),
        "pypy_version": list(sys.pypy_version_info)
        if hasattr(sys, "pypy_version_info")
        else None,
        "pypy_jit_defaults": jit_defaults,
        "pypy_translation": translation,
        "jit_overrides": os.environ.get("PYPYJIT"),
        "runtime_command": sys.orig_argv
        if hasattr(sys, "orig_argv")
        else sys.argv,
        "executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "visible_cpu_count": os.cpu_count(),
        "core_budget": 1,
        "git_head": commit,
        "working_tree": "uncommitted source identified by hashes below",
        "source_sha256": {
            str(path.relative_to(root)): hashlib.sha256(
                source_path(path).read_bytes()
            ).hexdigest()
            for path in paths
        },
    }


def _reference_primes(hi):
    """Return p < hi with an independent full-width Eratosthenes control."""
    flags = bytearray(b"\x01") * hi
    flags[:2] = b"\x00\x00"
    for prime in range(2, math.isqrt(hi - 1) + 1):
        if flags[prime]:
            start = prime * prime
            flags[start::prime] = b"\x00" * ((hi - 1 - start) // prime + 1)
    return [index for index in range(2, hi) if flags[index]]


def _measure_case(name, candidates, validate, repetitions, warmup_seconds):
    """Warm and time validated outputs; refuse invalid native candidates."""
    row = {"name": name, "candidates": {}}
    usable = {}

    for label, function in candidates.items():
        start = time.perf_counter()

        try:
            value = function()
            correct = bool(validate(value))
            warmup_calls = 1
            if correct:
                # Warm the actual workload, including JIT traces, not a sleep.
                while time.perf_counter() - start < warmup_seconds:
                    if not validate(function()):
                        raise AssertionError("invalid result during warmup")
                    warmup_calls += 1

            elapsed = time.perf_counter() - start
            row["candidates"][label] = {
                "correct": correct,
                "warmup_seconds": elapsed,
                "warmup_calls": warmup_calls,
                "samples_seconds": [],
            }
            if correct:
                usable[label] = function
        except Exception as error:
            row["candidates"][label] = {
                "correct": False,
                "error": repr(error),
                "samples_seconds": [],
            }

    if not row["candidates"]["v2_native"]["correct"]:
        raise AssertionError(f"invalid native benchmark: {name}")

    # Rotate order to reduce systematic thermal/cache-order bias.
    labels = list(usable)

    for round_index in range(repetitions):
        offset = round_index % len(labels) if labels else 0
        order = labels[offset:] + labels[:offset]

        for label in order:
            start = time.perf_counter()
            value = usable[label]()
            elapsed = time.perf_counter() - start

            # Check every timed output, keeping oracle cost outside timing.
            if not validate(value):
                raise AssertionError(
                    f"invalid measured result: {name}/{label}"
                )
            row["candidates"][label]["samples_seconds"].append(elapsed)

    for result in row["candidates"].values():
        samples = result["samples_seconds"]
        if samples:
            result["median_seconds"] = statistics.median(samples)
            result["min_seconds"] = min(samples)
            result["max_seconds"] = max(samples)

    if "v1_compat" in usable and "v2_native" in usable:
        previous = row["candidates"]["v1_compat"]["median_seconds"]
        current = row["candidates"]["v2_native"]["median_seconds"]
        row["v1_over_v2_time_ratio"] = previous / current
        row["v2_time_change_percent"] = (current / previous - 1) * 100
    return row


def run_benchmarks(repetitions=5, include_legacy=False, warmup_seconds=3.0):
    """Return seeded, matched batch timings and their acceptance evidence."""
    legacy = None
    if include_legacy:
        from ..support.legacy_loader import load_legacy

        legacy = load_legacy()
    generator = random.Random(CORPUS_SEED)
    results = []

    def measure(name, current, previous, validate):
        """Collect one case and emit compact progress after its warmup."""
        candidates = {"v2_native": current}
        if legacy is not None and previous is not None:
            candidates["v1_compat"] = previous
        row = _measure_case(
            name, candidates, validate, repetitions, warmup_seconds
        )
        results.append(row)
        print(
            name,
            json.dumps(
                {
                    label: round(item.get("median_seconds", 0) * 1000, 3)
                    if item["correct"]
                    else "invalid"
                    for label, item in row["candidates"].items()
                }
            ),
        )

    pairs = [
        (generator.getrandbits(256), generator.getrandbits(256))
        for _ in range(3000)
    ]
    expected = [math.gcd(a, b) for a, b in pairs]
    measure(
        "gcd_256bit_3000",
        lambda: [utils.gcd(a, b) for a, b in pairs],
        (lambda: [legacy["utils"].gcd(a, b) for a, b in pairs]),
        lambda value: value == expected,
    )

    modulus = 2**192 - 237
    residues = []
    while len(residues) < 256:
        value = generator.randrange(2, modulus)
        if math.gcd(value, modulus) == 1:
            residues.append(value)
    expected = [pow(value, -1, modulus) for value in residues]
    measure(
        "inverse_192bit_256",
        lambda: [utils.modular_inverse(value, modulus) for value in residues],
        lambda: [
            legacy["utils"].xgcd(modulus, value) % modulus
            for value in residues
        ],
        lambda value: value == expected,
    )

    array = list(range(0, 10000, 2))
    searches = [generator.randrange(-1, 10001) for _ in range(3000)]
    expected = [bisect_right(array, value) for value in searches]
    measure(
        "binary_search_distinct_3000",
        lambda: [utils.binary_search(value, array) for value in searches],
        lambda: [
            legacy["utils"].binary_search(value, array) for value in searches
        ],
        lambda value: value == expected,
    )

    primality_inputs = [
        1009,
        10007,
        65537,
        104729,
        2147483647,
        1000000007,
        999999937,
        2**61 - 1,
        1009 * 1013,
        341550071728321,
    ] * 30
    expected = [True] * 8 + [False] * 2
    expected *= 30
    measure(
        "primality_mixed_300",
        lambda: [utils.is_prime(n) for n in primality_inputs],
        lambda: [bool(legacy["utils"].is_prime(n)) for n in primality_inputs],
        lambda value: value == expected,
    )

    for hi in (100_000, 1_000_000, 3_500_001):
        expected = _reference_primes(hi)
        measure(
            f"prime_sieve_{hi}",
            lambda hi=hi: prime_sieve.prime_sieve(hi),
            lambda hi=hi: legacy["primeSieve"].prime_sieve(hi),
            lambda value, expected=expected: value == expected,
        )

    # Compare two correct v2 backends before retaining the default dispatcher.
    expected = _reference_primes(3_500_001)
    measure(
        "corrected_atkin_3500001",
        lambda: prime_sieve.sieve_of_atkin(3_500_001),
        None,
        lambda value: value == expected,
    )

    root_generator = random.Random(CORPUS_SEED + 1)
    root_inputs = [root_generator.getrandbits(256) for _ in range(500)]
    expected_roots = [math.isqrt(n) for n in root_inputs]
    measure(
        "exact_roots_256bit_500",
        lambda: [utils.isqrt(n) for n in root_inputs],
        lambda: [int(math.sqrt(n)) for n in root_inputs],
        lambda value: value == expected_roots,
    )

    for bits in (64, 192, 256, 512):
        modulus = (1 << bits) - 159
        values = [generator.randrange(modulus) for _ in range(7)]
        arguments = values[:2] + [modulus, values[6]]
        u_squared = (values[0] + values[1]) ** 2
        v_squared = (values[0] - values[1]) ** 2
        delta = u_squared - v_squared
        expected = (
            u_squared * v_squared % modulus,
            delta * (v_squared + values[6] * delta) % modulus,
        )
        measure(
            f"ecm_double_{bits}bit_3000",
            lambda: [ecm.point_double(*arguments) for _ in range(3000)],
            lambda: [
                legacy["ecm"].point_double(*arguments) for _ in range(3000)
            ],
            lambda value, expected=expected: all(v == expected for v in value),
        )

    inputs = [
        25013 * 25031,
        1000003 * 1000033,
        1009**3 * 1013**2,
        2**16 * 3**6 * 101,
        2147483647,
    ]

    def current_factorizations():
        """Run the same five inputs for each independently seeded v2 call."""
        return [factorize(n, seed=seed) for seed in range(5) for n in inputs]

    def previous_factorizations():
        """Run emulated v1 without its historical unsolicited output."""
        answers = []
        with redirect_stdout(StringIO()):
            for seed in range(5):
                for n in inputs:
                    random.seed(seed)
                    answers.append(legacy["factor"].factorize(n))

        return answers

    def validate_factorizations(answers):
        """Require 25 complete, reconstructible results from either API."""
        for n, answer in zip(inputs * 5, answers):
            if hasattr(answer, "complete"):
                if not answer.complete or answer.reconstruct() != n:
                    return False
            elif answer == -1 or math.prod(p**e for p, e in answer) != n:
                return False

        return len(answers) == 25

    measure(
        "complete_factorizations_5inputs_5seeds",
        current_factorizations,
        previous_factorizations,
        validate_factorizations,
    )

    correctness = {"native_p_minus_one_fixture": None}
    from ...pm1.core import factorize_pm1

    correctness["native_p_minus_one_fixture"] = factorize_pm1(
        607 * 1019, b1=10, b2=200, max_attempts=1
    )
    if legacy is not None:
        previous_bounds = legacy["pollardPm1"].compute_bounds
        legacy["pollardPm1"].compute_bounds = lambda n: (10, 200)

        try:
            correctness["legacy_p_minus_one_fixture"] = legacy[
                "pollardPm1"
            ].factorize_pm1(607 * 1019)
        finally:
            legacy["pollardPm1"].compute_bounds = previous_bounds
    return {
        "schema_version": 1,
        "milestone": "phase_1",
        "environment": environment(),
        "corpus_seed": CORPUS_SEED,
        "repetitions": repetitions,
        "warmup_min_seconds_per_candidate": warmup_seconds,
        "warmup_policy": "validated workload repeats before timed samples",
        "jit_policy": "default runtime settings; no JIT override in command",
        "cache_state": "warm modules/trial cache; fresh sieve buffers",
        "legacy_method": "lib2to3 + integer division emulation + math.gcd"
        if include_legacy
        else None,
        "limitations": [
            "Single-core warm microbenchmarks and 25 small complete runs",
            "Emulated v1 is not native Python 2; no 50-60 digit timing claim",
            "Phase 1 defaults differ; full-run ratios do not isolate one fix",
            "Legacy GCD already uses math.gcd through the adapter",
            "Invalid outputs are excluded from comparative speedup ratios",
            "No held-out parameter tuning or statistical promotion claim",
        ],
        "correctness_observations": correctness,
        "benchmarks": results,
    }


def main():
    """Validate timing options and save the chosen runtime's raw evidence."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-legacy", action="store_true")
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmup-seconds", type=float, default=3.0)
    args = parser.parse_args()
    if not math.isfinite(args.warmup_seconds) or args.warmup_seconds < 0:
        parser.error("warmup seconds must be finite and nonnegative")
    if args.repetitions < 1:
        parser.error("repetitions must be positive")
    data = run_benchmarks(
        args.repetitions, args.include_legacy, args.warmup_seconds
    )
    args.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
