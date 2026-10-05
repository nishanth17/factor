"""Measure M9 versus frozen M8 and optionally emulated v1 in one runtime."""

import argparse
import hashlib
import json
import math
import random
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

from v2 import ecm, prime_sieve, utils
from v2.factor import factorize

from .phase_one import _measure_case, _reference_primes, environment
from .snapshot_loader import load_snapshot

SEED = 20261003


def _cases():
    """Return the original cohort and a separately seeded broader control."""
    original = [
        (25013 * 25031, ((25013, 1), (25031, 1))),
        (1000003 * 1000033, ((1000003, 1), (1000033, 1))),
        (1009**3 * 1013**2, ((1009, 3), (1013, 2))),
        (2**16 * 3**6 * 101, ((2, 16), (3, 6), (101, 1))),
        (2147483647, ((2147483647, 1),)),
    ]
    generator = random.Random(20261005)
    primes = [p for p in _reference_primes(1_100_000) if p > 100_000]
    broader = []
    for _ in range(30):
        p, q = sorted(generator.sample(primes, 2))
        broader.append((p * q, ((p, 1), (q, 1))))
    for p in generator.sample(primes, 10):
        broader.append((1009 * p, ((1009, 1), (p, 1))))
    for p in generator.sample(primes, 10):
        broader.append((p, ((p, 1),)))
    for p in (2147483647, 2**61 - 1, 2**64 - 59):
        broader.append((p, ((p, 1),)))
    for p in (1009, 65537, 1000003):
        broader.append((p**2, ((p, 2),)))
    return original, broader


def _factor_calls(function, cases, legacy=False):
    """Pass integers/seeds to algorithms, keeping oracle data hidden."""
    answers = []
    with redirect_stdout(StringIO()):
        for seed in range(5):
            for n, _ in cases:
                if legacy:
                    random.seed(seed)
                    answer = function(n)
                else:
                    answer = function(n, seed=seed)

                answers.append(answer)

    return answers


def _valid_factorizations(answers, cases):
    """Require exact known prime multiplicities, proof, and reconstruction."""
    if len(answers) != len(cases) * 5:
        return False
    for answer, (n, expected) in zip(answers, cases * 5):
        if hasattr(answer, "complete"):
            if not answer.complete or not answer.proven:
                return False
            if answer.reconstruct() != n:
                return False
            actual = tuple((f.value, f.exponent) for f in answer.factors)
        elif isinstance(answer, list):
            actual = tuple(sorted(answer))
        else:
            return False

        if actual != expected:
            return False

    return True


def run(repetitions, warmup_seconds, include_legacy):
    """Warm all candidates and retain only ratios for validated outputs."""
    before = load_snapshot()
    legacy = None
    if include_legacy:
        from .legacy_loader import load_legacy

        legacy = load_legacy()
    rows = []

    def measure(name, current, previous, original, validate):
        """Rotate candidate order and attach exact before/after deltas."""
        candidates = {"v2_native": current, "v2_before": previous}
        if legacy is not None and original is not None:
            candidates["v1_compat"] = original
        row = _measure_case(
            name, candidates, validate, repetitions, warmup_seconds
        )

        for baseline in ("v2_before", "v1_compat"):
            result = row["candidates"].get(baseline)
            if result is not None and result["correct"]:
                old = result["median_seconds"]
                new = row["candidates"]["v2_native"]["median_seconds"]
                row[f"{baseline}_over_v2_time_ratio"] = old / new
                row[f"v2_time_change_vs_{baseline}_percent"] = (
                    new / old - 1
                ) * 100

        rows.append(row)
        print(
            name,
            {
                label: round(item.get("median_seconds", 0) * 1000, 3)
                if item["correct"]
                else "invalid"
                for label, item in row["candidates"].items()
            },
            flush=True,
        )

    generator = random.Random(SEED)
    pairs = [
        (generator.getrandbits(256), generator.getrandbits(256))
        for _ in range(3000)
    ]
    expected = [math.gcd(a, b) for a, b in pairs]
    measure(
        "gcd_256bit_3000",
        lambda: [utils.gcd(a, b) for a, b in pairs],
        lambda: [before.utils.gcd(a, b) for a, b in pairs],
        lambda: [legacy["utils"].gcd(a, b) for a, b in pairs],
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
        lambda: [utils.modular_inverse(a, modulus) for a in residues],
        lambda: [before.utils.modular_inverse(a, modulus) for a in residues],
        lambda: [legacy["utils"].xgcd(modulus, a) % modulus for a in residues],
        lambda value: value == expected,
    )

    array = list(range(0, 10000, 2))
    searches = [generator.randrange(-1, 10001) for _ in range(3000)]
    expected = [before.utils.binary_search(n, array) for n in searches]
    measure(
        "binary_search_distinct_3000",
        lambda: [utils.binary_search(n, array) for n in searches],
        lambda: [before.utils.binary_search(n, array) for n in searches],
        lambda: [legacy["utils"].binary_search(n, array) for n in searches],
        lambda value: value == expected,
    )

    inputs = [
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
    expected = ([True] * 8 + [False] * 2) * 30
    measure(
        "primality_mixed_300",
        lambda: [utils.is_prime(n) for n in inputs],
        lambda: [before.utils.is_prime(n) for n in inputs],
        lambda: [bool(legacy["utils"].is_prime(n)) for n in inputs],
        lambda value: value == expected,
    )

    for hi in (100_000, 1_000_000, 3_500_001):
        expected = _reference_primes(hi)
        measure(
            f"prime_sieve_{hi}",
            lambda hi=hi: prime_sieve.prime_sieve(hi),
            lambda hi=hi: before.prime_sieve.prime_sieve(hi),
            lambda hi=hi: legacy["primeSieve"].prime_sieve(hi),
            lambda value, expected=expected: value == expected,
        )

    root_rng = random.Random(SEED + 1)
    roots = [root_rng.getrandbits(256) for _ in range(500)]
    expected = [math.isqrt(n) for n in roots]
    measure(
        "exact_roots_256bit_500",
        lambda: [utils.isqrt(n) for n in roots],
        lambda: [before.utils.isqrt(n) for n in roots],
        lambda: [int(math.sqrt(n)) for n in roots],
        lambda value: value == expected,
    )

    for bits in (64, 192, 256, 512):
        modulus = (1 << bits) - 159
        values = [generator.randrange(modulus) for _ in range(7)]
        arguments = values[:2] + [modulus, values[6]]
        u2 = (values[0] + values[1]) ** 2
        v2 = (values[0] - values[1]) ** 2
        delta = u2 - v2
        expected = (
            u2 * v2 % modulus,
            delta * (v2 + values[6] * delta) % modulus,
        )
        measure(
            f"ecm_double_{bits}bit_3000",
            lambda: [ecm.point_double(*arguments) for _ in range(3000)],
            lambda: [before.ecm.point_double(*arguments) for _ in range(3000)],
            lambda: [
                legacy["ecm"].point_double(*arguments) for _ in range(3000)
            ],
            lambda value, expected=expected: all(
                item == expected for item in value
            ),
        )

    original, broader = _cases()

    for name, cases in (
        ("complete_original_5inputs_5seeds", original),
        ("complete_control_56inputs_5seeds", broader),
    ):
        measure(
            name,
            lambda cases=cases: _factor_calls(factorize, cases),
            lambda cases=cases: _factor_calls(before.factor.factorize, cases),
            lambda cases=cases: _factor_calls(
                legacy["factor"].factorize, cases, legacy=True
            ),
            lambda answers, cases=cases: _valid_factorizations(answers, cases),
        )

    snapshot = Path(__file__).resolve().parents[1] / "benchmarks"
    snapshot /= "inputs/baselines/m8_source_snapshot.json"
    return {
        "milestone": "m9",
        "environment": environment(),
        "baseline_snapshot_sha256": hashlib.sha256(
            snapshot.read_bytes()
        ).hexdigest(),
        "baseline_source_sha256": json.loads(snapshot.read_text())[
            "source_sha256"
        ],
        "warmup_seconds": warmup_seconds,
        "repetitions": repetitions,
        "seed": SEED,
        "broader_seed": 20261005,
        "factor_corpus": {"original": original, "broader": broader},
        "legacy_method": "lib2to3 + exact Py2 division adapter + math.gcd"
        if include_legacy
        else None,
        "limitations": [
            "Same-runtime warm before/after; invalid v1 has no speed ratio",
            "56-input control is small feasible validation, not broad tuning",
            "No deadline, aggregate RSS, or large-balanced completion claim",
            "Validated workload warmup per candidate then sample rotation",
        ],
        "benchmarks": rows,
    }


def main():
    """Save uniquely named evidence for the selected interpreter."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--legacy", action="store_true")
    parser.add_argument("--repetitions", type=int, default=15)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    args = parser.parse_args()
    if args.repetitions < 1 or not math.isfinite(args.warmup_seconds):
        parser.error("positive repetitions and finite warmup are required")
    if args.warmup_seconds < 0:
        parser.error("warmup must be nonnegative")
    data = run(args.repetitions, args.warmup_seconds, args.legacy)
    args.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
