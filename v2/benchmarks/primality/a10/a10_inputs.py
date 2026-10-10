"""Independent frozen prime proofs and composite divisors for A10.

These are benchmark/test oracles, not a production certificate API (B14).
Candidate filtering during generation is never accepted as expected truth.
"""

import argparse
import json
import random
from math import gcd, isqrt, prod
from pathlib import Path

from ...suites.build_phase_two_corpus import certified_prime
from ...support.paths import (
    BENCHMARK_ROOT,
    source_path,
)

CORPUS = BENCHMARK_ROOT / "inputs/corpora/a10_primality.json"
REPORTED_PRIME = 23220101083624511828731
REPORTED_INPUT = 38389398379837983789739873873
LIMITS = (
    9080191,
    4759123141,
    2**64,
    318665857834031151167461,
    3317044064679887385961981,
)


def verify(data):
    """Check finite proof nodes without Factor's primality implementation."""
    certificates = data["certificates"]
    verified = set()

    def prove(n, active):
        if n in verified:
            return
        if n in active or n < 2:
            raise ValueError("invalid oracle proof graph")
        node = certificates[str(n)]
        if node["kind"] == "trial":
            if n > 2**32 or any(n % d == 0 for d in range(2, isqrt(n) + 1)):
                raise ValueError("invalid trial proof")
        elif node["kind"] == "pocklington":
            q, a = node["q"], node["witness"]
            prove(q, active | {n})
            if (
                (n - 1) % q
                or q * q <= n
                or pow(a, n - 1, n) != 1
                or gcd(pow(a, (n - 1) // q, n) - 1, n) != 1
            ):
                raise ValueError("invalid Pocklington proof")
        elif node["kind"] == "lucas":
            factors, a = node["factors"], node["witness"]
            if prod(q**e for q, e in factors) != n - 1:
                raise ValueError("incomplete n-1 factorization")
            if pow(a, n - 1, n) != 1:
                raise ValueError("invalid Lucas Fermat congruence")
            for q, e in factors:
                if e < 1:
                    raise ValueError("invalid oracle exponent")
                prove(q, active | {n})
                if gcd(pow(a, (n - 1) // q, n) - 1, n) != 1:
                    raise ValueError("invalid Lucas order congruence")
        else:
            raise ValueError("unknown oracle proof")
        verified.add(n)

    for value in certificates:
        prove(int(value), set())
    for fixture in data["fixtures"]:
        n = fixture["n"]
        if fixture["prime"]:
            prove(n, set())
        elif n >= 2 and not (
            1 < fixture["divisor"] < n and n % fixture["divisor"] == 0
        ):
            raise ValueError("invalid composite oracle")
    for fixture in data["factoring"]:
        if prod(p**e for p, e in fixture["factors"]) != fixture["n"]:
            raise ValueError("invalid factoring oracle")
        for p, _ in fixture["factors"]:
            prove(p, set())


def load_corpus():
    data = json.loads(source_path(CORPUS).read_text())
    verify(data)
    return data


def build():
    rng = random.Random(20261009)
    certs, fixtures, factoring = {}, [], []
    for p in (2, 3, 5, 61, 27103, 562711, 581557, 2365185233):
        certs[str(p)] = {"kind": "trial"}
    certs[str(REPORTED_PRIME)] = {
        "kind": "lucas",
        "witness": 2,
        "factors": [[p, 1] for p in (2, 3, 5, 562711, 581557, 2365185233)],
    }

    def add(n, prime, cohort, split, divisor=None):
        fixtures.append(
            dict(n=n, prime=prime, cohort=cohort, split=split, divisor=divisor)
        )

    add(REPORTED_PRIME, True, "wide_primes", "confirmation")
    factoring.append(
        dict(
            id="reported",
            split="regression",
            n=REPORTED_INPUT,
            factors=[[61, 1], [27103, 1], [REPORTED_PRIME, 1]],
        )
    )
    for _ in range(128):
        upper_prime = certified_prime(82, rng, certs)
        if LIMITS[-2] <= upper_prime < LIMITS[-1]:
            add(upper_prime, True, "wide_primes", "confirmation")
            factoring.append(
                dict(
                    id="upper_13_bases",
                    split="confirmation",
                    n=upper_prime,
                    factors=[[upper_prime, 1]],
                )
            )
            break
    else:
        raise ValueError("failed to generate the bounded upper-range oracle")
    for bits in (20, 24, 32, 48, 63, 64, 70, 78, 81, 82, 96, 127):
        for index in range(4):
            p = certified_prime(bits, rng, certs)
            split = "training" if index < 2 else "confirmation"
            cohort = "wide_primes" if 64 < bits <= 81 else "other_primes"
            add(p, True, cohort, split)
            q = certified_prime(min(bits, 17), rng, certs)
            add(p * q, False, "rough_composites", split, q)
            add(p * 43, False, "filtered_composites", split, 43)
            add(p * p, False, "powers", split, p)
            for divisor in (2, 3, 37):
                add(p * divisor, False, "small_divisibles", split, divisor)
            if index in (0, 2) and bits <= 81:
                factoring.append(
                    dict(
                        id=f"prime_{bits}_{index}",
                        split=split,
                        n=p,
                        factors=[[p, 1]],
                    )
                )
                if bits in (20, 70, 81):
                    factoring.append(
                        dict(
                            id=f"power_{bits}_{index}",
                            split=split,
                            n=p * p,
                            factors=[[p, 2]],
                        )
                    )
    for n, d in (
        (9080191, 2131),
        (4759123141, 48781),
        (341550071728321, 10670053),
        (3825123056546413051, 149491),
        (318665857834031151167461, 399165290221),
        (3317044064679887385961981, 1287836182261),
        (561, 3),
        (1105, 5),
        (1729, 7),
        (2465, 5),
        (2821, 7),
        (6601, 7),
        (3215031751, 151),
    ):
        add(n, False, "pseudoprimes", "adversarial", d)
    for boundary in LIMITS:
        for n in (boundary - 1, boundary + 1):
            # The word boundary is even; its odd neighbors have known factors.
            d = 3 if n == 2**64 - 1 else 274177 if n == 2**64 + 1 else 2
            add(n, False, "boundaries", "adversarial", d)
    add(2**64, False, "boundaries", "adversarial", 2)
    for n in (-1, 0, 1):
        add(n, False, "trivial", "adversarial")
    data = dict(
        schema=1,
        seed=20261009,
        seeds=[7, 104729, 130363],
        work_limit=1000000,
        certificates=certs,
        fixtures=fixtures,
        factoring=factoring,
    )
    verify(data)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CORPUS)
    args = parser.parse_args()
    args.output.write_text(json.dumps(build(), indent=2) + "\n")


if __name__ == "__main__":
    main()
