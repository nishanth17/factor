"""Freeze generated inputs with independent recursive Pocklington proofs."""

import argparse
import json
import random
from math import gcd, isqrt
from pathlib import Path

from ...common import utils

GENERATION_SEED = 20261003


def certified_prime(bits, generator, certificates):
    """Generate p=k*q+1 with a certified q > sqrt(p), then prove p.

    The Factor classifier is only a candidate filter. The certificate proves
    primality independently via Pocklington's theorem; small leaves use exact
    trial division. Known factors/certificates are never given to algorithms.
    """
    if bits <= 16:
        while True:
            candidate = generator.randrange(2 ** (bits - 1), 2**bits) | 1
            if all(candidate % d for d in range(2, isqrt(candidate) + 1)):
                certificates[str(candidate)] = {"kind": "trial"}
                return candidate

    q = certified_prime(bits // 2 + 2, generator, certificates)
    lower = ((2 ** (bits - 1) - 1) // q + 2) // 2
    upper = (2**bits - 2) // (2 * q)

    while True:
        k = 2 * generator.randint(lower, upper)
        candidate = k * q + 1
        if not utils.is_prime(candidate, rng=generator):
            continue
        for witness in range(2, 100):
            if (
                pow(witness, candidate - 1, candidate) == 1
                and gcd(pow(witness, k, candidate) - 1, candidate) == 1
                and q * q > candidate
            ):
                certificates[str(candidate)] = {
                    "kind": "pocklington",
                    "q": q,
                    "witness": witness,
                }
                return candidate


def verify_certificates(certificates):
    """Verify all proof nodes without calling Factor's primality test."""
    verified = set()

    def verify(n):
        if n in verified:
            return
        proof = certificates[str(n)]
        if proof["kind"] == "trial":
            if (
                n < 2
                or n >= 2**16
                or any(n % divisor == 0 for divisor in range(2, isqrt(n) + 1))
            ):
                raise AssertionError("invalid trial certificate")
        else:
            q, witness = proof["q"], proof["witness"]
            if not 2 <= q < n or (n - 1) % q:
                raise AssertionError("invalid certificate dependency")
            verify(q)
            if not (
                q * q > n
                and pow(witness, n - 1, n) == 1
                and gcd(pow(witness, (n - 1) // q, n) - 1, n) == 1
            ):
                raise AssertionError("invalid Pocklington certificate")

        verified.add(n)

    for n in certificates:
        verify(int(n))


def build(count=40):
    """Create disjoint cohorts; large bands are capped exploration."""
    generator = random.Random(GENERATION_SEED)
    certificates = {}
    fixtures = []

    def add(band, factors, index):
        n = 1
        for prime, exponent in factors:
            n *= prime**exponent
        fixtures.append(
            {
                "id": f"{band}_{index:03d}",
                "band": band,
                "n": n,
                "factors": sorted(factors),
                "split": "training" if index < count // 2 else "held_out",
                "digits": len(str(n)),
            }
        )

    for digits in range(20, 81, 10):
        for index in range(count):
            # Choose exact decimal bands, not a float approximation to n.
            while True:
                bits = (digits * 3322 // 1000 + 1) // 2
                p = certified_prime(bits, generator, certificates)
                q = certified_prime(bits, generator, certificates)
                if p != q and len(str(p * q)) == digits:
                    break

            add(f"balanced_{digits}d", [(p, 1), (q, 1)], index)

        print("generated", digits, "digit balanced band", flush=True)

    for digits in (60, 100):
        for small_digits in (5, 10, 20, 30):
            for index in range(count):
                p = certified_prime(
                    small_digits * 3322 // 1000, generator, certificates
                )

                while True:
                    q = certified_prime(
                        (digits - small_digits) * 3322 // 1000 + 1,
                        generator,
                        certificates,
                    )
                    if len(str(p * q)) == digits:
                        break

                add(
                    f"unbalanced_{digits}d_small_{small_digits}d",
                    [(p, 1), (q, 1)],
                    index,
                )

    for index in range(count):
        p = certified_prime(24 + index % 8, generator, certificates)
        q = certified_prime(28 + index % 8, generator, certificates)
        add("powers", [(p, (2, 3, 5, 7)[index % 4])], index)
        add("random_small", [(p, 1), (q, 1)], index)
        large = certified_prime(128 + index % 16, generator, certificates)
        add("primes", [(large, 1)], index)
        # Exact small leaves make a close-factor/Fermat cohort independent.
        p = certified_prime(15, generator, certificates)
        q = p + 2
        while any(q % d == 0 for d in range(2, isqrt(q) + 1)):
            q += 2
        certificates[str(q)] = {"kind": "trial"}
        if q >= 2**16:
            continue
        add("close_small", [(p, 1), (q, 1)], index)
        # Stage-one, first/last stage-two primes, and mixed saturation cases.
        smooth = (13, 19, 31, 61, 211, 607)[index % 6]
        other = (1019, 1009, 1013)[index % 3]
        for value in (smooth, other):
            certificates[str(value)] = {"kind": "trial"}
        add("pm1_boundaries", [(smooth, 1), (other, 1)], index)
        plus = (11, 23, 59, 107, 179)[index % 5]
        certificates[str(plus)] = {"kind": "trial"}
        add("pp1_boundaries", [(plus, 1), (other, 1)], index)

    for index, factors in enumerate(
        (
            [(151, 1), (751, 1), (28351, 1)],
            [(10670053, 1), (32010157, 1)],
        )
    ):
        # Independent exact trial proofs suffice for these published fixtures.
        for prime, _ in factors:
            if prime >= 2**16:
                continue
            certificates[str(prime)] = {"kind": "trial"}
        n = 1
        for prime, exponent in factors:
            n *= prime**exponent
        fixtures.append(
            {
                "id": f"pseudoprime_{index}",
                "band": "pseudoprimes",
                "n": n,
                "factors": factors,
                "split": "held_out",
                "digits": len(str(n)),
            }
        )

    verify_certificates(certificates)
    return {
        "schema": 1,
        "generation_seed": GENERATION_SEED,
        "inputs_per_main_band": count,
        "certificates": certificates,
        "seeds": [104729, 130363, 155921, 181081, 206369],
        "fixtures": fixtures,
        "notes": [
            "Balanced 20-80 digit bands are capped exploration until measured",
            "No p+1 engine yet; p+1 fixtures measure portfolio coverage",
            "Verify larger pseudoprime factors by independent trial division",
            "Training and held-out fixture IDs are disjoint",
        ],
    }


def main():
    """Generate once; refuse to replace a frozen corpus accidentally."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--count", type=int, default=40)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output exists; select a new corpus filename")
    if args.count < 40 or args.count % 2:
        parser.error("use at least 40 inputs per band and equal splits")
    args.output.write_text(json.dumps(build(args.count), indent=2) + "\n")


if __name__ == "__main__":
    main()
