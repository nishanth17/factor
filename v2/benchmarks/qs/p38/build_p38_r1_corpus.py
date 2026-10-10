"""Certified sub-100-digit training and post-freeze confirmation cohorts."""

import argparse
import hashlib
import json
import random
import sys
from math import gcd, isqrt, prod
from pathlib import Path

from ....common import utils
from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from ...support.paths import (
    source_path,
)

BANDS = (30, 40, 60, 70, 80, 90, 99)


def build(seed, *, split="training", count=2, frozen_sha256=None):
    """Generate independent proofs; methods receive only n and their seed.

    Pocklington generation conditions primes on a large factor of p-1.
    This sampling bias is disclosed; no RSA distribution is claimed.
    """
    generator = random.Random(seed)
    certificates, fixtures, seen = {}, [], set()

    def prime(digits):
        bits = (10**digits - 1).bit_length()
        for _ in range(10000):
            value = certified_prime(bits, generator, certificates)
            if len(str(value)) == digits:
                return value
        raise RuntimeError("prime generation allowance exhausted")

    def pair(digits, smaller=None):
        small_digits = smaller or (digits + 1) // 2
        # Odd balanced bands need two ceil(d/2)-digit factors.
        large_digits = digits - smaller if smaller else (digits + 1) // 2

        for _ in range(100000):
            # Vary bit length to cover the low end of a decimal band.
            bits = (digits * 3322 // 1000 + 1) // 2
            p = (
                prime(small_digits)
                if smaller
                else certified_prime(bits, generator, certificates)
            )
            q = (
                prime(large_digits)
                if smaller
                else certified_prime(bits, generator, certificates)
            )
            if (
                p != q
                and len(str(p * q)) == digits
                and p * q not in seen
                and (smaller or 100 * abs(p - q) > min(p, q))
            ):
                return p, q

        raise RuntimeError("pair generation allowance exhausted")

    def add(kind, digits, index, factors, **metadata):
        n = prod(factors)

        if n in seen or len(str(n)) != digits:
            raise AssertionError("duplicate or incorrectly sized fixture")
        seen.add(n)
        fixtures.append(
            dict(
                id=f"{kind}_{digits}d_{split}_{index}",
                kind=kind,
                digits=digits,
                split=split,
                n=n,
                factors=sorted(factors),
                factor_digits=sorted(len(str(p)) for p in factors),
                smaller_factor_digits=min(len(str(p)) for p in factors),
                **metadata,
            )
        )

    for digits in BANDS:
        for index in range(count):
            add("balanced", digits, index, pair(digits))
        for smaller in (5, 10, 20, 30):
            if smaller >= digits // 2:
                continue
            for index in range(count):
                add(f"uneven_{smaller}", digits, index, pair(digits, smaller))

        for kind, p, smooth in (
            ("pm1_smooth", 65521, [[2, 4], [3, 2], [5, 1], [7, 1], [13, 1]]),
            ("pp1_smooth", 65519, [[2, 4], [3, 2], [5, 1], [7, 1], [13, 1]]),
        ):
            certificates[str(p)] = {"kind": "trial"}
            assert prod(p**e for p, e in smooth) == (65520)
            for index in range(count):
                for _ in range(10000):
                    q = prime(digits - 5)
                    if len(str(p * q)) == digits and p * q not in seen:
                        break
                else:
                    raise RuntimeError(
                        "smooth-control generation allowance exhausted"
                    )

                add(kind, digits, index, (p, q), smooth_neighbor=smooth)

        for index in range(count):
            bits = (digits * 3322 // 1000 + 1) // 2

            for _ in range(10000):
                p = certified_prime(bits, generator, certificates)
                if len(str(p * p)) == digits and p * p not in seen:
                    add("power", digits, index, (p, p))
                    break
            else:
                raise RuntimeError("power generation allowance exhausted")

            for _ in range(10000):
                p = certified_prime(bits, generator, certificates)
                common = certificates[str(p)]["q"]
                found = None

                for step in range(1, 257):
                    q = p + 2 * common * step
                    if (
                        common * common <= q
                        or len(str(p * q)) != digits
                        or p * q in seen
                        or not utils.is_prime(q, rng=generator)
                    ):
                        continue

                    for witness in range(2, 100):
                        if (
                            pow(witness, q - 1, q) == 1
                            and gcd(pow(witness, (q - 1) // common, q) - 1, q)
                            == 1
                        ):
                            certificates[str(q)] = dict(
                                kind="pocklington", q=common, witness=witness
                            )
                            found = q
                            break

                    if found:
                        break

                if found:
                    start = isqrt(p * found)
                    start += start * start < p * found
                    add(
                        "close",
                        digits,
                        index,
                        (p, found),
                        fermat_iterations=(p + found) // 2 - start,
                    )
                    break
            else:
                raise RuntimeError(
                    "close-control generation allowance exhausted"
                )

        print("certified", split, digits, flush=True)

    verify_certificates(certificates)
    return dict(
        schema=1,
        seed=seed,
        split=split,
        frozen_sha256=frozen_sha256,
        bands=BANDS,
        seeds=[7, 29],
        fixtures=fixtures,
        certificates=certificates,
        sampling="Pocklington primes have a large certified factor of p-1. "
        "Close pairs share that proof factor. Smooth controls have small "
        "proven factors. Populations are labelled separately. No known "
        "factors enter configuration selection.",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--count", type=int, default=2)
    parser.add_argument("--frozen", type=Path)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    if not 1 <= args.count <= 32:
        parser.error("count must be between 1 and 32")
    digest = (
        hashlib.sha256(source_path(args.frozen).read_bytes()).hexdigest()
        if args.frozen
        else None
    )
    data = build(
        args.seed,
        split="confirmation" if digest else "training",
        count=args.count,
        frozen_sha256=digest,
    )
    with args.output.open("x") as stream:
        json.dump(data, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
