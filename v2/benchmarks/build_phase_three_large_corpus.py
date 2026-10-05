"""Fresh independently certified balanced and adversarial P3.4 inputs."""

import argparse
import json
import random
from math import gcd, isqrt, prod
from pathlib import Path

from .. import utils
from .build_phase_two_corpus import certified_prime, verify_certificates

SEED = 31042026


def build():
    certificates, fixtures, seen = {}, [], set()
    generator = random.Random(SEED)

    def add(kind, digits, split, index, factors, **metadata):
        n = prod(factors)

        if n in seen:
            raise AssertionError("duplicate independent fixture")
        seen.add(n)
        fixtures.append(
            dict(
                id=f"{kind}_{digits}d_{split}_{index}",
                kind=kind,
                digits=len(str(n)),
                band=f"{kind}_{digits}d",
                split=split,
                n=n,
                factors=sorted(factors),
                residue8=n % 8,
                factor_digits=sorted(len(str(p)) for p in factors),
                **metadata,
            )
        )

    def pair(digits, smaller_digits=None, residue=None):
        while True:
            small_bits = (
                (digits * 3322 // 1000 + 1) // 2
                if smaller_digits is None
                else smaller_digits * 3322 // 1000 + 1
            )
            large_bits = (
                small_bits
                if smaller_digits is None
                else (digits - smaller_digits) * 3322 // 1000 + 1
            )
            p = certified_prime(small_bits, generator, certificates)
            q = certified_prime(large_bits, generator, certificates)
            n = p * q
            if smaller_digits is not None and (
                len(str(p)) != smaller_digits
                or len(str(q)) != digits - smaller_digits
            ):
                continue

            if (
                p == q
                or len(str(n)) != digits
                or (residue is not None and n % 8 != residue)
                or n in seen
            ):
                continue

            if smaller_digits is None and 100 * abs(p - q) <= min(p, q):
                continue
            return p, q

    for digits in (30, 40, 50, 60, 70, 80):
        for split, residues in (
            ("training", (1, 5)),
            ("held_out", (1, 3, 5, 7)),
        ):
            for index, residue in enumerate(residues):
                add(
                    "balanced",
                    digits,
                    split,
                    index,
                    pair(digits, residue=residue),
                    performance_representative=split == "held_out"
                    and index == 0,
                )

        print("certified balanced", digits, flush=True)

    for digits in (60, 80):
        for smaller in (5, 10, 20):
            for split, count in (("training", 1), ("held_out", 2)):
                for index in range(count):
                    add(
                        f"unbalanced_small{smaller}",
                        digits,
                        split,
                        index,
                        pair(digits, smaller_digits=smaller),
                    )

    for digits in (50, 60, 80):
        bits = (digits * 3322 // 1000 + 1) // 2

        for index in range(2):
            while True:
                p = certified_prime(bits, generator, certificates)
                proof = certificates[str(p)]
                common = proof["q"]
                found = None

                for step in range(1, 257):
                    q = p + 2 * common * step
                    if (
                        len(str(p * q)) != digits
                        or common * common <= q
                        or not utils.is_prime(q, rng=generator)
                    ):
                        continue

                    for witness in range(2, 100):
                        if (
                            pow(witness, q - 1, q) == 1
                            and gcd(pow(witness, (q - 1) // common, q) - 1, q)
                            == 1
                        ):
                            certificates[str(q)] = {
                                "kind": "pocklington",
                                "q": common,
                                "witness": witness,
                            }
                            start = isqrt(p * q)
                            start += start * start < p * q
                            distance = (p + q) // 2 - start
                            if distance <= 100000:
                                found = (q, distance)
                            break

                    if found:
                        break

                if found and p * found[0] not in seen:
                    break

            add(
                "close",
                digits,
                "held_out",
                index,
                (p, found[0]),
                fermat_iterations=found[1],
            )
            while True:
                p = certified_prime(bits, generator, certificates)
                if len(str(p * p)) == digits and p * p not in seen:
                    break
            add("square", digits, "held_out", index, (p, p))
            while True:
                p = certified_prime(
                    digits * 3322 // 1000 + 1, generator, certificates
                )
                if len(str(p)) == digits:
                    break

            add("prime", digits, "held_out", index, (p,))

    for index in range(4):
        factors = tuple(
            certified_prime(55 + index % 2, generator, certificates)
            for _ in range(3)
        )
        add("three_prime", len(str(prod(factors))), "held_out", index, factors)

    for kind, p in (("smooth_pm1", 65521), ("smooth_pp1", 65519)):
        certificates[str(p)] = {"kind": "trial"}
        for index in range(2):
            _, q = pair(60, smaller_digits=5)
            add(kind, len(str(p * q)), "held_out", index, (p, q))

    verify_certificates(certificates)
    return dict(
        schema=1,
        seed=SEED,
        seeds=[7, 29],
        fixtures=fixtures,
        certificates=certificates,
        notes=[
            (
                "All prime factors have independently checked "
                "Pocklington or exact trial proofs."
            ),
            (
                "Balanced factors are comparable in size and "
                "exclude extremely close pairs."
            ),
            (
                "Close pairs are labelled separately; their "
                "proof-generating structure is disclosed."
            ),
            "Representatives for controlled probes are preselected.",
            (
                "Factors/certificates/shape metadata are withheld "
                "from factoring methods."
            ),
        ],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("preserve existing corpus bytes")
    args.output.write_text(json.dumps(build(), indent=2) + "\n")


if __name__ == "__main__":
    main()
