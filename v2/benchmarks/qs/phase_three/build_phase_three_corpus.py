"""Independent Pocklington inputs, frozen before P3.4 parameter experiments."""

import argparse
import json
import random
from pathlib import Path

from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)

SEED = 30342026


def build():
    certificates, fixtures, seen = {}, [], set()

    for split, offset, count in (("training", 0, 12), ("held_out", 1, 32)):
        generator = random.Random(SEED + offset)

        for i in range(count):
            while True:
                p = certified_prime(
                    generator.choice((12, 13)), generator, certificates
                )
                q = certified_prime(13, generator, certificates)
                if p != q and p * q not in seen:
                    break

            seen.add(p * q)
            fixtures.append(
                dict(
                    id=f"small_{split}_{i}",
                    split=split,
                    band="small",
                    n=p * q,
                    factors=sorted((p, q)),
                    residue8=p * q % 8,
                )
            )

    for split, offset in (("training", 2), ("held_out", 3)):
        generator = random.Random(SEED + offset)

        for digits in (30, 40, 50, 60, 70, 80):
            residues = (
                (1,)
                if split == "training"
                else ((1, 3, 5, 7) if digits <= 60 else (1, 5))
            )

            for residue in residues:
                while True:
                    bits = (digits * 3322 // 1000 + 1) // 2
                    p = certified_prime(bits, generator, certificates)
                    q = certified_prime(bits, generator, certificates)
                    n = p * q
                    if (
                        p != q
                        and n not in seen
                        and len(str(n)) == digits
                        and n % 8 == residue
                    ):
                        break

                seen.add(n)
                fixtures.append(
                    dict(
                        id=f"balanced{digits}_{split}_r{residue}",
                        split=split,
                        band=f"balanced_{digits}d",
                        n=n,
                        factors=sorted((p, q)),
                        residue8=residue,
                    )
                )
                print("certified", split, digits, residue, flush=True)

    verify_certificates(certificates)
    return dict(
        schema=1,
        seed=SEED,
        fixtures=fixtures,
        certificates=certificates,
        seeds=[7, 29],
        notes=[
            "Pocklington proofs; exact trial leaves.",
            "No known factors or certificates are supplied to algorithms.",
            "Main bands cover four residues; 70/80 are capped exploration.",
        ],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("preserve existing inputs; choose a new path")
    args.output.write_text(json.dumps(build(), indent=2) + "\n")


if __name__ == "__main__":
    main()
