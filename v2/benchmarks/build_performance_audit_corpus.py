"""Fresh independent certificates for the P2/P3.6.1 performance audit."""

import argparse
import json
import random
from math import prod
from pathlib import Path

from .build_phase_two_corpus import certified_prime, verify_certificates


def build(seed=361042026, selected_split="all", p2_shapes=False):
    """Freeze training and confirmation inputs before parameter selection."""
    certificates, fixtures, seen = {}, [], set()
    for split, offset, count in (("training", 0, 3), ("held_out", 1, 4)):
        if selected_split not in ("all", split):
            continue
        generator = random.Random(seed + offset)
        for band, bits, digits in (
            ("small", 13, None),
            ("medium", 21, None),
            ("30d", 50, 30),
            ("40d", 67, 40),
            ("60d", 100, 60),
            ("80d", 133, 80),
        ):
            for index in range(count if digits in (None, 30) else 1):
                while True:
                    p = certified_prime(bits, generator, certificates)
                    q = certified_prime(bits, generator, certificates)
                    n = p * q
                    if (
                        p != q
                        and n not in seen
                        and (digits is None or len(str(n)) == digits)
                    ):
                        break
                seen.add(n)
                fixtures.append(
                    dict(
                        id=f"p361_{band}_{split}_{index}",
                        band=band,
                        split=split,
                        n=n,
                        factors=sorted((p, q)),
                    )
                )
        if p2_shapes:
            # Deliberate regression classes, separate from balanced search.
            # 65521-1 = 2**4 * 3**2 * 5 * 7 * 13; trial proof is independent.
            certificates["65521"] = {"kind": "trial"}
            for index in range(count):
                small = certified_prime(21, generator, certificates)
                large = certified_prime(100, generator, certificates)
                third = certified_prime(23, generator, certificates)
                for shape, factors in (
                    ("power", [small] * (2, 3, 5, 7)[index]),
                    ("uneven", [small, large]),
                    ("many", [small, large, third]),
                    ("prime", [large]),
                    ("smooth", [65521, large]),
                ):
                    n = prod(factors)
                    if n in seen:
                        raise AssertionError("repeated constructed input")
                    seen.add(n)
                    fixtures.append(
                        dict(
                            id=f"p361_p2_{shape}_{split}_{index}",
                            band="p2_" + shape,
                            split=split,
                            n=n,
                            factors=sorted(factors),
                        )
                    )
    verify_certificates(certificates)
    return dict(
        schema=1,
        seed=seed,
        seeds=[7, 29],
        fixtures=fixtures,
        certificates=certificates,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=361042026)
    parser.add_argument("--p2-shapes", action="store_true")
    parser.add_argument(
        "--split", choices=("all", "training", "held_out"), default="all"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("preserve existing inputs; choose a new path")
    args.output.write_text(
        json.dumps(build(args.seed, args.split, args.p2_shapes), indent=2)
        + "\n"
    )


if __name__ == "__main__":
    main()
