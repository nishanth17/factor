"""Build independent certified P3.6 training and held-out workloads."""

import argparse
import json
import random
import sys
from math import isqrt
from pathlib import Path

from .. import utils

OUTPUT = Path(__file__).parent / "inputs/phase_three_p36_corpus.json"


def build():
    generator = random.Random(36042026)
    fixtures = []
    for band, low, high in (
        ("small", 4000, 8000),
        ("medium", 1_000_000, 4_000_000),
    ):
        for split, count in (("training", 4), ("held_out", 4)):
            for index in range(count):
                primes = []
                while len(primes) < 2:
                    value = generator.randrange(low, high) | 1
                    while (
                        utils.classify_prime(value) != utils.Primality.PROVEN
                    ):
                        value += 2
                    if value not in primes:
                        primes.append(value)
                primes.sort()
                for prime in primes:
                    if any(
                        prime % divisor == 0
                        for divisor in range(2, isqrt(prime) + 1)
                    ):
                        raise AssertionError(
                            "independent prime certification failed"
                        )
                fixtures.append(
                    dict(
                        id=f"{band}_{split}_{index}",
                        band=band,
                        split=split,
                        n=primes[0] * primes[1],
                        factors=primes,
                    )
                )
    return dict(schema=1, seed=36042026, fixtures=fixtures)


def main():
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("use PyPy Python 3.11")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    args.output.write_text(json.dumps(build(), indent=2) + "\n")


if __name__ == "__main__":
    main()
