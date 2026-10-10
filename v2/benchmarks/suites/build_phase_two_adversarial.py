"""Add independently verified Carmichael cases without replacing the corpus."""

import argparse
import json
import random
from math import isqrt, prod
from pathlib import Path

from ..support.paths import (
    BENCHMARK_ROOT,
    source_path,
)
from .build_phase_two_corpus import GENERATION_SEED, verify_certificates


def supplement(source):
    """Append forty Chernick triples with shuffled train/test allocation."""
    corpus = json.loads(source_path(source).read_text())
    additions = []

    for parameter in range(1, 3641):
        primes = [6 * parameter + 1, 12 * parameter + 1, 18 * parameter + 1]
        if not all(
            all(prime % d for d in range(2, isqrt(prime) + 1))
            for prime in primes
        ):
            continue

        n = prod(primes)

        # Korselt's criterion: odd, square-free, and each p-1 divides n-1.
        # Also verify the actual base-2 pseudoprime relation independently.
        if (
            any((n - 1) % (prime - 1) for prime in primes)
            or pow(2, n - 1, n) != 1
        ):
            raise AssertionError("invalid Carmichael construction")

        for prime in primes:
            corpus["certificates"].setdefault(str(prime), {"kind": "trial"})
        additions.append(
            {
                "id": f"pseudoprime_carmichael_m{parameter}",
                "band": "pseudoprimes",
                "n": n,
                "factors": [[prime, 1] for prime in primes],
                "digits": len(str(n)),
                "construction": "Chernick",
            }
        )
        if len(additions) == 40:
            break

    if len(additions) != 40:
        raise AssertionError("insufficient independently verified fixtures")
    random.Random(GENERATION_SEED + 1).shuffle(additions)
    for index, fixture in enumerate(additions):
        fixture["split"] = "training" if index < 20 else "held_out"
    corpus["fixtures"].extend(additions)
    corpus["notes"].append(
        "Forty independently verified Carmichael cases supplement "
        "two MR regressions"
    )
    corpus["supplement_seed"] = GENERATION_SEED + 1
    verify_certificates(corpus["certificates"])
    return corpus


def main():
    """Save a new combined fixture file and preserve the original capture."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=(BENCHMARK_ROOT / "inputs/corpora/phase_two_corpus.json"),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output exists; select a new combined corpus filename")
    args.output.write_text(
        json.dumps(supplement(args.source), indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
