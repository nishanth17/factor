"""Generate and verify the disjoint G4 corpus; never tune on these inputs."""

import json
from math import gcd, prod

from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from ..c9a.inputs import GenerationRandom
from .combined import fixtures
from .protocol import INPUTS

BANDS = (40, 50, 55, 60, 65, 70)
SEED = 202611071442
PATH = INPUTS / "confirmation_g4.json"


def validate(data, previous):
    """Verify proofs, balanced sizes and disjointness."""
    verify_certificates(data["certificates"])
    expected = {f"g4_{d}_{i}" for d in BANDS for i in range(9)}
    if {row["id"] for row in data["fixtures"]} != expected:
        raise ValueError("confirmation coverage mismatch")
    if len(data["fixtures"]) != len(expected) or data["seed"] != SEED:
        raise ValueError("confirmation identity mismatch")
    seen = list(previous)
    for row in data["fixtures"]:
        n, (p, q) = row["n"], row["factors"]
        if (
            p == q
            or prod((p, q)) != n
            or len(str(n)) != row["digits"]
            or abs(p.bit_length() - q.bit_length()) > 1
            or any(str(x) not in data["certificates"] for x in (p, q))
            or any(gcd(n, old) != 1 for old in seen)
        ):
            raise ValueError("invalid or overlapping confirmation modulus")
        seen.append(n)


def generate():
    """Create the corpus exclusively, with bounded prime generation."""
    if PATH.exists():
        raise FileExistsError(PATH)
    previous = [r["n"] for r in fixtures("regimes").values() if r["n"] > 1]
    rng = GenerationRandom(SEED, seconds=120)
    proofs, rows, seen = {}, [], list(previous)
    for digits in BANDS:
        bits = (digits - 1) * 3322 // 1000 + 3
        for index in range(9):
            for _ in range(10000):
                p = certified_prime(bits // 2, rng, proofs)
                q = certified_prime(bits - bits // 2, rng, proofs)
                n = p * q
                if (
                    p != q
                    and len(str(n)) == digits
                    and all(gcd(n, old) == 1 for old in seen)
                ):
                    break
            else:
                raise RuntimeError("finite generation allowance exhausted")
            seen.append(n)
            rows.append(
                dict(
                    id=f"g4_{digits}_{index}",
                    digits=digits,
                    n=n,
                    factors=[p, q],
                    purpose="confirmation",
                )
            )
    data = dict(
        schema=1,
        seed=SEED,
        fixtures=rows,
        certificates=proofs,
        conditioning="Exact size and balanced certified primes only",
    )
    validate(data, previous)
    with PATH.open("x") as stream:
        json.dump(data, stream, indent=2)
        stream.write("\n")
    return len(rows)


if __name__ == "__main__":
    print(generate())
