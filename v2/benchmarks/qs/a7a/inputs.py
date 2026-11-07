"""Prospective disjoint confirmation inputs; no generation on import."""

import json
from math import gcd
from pathlib import Path

from .protocol import corpus


def generate_confirmation(output):
    """Call only in the later locked, reserved preparation lane (600 s)."""
    from ...suites.build_phase_two_corpus import (
        certified_prime,
        verify_certificates,
    )
    from ..c9a.inputs import GenerationRandom

    output = Path(output)
    if output.exists():
        raise ValueError("preserve prior confirmation generation")
    training = {r["n"] for r in corpus()["fixtures"]}
    rng = GenerationRandom(202611040799, seconds=590)
    certificates, fixtures, used = {}, [], set()
    for digits in (60, 70):
        bits = (digits - 1) * 3322 // 1000 + 3
        for index in range(9):
            for attempt in range(10000):
                p = certified_prime(bits // 2, rng, certificates)
                q = certified_prime(bits - bits // 2, rng, certificates)
                n = p * q
                if (
                    p != q
                    and p not in used
                    and q not in used
                    and len(str(n)) == digits
                    and all(gcd(n, old) == 1 for old in training)
                ):
                    break
            else:
                raise RuntimeError("finite confirmation generation exhausted")
            used.update((p, q))
            training.add(n)
            fixtures.append(
                dict(
                    id=f"fresh_{digits}_{index:02}",
                    n=n,
                    digits=digits,
                    factors=[p, q],
                    attempts=attempt + 1,
                )
            )
    verify_certificates(certificates)
    result = dict(
        schema=1,
        seed=202611040799,
        fixtures=fixtures,
        certificates=certificates,
        draws=rng.draws,
        conditioning=(
            "Exact decimal size and recursive Pocklington primes; "
            "no timing/outcome conditioning."
        ),
    )
    with output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    return result


def validate_committed():
    """Check source hashes, products and certificates during Stage Q."""
    import hashlib
    from math import prod

    from ...suites.build_phase_two_corpus import verify_certificates
    from .protocol import ROOT

    data = corpus()
    for name, expected in data["sources"].items():
        raw = (ROOT / name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != expected:
            raise ValueError("source corpus changed")
        verify_certificates(json.loads(raw)["certificates"])
    seen = set()
    for row in data["fixtures"]:
        if row["id"] in seen or len(str(row["n"])) != row["digits"]:
            raise ValueError("corpus identity/size changed")
        seen.add(row["id"])
        if "factors" in row and prod(row["factors"]) != row["n"]:
            raise ValueError("corpus product does not reconstruct")
