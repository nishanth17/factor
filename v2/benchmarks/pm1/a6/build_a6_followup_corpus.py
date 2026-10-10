"""Freeze fresh balanced inputs before A6 follow-up confirmation timing."""

import json
import random

from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from ...support.paths import (
    BENCHMARK_ROOT,
)

DESTINATION = BENCHMARK_ROOT / "inputs/corpora/a6_followup_corpus.json"


def build():
    generator = random.Random(2026100919)
    certificates, fixtures = {}, []
    for digits in (20, 50, 100):
        bits = digits * 3322 // 1000 // 2
        for _ in range(64):
            p = certified_prime(bits, generator, certificates)
            q = certified_prime(bits, generator, certificates)
            if p != q and len(str(p * q)) == digits:
                fixtures.append(
                    {
                        "id": f"a6_fresh_balanced_{digits}d",
                        "digits": digits,
                        "n": p * q,
                        "factors": [p, q],
                    }
                )
                break
        else:
            raise RuntimeError("finite fixture construction exhausted")
    verify_certificates(certificates)
    return {
        "schema": 1,
        "seed": 2026100919,
        "construction": (
            "Recursive Pocklington; biased p-1 structure, "
            "not an RSA population model. No engine or timing selects inputs."
        ),
        "certificates": certificates,
        "fixtures": fixtures,
    }


if __name__ == "__main__":
    DESTINATION.write_text(json.dumps(build(), indent=2) + "\n")
    print(DESTINATION)
