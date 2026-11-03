"""Finite, auditable corpus proposals; no generation during preparation."""

import hashlib
import json
import random
import time
from math import prod
from pathlib import Path

TOTAL_DIGITS = (60, 80, 100)
FACTOR_BANDS = ((10, 14), (15, 19), (20, 24), (25, 29), (30, 35))


def random_integers(seed, count_per_band):
    """Uniform exact-digit draws, including primes and even numbers."""
    if type(count_per_band) is not int or not 1 <= count_per_band <= 1000:
        raise ValueError("finite panel count required")
    generator = random.Random(seed)
    return [
        dict(
            id=f"uniform-{digits}-{index}",
            n=generator.randrange(10 ** (digits - 1), 10**digits),
            digits=digits,
            population="uniform",
        )
        for digits in TOTAL_DIGITS
        for index in range(count_per_band)
    ]


def enriched_proposals(seed, *, count=4, seconds=120, draws=100_000):
    """Future certified enriched panel, under a finite cumulative draw clock.

    These are conditioned composites, separate from primary random panels.
    Existing Pocklington generator and its certificate verifier are reused.
    """
    from ...suites.build_phase_two_corpus import certified_prime

    if (
        not 1 <= count <= 100
        or not 0 < seconds <= 120
        or not 1 <= draws <= 100_000
    ):
        raise ValueError("generation request exceeds bootstrap reservation")
    started = time.monotonic()

    class LimitedRandom(random.Random):
        def randrange(self, *args):
            self.draws += 1
            if self.draws > draws or time.monotonic() - started >= seconds:
                raise TimeoutError("generation allowance exhausted")
            return super().randrange(*args)

    generator = LimitedRandom(seed)
    generator.draws = 0
    certificates, fixtures = {}, []
    for digits in TOTAL_DIGITS:
        for low, high in FACTOR_BANDS:
            # A smaller factor cannot exceed half the total logarithmic size.
            eligible = [d for d in range(low, high + 1) if 2 * d <= digits]
            if not eligible:
                continue
            for index in range(count):
                small = generator.choice(eligible)
                while True:
                    p = certified_prime(
                        (10**small - 1).bit_length(), generator, certificates
                    )
                    q = certified_prime(
                        (10 ** (digits - small) - 1).bit_length(),
                        generator,
                        certificates,
                    )
                    if len(str(p)) != small or len(str(q)) != digits - small:
                        continue
                    n = p * q
                    if len(str(n)) == digits and p <= q:
                        break
                fixtures.append(
                    dict(
                        id=f"ecm-{digits}-{low}-{index}",
                        n=n,
                        factors=[[p, 1], [q, 1]],
                        digits=digits,
                        smaller_digits=len(str(p)),
                        population="enriched",
                    )
                )
    return dict(
        fixtures=fixtures,
        certificates=certificates,
        seed=seed,
        draws=generator.draws,
    )


def validate_fixture(fixture):
    """Use certificates/factors for independent validation, never dispatch."""
    factors = fixture.get("factors")
    if factors is not None:
        pairs = (
            factors
            if isinstance(factors[0], list)
            else [[p, 1] for p in factors]
        )
        if prod(p**e for p, e in pairs) != fixture["n"]:
            raise ValueError("fixture factorization does not reconstruct")
    if len(str(fixture["n"])) != fixture.get("digits", len(str(fixture["n"]))):
        raise ValueError("fixture size label differs")


def retained(root, name, fixture_id, *, split="training"):
    """Require unambiguous exposure/split selection from committed bytes."""
    path = Path(root) / "v2/benchmarks/inputs/corpora" / name
    raw = path.read_bytes()
    data = json.loads(raw)
    matches = [
        f
        for f in data["fixtures"]
        if f["id"] == fixture_id and f.get("split", "training") == split
    ]
    if len(matches) != 1:
        raise ValueError("ambiguous or missing retained fixture")
    fixture = matches[0]
    validate_fixture(fixture)
    return fixture, hashlib.sha256(raw).hexdigest()
