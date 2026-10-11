"""Disjoint, exactly labelled C3 training and later confirmation inputs."""

import json
import math
import random
from collections import Counter

from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from .c3_study import INPUTS, save, training

TRAINING = INPUTS / "corpora/c3_round2_training.json"


class GenerationRandom(random.Random):
    """Bound rejection sampling, including recursive proof generation."""

    def __init__(self, seed, draws=1_000_000):
        super().__init__(seed)
        self.remaining = draws

    def getrandbits(self, count):
        self.remaining -= 1
        if self.remaining < 0:
            raise RuntimeError("C3 round-two generation allowance exhausted")
        return super().getrandbits(count)


def generate_semiprimes(seed, *, balanced_count=4, uneven_count=3):
    """Cover actual decimal factor bands; proofs never enter dispatch."""
    generator = GenerationRandom(seed)
    certificates, fixtures = {}, []

    def prime(digits):
        lower, upper = 10 ** (digits - 1), 10**digits - 1
        for _ in range(1000):
            bits = generator.randint(lower.bit_length(), upper.bit_length())
            value = certified_prime(bits, generator, certificates)
            if lower <= value <= upper:
                return value
        raise RuntimeError("prime decimal-band generation exhausted")

    for digits in (30, 40):
        strata = [(digits // 2, balanced_count)]
        strata += [(small, uneven_count) for small in (8, 10, 12)]
        strata.append((14 if digits == 30 else 16, uneven_count))
        for small_digits, count in strata:
            for index in range(count):
                for _ in range(1000):
                    p, q = prime(small_digits), prime(digits - small_digits)
                    n = p * q
                    if p != q and len(str(n)) == digits:
                        break
                else:
                    raise RuntimeError(
                        "composite decimal-band generation exhausted"
                    )
                kind = (
                    "balanced"
                    if small_digits == digits // 2
                    else f"uneven_{small_digits}"
                )
                fixtures.append(
                    dict(
                        id=f"r2_{digits}_{kind}_{index}",
                        kind=kind,
                        n=n,
                        factors=sorted([(p, 1), (q, 1)]),
                        digits=digits,
                        small_digits=small_digits,
                    )
                )
    corpus = dict(
        seed=seed,
        fixtures=fixtures,
        certificates=certificates,
        bias="Nonuniform Pocklington p=kq+1; exact decimal rejection",
        generation_draws=1_000_000 - generator.remaining,
    )
    validate_corpus(corpus)
    return corpus


def validate_corpus(corpus):
    """Verify proofs, actual strata, uniqueness and exact reconstruction."""
    verify_certificates(corpus["certificates"])
    identities, numbers = set(), set()
    for fixture in corpus["fixtures"]:
        if fixture["id"] in identities or fixture["n"] in numbers:
            raise ValueError("duplicate C3 round-two fixture")
        identities.add(fixture["id"])
        numbers.add(fixture["n"])
        factors = fixture["factors"]
        if (
            any(str(p) not in corpus["certificates"] for p, _ in factors)
            or math.prod(p**e for p, e in factors) != fixture["n"]
            or len(str(fixture["n"])) != fixture["digits"]
        ):
            raise ValueError("uncertified or incorrectly labelled fixture")
        if "small_digits" in fixture and (
            len(str(min(p for p, _ in factors))) != fixture["small_digits"]
        ):
            raise ValueError("incorrect smaller-factor decimal stratum")


def load_training():
    """Previously revealed confirmation is regression/training evidence."""
    generated = json.loads(TRAINING.read_text())
    validate_corpus(generated)
    historical = json.loads(
        (INPUTS / "corpora/c3_confirmation.json").read_text()
    )
    validate_corpus(historical)
    fixtures = training() + historical["fixtures"] + generated["fixtures"]
    numbers = [fixture["n"] for fixture in fixtures]
    if len(numbers) != len(set(numbers)):
        raise ValueError("overlapping C3 training inputs")
    return fixtures


def write_training():
    corpus = generate_semiprimes(2026101107)
    old_numbers = {f["n"] for f in training()}
    old = json.loads((INPUTS / "corpora/c3_confirmation.json").read_text())
    old_numbers.update(f["n"] for f in old["fixtures"])
    if any(f["n"] in old_numbers for f in corpus["fixtures"]):
        raise ValueError("generated training overlaps historical inputs")
    save(TRAINING, corpus)
    return Counter((f["digits"], f["kind"]) for f in corpus["fixtures"])
