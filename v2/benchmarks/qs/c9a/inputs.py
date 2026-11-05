"""Bounded exposed sampling and independent recursive prime certificates."""

import random
import time
from math import gcd

from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from .configuration import SEEDS, SIZES

RESIDUES = dict(
    zip(
        SIZES, ((1, 3), (5, 7), (3, 1), (7, 5), (1, 3), (5, 7), (3, 1), (7, 5))
    )
)


class GenerationRandom(random.Random):
    """Enforce finite draw, wall and CPU limits even inside primality loops."""

    def __init__(self, seed, seconds=590):
        super().__init__(seed)
        self.draws = 0
        self.seconds = seconds
        self.wall = time.monotonic()
        self.cpu = time.process_time()

    def randrange(self, *args, **kwargs):
        self.draws += 1
        if (
            self.draws > 1000000
            or time.monotonic() - self.wall >= self.seconds
            or time.process_time() - self.cpu >= self.seconds
        ):
            raise RuntimeError("finite generation allowance exhausted")
        return super().randrange(*args, **kwargs)


def specification(band=None):
    """Confirmation is generated only after its one band is selected."""
    if band is not None:
        if band not in (60, 62, 65, 70):
            raise ValueError("unfrozen confirmation band")
        return [
            (band, i >= 9, (1, 3, 5, 7)[i % 4], f"fresh_{i:02}")
            for i in range(12)
        ]
    cells = [
        (d, bool(i), r, f"exposed_{d}_{i}")
        for d, residues in RESIDUES.items()
        for i, r in enumerate(residues)
    ]
    cells += [
        (d, False, r, f"regression_{d}_{i}")
        for d in (50, 55)
        for i, r in enumerate((1, 3, 5))
    ]
    cells += [(d, False, 7, f"routing_{d}") for d in (40, 51, 52, 56, 59, 61)]
    return cells


def generate(exclusions, *, band=None, seconds=590):
    """No factoring outcomes condition input selection; no smooth p−1 sieve."""
    generator = GenerationRandom(
        202611030901 if band is None else 202611030902 + band, seconds
    )
    certificates, fixtures, used = {}, [], set()
    previous = tuple(int(n) for n in exclusions)
    for digits, uneven, residue, identity in specification(band):
        total_bits = (digits - 1) * 3322 // 1000 + 3
        pbits = total_bits * 35 // 100 if uneven else total_bits // 2
        qbits = total_bits - pbits
        for attempt in range(10000):
            p = certified_prime(pbits, generator, certificates)
            q = certified_prime(qbits, generator, certificates)
            n = p * q
            if (
                p != q
                and p not in used
                and q not in used
                and len(str(n)) == digits
                and n % 8 == residue
                and all(gcd(n, old) == 1 for old in previous)
            ):
                break
        else:
            raise RuntimeError("finite input-cell attempts exhausted")
        used.update((p, q))
        fixtures.append(
            dict(
                id=identity,
                digits=digits,
                uneven=uneven,
                n=n,
                factors=[p, q],
                factor_bits=[pbits, qbits],
                n_mod8=residue,
                attempts=attempt + 1,
            )
        )
    corpus = dict(
        schema=1,
        fixtures=fixtures,
        certificates=certificates,
        seeds=list(SEEDS),
        draws=generator.draws,
        conditioning=(
            "Recursive Pocklington p=kq+1 with q>sqrt(p), decimal "
            "size, bit split and residue. Not uniform RSA primes; "
            "p-1 has a large certified factor, not forced smooth."
        ),
        overlap_checked_integers=len(previous),
        work=generator.draws * 100000,
    )
    validate(corpus, band=band)
    return corpus


def validate(corpus, *, band=None):
    verify_certificates(corpus["certificates"])
    expected = specification(band)
    if len(corpus["fixtures"]) != len(expected):
        raise ValueError("corpus cell count changed")
    seen = set()
    for row, (digits, uneven, residue, identity) in zip(
        corpus["fixtures"], expected
    ):
        p, q = row["factors"]
        ratio = min(p.bit_length(), q.bit_length()) / (
            p.bit_length() + q.bit_length()
        )
        if not (
            row["id"] == identity
            and p * q == row["n"]
            and p != q
            and p not in seen
            and q not in seen
            and len(str(row["n"])) == digits
            and row["n"] % 8 == residue
            and row["digits"] == digits
            and row["uneven"] == uneven
            and row["factor_bits"] == [p.bit_length(), q.bit_length()]
            and all(str(v) in corpus["certificates"] for v in (p, q))
            and abs(ratio - (0.35 if uneven else 0.5)) < 0.01
        ):
            raise ValueError("corpus arithmetic, class or independence failed")
        seen.update((p, q))
