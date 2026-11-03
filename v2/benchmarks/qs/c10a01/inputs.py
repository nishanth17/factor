"""One-shot exposed corpus with independent recursive prime certificates."""

import json
import random
import time
from pathlib import Path

from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)

CORPUS = (
    Path(__file__).resolve().parents[2] / "inputs/corpora/c10a01_exposed.json"
)
RESIDUES = {60: (1, 3), 62: (5, 7), 65: (3, 1), 70: (7, 5)}


class GenerationRandom(random.Random):
    """Bound candidate and classifier draws in the reused sampler."""

    def __init__(self):
        super().__init__(202611020100)
        self.draws = 0
        self.wall = time.monotonic()
        self.cpu = time.process_time()

    def randrange(self, *args, **kwargs):
        self.draws += 1
        if (
            self.draws > 1000000
            or time.monotonic() - self.wall >= 590
            or time.process_time() - self.cpu >= 590
        ):
            raise RuntimeError("generation allowance exhausted")
        return super().randrange(*args, **kwargs)


def validate(corpus):
    """Check sizes, splits, unique factors, residues and prime proofs."""
    verify_certificates(corpus["certificates"])
    seen = set()
    assert len(corpus["fixtures"]) == 8
    for row in corpus["fixtures"]:
        p, q = row["factors"]
        assert p != q and p not in seen and q not in seen
        seen.update((p, q))
        assert p * q == row["n"] and len(str(row["n"])) == row["digits"]
        assert row["n"] % 8 == RESIDUES[row["digits"]][row["balance_index"]]
        assert row["factor_bits"] == [p.bit_length(), q.bit_length()]
        assert all(str(v) in corpus["certificates"] for v in (p, q))
        ratio = p.bit_length() / (p.bit_length() + q.bit_length())
        assert abs(ratio - (0.5 if row["balance_index"] == 0 else 0.35)) < 0.01


def generate():
    """Refuse replacement; never condition an input on factoring outcomes."""
    if CORPUS.exists():
        raise ValueError("exposed corpus already exists")
    generator, certificates, rows, used = GenerationRandom(), {}, [], set()
    for digits, residues in RESIDUES.items():
        total_bits = (digits - 1) * 3322 // 1000 + 3
        for index, residue in enumerate(residues):
            pbits = total_bits // 2 if index == 0 else total_bits * 35 // 100
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
                ):
                    break
            else:
                raise RuntimeError("finite input-cell attempts exhausted")
            used.update((p, q))
            balance = "balanced" if index == 0 else "35_65"
            rows.append(
                dict(
                    id=f"c10a01_{digits}_{balance}",
                    digits=digits,
                    balance_index=index,
                    n=n,
                    factors=[p, q],
                    factor_bits=[pbits, qbits],
                    n_mod8=residue,
                    generation_attempts=attempt + 1,
                )
            )
    corpus = dict(
        schema=1,
        seed=202611020100,
        fixtures=rows,
        certificates=certificates,
        draws=generator.draws,
        conditioning="Pocklington shape, decimal size, residue and split; "
        "exposed exploratory inputs, not uniform RSA or fresh confirmation",
    )
    validate(corpus)
    CORPUS.write_text(json.dumps(corpus, indent=2, sort_keys=True) + "\n")
    return dict(
        work=generator.draws * 100000,
        draws=generator.draws,
        fixtures=len(rows),
        wall=time.monotonic() - generator.wall,
        cpu=time.process_time() - generator.cpu,
    )
