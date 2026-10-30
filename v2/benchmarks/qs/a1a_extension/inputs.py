"""Untouched certified populations for the two fixed integration bands."""

import json
from pathlib import Path

from ...infrastructure.performance.performance_audit import verify_corpus
from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from ..c10a04.accounting import receipt
from ..c10a05_plus.selection import GenerationRandom

ROOT = Path(__file__).resolve().parents[4]
INPUTS = ROOT / "v2/benchmarks/inputs"
CORPUS = INPUTS / "corpora/a1a_extension_fresh.json"
EXCLUSIONS = INPUTS / "controls/a1a_extension_exclusions.json"
SEED = 20261030025055
STRATA = tuple(
    (r, s) for s in ("balanced", "uneven") for r in (1, 3, 5, 7)
) + ((1, "balanced"),)


def check_warmup(corpus):
    """Check the revealed fixtures independently before calibration."""
    verify_certificates(corpus["certificates"])
    verify_corpus(corpus)
    if sorted(row["digits"] for row in corpus["fixtures"]) != [50, 55]:
        raise ValueError("warmup needs exactly one revealed input per band")


def check(corpus, binding):
    """Verify independent proofs, shapes and cross-study prime exclusions."""
    if corpus["binding"] != binding or corpus["seed"] != SEED:
        raise ValueError("fresh corpus binding changed")
    verify_certificates(corpus["certificates"])
    verify_corpus(corpus)
    exclusions = json.loads(EXCLUSIONS.read_text())
    used, products = set(exclusions["primes"]), set(exclusions["products"])
    rows = corpus["fixtures"]
    if len(rows) != 18:
        raise ValueError("expected nine independent inputs per band")
    for digits in (50, 55):
        group = [r for r in rows if r["digits"] == digits]
        if tuple((r["n"] % 8, r["kind"]) for r in group) != STRATA:
            raise ValueError("fresh strata changed")
        bits = digits * 3322 // 1000 + 1
        for row in group:
            p, q = row["factors"]
            small = bits // 2 if row["kind"] == "balanced" else bits * 2 // 5
            if (
                {p.bit_length(), q.bit_length()} != {small, bits - small}
                or p == q
                or len(str(row["n"])) != digits
                or {p, q} & used
                or row["n"] in products
                or row["residue_mod8"] != row["n"] % 8
            ):
                raise ValueError("invalid or reused fresh input")
            used.update((p, q))
            products.add(row["n"])


def generate(binding):
    """Generate once after the source/protocol freeze, before timing."""
    if CORPUS.exists():
        raise ValueError("fresh population already exists; never regenerate")
    exclusions = json.loads(EXCLUSIONS.read_text())
    used, products = set(exclusions["primes"]), set(exclusions["products"])
    generator = GenerationRandom(SEED)
    certificates, rows = {}, []
    for digits in (50, 55):
        bits = digits * 3322 // 1000 + 1
        for index, (residue, shape) in enumerate(STRATA):
            small = bits // 2 if shape == "balanced" else bits * 2 // 5
            for _ in range(10000):
                p = certified_prime(small, generator, certificates)
                q = certified_prime(bits - small, generator, certificates)
                n = p * q
                if len(certificates) > 8192:
                    raise RuntimeError("finite certificate-entry cap")
                if (
                    p != q
                    and not {p, q} & used
                    and n not in products
                    and len(str(n)) == digits
                    and n % 8 == residue
                ):
                    break
            else:
                raise RuntimeError("finite per-stratum generation cap")
            rows.append(
                dict(
                    id=f"a1a_extension_{digits}_{index}",
                    digits=digits,
                    n=n,
                    factors=sorted((p, q)),
                    kind=shape,
                    residue_mod8=residue,
                )
            )
            used.update((p, q))
            products.add(n)
    corpus = dict(
        schema=1,
        binding=binding,
        seed=SEED,
        fixtures=rows,
        certificates=certificates,
        sampling=(
            "Certified products conditioned on a large p-1 factor; four odd "
            "residues crossed with balanced/40:60 bit shapes, plus one "
            "balanced residue-1 input."
        ),
    )
    check(corpus, binding)
    receipt(CORPUS, corpus)
    return {"fixtures": len(rows), "draws": generator.calls, "work": 0}
