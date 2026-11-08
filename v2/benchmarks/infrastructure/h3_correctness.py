"""Fixed arithmetic prefixes for H3 correctness; no performance measurement."""

import argparse
import inspect
import json
from dataclasses import asdict
from math import gcd, prod
from pathlib import Path

from ...common.prime_sieve import prime_sieve
from ...ecm.core import EcmStats, factorize_ecm
from ...ecm.workers import arithmetic_options
from ...execution.budget import Budget
from ...execution.portfolio_config import PortfolioConfig
from ...execution.worker_contract import Manifest, WorkerSettings
from ...qs.factor_base import build_factor_base
from ...qs.polynomial import qs_polynomial
from ...qs.sieve_collector import SieveCollector, SieveConfig


def budget():
    """Use the same non-refusing work grant on the committed reference."""
    if "work_limit" in inspect.signature(Budget).parameters:
        return Budget(work_limit=10**13, seconds=None, cpu_seconds=None)
    return Budget(seconds=None, cpu_seconds=None)


def cases():
    directory = Path(__file__).resolve().parents[1] / "inputs" / "corpora"
    selected = []
    for filename, digits in (
        ("c9_fresh_corpus.json", 50),
        ("c10a04_training.json", 65),
    ):
        data = json.loads((directory / filename).read_text())
        fixtures = [row for row in data["fixtures"] if row["digits"] == digits]
        if len(fixtures) < 3:
            raise ValueError("three fixed fixtures are required per band")
        selected.extend((row["id"], row["n"]) for row in fixtures[:3])

    # Separate exact smooth witnesses ensure nonempty accepted relation
    # signatures even if the six representative short prefixes find none.
    radical = prod(prime_sieve(1000))
    for digits, start in ((50, 4 * 10**24), (65, 2 * 10**32)):
        for offset in range(5000):
            root = start + 2 * offset
            n = root * root - 1
            if gcd(n, radical) == 1:
                selected.append(("smooth-witness-" + str(digits), n))
                break
        else:
            raise ValueError("finite witness construction failed")
    return selected


def qs_prefix(identity, n):
    owner = budget()
    base = build_factor_base(n, bound=1000, budget=owner).factor_base
    if base is None:
        raise ValueError("fixture has an early base divisor")
    options = dict(block_width=17, residual_bound=500, score_policy="powers")
    if "memory_bytes" in inspect.signature(SieveConfig).parameters:
        options["memory_bytes"] = 256 * 2**20
    collector = SieveCollector(
        qs_polynomial(base), base, config=SieveConfig(**options), budget=owner
    )
    result = collector.collect(-128, 129)
    atoms = [
        (
            atom.relation_id,
            atom.polynomial.a,
            atom.polynomial.b,
            atom.position,
            atom.sign,
            tuple(atom.exponents),
            atom.residual,
            atom.large_primes,
        )
        for atom in collector._atoms.values()
    ]
    rows = [
        (
            getattr(row, "relation_id", None),
            getattr(row, "atom_ids", ()),
            tuple(row.exponents),
            getattr(row, "u", None),
            getattr(row, "square_correction", None),
        )
        for row in collector._rows
    ]
    return dict(
        id=identity,
        n=n,
        next_position=result.next_position,
        scanned=result.stats["scanned"],
        reason=result.reason,
        atoms=atoms,
        rows=rows,
    )


def ecm_cases():
    results = []
    for n in (1009 * 1013, 10007 * 10009, 1000003 * 1000033):
        for seed in (7, 7919, 104729):
            stats = EcmStats()
            divisor = factorize_ecm(
                n,
                seed=seed,
                max_curves=4,
                b1=30,
                b2=100,
                stats=stats,
                _known_composite=True,
            )
            if divisor is not None and (not 1 < divisor < n or n % divisor):
                raise ArithmeticError("invalid ECM divisor")
            results.append(
                dict(
                    n=n,
                    seed=seed,
                    curves=4,
                    b1=30,
                    b2=100,
                    divisor=divisor,
                    transitions=asdict(stats),
                )
            )
    return results


def manifest_cases():
    config = PortfolioConfig(
        trial_bound=2,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=((30, 100, 4),),
        ecm_chain_mode="off",
        max_input_bits=512,
        segment_size=16,
    )
    manifest = Manifest(
        1009 * 1013,
        7,
        "fixed-h3",
        arithmetic_options(config, WorkerSettings()),
    )
    return [manifest.assignment(position) for position in range(4)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = dict(
        qs=[qs_prefix(identity, n) for identity, n in cases()], ecm=ecm_cases()
    )
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")


if __name__ == "__main__":
    main()
