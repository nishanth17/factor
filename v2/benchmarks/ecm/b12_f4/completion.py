"""Explicitly funded balanced completion configs and exploratory comparisons.

Run --prepare in each isolated source checkout before freezing the study.
This changes resource/fallback configuration, never the incumbent source.
"""

import argparse
import json
import math
import sys
from dataclasses import asdict, replace


def configuration(n, source_root, *, funded=True):
    """Retain each source's ECM defaults and use a common funded QS child."""
    sys.path.insert(0, str(source_root))
    from v2.portfolio import PortfolioConfig

    digits = len(str(n))
    if digits not in (60, 70):
        raise ValueError("balanced completion preparation supports 60/70 only")
    config = PortfolioConfig().for_input(n)
    if not funded:
        return config
    if digits == 60:
        if config.siqs is None:
            raise ValueError(
                "the incumbent calibrated 60-digit fallback is missing"
            )
        # Select the calibrated numerical policies before increasing capacity:
        # rebuilding defaults under a different grant can change solve_policy.
        fallback = replace(config.siqs)
        2 * 2**30 + 16 * 2**20
    else:
        from v2.benchmarks.qs.c10a01.configuration_v4 import (
            configuration as qs_configuration,
        )

        # This is an explicit common exploratory fallback, not stock 70-digit
        # automatic routing. These are the retained C10a completion policies.
        fallback = replace(
            qs_configuration(70, "coupled", multiplier_policy="parity2_wide")
        )
        4 * 2**30
    return replace(config, siqs=fallback)


def validate_known_factorization(row, factors):
    """Validate with certified factors held outside the factoring process."""
    from .selection_child import validate_reconstruction

    validate_reconstruction(row)
    expected = dict(factors)
    if math.prod(p**e for p, e in expected.items()) != row["n"]:
        raise ValueError("oracle belongs to a different input")
    observed = {}
    for prime, exponent, _ in row["factors"]:
        observed[prime] = observed.get(prime, 0) + exponent
    if any(p not in expected or e > expected[p] for p, e in observed.items()):
        raise ValueError("output contradicts the independent factor oracle")
    if row["factorization_complete"] and observed != expected:
        raise ValueError("complete output contradicts the factor oracle")


def point_comparison(control, candidate):
    """A censored outcome never produces a completed-call speedup ratio."""
    rows = (control, candidate)
    for row in rows:
        if type(row["complete"]) is not bool:
            raise ValueError("invalid completion indicator")
        if not math.isfinite(row["seconds"]) or row["seconds"] <= 0:
            raise ValueError("invalid completion time")
    both = all(row["complete"] for row in rows)
    return dict(
        both_complete=both,
        speedup=control["seconds"] / candidate["seconds"] if both else None,
        candidate_fraction=candidate["seconds"] / control["seconds"]
        if both
        else None,
        inference="exposed single-input point; no population claim",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument(
        "--resources", choices=("funded", "stock"), default="funded"
    )
    arguments = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 is required")
    print(
        json.dumps(
            asdict(
                configuration(
                    arguments.n,
                    arguments.source_root,
                    funded=arguments.resources == "funded",
                )
            )
        )
    )


if __name__ == "__main__":
    main()
