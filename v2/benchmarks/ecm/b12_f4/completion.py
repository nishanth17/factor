"""Explicitly funded balanced completion configs and exploratory comparisons.

Run --prepare in each isolated source checkout before freezing the study.
This changes resource/fallback configuration, never the incumbent source.
"""

import argparse
import json
import math
import sys
from dataclasses import asdict, replace


def configuration(n, source_root):
    """Retain each source's ECM defaults and use a common funded QS child."""
    sys.path.insert(0, str(source_root))
    from v2.portfolio import PortfolioConfig

    digits = len(str(n))
    if digits not in (60, 70):
        raise ValueError("balanced completion preparation supports 60/70 only")
    config = PortfolioConfig().for_input(n)
    if digits == 60:
        if config.siqs is None:
            raise ValueError("the incumbent calibrated 60-digit fallback is missing")
        # Select the calibrated numerical policies before increasing capacity:
        # rebuilding defaults under a different grant can change solve_policy.
        fallback = replace(
            config.siqs, memory_bytes=2 * 2**30, checkpoint_bytes=64 * 2**20
        )
        parent_memory = 2 * 2**30 + 16 * 2**20
    else:
        from v2.benchmarks.qs.c10a01.configuration_v4 import (
            configuration as qs_configuration,
        )

        # This is an explicit common exploratory fallback, not stock 70-digit
        # automatic routing. These are the retained C10a completion policies.
        fallback = replace(
            qs_configuration(70, "coupled", multiplier_policy="parity2_wide"),
            checkpoint_bytes=64 * 2**20,
        )
        parent_memory = 4 * 2**30
    return replace(config, memory_bytes=parent_memory, siqs=fallback)


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
        if both else None,
        inference="exposed single-input exploratory point; no population claim",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--n", type=int, required=True)
    arguments = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 is required")
    print(json.dumps(asdict(configuration(arguments.n, arguments.source_root))))


if __name__ == "__main__":
    main()
