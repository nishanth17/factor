"""P3.4 family/root foundation costs and small verified QS interoperability."""

import argparse
import json
import math
import resource
import subprocess
import sys
import time
from functools import partial
from math import prod
from pathlib import Path

from .. import utils
from ..budget import Budget
from ..qs import (
    Polynomial,
    PolynomialFamily,
    QSJob,
    SieveCollector,
    SieveConfig,
    build_factor_base,
    family_assignments,
    polynomial_roots,
)
from .phase_one import environment
from .phase_three_pipeline import _corpus, _record
from .phase_three_reference import _rss_bytes

SEED = 2936
MEMORY_BYTES = 32 * 1024 * 1024


def _budget():
    """Declare one setup-to-extraction allowance per input in both arms."""
    return Budget(work_limit=200_000_000, seconds=10, cpu_seconds=10)


def full_steps(base, primes, budget):
    """Direct CRT sign sums and full inverse/root recomputation control."""
    a = prod(primes)
    entries = {entry.prime: entry for entry in base.entries}
    terms = []

    for prime in primes:
        budget.consume(a.bit_length() + prime.bit_length() ** 2)
        quotient = a // prime
        terms.append(
            quotient
            * utils.modular_inverse(quotient, prime)
            * entries[prime].square_roots[0]
            % a
        )

    for index in range(1 << (len(primes) - 1)):
        budget.consume(
            (len(base.entries) + len(terms) + 1)
            * (a.bit_length() + base.n.bit_length() + 1)
        )
        gray = index ^ (index >> 1)
        raw = terms[0] + sum(
            -term if gray & (1 << bit) else term
            for bit, term in enumerate(terms[1:])
        )
        b = (raw + a // 2) % a - a // 2
        polynomial = Polynomial(base.n, base.multiplier, a, b)
        roots = tuple(
            polynomial_roots(polynomial, base, entry, budget=budget)
            for entry in base.entries
        )
        yield polynomial, roots


def cached_steps(base, primes, budget):
    """Consume the entire bounded stream of exact Gray/cache steps."""
    family = PolynomialFamily(
        base, primes, budget=budget, memory_bytes=MEMORY_BYTES
    )
    while (step := family.next()) is not None:
        yield step.polynomial, step.roots


def root_utility(cached):
    """Include base, CRT/cache setup and all eight polynomial root sets."""
    budget = _budget()
    base = build_factor_base(
        4001 * 5003, bound=2000, budget=budget
    ).factor_base
    primes = family_assignments(
        base,
        256,
        factor_count=4,
        family_count=1,
        budget=budget,
        memory_bytes=MEMORY_BYTES,
    )[0]
    steps = cached_steps if cached else full_steps
    return [
        (
            polynomial.a,
            polynomial.b,
            polynomial.c,
            tuple(
                (item.prime, item.roots, item.all_positions) for item in roots
            ),
        )
        for polynomial, roots in steps(base, primes, budget)
    ]


def attempt_cohort(cached):
    """Use identical families and fresh rows under one allowance per input.

    This interoperability probe does not yet share relation stores across
    polynomials, serialize full jobs, or integrate the production dispatcher.
    """
    output = []

    for fixture in _corpus(SEED, 16):
        budget = _budget()
        n = fixture["n"]
        base = build_factor_base(n, bound=200, budget=budget).factor_base
        assignments = family_assignments(
            base,
            256,
            factor_count=3,
            family_count=8,
            seed=2934,
            budget=budget,
            memory_bytes=MEMORY_BYTES,
        )
        found = None

        for primes in assignments:
            steps = cached_steps if cached else full_steps

            for polynomial, roots in steps(base, primes, budget):
                collector = partial(SieveCollector, precomputed_roots=roots)
                # Both arms supply independently validated roots to the same
                # collector, avoiding an output-contract or validation bias.
                family_reserve = 32768 + 640 * len(base.entries) + 8192
                config = SieveConfig(
                    division="bucket",
                    residual_bound=1,
                    max_atoms=4096,
                    max_relations=4096,
                    max_partials=4096,
                    memory_bytes=MEMORY_BYTES - family_reserve,
                )
                job = QSJob(
                    polynomial,
                    base,
                    -256,
                    257,
                    budget=budget,
                    config=config,
                    collector_class=collector,
                    weight_two=True,
                )

                result = job.run()
                if result.divisor is not None:
                    found = sorted((result.divisor, result.cofactor))
                    assert found == fixture["factors"] and prod(found) == n
                    assert all(
                        utils.classify_prime(p) == utils.Primality.PROVEN
                        for p in found
                    )
                    break

            if found is not None:
                break

        output.append(
            {
                "n": n,
                "factors": found or [],
                "remaining": 1 if found is not None else n,
                "reason": (
                    "factor_found"
                    if found is not None
                    else "families_exhausted"
                ),
            }
        )
        assert prod(output[-1]["factors"]) * output[-1]["remaining"] == n

    return output


def main():
    """Record utility and all-outcome costs; P3.4 remains open."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--cold", choices=("full", "cached"))
    args = parser.parse_args()
    if args.cold:
        output = attempt_cohort(args.cold == "cached")
        usage = resource.getrusage(resource.RUSAGE_SELF)
        print(
            json.dumps(
                {
                    "output": output,
                    "rss_bytes": _rss_bytes(),
                    "cpu_seconds": usage.ru_utime + usage.ru_stime,
                }
            )
        )
        return

    if (
        args.output.exists()
        or not math.isfinite(args.warmup_seconds)
        or args.warmup_seconds < 3
        or args.repetitions < 9
    ):
        parser.error("use a new capture, three-second warmup and nine samples")

    measured_environment = environment()
    expected_roots = root_utility(False)
    assert root_utility(True) == expected_roots
    expected_attempts = attempt_cohort(False)
    assert attempt_cohort(True) == expected_attempts
    results = []

    for cached in (False, True):
        name = "cached" if cached else "full"
        results.append(
            _record(
                name + "_roots",
                lambda: root_utility(cached),
                expected_roots,
                args,
            )
        )
        results.append(
            _record(
                name + "_attempts",
                lambda: attempt_cohort(cached),
                expected_attempts,
                args,
            )
        )

    cold = []

    for name in ("full", "cached"):
        for _ in range(9):
            start = time.perf_counter()
            prior = resource.getrusage(resource.RUSAGE_CHILDREN)

            child = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "v2.benchmarks.phase_three_families",
                    "--output",
                    str(args.output),
                    "--cold",
                    name,
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=30,
            )
            value = json.loads(child.stdout)
            assert value["output"] == expected_attempts
            usage = resource.getrusage(resource.RUSAGE_CHILDREN)
            value.update(
                arm=name,
                lifecycle_seconds=time.perf_counter() - start,
                lifecycle_cpu_seconds=usage.ru_utime
                + usage.ru_stime
                - prior.ru_utime
                - prior.ru_stime,
            )
            cold.append(value)

    assert (
        measured_environment["source_sha256"] == environment()["source_sha256"]
    )
    args.output.write_text(
        json.dumps(
            {
                "milestone": "M29 / P3.4 family and root foundation",
                "environment": measured_environment,
                "command": sys.orig_argv,
                "results": results,
                "cold_samples": cold,
                "corpus": _corpus(SEED, 16),
                "validated_outcomes": expected_attempts,
                "completed_count": sum(
                    row["reason"] == "factor_found"
                    for row in expected_attempts
                ),
                "outcome_count": len(expected_attempts),
                "timing_scope": "All outcomes, including unfinished attempts",
                "held_out_seed": SEED,
                "scope": "Bounded families and small QS interoperability",
                "P3.4_complete": False,
                "remaining": "Shared jobs, dispatch and large bands",
                "rss_scope": "Warm RSS is shared; cold RSS is per child",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
