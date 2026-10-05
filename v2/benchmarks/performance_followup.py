"""Validated follow-up budget, coordination and setup controls on PyPy 3.11.

Invoke this file directly so each immutable runtime is imported in isolation.
The shared audit runner enforces warmup, repetition and stability gates.
Use one band and case per process for timing confirmation to avoid preceding
cases influencing PyPy's trace history. Capacity outcomes are not speedups.
"""

import argparse
import hashlib
import importlib.util
import json
import sys
import tempfile
import time
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[2]
    )
    parser.add_argument("--baseline-json", type=Path)
    parser.add_argument("--variant", default="none")
    parser.add_argument("--experiments", type=Path)
    parser.add_argument("--corpus", type=Path)
    parser.add_argument(
        "--split", choices=("training", "held_out"), default="training"
    )
    parser.add_argument("--bands", default="small")
    parser.add_argument(
        "--cases",
        default="parallel_lease1m_serial",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--repetitions", type=int, default=9)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 is required")
    if args.output.exists():
        parser.error("choose a new capture path")
    bench = args.root / "v2/benchmarks"
    baseline = (
        args.baseline_json
        or bench / "inputs/baselines/performance_followup_baseline.json"
    )
    corpus = (
        args.corpus or bench / "inputs/corpora/performance_audit_corpus.json"
    )
    specification = importlib.util.spec_from_file_location(
        "followup_audit", bench / "performance_audit.py"
    )
    audit = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(audit)
    original_config, original_run = audit.make_config, audit.run_one

    def configuration(case, band, modules):
        if case == "resieve_capacity":
            return dict(
                base_bound=10000, block_width=4096, memory_bytes=64 * 2**20
            )
        if case == "family_setup":
            return dict(base_bound=10000, family_count=4, repetitions=16)
        if case == "p2_disabled":
            return dict(
                pm1_attempts=0, pm1_b2=10**12, ecm_tiers=((10**6, 10**12, 0),)
            )
        config = original_config(case, band, modules)
        if "lease" in case:
            config = replace(config, assignment_work=10_000_000)
        if "mem10" in case:
            config = replace(
                config, base_bound=10000, worker_memory_bytes=10 * 2**20
            )
        return config

    def run(case, fixture, seed, config, modules, pool=None):
        if case == "resieve_capacity":
            return resieve_capacity(fixture, seed, config, modules)
        if case == "family_setup":
            return family_setup(fixture, seed, config, modules)
        if "lease" in case:
            original_budget = modules["budget"]
            limit = 1_000_000 if "lease1m" in case else 15_000_000

            def allowance(**options):
                options["work_limit"] = limit
                return original_budget.Budget(**options)

            modules = dict(
                modules,
                budget=SimpleNamespace(
                    Budget=allowance,
                    BudgetExhaustedError=original_budget.BudgetExhaustedError,
                ),
            )
        return original_run(case, fixture, seed, config, modules, pool)

    audit.make_config, audit.run_one = configuration, run
    with tempfile.TemporaryDirectory(prefix="factor-followup-") as directory:
        runtime = Path(directory)
        audit.materialize_baseline(baseline, runtime)
        if args.variant != "none":
            experiments = json.loads(
                (
                    args.experiments
                    or bench / "history/performance_followup_experiments.json"
                ).read_text()
            )
            if (
                hashlib.sha256(baseline.read_bytes()).hexdigest()
                != experiments["baseline_sha256"]
            ):
                raise ValueError(
                    "variant belongs to a different frozen baseline"
                )

            variant = experiments["variants"][args.variant]
            original = audit.checked_baseline(baseline)["source"]

            for name, source in variant["source"].items():
                if (
                    name not in original
                    or hashlib.sha256(source.encode()).hexdigest()
                    != variant["source_sha256"][name]
                ):
                    raise ValueError("invalid variant source")

                (runtime / name).write_text(source)

        audit.run(
            SimpleNamespace(
                runtime_root=runtime,
                baseline_json=None,
                corpus=corpus,
                split=args.split,
                bands=args.bands,
                cases=args.cases,
                output=args.output,
                repetitions=args.repetitions,
                warmup_seconds=args.warmup_seconds,
                profile=False,
                revert="none",
                polynomial_order="generated",
            )
        )

    data = json.loads(args.output.read_text())
    data["followup"] = dict(
        variant=args.variant,
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        lease1m_work=1_000_000,
        lease15m_work=15_000_000,
        note="Lease work overrides apply only to named lease cases; "
        "wall/CPU and memory limits stay matched.",
    )
    args.output.write_text(json.dumps(data, indent=2) + "\n")


def family_setup(fixture, seed, config, modules):
    """Check finite setup CRT/inverse outputs, including repeated bases."""
    started = time.perf_counter()
    budget = modules["budget"].Budget(
        work_limit=2_000_000_000, seconds=5, cpu_seconds=5
    )
    parallel = modules["parallel"]
    base = parallel.build_factor_base(
        fixture["n"],
        bound=config["base_bound"],
        budget=budget,
        memory_bytes=32 * 2**20,
    ).factor_base
    if base is None:
        raise ValueError("family fixture must not split in base setup")
    assignments = parallel.family_assignments(
        base,
        512,
        factor_count=3,
        family_count=config["family_count"],
        pool_size=16,
        seed=seed,
        budget=budget,
    )
    identities = []

    for _ in range(config["repetitions"]):
        for primes in assignments:
            family = parallel.PolynomialFamily(
                base, primes, budget=budget, memory_bytes=32 * 2**20
            )

            for entry, inverse in zip(base.entries, family.inverses):
                if (
                    inverse is not None
                    and family.a * inverse % entry.prime != 1
                ):
                    raise AssertionError("invalid cached family inverse")

            b = sum(family.terms) % family.a
            if (b * b - base.n_prime) % family.a:
                raise AssertionError("invalid family CRT output")
            identities.append(family.identity)

    return dict(
        id=fixture["id"],
        seed=seed,
        seconds=time.perf_counter() - started,
        cpu_seconds=budget.cpu_used,
        work_used=budget.used,
        completed=False,
        factors=[],
        remaining=[fixture["n"]],
        certainty=[],
        reason="setup_complete",
        stats=dict(families=len(identities)),
        signature=hashlib.sha256(repr(identities).encode()).hexdigest(),
    )


def resieve_capacity(fixture, seed, config, modules):
    """Validate each returned atom at a matched finite storage capacity."""
    started = time.perf_counter()
    budget = modules["budget"].Budget(
        work_limit=2_000_000_000, seconds=5, cpu_seconds=5
    )
    stats, signature, divisor = {}, None, None
    base = modules["factor_base"].build_factor_base(
        fixture["n"],
        bound=config["base_bound"],
        budget=budget,
        memory_bytes=config["memory_bytes"],
    )
    if base.divisor:
        raise ValueError("capacity fixture must not split during base setup")
    polynomial = modules["polynomial"].qs_polynomial(base.factor_base)

    try:
        collector = modules["sieve_collector"].SieveCollector(
            polynomial,
            base.factor_base,
            budget=budget,
            config=modules["sieve_collector"].SieveConfig(
                score_policy="powers",
                division="resieve",
                block_width=config["block_width"],
                memory_bytes=config["memory_bytes"],
                residual_bound=10000,
            ),
        )

        result = collector.collect(-2048, 2048)
        reason, stats, divisor = result.reason, result.stats, result.divisor
        stats = dict(stats, workspace_bytes=result.workspace_bytes)

        for atom in result.atoms:
            expected = abs(
                atom.polynomial.a * atom.polynomial.value(atom.position)
            )
            actual = atom.residual
            for prime, exponent in atom.exponents:
                actual *= prime**exponent
            if actual != expected:
                raise AssertionError("invalid recovered exact exponents")

        signature = hashlib.sha256(
            repr(
                sorted(
                    (a.position, a.exponents, a.residual) for a in result.atoms
                )
            ).encode()
        ).hexdigest()
    except MemoryError:
        reason = "memory_limit"

    if divisor is not None and not modules["utils"].valid_divisor(
        divisor, fixture["n"]
    ):
        raise AssertionError("invalid split")
    return dict(
        id=fixture["id"],
        seed=seed,
        seconds=time.perf_counter() - started,
        cpu_seconds=budget.cpu_used,
        work_used=budget.used,
        completed=False,
        factors=[],
        remaining=[fixture["n"]],
        certainty=[],
        reason=reason,
        stats=stats,
        signature=signature,
    )


if __name__ == "__main__":
    main()
