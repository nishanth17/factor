"""Diagnostic profiles of matched APIs; these are not speedup measurements."""

import argparse
import cProfile
import gzip
import hashlib
import json
import pstats
import time
from pathlib import Path

from v2.benchmarks.phase_one import environment
from v2.budget import Budget
from v2.factor import factorize
from v2.portfolio import PortfolioConfig, factorize_bounded


def run(args):
    """Profile fixed held-out inputs with finite work but no clock caps."""
    root = Path(__file__).resolve().parents[2]
    control_path = root / "v2/benchmarks/phase_two_m12_feasible_pypy.json.gz"
    # Storage changed; the manifest preserves the original capture's hash.
    with gzip.open(control_path, "rt", encoding="utf-8") as source:
        control = json.load(source)
    corpus_path = root / control["corpus"]
    corpus = json.loads(corpus_path.read_text())
    bands = ("balanced_20d", "random_small", "close_small", "primes")
    fixtures = []
    for band in bands:
        fixtures.extend(
            [
                item
                for item in corpus["fixtures"]
                if item["split"] == "held_out" and item["band"] == band
            ][:5]
        )
    config = PortfolioConfig(**control["config"])
    seeds = control["seeds"]
    paths = {
        engine: args.output.with_name(f"{args.output.stem}_{engine}.prof")
        for engine in ("phase_one", "bounded")
    }
    if args.output.exists() or any(path.exists() for path in paths.values()):
        raise FileExistsError("select fresh profile output filenames")

    def invoke(engine, fixture, seed):
        """Validate against the oracle after the API has returned."""
        n = fixture["n"]
        if engine == "bounded":
            run = factorize_bounded(
                n,
                seed=seed,
                config=config,
                budget=Budget(
                    work_limit=control["work_cap_bounded_only"],
                    seconds=None,
                    cpu_seconds=None,
                ),
            )
            answer, reason = run.result, run.reason
        else:
            b1, b2, curves = config.ecm_tiers[0]
            answer = factorize(
                n,
                seed=seed,
                trial_bound=config.trial_bound,
                rho_attempts=config.rho_attempts,
                rho_evaluations=config.rho_evaluations,
                ecm_b1=b1,
                ecm_b2=b2,
                ecm_curves=curves,
            )
            reason = "complete" if answer.complete else "exhausted"
        actual = {item.value: item.exponent for item in answer.factors}
        expected = dict(fixture["factors"])
        if answer.reconstruct() != n or any(
            prime not in expected or exponent > expected[prime]
            for prime, exponent in actual.items()
        ):
            raise AssertionError("invalid profiled factorization")
        if answer.complete and actual != expected:
            raise AssertionError("profiled completion disagrees with oracle")
        return {"id": fixture["id"], "seed": seed, "reason": reason}

    rows = []
    for engine in paths:
        started = time.perf_counter()
        calls = 0
        while time.perf_counter() - started < args.warmup_seconds:
            for fixture in fixtures:
                invoke(engine, fixture, seeds[0])
                calls += 1
        warmup = {"seconds": time.perf_counter() - started, "calls": calls}
        profiler = cProfile.Profile()
        outcomes = [
            profiler.runcall(invoke, engine, fixture, seed)
            for fixture in fixtures
            for seed in seeds
        ]
        profiler.dump_stats(str(paths[engine]))
        stats = pstats.Stats(profiler)
        functions = []
        for (filename, line, name), counts in stats.stats.items():
            primitive, total, exclusive, cumulative, _ = counts
            functions.append(
                {
                    "file": filename,
                    "line": line,
                    "function": name,
                    "primitive_calls": primitive,
                    "total_calls": total,
                    "exclusive_seconds": exclusive,
                    "cumulative_seconds": cumulative,
                }
            )
        rows.append(
            {
                "engine": engine,
                "warmup": warmup,
                "outcomes": outcomes,
                "profile": str(paths[engine]),
                "functions": sorted(
                    functions,
                    key=lambda item: item["exclusive_seconds"],
                    reverse=True,
                ),
            }
        )
        print(engine, "profile done", flush=True)
    return {
        "milestone": "phase_two_m12_diagnostic_profile",
        "environment": environment(),
        "diagnostic_source_sha256": hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        "control": str(control_path),
        "config": control["config"],
        "fixtures": [item["id"] for item in fixtures],
        "seeds": seeds,
        "work_limit": control["work_cap_bounded_only"],
        "profiles": rows,
        "limitations": [
            "Clock caps disabled to inspect complete local "
            "candidate schedules",
            "Profiling distorts JIT behavior and supplies no speed ratios",
            "Cumulative function times overlap and must not be added",
            "Five diagnostic inputs per band are not a promotion cohort",
        ],
    }


def main():
    """Save profiles and attribution metadata without replacing evidence."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    args = parser.parse_args()
    if args.warmup_seconds < 3:
        parser.error("use at least three seconds of workload warmup")
    args.output.write_text(json.dumps(run(args), indent=2) + "\n")


if __name__ == "__main__":
    main()
