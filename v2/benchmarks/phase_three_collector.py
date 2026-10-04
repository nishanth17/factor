"""P3.2 matched collection costs, tuning sweeps and held-out coverage."""

import argparse
import cProfile
import hashlib
import json
import math
import statistics
import subprocess
import sys
import time
from collections import Counter
from dataclasses import asdict, replace
from pathlib import Path

from ..qs import (
    Polynomial,
    SieveCollector,
    SieveConfig,
    SieveResult,
    build_factor_base,
    collect_block,
    verify_atomic,
    verify_combined,
)
from ..tests.test_qs import reference_positions
from .phase_one import environment
from .phase_three_reference import (
    BASELINE,
    ROOT,
    _budget,
    _cases,
    _measure,
    _rss_bytes,
)

# Training only: coefficients and windows are supplied, not SIQS-selected.
TRAINING = [
    dict(name="qs", n=10403, h=1, a=1, b=102, bound=40, residual=97),
    dict(name="mpqs", n=10403, h=1, a=49, b=8, bound=100, residual=97),
    dict(
        name="multiplier", n=1022117, h=3, a=1, b=1752, bound=100, residual=500
    ),
]
HELD_OUT = [
    dict(
        name="balanced_30bit",
        n=25013 * 25031,
        h=1,
        a=1,
        b=25022,
        bound=100,
        residual=500,
    ),
    dict(
        name="balanced_24bit",
        n=4001 * 4003,
        h=1,
        a=1,
        b=4002,
        bound=100,
        residual=500,
    ),
]
LO, HI = -128, 129


def _signature(atoms):
    """Consume exact signs, multiplicities and residuals."""
    return {
        atom.position: (atom.sign, atom.exponents, atom.residual)
        for atom in atoms
    }


def _setup(fixture, budget):
    """Build a checked base and supplied polynomial under a shared budget."""
    base = build_factor_base(
        fixture["n"],
        multiplier=fixture["h"],
        bound=fixture["bound"],
        budget=budget,
    ).factor_base
    polynomial = Polynomial(
        fixture["n"], fixture["h"], fixture["a"], fixture["b"]
    )
    return base, polynomial


def _run(fixture, config, *, reuse=False, reference=False):
    """Create a validated run including output consumption."""
    shared = _budget()
    base, polynomial = _setup(fixture, shared)
    worker = None
    if reuse:
        worker = SieveCollector(polynomial, base, config=config, budget=shared)

    def call():
        """Validate full/combined provenance and expose diagnostics."""
        nonlocal worker
        budget = _budget()
        if not reuse:
            local_base, local_polynomial = _setup(fixture, budget)
        else:
            local_base, local_polynomial = base, polynomial
        if reference:
            exhaustive = collect_block(
                local_polynomial,
                local_base,
                LO,
                HI,
                residual_bound=fixture["residual"],
                budget=budget,
            )
            if (
                exhaustive.reason != "complete"
                or exhaustive.divisor is not None
            ):
                raise AssertionError("reference failed collection")
            worker = SieveCollector(
                local_polynomial, local_base, config=config, budget=budget
            )
            admission_stats = dict(
                duplicates=0,
                dropped_partials=0,
                evictions=0,
                admitted_atoms=0,
                matches=0,
            )
            for atom in exhaustive.relations:
                refusal = worker._admit(atom, admission_stats)
                if refusal is not None:
                    raise AssertionError(
                        "matched reference store refused atom"
                    )
            # The baseline keeps the exhaustive output alive during transfer;
            # charge this extra reservation instead of hiding peak workspace.
            workspace = worker._workspace + exhaustive.workspace_bytes
            if workspace > config.memory_bytes:
                raise AssertionError(
                    "matched baseline exceeds memory allowance"
                )
            result = SieveResult(
                tuple(worker._atoms.values()),
                tuple(worker._full),
                tuple(worker._combined),
                tuple(worker._pending.values()),
                None,
                HI,
                "complete",
                {"scanned": exhaustive.scanned, **admission_stats},
                workspace,
            )
        elif not reuse:
            worker = SieveCollector(
                local_polynomial, local_base, config=config, budget=budget
            )
        else:
            # Reuse metadata/scores, not prior output or duplicate suppression.
            # This internal benchmark reset is not a public collector API.
            worker._workspace = initial_workspace
            worker._atoms.clear()
            worker._pending.clear()
            worker._atom_bytes.clear()
            worker._full.clear()
            worker._combined.clear()
            worker.budget = budget
        if not reference:
            result = worker.collect(LO, HI)
        if result.reason != "complete" or result.divisor is not None:
            raise AssertionError(f"collector failed: {result.reason}")
        store = {atom.relation_id: atom for atom in result.atoms}
        for atom in result.full_relations:
            verify_atomic(atom, local_base, budget=budget)
        for relation in result.combined_relations:
            verify_combined(relation, local_base, store, budget=budget)
        counts = Counter(atom.residual for atom in result.atoms)
        if len(result.full_relations) != counts[1] or (
            len(result.combined_relations)
            != sum(
                count // 2
                for residual, count in counts.items()
                if residual != 1
            )
        ):
            raise AssertionError("bounded store lost a full relation or match")
        info = {
            **result.stats,
            "full": len(result.full_relations),
            "combined": len(result.combined_relations),
            "unmatched": len(result.partial_ids),
            "workspace_bytes": result.workspace_bytes,
            "work_used": budget.used,
            "reason": result.reason,
        }
        return _signature(result.atoms), info

    initial_workspace = worker._workspace if worker is not None else 0
    return call


def _coverage(fixture, config):
    """Compare candidate stages to independent complete norm factorization."""
    base, polynomial = _setup(fixture, _budget())
    expected = reference_positions(
        polynomial, base, LO, HI, fixture["residual"]
    )
    call = _run(fixture, config)
    signature, info = call()
    if any(expected.get(x) != payload for x, payload in signature.items()):
        raise AssertionError("collector admitted invalid payload")
    missed = sorted(set(expected) - set(signature))
    if missed and config.threshold_extra == 0:
        raise AssertionError("safe score missed admissible positions")
    # Independently enumerate which values at each candidate stage are useful;
    # Factorizations come from the generic norm oracle.
    worker = SieveCollector(polynomial, base, config=config, budget=_budget())
    selected = set()
    for lo in range(LO, HI, config.block_width):
        hi = min(HI, lo + config.block_width)
        threshold = worker._sieve(
            lo,
            hi,
            dict(
                root_hits=0,
                slice_bytes=0,
                slice_allocations=0,
                blocks=0,
                threshold_min=2**63,
                threshold_max=0,
            ),
        )
        selected.update(
            lo + offset
            for offset in range(hi - lo)
            if worker._scores[offset] >= threshold
        )
    return (
        call,
        signature,
        {
            "fixture": fixture["name"],
            "expected_atoms": len(expected),
            "verified_atoms": len(signature),
            "missed_positions": missed,
            "missed_candidate_positions": sorted(set(expected) - selected),
            "false_candidates": len(selected - set(expected)),
            "candidate_positions": len(selected),
            **info,
        },
    )


def _config(fixture, **options):
    """Keep eviction out of score coverage comparisons."""
    return SieveConfig(
        residual_bound=fixture["residual"],
        max_partials=4096,
        max_relations=4096,
        max_atoms=4096,
        **options,
    )


def _capture(
    name, stage, fixture, config, args, *, reference=False, reuse=False
):
    """Keep raw costs and independent candidate coverage counts."""
    call, expected, coverage = _coverage(fixture, config)
    if reference or reuse:
        call = _run(fixture, config, reference=reference, reuse=reuse)
    if reference:
        payload, reference_info = call()
        if payload != expected:
            raise AssertionError("matched baseline differs from exact oracle")
        coverage = {
            "fixture": fixture["name"],
            "expected_atoms": len(expected),
            "verified_atoms": len(payload),
            **reference_info,
        }
    measured = _measure(
        lambda: call()[0], expected, 10, args.repetitions, args.warmup_seconds
    )
    print(
        name,
        f"{measured['median_seconds'] * 1000:.3f} ms",
        f"stable={measured['stable']}",
        flush=True,
    )
    return {
        "name": name,
        "stage": stage,
        "fixture": fixture,
        "config": asdict(config),
        "reference": reference,
        "reuse": reuse,
        "coverage": coverage,
        "measurement": measured,
    }


def _cold_worker():
    """Include imports, setup, division, matching and checked output."""
    fixture = TRAINING[2]
    config = _config(fixture)
    signature, info = _run(fixture, config)()
    return {
        "signature": signature,
        "stats": info,
        "peak_rss_bytes": _rss_bytes(),
    }


def _division_utilities(args):
    """Time exact recovery on identical candidates with cached root hits."""
    fixture = TRAINING[2]
    base, polynomial = _setup(fixture, _budget())
    expected = reference_positions(
        polynomial, base, LO, HI, fixture["residual"]
    )
    rows = []
    for mode in ("full", "roots", "bucket"):
        config = _config(fixture, block_width=512, division=mode)
        worker = SieveCollector(
            polynomial, base, config=config, budget=_budget()
        )
        stats = dict(
            root_hits=0,
            slice_bytes=0,
            slice_allocations=0,
            blocks=0,
            threshold_min=2**63,
            threshold_max=0,
        )
        worker._sieve(LO, HI, stats)

        def divide():
            """Include full exponent recovery, primality, GCD and verifier."""
            worker.budget = _budget()
            atoms = []
            stats = dict(
                candidates=0,
                zeros=0,
                division_primes=0,
                division_steps=0,
                composite_residuals=0,
            )
            for position in range(LO, HI):
                atom, divisor = worker._divide(position, position - LO, stats)
                if divisor is not None:
                    raise AssertionError("unexpected divisor in utility")
                if atom is not None:
                    atoms.append(atom)
            return _signature(atoms)

        measured = _measure(
            divide, expected, 10, args.repetitions, args.warmup_seconds
        )
        print(
            "division_" + mode,
            f"{measured['median_seconds'] * 1000:.3f} ms",
            f"stable={measured['stable']}",
            flush=True,
        )
        rows.append(
            {
                "name": "division_" + mode,
                "measurement": measured,
                "candidate_assignment": [LO, HI],
                "expected_atoms": len(expected),
            }
        )
    return rows


def main():
    """Freeze tuning before held-out runs; save runtime/source/raw samples."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--cold-worker", action="store_true")
    args = parser.parse_args()
    if args.cold_worker:
        print(json.dumps(_cold_worker()))
        return
    if args.output is None:
        parser.error("--output is required")
    if (
        not math.isfinite(args.warmup_seconds)
        or args.warmup_seconds < 3
        or args.repetitions < 9
    ):
        parser.error("require at least three seconds warmup and nine samples")
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("supported runtime is PyPy Python 3.11")
    measured_environment = environment()
    frozen = json.loads(BASELINE.read_text())
    unchanged = all(
        hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest
        for path, digest in frozen["production_source_sha256"].items()
    )
    if not unchanged:
        raise AssertionError("M17 production control changed")
    rows, tuning, coverage_sweeps = [], [], []
    fixture = TRAINING[2]
    options = [("base", {})]
    for key, values in (
        ("score_backend", ("bytearray", "array")),
        ("marking", ("dense", "bucket")),
        ("division", ("full", "bucket")),
        ("block_width", (16, 64, 1024)),
        ("metadata_chunk", (1, 16, 1024)),
        ("small_prime_cutoff", (3, 8, 32)),
        ("threshold_extra", (4, 16, 64)),
    ):
        options.extend((f"{key}_{v}", {key: v}) for v in values)
    for name, overrides in options:
        for coverage_fixture in TRAINING:
            _, _, coverage = _coverage(
                coverage_fixture, _config(coverage_fixture, **overrides)
            )
            coverage_sweeps.append({"configuration": name, **coverage})
        row = _capture(
            "training_" + name,
            "training",
            fixture,
            _config(fixture, **overrides),
            args,
        )
        rows.append(row)
        if row["measurement"]["stable"] and not overrides.get(
            "threshold_extra"
        ):
            tuning.append(row)
    # Choose a measured complete-coverage configuration from training only.
    winner = min(tuning, key=lambda row: row["measurement"]["median_seconds"])
    selected = SieveConfig(**winner["config"])
    freeze = {
        "training_winner": winner["name"],
        "config": asdict(selected),
        "selection": "lowest stable training median among safe-score arms",
        "timestamp_unix": time.time(),
    }
    # Persist the choice before any held-out measurements or source adaptation.
    freeze_path = args.output.with_suffix(".frozen.json")
    freeze_path.write_text(json.dumps(freeze, indent=2) + "\n")
    for fixture in TRAINING + HELD_OUT:
        stage = "held_out" if fixture in HELD_OUT else "matched_training"
        config = replace(selected, residual_bound=fixture["residual"])
        rows.append(
            _capture(fixture["name"] + "_sieve", stage, fixture, config, args)
        )
        rows.append(
            _capture(
                fixture["name"] + "_exhaustive",
                stage,
                fixture,
                config,
                args,
                reference=True,
            )
        )
    rows.append(
        _capture(
            "cached_metadata_and_buffers",
            "utility",
            TRAINING[2],
            selected,
            args,
            reuse=True,
        )
    )
    utilities = _division_utilities(args)
    # Measure the production portfolio with unchanged code and identical seed.
    cases, control_corpus = _cases()
    name, candidates, expected, iterations = cases[-1]
    control = _measure(
        next(iter(candidates.values())),
        expected,
        iterations,
        args.repetitions,
        args.warmup_seconds,
    )
    print(name, f"{control['median_seconds'] * 1000:.3f} ms", flush=True)
    # Diagnostic profiling is separate from timing whenever sieve cost exceeds
    # the exhaustive cost; profiler overhead is never a speed ratio.
    profiles = []
    for fixture in TRAINING + HELD_OUT:
        pair = [
            row
            for row in rows
            if row["name"]
            in (fixture["name"] + "_sieve", fixture["name"] + "_exhaustive")
        ]
        if (
            pair[0]["measurement"]["median_seconds"]
            > (pair[1]["measurement"]["median_seconds"])
        ):
            profile = cProfile.Profile()
            call = _run(
                fixture, replace(selected, residual_bound=fixture["residual"])
            )
            profile.enable()
            for _ in range(200):
                call()
            profile.disable()
            profile_path = args.output.with_suffix(
                "." + fixture["name"] + ".prof"
            )
            profile.dump_stats(str(profile_path))
            profiles.append(str(profile_path))
    cold = []
    expected_cold = _cold_worker()["signature"]
    for _ in range(args.repetitions):
        started = time.perf_counter()
        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "v2.benchmarks.phase_three_collector",
                "--cold-worker",
            ],
            capture_output=True,
            text=True,
            check=True,
            timeout=30,
        )
        elapsed = time.perf_counter() - started
        response = json.loads(completed.stdout)
        if response["signature"] != json.loads(json.dumps(expected_cold)):
            raise AssertionError("invalid cold collection")
        # Payloads were consumed and verified. Avoid copying duplicate outputs
        # into every cold sample; preserve their common signature hash instead.
        del response["signature"]
        cold.append(
            {"lifecycle_seconds": elapsed, "correct": True, **response}
        )
    if environment()["source_sha256"] != measured_environment["source_sha256"]:
        raise AssertionError("measured source changed during capture")
    args.output.write_text(
        json.dumps(
            {
                "milestone": "M22 / P3.2",
                "environment": measured_environment,
                "command": sys.orig_argv,
                "production_m17_hashes_unchanged": unchanged,
                "training": TRAINING,
                "held_out": HELD_OUT,
                "window": [LO, HI],
                "freeze": freeze,
                "results": rows,
                "coverage_sweeps": coverage_sweeps,
                "division_utilities": utilities,
                "complete_factorization_control": control,
                "control_inputs": control_corpus["complete_inputs"],
                "control_seed": 7,
                "profiles": profiles,
                "cold_samples": cold,
                "cold_median_seconds": statistics.median(
                    sample["lifecycle_seconds"] for sample in cold
                ),
                "cold_signature_sha256": hashlib.sha256(
                    json.dumps(expected_cold, sort_keys=True).encode()
                ).hexdigest(),
                "budgets": {
                    "work": 20_000_000,
                    "wall_seconds": 5,
                    "cpu_seconds": 5,
                    "owned_bytes": 8_388_608,
                },
                "scope": "collection completion; QS extraction awaits P3.3",
                "timing": "setup through checked output consumption",
                "jit": "default JIT; >=3s warmup per arm; adaptive sampling",
                "rss": "high-water includes JIT, oracles and prior work",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
