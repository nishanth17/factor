"""Matched complete two-stage ECM attempts on certified 40--80 digit inputs.

The control calls the real factorize_ecm API. The experimental arm changes
stage one only: bound-owned verified prime-power records, checked execution
and bounded prime-unit recovery. This is an experiment, not a production
stage_jobs integration. Known factors are only used by output validators.
"""

import argparse
import hashlib
import json
import platform
import signal
import statistics
import subprocess
import sys
import time
from dataclasses import asdict
from functools import lru_cache
from math import gcd, prod
from pathlib import Path
from types import SimpleNamespace

from .... import constants
from ....common import prime_sieve, utils
from ....ecm import core as ecm
from ....ecm import prac
from ...suites.build_phase_two_corpus import verify_certificates
from ...suites.phase_one import environment
from ...support.paths import (
    BENCHMARK_ROOT,
    source_path,
)

CORPUS = BENCHMARK_ROOT / "inputs/corpora/p41_ecm_40_80_corpus.json"
ARMS = ("production_ladder", "checked_prac", "gmp_ladder", "gmp_prac")
PYTHON_BACKEND = SimpleNamespace(ecm=ecm, prac=prac, gcd=gcd, integer=int)
TIERS = {
    "current": (2000, 147396, 32),
    "factor20": (11000, 1873422, 32),
}
MAX_PROGRAM_RECORDS = 4096


class CampaignTimeoutError(Exception):
    """The shared wall or CPU allowance was exhausted."""


def load_corpus():
    data = json.loads(source_path(CORPUS).read_text())
    verify_certificates(data["certificates"])
    for case in data["fixtures"]:
        factors = case["factors"]
        assert prod(factors) == case["n"]
        assert len(str(case["n"])) == case["digits"]
        assert [len(str(p)) for p in factors] == case["factor_digits"]
        assert all(str(p) in data["certificates"] for p in factors)
    return data


def build_program(b1):
    """Build once per attempt, with an explicit bound on owned records.

    External immutable records avoid LRU thrashing when B1 exceeds the
    scalar cache's capacity. Public execution re-verifies each record; that
    cost, construction, and prime-unit recovery are all inside the timer.
    """
    utils.require_integer(b1, "b1", 2)
    if b1 > 11000:
        raise ValueError("this experiment supports B1 <= 11000")
    primes = tuple(prime_sieve.prime_sieve(b1 + 1))
    if 2 * len(primes) > MAX_PROGRAM_RECORDS:
        raise ValueError("campaign program exceeds its finite record cap")
    prac.clear_cache()
    return tuple(
        (
            p,
            utils.prime_power(p, b1),
            prac.get_chain(utils.prime_power(p, b1)),
            prac.get_chain(p),
        )
        for p in primes
    )


def apply_record(scalar, record, point, n, a24, backend):
    return backend.prac.multiply(
        scalar,
        point,
        n,
        a24,
        backend.ecm.point_add,
        backend.ecm.point_double,
        chain=record,
    )


def candidate_stage_one(point, n, a24, program, extra, backend=PYTHON_BACKEND):
    for prime, power, record, prime_record in program:
        original = point
        try:
            point = apply_record(power, record, point, n, a24, backend)
            if backend.gcd(point[1], n) != n:
                continue
        except prac.NonunitPointError as result:
            if result.factor is not None:
                return None, result.factor

        # A prime-power collapse is retried only from its own starting point,
        # with at most log_prime(power) units. Never replay a whole campaign.
        extra["prime_power_replays"] += 1
        point, remaining = original, power
        while remaining > 1:
            extra["prime_units_replayed"] += 1
            try:
                point = apply_record(
                    prime, prime_record, point, n, a24, backend
                )
            except prac.NonunitPointError as result:
                return None, result.factor
            if backend.gcd(point[1], n) == n:
                return None, None
            remaining //= prime
    return point, None


def candidate(
    n, *, b1, b2, seed, max_curves, stats, extra, backend=PYTHON_BACKEND
):
    """Mirror the standalone production campaign with changed stage one."""
    program = build_program(b1)
    extra["program_entries"] = len(program)
    extra["program_instructions"] = sum(
        len(record.instructions) + len(unit.instructions)
        for _, _, record, unit in program
    )
    generator = utils.resolve_rng(seed, None)
    stage_two_primes = None
    for _ in range(max_curves):
        stats.curves += 1
        setup = backend.ecm.setup_curve(
            n, generator.randint(6, constants.MAX_RANDOM_ECM)
        )
        if setup.factor is not None:
            return setup.factor
        if setup.retry:
            stats.setup_retries += 1
            continue

        stats.stage_one_calls += 1
        point, factor = candidate_stage_one(
            setup.point, n, setup.a24, program, extra, backend
        )
        if factor is not None:
            return factor
        if point is None:
            stats.stage_one_saturations += 1
            continue

        if stage_two_primes is None:
            stage_two_primes = prime_sieve.segmented_sieve(b1 + 1, b2 + 1)
        stats.stage_two_calls += 1
        factor, saturated = backend.ecm.stage_two(
            point,
            n,
            setup.a24,
            b1,
            stage_two_primes,
            constants.GCD_BATCH_SIZE,
        )
        if saturated:
            stats.stage_two_saturations += 1
        if backend.ecm.utils.valid_divisor(factor, n):
            return factor
    return None


def validate_result(result, case):
    factor = result["factor"]
    if factor is None:
        assert result["unresolved"] == case["n"]
        assert result["cofactor"] is None
    else:
        assert utils.valid_divisor(factor, case["n"])
        assert factor * result["cofactor"] == case["n"]
        assert sorted((factor, result["cofactor"])) == sorted(case["factors"])
        assert result["unresolved"] is None


def validate_backend_pairs(samples):
    """Representation changes must preserve every uncensored transition."""
    indexed = {(r["seed"], r["arm"]): r for r in samples}
    checks = 0
    for seed, arm in indexed:
        if arm not in ("production_ladder", "checked_prac"):
            continue
        partner = "gmp_ladder" if arm.endswith("ladder") else "gmp_prac"
        if (seed, partner) not in indexed:
            continue
        left, right = indexed[seed, arm], indexed[seed, partner]
        if left["timed_out"] or right["timed_out"]:
            continue
        for field in ("factor", "cofactor", "unresolved", "stats", "extra"):
            assert left[field] == right[field], (seed, arm, field)
        checks += 1
    return checks


@lru_cache(maxsize=1)
def gmp_backend():
    from .p41_gmp import load

    return load()


def attempt(arm, n, seed, tier, seconds):
    backend = gmp_backend() if arm.startswith("gmp_") else PYTHON_BACKEND
    b1, b2, curves = TIERS[tier]
    stats = ecm.EcmStats()
    extra = {"prime_power_replays": 0, "prime_units_replayed": 0}

    def expired(signum, frame):
        raise CampaignTimeoutError

    old_wall = signal.signal(signal.SIGALRM, expired)
    old_cpu = signal.signal(signal.SIGPROF, expired)
    started_wall, started_cpu = time.perf_counter(), time.process_time()
    factor, timed_out = None, False
    signal.setitimer(signal.ITIMER_REAL, seconds)
    signal.setitimer(signal.ITIMER_PROF, seconds)
    try:
        native_n = backend.integer(n)
        if arm.endswith("ladder"):
            factor = backend.ecm.factorize_ecm(
                native_n,
                b1=b1,
                b2=b2,
                seed=seed,
                max_curves=curves,
                stats=stats,
                _known_composite=True,
            )
        else:
            factor = candidate(
                native_n,
                b1=b1,
                b2=b2,
                seed=seed,
                max_curves=curves,
                stats=stats,
                extra=extra,
                backend=backend,
            )
    except CampaignTimeoutError:
        timed_out = True
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.setitimer(signal.ITIMER_PROF, 0)
        elapsed = time.perf_counter() - started_wall
        cpu = time.process_time() - started_cpu
        signal.signal(signal.SIGALRM, old_wall)
        signal.signal(signal.SIGPROF, old_cpu)
    factor = int(factor) if factor is not None else None
    return {
        "arm": arm,
        "seed": seed,
        "wall_seconds": elapsed,
        "cpu_seconds": cpu,
        "factor": factor,
        "cofactor": n // factor if factor else None,
        "unresolved": None if factor else n,
        "timed_out": timed_out,
        "stats": asdict(stats),
        "extra": extra,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tier", choices=TIERS, default="current")
    parser.add_argument("--case", action="append")
    parser.add_argument("--cold-case")
    parser.add_argument("--seconds", type=float, default=60)
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--gmp", action="store_true")
    parser.add_argument("--cold-arm", choices=ARMS)
    args = parser.parse_args()
    if platform.python_implementation() != "PyPy" or sys.version_info[:2] != (
        3,
        11,
    ):
        parser.error("requires PyPy implementing Python 3.11")
    if args.output.exists():
        parser.error("refusing to overwrite evidence")
    if args.seconds <= 0:
        parser.error("seconds must be positive")
    arms = ARMS if args.gmp else ARMS[:2]
    backend_identity = gmp_backend().identity if args.gmp else None
    corpus = load_corpus()
    cases = [
        c for c in corpus["fixtures"] if not args.case or c["id"] in args.case
    ]
    if args.case and {c["id"] for c in cases} != set(args.case):
        parser.error("unknown case")
    if args.tier == "factor20" and not args.case:
        cases = [c for c in cases if min(c["factor_digits"]) == 20]
    cold_cases = [
        c
        for c in corpus["fixtures"]
        if c["id"] == (args.cold_case or cases[0]["id"])
    ]
    if not cold_cases:
        parser.error("unknown cold case")
    cold_case = cold_cases[0]
    if args.cold_arm:
        result = attempt(
            args.cold_arm,
            cases[0]["n"],
            corpus["seeds"][0],
            args.tier,
            args.seconds,
        )
        validate_result(result, cases[0])
        args.output.write_text(json.dumps(result))
        return

    report = {
        "environment": environment(),
        "tier": args.tier,
        "gmp_backend": backend_identity,
        "arms": arms,
        "b1_b2_curves": TIERS[args.tier],
        "seconds_per_arm": args.seconds,
        "pilot": args.pilot,
        "corpus_sha256": hashlib.sha256(
            source_path(CORPUS).read_bytes()
        ).hexdigest(),
        "source_sha256": {},
        "cases": [],
        "cold": [],
        "cold_case": cold_case["id"],
        "limitations": [
            "fixed inspected fixtures, not a population-level speedup claim",
            "32-curve attempts; unresolved composites remain explicit",
            "balanced 25--40 digit factors exceed the factor20 target tier",
            "candidate program and replay are benchmark-only, outside B3",
            "public record verification and construction charged per attempt",
        ],
    }
    for name in (
        "v2/ecm.py",
        "v2/prac.py",
        __file__,
        str(source_path(BENCHMARK_ROOT / "p41_gmp.py")),
    ):
        path = Path(name)
        report["source_sha256"][path.name] = hashlib.sha256(
            source_path(path).read_bytes()
        ).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    seeds = corpus["seeds"][:1] if args.pilot else corpus["seeds"]
    for case in cases:
        row = {
            "id": case["id"],
            "digits": case["digits"],
            "factor_digits": case["factor_digits"],
            "warmup": [],
            "samples": [],
            "repeat_samples": [],
        }
        if not args.pilot:
            for arm in arms:
                warm_start, runs = time.perf_counter(), 0
                while time.perf_counter() - warm_start < 3:
                    result = attempt(
                        arm, case["n"], seeds[0], args.tier, args.seconds
                    )
                    validate_result(result, case)
                    runs += 1
                row["warmup"].append(
                    {
                        "arm": arm,
                        "seconds": time.perf_counter() - warm_start,
                        "validated_runs": runs,
                    }
                )

        for index, seed in enumerate(seeds):
            for arm in arms[:: 1 if index % 2 == 0 else -1]:
                result = attempt(arm, case["n"], seed, args.tier, args.seconds)
                validate_result(result, case)
                row["samples"].append(result)
        # Seed-dependent time-to-factor spread is expected. Repeat the same
        # nine paired seeds when the group is variable; retain both sets.
        if not args.pilot and any(
            (max(times) - min(times)) / statistics.median(times) > 0.25
            for arm in arms
            for times in [
                [r["wall_seconds"] for r in row["samples"] if r["arm"] == arm]
            ]
        ):
            for index, seed in enumerate(seeds):
                for arm in arms[:: -1 if index % 2 == 0 else 1]:
                    result = attempt(
                        arm, case["n"], seed, args.tier, args.seconds
                    )
                    validate_result(result, case)
                    row["repeat_samples"].append(result)
        row["backend_pairs_verified"] = sum(
            validate_backend_pairs(row[key])
            for key in ("samples", "repeat_samples")
        )
        report["cases"].append(row)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(
            case["id"],
            {
                arm: round(
                    statistics.median(
                        r["wall_seconds"]
                        for r in row["samples"]
                        if r["arm"] == arm
                    ),
                    4,
                )
                for arm in arms
            },
            flush=True,
        )

    if not args.pilot:
        # Cold process startup, imports, certificate validation and first
        # complete attempt are deliberately separate from warmed samples.
        cold_path = args.output.with_suffix(".cold.json")
        for arm in arms:
            for _ in seeds:
                started = time.perf_counter()
                subprocess.run(
                    [
                        sys.executable,
                        "-B",
                        "-m",
                        __spec__.name,
                        "--cold-arm",
                        arm,
                        "--case",
                        cold_case["id"],
                        "--tier",
                        args.tier,
                        "--seconds",
                        str(args.seconds),
                        "--output",
                        str(cold_path),
                    ],
                    check=True,
                )
                total = time.perf_counter() - started
                result = json.loads(source_path(cold_path).read_text())
                validate_result(result, cold_case)
                cold_path.unlink()
                report["cold"].append({"total_seconds": total, **result})
        args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
