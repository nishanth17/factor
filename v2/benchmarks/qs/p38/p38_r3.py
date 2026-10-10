"""R3 preparation, matrix and complete-factor experiments on PyPy."""

import argparse
import ast
import cProfile
import hashlib
import importlib
import json
import platform
import random
import statistics
import subprocess
import sys
import tempfile
import time
import types
from dataclasses import asdict, replace
from functools import lru_cache
from pathlib import Path

from ...infrastructure.performance.performance_audit import (
    checked_baseline,
    fingerprint,
    verify_corpus,
)
from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from ...support.paths import (
    BENCHMARK_ROOT,
    REPOSITORY_ROOT,
    runtime_module,
    source_path,
)
from .p38_r3_experiments import history_filter, live_compaction

ROOT = REPOSITORY_ROOT
HERE = BENCHMARK_ROOT
MEMORY = 128 * 2**20
WORK = 10**10
SECONDS = 0.2


def allowance(module, *, timed=False):
    return module.Budget(
        work_limit=WORK,
        seconds=SECONDS if timed else 30,
        cpu_seconds=SECONDS if timed else 30,
    )


def build_corpus(split):
    """Generate proof-backed inputs independently after freezing policies."""
    generator = random.Random(3803042026 + (split == "held_out"))
    certificates, fixtures, seen = {}, [], set()

    for band, bits in (
        ("small", 13),
        ("medium", 21),
        ("20d", 34),
        ("30d", 50),
    ):
        for index in range(2 if split == "training" else 3):
            for _ in range(10000):
                p = certified_prime(bits, generator, certificates)
                q = certified_prime(bits, generator, certificates)
                if p != q and p * q not in seen:
                    break
            else:
                raise RuntimeError("fixture generation cap exceeded")

            seen.add(p * q)
            fixtures.append(
                dict(
                    id=f"r3_{band}_{split}_{index}",
                    band=band,
                    split=split,
                    n=p * q,
                    factors=sorted((p, q)),
                )
            )

    verify_certificates(certificates)
    return dict(
        schema=1,
        split=split,
        fixtures=fixtures,
        certificates=certificates,
        seeds=[7, 29],
        sampling="Certified Pocklington primes have a large p-1 "
        "factor; these bounded cohorts are not RSA samples.",
    )


def load_control(directory):
    """Materialize only hash-checked owned source bytes in a temporary root."""
    data = checked_baseline(HERE / "inputs/baselines/p38_r3_baseline.json")
    for name, source in data["source"].items():
        path = Path(directory) / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(source)
    package = types.ModuleType("_p38_r3_control")
    package.__path__ = [str(Path(directory) / "v2")]
    package.__package__ = package.__name__
    sys.modules[package.__name__] = package
    importlib.import_module(package.__name__ + ".qs")
    return package


def modules(package):
    return {
        name: runtime_module(name, package.__name__)
        for name in (
            "budget",
            "qs.pipeline",
            "qs.extraction",
            "qs.linear_algebra",
            "qs.siqs",
            "qs.factor_base",
            "qs.polynomial",
            "qs.sieve_collector",
        )
    }


@lru_cache(maxsize=2)
def rebuilt_filter(module):
    """Load only the old repeated-rebuild function, retaining exact lifting."""
    path = HERE / "inputs/baselines/qs_m26_baseline.json"
    data = json.loads(source_path(path).read_text())
    key = "v2/qs/linear_algebra.py"
    source = data["source"][key]
    if (
        hashlib.sha256(source.encode()).hexdigest()
        != data["source_sha256"][key]
    ):
        raise ValueError("corrupt repeated-rebuild control")

    node = next(
        n
        for n in ast.parse(source).body
        if isinstance(n, ast.FunctionDef) and n.name == "filter_matrix"
    )
    namespace = {}
    exec(
        compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"),
        module.__dict__,
        namespace,
    )
    return namespace["filter_matrix"]


def row_fixture(kind, count, seed):
    """Freeze graph/rank shapes without deriving rows from known factors."""
    if kind in ("cycle", "wide_cycle"):
        stride = 97 if kind == "wide_cycle" else 1
        return tuple(
            (1 << (i * stride)) | (1 << (((i + 1) % count) * stride))
            for i in range(count)
        )

    if kind == "pivots":
        return tuple(1 << i for i in range(count) for _ in range(3))
    if kind == "cascade":
        return tuple((1 << i) | (1 << (i + 1)) for i in range(count - 1)) + (
            1 << (count - 1),
        )
    generator = random.Random(seed)
    if kind == "live_gaps":
        live = tuple(
            sum(1 << (32 * i) for i in generator.sample(range(16), 4))
            for _ in range(count // 2)
        )
        # The initial union is dense; singleton elimination leaves wide gaps.
        return live + tuple(1 << i for i in range(512) if i % 32)
    return tuple(
        sum(1 << i for i in generator.sample(range(count // 2), 8))
        for _ in range(count)
    )


def rank_oracle(rows):
    """Independent set-valued column elimination gives expected nullity."""
    pivots = {}

    for row in rows:
        columns = set()
        bits = row
        while bits:
            bit = bits & -bits
            columns.add(bit.bit_length() - 1)
            bits ^= bit
        while columns:
            column = min(columns)
            if column not in pivots:
                pivots[column] = columns
                break
            columns = columns.symmetric_difference(pivots[column])

    return len(pivots)


def stable(samples):
    times = [s["seconds"] for s in samples]
    median = statistics.median(times)
    q1, _, q3 = statistics.quantiles(times, n=4)
    drift = abs(
        statistics.median(times[:3]) / statistics.median(times[-3:]) - 1
    )
    spread = (q3 - q1) / median
    return dict(
        median_seconds=median,
        drift=drift,
        relative_iqr=spread,
        stable=drift <= 0.15 and spread <= 0.2,
        samples=samples,
    )


def paired_measure(calls):
    """Interleave arms, validate every warmup/sample and extend noisy sets."""
    attempts = []

    for attempt in range(3):
        warmups = {}

        for name, call in calls.items():
            started, count = time.perf_counter(), 0
            while time.perf_counter() - started < (3 if not attempt else 5):
                call()
                count += 1
            warmups[name] = dict(
                seconds=time.perf_counter() - started, validated_calls=count
            )

        samples = {name: [] for name in calls}

        for index in range(9 if not attempt else 15):
            names = list(calls)
            if index % 2:
                names.reverse()
            for name in names:
                started = time.perf_counter()
                result = calls[name]()
                samples[name].append(
                    dict(seconds=time.perf_counter() - started, result=result)
                )

        summaries = {name: stable(values) for name, values in samples.items()}
        attempts.append(dict(warmups=warmups, arms=summaries))
        print(
            "sample gate",
            attempt + 1,
            {
                name: (round(s["median_seconds"], 6), s["stable"])
                for name, s in summaries.items()
            },
            flush=True,
        )
        if all(s["stable"] for s in summaries.values()):
            break

    return dict(
        stable=all(s["stable"] for s in summaries.values()), attempts=attempts
    )


def comparison(measurement, control, candidate):
    arms = measurement["attempts"][-1]["arms"]
    a, b = arms[control]["samples"], arms[candidate]["samples"]
    generator = random.Random(380303)
    gains = []

    for _ in range(999):
        indices = [generator.randrange(len(a)) for _ in a]
        gains.append(
            1
            - statistics.median(b[i]["seconds"] for i in indices)
            / statistics.median(a[i]["seconds"] for i in indices)
        )

    gains.sort()
    deltas = [
        sum(row["complete"] for row in second["result"])
        / len(second["result"])
        - sum(row["complete"] for row in first["result"])
        / len(first["result"])
        for first, second in zip(a, b)
    ]
    completion = sorted(
        statistics.mean(generator.choice(deltas) for _ in deltas)
        for _ in range(999)
    )
    return dict(
        completion_difference=statistics.mean(deltas),
        completion_ci95=[completion[24], completion[974]],
        uncertainty_scope="Paired repeats conditional on this fixed cohort.",
        median_reduction=1
        - arms[candidate]["median_seconds"] / arms[control]["median_seconds"],
        ci95=[gains[24], gains[974]],
        stable=arms[control]["stable"] and arms[candidate]["stable"],
    )


def matrix_run(module, rows, variant, expected_rank, prepared=None):
    budget = allowance(
        runtime_module("budget", module.__package__.rsplit(".", 1)[0])
    )
    retained = prepared.workspace_bytes if prepared is not None else 0
    remaining = MEMORY - retained
    if variant == "rebuild":
        matrix = rebuilt_filter(module)(
            rows, weight_two=True, budget=budget, memory_bytes=remaining
        )
    elif variant.startswith(("history", "dense_batch")):
        matrix = history_filter(
            rows,
            module,
            budget,
            remaining,
            batch_size=32 if variant.endswith("32") else 1,
            use_history=variant.startswith("history"),
        )
    else:
        matrix = module.filter_matrix(
            rows, weight_two=True, budget=budget, memory_bytes=remaining
        )
        if variant == "live":
            matrix = live_compaction(matrix, budget, remaining)

    solver = module.DependencySolver(matrix, budget=budget)

    dependencies = solver.run()
    for mask in dependencies:
        module.verify_dependency(mask, rows)
    if len(dependencies) != len(rows) - expected_rank:
        raise AssertionError("filtered kernel dimension differs from oracle")
    congruences = []
    if prepared is not None:
        from ....qs.extraction import extract_dependency

        congruences = [
            extract_dependency(prepared, mask, budget=budget)
            for mask in dependencies
        ]

    if budget.used > WORK or retained + matrix.workspace_bytes > MEMORY:
        raise AssertionError("matrix allowance exceeded")
    return dict(
        dependencies=len(dependencies),
        work=budget.used,
        workspace=matrix.workspace_bytes,
        simultaneous_workspace=retained + matrix.workspace_bytes,
        stats=matrix.stats,
        xors=solver.xors,
        pivot_nonzeros=solver.peak_nonzeros,
        extracted_dependencies=len(congruences),
    )


def config(module, band):
    if band == "30d":
        return module.SIQSConfig(
            base_bound=10000,
            half_width=8192,
            max_half_width=8192,
            factor_count=4,
            family_count=64,
            pool_size=32,
            max_stalled=8192,
            max_trivial=4096,
            row_excess=32,
            batch_width=4096,
            memory_bytes=MEMORY,
            checkpoint_bytes=2**20,
            collector=module.SieveConfig(
                block_width=4096,
                score_policy="powers",
                division="bucket",
                residual_bound=10000000,
                max_atoms=65536,
                max_relations=32768,
                max_partials=32768,
            ),
        )

    bound, width, count = {
        "small": (200, 256, 1),
        "medium": (1000, 512, 3),
        "20d": (3000, 2048, 3),
        "30d": (3000, 4096, 5),
    }[band]
    return module.SIQSConfig(
        base_bound=bound,
        half_width=width,
        max_half_width=width,
        factor_count=count,
        family_count=64,
        pool_size=32,
        max_stalled=4096,
        max_trivial=4096,
        memory_bytes=MEMORY,
        collector=module.SieveConfig(
            score_policy="powers",
            division="bucket",
            residual_bound=bound**2,
            max_atoms=8192,
            max_relations=4096,
            max_partials=2048,
        ),
    )


def factor_run(mods, fixture, seed, variant):
    module = mods["qs.siqs"]
    seconds = 2 if fixture["band"] == "30d" else 0.2
    budget = mods["budget"].Budget(
        work_limit=WORK, seconds=seconds, cpu_seconds=seconds
    )
    selected = config(module, fixture["band"])
    if variant.startswith("cadence"):
        selected = replace(selected, filter_row_growth=int(variant[7:]))
    if variant == "dependencies":
        selected = replace(selected, tested_dependencies=True)
    job = module.SIQSJob(
        fixture["n"], seed=seed, config=selected, budget=budget
    )

    result = job.run()

    if (result.divisor or 1) * result.cofactor != fixture["n"]:
        raise AssertionError("factor/cofactor reconstruction failed")
    if result.divisor is not None and (
        result.divisor not in fixture["factors"]
        or not 1 < result.divisor < fixture["n"]
    ):
        raise AssertionError("improper divisor or independent proof mismatch")

    if budget.used > WORK or result.stats.get("workspace_bytes", 0) > MEMORY:
        raise AssertionError("factor allowance exceeded")
    labels, reason = [], result.reason
    if result.divisor:
        try:
            budget.consume(
                sum(
                    p.bit_length() ** 2
                    for p in (result.divisor, result.cofactor)
                )
            )
            labels = [
                module.utils.classify_prime(p).value
                for p in (result.divisor, result.cofactor)
            ]
            budget.consume(0)
        except mods["budget"].BudgetExhaustedError:
            labels = []
            reason = "classification_" + budget.reason
    return dict(
        id=fixture["id"],
        seed=seed,
        complete=bool(labels),
        reason=reason,
        divisor=result.divisor,
        cofactor=result.cofactor,
        certainty=labels,
        work=budget.used,
        wall=budget.wall_used,
        cpu=budget.cpu_used,
        stats=result.stats,
    )


def prepare_fixture(mods, fixture):
    budget = allowance(mods["budget"])
    base = (
        mods["qs.factor_base"]
        .build_factor_base(fixture["n"], bound=400, budget=budget)
        .factor_base
    )
    if base is None:
        raise ValueError("profile input split during base setup")
    polynomial = mods["qs.polynomial"].qs_polynomial(base)
    collector = mods["qs.sieve_collector"].SieveCollector(
        polynomial,
        base,
        budget=budget,
        config=mods["qs.sieve_collector"].SieveConfig(
            residual_bound=1000,
            max_atoms=8192,
            max_relations=4096,
            max_partials=2048,
            memory_bytes=MEMORY,
        ),
    )

    result = collector.collect(-2048, 2048)
    return (
        base,
        result.full_relations + result.combined_relations,
        {atom.relation_id: atom for atom in result.atoms},
    )


def prepare_run(mods, fixture_data, cached):
    base, relations, atoms = fixture_data
    extraction = mods["qs.extraction"]
    cache = extraction._VerificationCache(2 * 2**20) if cached else None
    costs, last = [], None
    budget = allowance(mods["budget"])
    cache_reserve = cache.memory_bytes if cache else 0

    for count in sorted(
        set(min(len(relations), n) for n in (8, 16, 32, 64, 128))
    ):
        last = None
        last = extraction._prepare_relations(
            relations[:count],
            base,
            atoms,
            budget=budget,
            memory_bytes=MEMORY - cache_reserve,
            verification_cache=cache,
        )
        costs.append(
            dict(count=count, work=budget.used, workspace=last.workspace_bytes)
        )

    matrix = mods["qs.linear_algebra"].filter_matrix(
        last.rows,
        budget=budget,
        memory_bytes=MEMORY - last.workspace_bytes - cache_reserve,
    )

    dependencies = (
        mods["qs.linear_algebra"].DependencySolver(matrix, budget=budget).run()
    )
    for mask in dependencies:
        extraction.extract_dependency(last, mask, budget=budget)
    if budget.used > WORK:
        raise AssertionError("preparation pipeline work cap exceeded")
    return dict(
        total_work=budget.used,
        simultaneous_workspace=(
            last.workspace_bytes + cache_reserve + matrix.workspace_bytes
        ),
        costs=costs,
        dependencies=len(dependencies),
        cache_hits=cache.hits if cache else 0,
        cache_used=cache.used if cache else 0,
    )


def cold_run(root, mods, fixture):
    """Measure new PyPy startup, imports, setup and factoring separately."""
    code = """
import json, sys
sys.path.insert(0, sys.argv[1])
from v2.execution.budget import Budget
from v2.qs.siqs import SIQSConfig, SIQSJob
from v2.qs.sieve_collector import SieveConfig
values=json.loads(sys.argv[2])
values['collector']=SieveConfig(**values['collector'])
budget=Budget(work_limit=10**10, seconds=.2, cpu_seconds=.2)
job=SIQSJob(int(sys.argv[3]), seed=7,
            config=SIQSConfig(**values), budget=budget)
result=job.run()
print(json.dumps(dict(divisor=result.divisor, cofactor=result.cofactor,
                      reason=result.reason, work=budget.used)))
"""
    samples = []

    for _ in range(9):
        started = time.perf_counter()

        process = subprocess.run(
            [
                sys.executable,
                "-B",
                "-c",
                code,
                str(root),
                json.dumps(asdict(config(mods["qs.siqs"], fixture["band"]))),
                str(fixture["n"]),
            ],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        row = json.loads(process.stdout)

        if (row["divisor"] or 1) * row["cofactor"] != fixture["n"]:
            raise AssertionError("cold result reconstruction")
        if row["divisor"] and row["divisor"] not in fixture["factors"]:
            raise AssertionError("cold split disagrees with proof")
        samples.append(dict(seconds=time.perf_counter() - started, result=row))

    return dict(
        scope="Cold first-split process/import/setup, without a warm claim.",
        samples=samples,
        median_seconds=statistics.median(row["seconds"] for row in samples),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase",
        choices=(
            "build",
            "profile",
            "matrix",
            "prepare",
            "cold",
            "training",
            "confirmation",
        ),
        required=True,
    )
    parser.add_argument(
        "--split", choices=("training", "held_out"), default="training"
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--variants", default="control,current,dependencies,cadence8,cadence32"
    )
    parser.add_argument("--bands", default="small,medium,20d,30d")
    parser.add_argument("--short-repeats", type=int, default=5)
    args = parser.parse_args()
    if not 1 <= args.short_repeats <= 20:
        parser.error("short repeats must be between 1 and 20")
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    if args.output.exists():
        parser.error("preserve captures; choose a new output path")
    if args.phase == "build":
        corpus = build_corpus(args.split)
        if args.split == "held_out":
            corpus["frozen_sha256"] = hashlib.sha256(
                (
                    source_path(HERE / "inputs/controls/p38_r3_frozen.json")
                ).read_bytes()
            ).hexdigest()
        args.output.write_text(json.dumps(corpus, indent=2) + "\n")
        return

    corpus_path = (
        HERE / "inputs/corpora" / ("p38_r3_" + args.split + "_corpus.json")
    )
    corpus = json.loads(source_path(corpus_path).read_text())
    verify_corpus(corpus)
    if (
        args.split == "held_out"
        and corpus["frozen_sha256"]
        != hashlib.sha256(
            (
                source_path(HERE / "inputs/controls/p38_r3_frozen.json")
            ).read_bytes()
        ).hexdigest()
    ):
        raise ValueError("held-out corpus does not match frozen policies")

    before = fingerprint(ROOT)
    with tempfile.TemporaryDirectory(prefix="factor-r3-") as directory:
        control = modules(load_control(directory))
        current = modules(importlib.import_module("v2"))
        results = {}
        if args.phase == "profile":
            profile_fixtures = corpus["fixtures"][:4] + [
                next(f for f in corpus["fixtures"] if f["band"] == band)
                for band in ("20d", "30d")
            ]

            for fixture in profile_fixtures:
                data = prepare_fixture(control, fixture)
                profiler = cProfile.Profile()
                profiler.runcall(prepare_run, control, data, True)
                profiler.dump_stats(
                    str(args.output) + "." + fixture["id"] + ".prof"
                )
                results[fixture["id"] + "/prepare"] = prepare_run(
                    control, data, True
                )
                profiler = cProfile.Profile()
                results[fixture["id"] + "/factor"] = profiler.runcall(
                    factor_run, control, fixture, 7, "current"
                )
                profiler.dump_stats(
                    str(args.output) + "." + fixture["id"] + ".factor.prof"
                )
        elif args.phase == "cold":
            for fixture in corpus["fixtures"][:2]:
                results[fixture["id"]] = {
                    "control": cold_run(directory, control, fixture),
                    "current": cold_run(ROOT, current, fixture),
                }
        elif args.phase == "matrix":
            kinds = (
                "cycle",
                "cascade",
                "dense",
                "pivots",
                "live_gaps",
                "wide_cycle",
                "verified_store",
            )

            for kind in kinds:
                prepared = None
                if kind == "verified_store":
                    base, relations, atoms = prepare_fixture(
                        current, corpus["fixtures"][0]
                    )
                    prepared = current["qs.extraction"].prepare_relations(
                        relations[:128],
                        base,
                        atoms,
                        budget=allowance(current["budget"]),
                        memory_bytes=MEMORY,
                    )
                    rows = prepared.rows
                else:
                    rows = row_fixture(kind, 512, 3803)

                expected_rank = rank_oracle(rows)
                variants = ["rebuild", "control", "current", "live"]
                if kind != "wide_cycle":
                    variants += [
                        "dense_batch1",
                        "dense_batch32",
                        "history1",
                        "history32",
                    ]

                calls = {}

                for variant in variants:
                    module = (
                        control["qs.linear_algebra"]
                        if variant == "control"
                        else current["qs.linear_algebra"]
                    )
                    calls[variant] = (
                        lambda variant=variant, module=module: matrix_run(
                            module, rows, variant, expected_rank, prepared
                        )
                    )

                results[kind] = paired_measure(calls)
        elif args.phase == "prepare":
            for fixture in corpus["fixtures"][:2]:
                fixtures = {
                    "control": prepare_fixture(control, fixture),
                    "current": prepare_fixture(current, fixture),
                }
                calls = {
                    "uncached": lambda: prepare_run(
                        control, fixtures["control"], False
                    ),
                    "control": lambda: prepare_run(
                        control, fixtures["control"], True
                    ),
                    "current": lambda: prepare_run(
                        current, fixtures["current"], True
                    ),
                }
                results[fixture["id"]] = paired_measure(calls)
        else:
            for band in args.bands.split(","):
                cohort = [f for f in corpus["fixtures"] if f["band"] == band]
                calls = {}

                for variant in args.variants.split(","):
                    mods = control if variant == "control" else current
                    calls[variant] = lambda variant=variant, mods=mods: [
                        factor_run(mods, fixture, seed, variant)
                        for _ in range(
                            args.short_repeats if band == "small" else 1
                        )
                        for fixture in cohort
                        for seed in corpus["seeds"]
                    ]

                results[band] = paired_measure(calls)
                results[band]["comparisons"] = {
                    candidate: comparison(results[band], "control", candidate)
                    for candidate in calls
                    if candidate != "control"
                }

        if fingerprint(ROOT) != before:
            raise AssertionError("runtime changed during capture")
        args.output.write_text(
            json.dumps(
                dict(
                    phase=args.phase,
                    split=args.split,
                    python=platform.python_version(),
                    implementation=platform.python_implementation(),
                    source_sha256=before,
                    driver_sha256={
                        name: hashlib.sha256(
                            (source_path(HERE / name)).read_bytes()
                        ).hexdigest()
                        for name in ("p38_r3.py", "p38_r3_experiments.py")
                    },
                    control_sha256=hashlib.sha256(
                        (
                            source_path(
                                HERE / "inputs/baselines/p38_r3_baseline.json"
                            )
                        ).read_bytes()
                    ).hexdigest(),
                    corpus_sha256=hashlib.sha256(
                        source_path(corpus_path).read_bytes()
                    ).hexdigest(),
                    budgets=dict(
                        work=WORK,
                        wall=SECONDS,
                        cpu=SECONDS,
                        larger_wall_cpu=2,
                        matrix_preparation_wall_cpu=30,
                        owned_bytes=MEMORY,
                    ),
                    short_repeats=args.short_repeats,
                    scope="Profiles diagnose costs. Warmed samples include "
                    "complete calls and validation. Censored times are not "
                    "successful-factorization speedups. No worker overlap.",
                    results=results,
                ),
                indent=2,
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
