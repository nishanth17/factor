"""Frozen, validated R2 candidate-cascade comparisons on PyPy 3.11."""

import argparse
import cProfile
import hashlib
import importlib
import json
import platform
import random
import resource
import subprocess
import sys
import tempfile
import time
import types
from dataclasses import asdict, replace
from pathlib import Path

from .build_phase_two_corpus import certified_prime, verify_certificates
from .p38_r2_experiments import experiment_collector
from .p38_r3 import comparison, config, modules, paired_measure
from .performance_audit import checked_baseline, fingerprint, verify_corpus

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
MEMORY = 128 * 2**20
WORK = 10**10
VARIANTS = (
    "control",
    "current",
    "conservative",
    "fixed",
    "plans",
    "fixed_plans",
    "cutoff",
    "resieve",
    "tiny",
    "batch",
    "chunks",
)


def build_corpus(split):
    """Generate independent proof-backed cohorts, after freeze for held-out."""
    generator = random.Random(3802042026 + (split == "held_out"))
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
                    id=f"r2_{band}_{split}_{index}",
                    band=band,
                    n=p * q,
                    factors=sorted((p, q)),
                    split=split,
                )
            )
    verify_certificates(certificates)
    corpus = dict(
        schema=1,
        split=split,
        fixtures=fixtures,
        certificates=certificates,
        seeds=[7, 29],
        sampling="Certified Pocklington primes have a large "
        "p-1 factor; these bounded balanced cohorts are not RSA samples.",
    )
    if split == "held_out":
        corpus["frozen_sha256"] = digest(HERE / "p38_r2_frozen.json")
    return corpus


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_control(directory):
    data = checked_baseline(HERE / "p38_r2_baseline.json")
    for name, source in data["source"].items():
        path = Path(directory) / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(source)
    package = types.ModuleType("_p38_r2_control")
    package.__path__ = [str(Path(directory) / "v2")]
    package.__package__ = package.__name__
    sys.modules[package.__name__] = package
    importlib.import_module(package.__name__ + ".qs")
    return package


def selected_config(module, band, variant):
    selected = config(module, band)
    # R3's accepted scoped cadence avoids repeatedly measuring preparation
    # in the feasible nominal 30-digit collector comparison; same in all arms.
    if band == "30d":
        selected = replace(selected, filter_row_growth=32)
    changes = {"memory_bytes": MEMORY}
    if variant == "conservative":
        changes["score_policy"] = "candidate"
    if variant in ("fixed", "fixed_plans"):
        changes["score_policy"] = "fixed"
    if variant in ("plans", "fixed_plans"):
        changes["power_plan_bytes"] = 2**20
    if variant == "cutoff":
        changes["small_prime_cutoff"] = 5
    if variant == "resieve":
        changes["division"] = "resieve"
    return replace(selected, collector=replace(selected.collector, **changes))


def factor_run(mods, fixture, seed, variant):
    module = mods["qs.siqs"]
    budget = mods["budget"].Budget(
        work_limit=WORK,
        seconds=2 if fixture["band"] == "30d" else 0.2,
        cpu_seconds=2 if fixture["band"] == "30d" else 0.2,
    )
    collector = module.SieveCollector
    if variant in ("tiny", "batch", "chunks"):
        module.SieveCollector = experiment_collector(
            mods["qs.sieve_collector"], variant
        )
    try:
        job = module.SIQSJob(
            fixture["n"],
            seed=seed,
            config=selected_config(module, fixture["band"], variant),
            budget=budget,
        )
        result = job.run()
    finally:
        module.SieveCollector = collector
    if (result.divisor or 1) * result.cofactor != fixture["n"]:
        raise AssertionError("unresolved input reconstruction failed")
    if result.divisor is not None and (
        result.divisor not in fixture["factors"]
        or not 1 < result.divisor < fixture["n"]
    ):
        raise AssertionError("proper divisor/proof mismatch")
    if budget.used > WORK or result.stats.get("workspace_bytes", 0) > MEMORY:
        raise AssertionError("shared finite allowance exceeded")
    labels = []
    reason = result.reason
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
        divisor=result.divisor,
        cofactor=result.cofactor,
        certainty=labels,
        reason=reason,
        work=budget.used,
        wall=budget.wall_used,
        cpu=budget.cpu_used,
        stats=result.stats,
    )


def collector_run(mods, fixture, seed, variant):
    module = mods["qs.siqs"]
    selected = selected_config(module, fixture["band"], variant)
    budget = mods["budget"].Budget(work_limit=WORK, seconds=30, cpu_seconds=30)
    job = module.SIQSJob(
        fixture["n"], seed=seed, config=selected, budget=budget
    )
    job._setup()
    if job.divisor is not None:
        raise ValueError("collector fixture split during setup")
    polynomial, roots = job._next_step()
    cls = module.SieveCollector
    if variant in ("tiny", "batch", "chunks"):
        cls = experiment_collector(mods["qs.sieve_collector"], variant)
    try:
        worker = cls(
            polynomial,
            job.base,
            config=selected.collector,
            precomputed_roots=roots,
            budget=budget,
        )
    except MemoryError:
        return dict(
            id=fixture["id"],
            seed=seed,
            complete=False,
            reason="memory_limit",
            cofactor=fixture["n"],
            work=budget.used,
            stats={},
            rows=0,
            setup_refused=True,
        )
    # Several blocks exercise reuse; complete polynomial interval decides
    # eligibility for a future family-wide CRT arm, never block width.
    lo, hi = -selected.half_width, selected.half_width + 1
    result = worker.collect(lo, hi)
    if result.divisor is not None and result.divisor not in fixture["factors"]:
        raise AssertionError("collector returned an improper divisor")
    atoms = {a.relation_id: a for a in result.atoms}
    extraction = mods["qs.extraction"]
    prepared = extraction.prepare_relations(
        worker.matrix_relations,
        job.base,
        atoms,
        budget=budget,
        memory_bytes=MEMORY,
    )
    matrix = mods["qs.linear_algebra"].filter_matrix(
        prepared.rows,
        budget=budget,
        memory_bytes=MEMORY,
    )
    large_primes = sum(e.prime > hi - lo for e in job.base.entries)
    return dict(
        id=fixture["id"],
        seed=seed,
        complete=result.reason in ("complete", "factor_found"),
        reason=result.reason,
        stats=result.stats,
        atoms=len(result.atoms),
        rows=len(worker.matrix_relations),
        pending=len(result.partial_ids),
        post_filter=matrix.stats,
        work=budget.used,
        workspace=result.workspace_bytes,
        divisor=result.divisor,
        cofactor=fixture["n"] // (result.divisor or 1),
        filter_excess=len(matrix.rows) - matrix.stats["output_columns"],
        filter_kernels=len(matrix.zero_dependencies),
        full_interval=hi - lo,
        crt_eligible_primes=large_primes,
        signature=hashlib.sha256(
            repr(
                sorted(
                    (a.position, a.sign, a.exponents, a.residual)
                    for a in result.atoms
                )
            ).encode()
        ).hexdigest(),
    )


def cold_run(root, mods, fixture, variant):
    """Measure nine fresh PyPy processes including startup and validation."""
    options = asdict(
        selected_config(mods["qs.siqs"], fixture["band"], variant)
    )
    code = """import json,sys,time
started=time.perf_counter()
from v2.budget import Budget
from v2.qs.siqs import SIQSConfig,SIQSJob,SieveConfig
options=json.loads(sys.argv[1])
options['collector']=SieveConfig(**options['collector'])
n=int(sys.argv[2]);b=Budget(work_limit=10**10,seconds=2,cpu_seconds=2)
r=SIQSJob(n,seed=7,config=SIQSConfig(**options),budget=b).run()
print(json.dumps(dict(seconds=time.perf_counter()-started,reason=r.reason,
divisor=r.divisor,cofactor=r.cofactor)))
"""
    samples = []
    for _ in range(9):
        started = time.perf_counter()
        output = subprocess.run(
            [
                sys.executable,
                "-c",
                code,
                json.dumps(options),
                str(fixture["n"]),
            ],
            cwd=root,
            text=True,
            capture_output=True,
            check=True,
            timeout=30,
        )
        result = json.loads(output.stdout)
        result["runtime_seconds"] = result["seconds"]
        result["seconds"] = time.perf_counter() - started
        if (result["divisor"] or 1) * result["cofactor"] != fixture["n"]:
            raise AssertionError("cold result failed reconstruction")
        if result["divisor"] and result["divisor"] not in fixture["factors"]:
            raise AssertionError("cold divisor proof mismatch")
        samples.append(result)
    return samples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase",
        required=True,
        choices=("build", "collector", "profile", "factor", "cold"),
    )
    parser.add_argument(
        "--split", default="training", choices=("training", "held_out")
    )
    parser.add_argument("--bands", default="small,medium,20d,30d")
    parser.add_argument("--variants", default=",".join(VARIANTS))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    if args.output.exists():
        parser.error("choose a new capture path; overwrites are refused")
    if args.output.with_suffix(".partial.json").exists():
        parser.error(
            "choose a new capture path; partial captures are retained"
        )
    if args.phase == "build":
        args.output.write_text(
            json.dumps(build_corpus(args.split), indent=2) + "\n"
        )
        return
    variants = args.variants.split(",")
    if set(variants) - set(VARIANTS):
        parser.error("unknown variant")
    if args.phase == "cold" and set(variants) & {"tiny", "batch", "chunks"}:
        parser.error("cold mode accepts runtime policies only")
    corpus_path = HERE / f"p38_r2_{args.split}_corpus.json"
    corpus = json.loads(corpus_path.read_text())
    verify_corpus(corpus)
    if args.split == "held_out" and corpus["frozen_sha256"] != digest(
        HERE / "p38_r2_frozen.json"
    ):
        raise ValueError("held-out population does not match frozen policy")
    before = fingerprint(ROOT)
    driver_hashes = {
        n: digest(HERE / n)
        for n in (
            "p38_r2.py",
            "p38_r2_experiments.py",
            "p38_r3.py",
        )
    }
    if args.split == "held_out":
        frozen = json.loads((HERE / "p38_r2_frozen.json").read_text())
        if frozen["source_sha256"] != before:
            raise ValueError("confirmation runtime differs from frozen source")
        if frozen["driver_sha256"] != driver_hashes:
            raise ValueError("confirmation driver differs from frozen source")
    results = {}
    with tempfile.TemporaryDirectory(prefix="factor-r2-") as directory:
        control = modules(load_control(directory))
        current = modules(importlib.import_module("v2"))
        for band in args.bands.split(","):
            fixtures = [f for f in corpus["fixtures"] if f["band"] == band]
            if not fixtures:
                parser.error("unknown or empty band")
            calls = {}
            for variant in variants:
                mods = control if variant == "control" else current
                runner = (
                    collector_run if args.phase == "collector" else factor_run
                )
                calls[variant] = lambda variant=variant, mods=mods: [
                    runner(mods, fixture, seed, variant)
                    for fixture in fixtures
                    for seed in corpus["seeds"]
                ]
            print("START", band, variants, flush=True)
            if args.phase == "cold":
                results[band] = {
                    variant: cold_run(
                        directory if variant == "control" else ROOT,
                        control if variant == "control" else current,
                        fixtures[0],
                        variant,
                    )
                    for variant in variants
                }
            elif args.phase == "profile":
                profiler = cProfile.Profile()
                results[band] = profiler.runcall(calls[variants[0]])
                profiler.dump_stats(str(args.output) + "." + band + ".prof")
            else:
                results[band] = paired_measure(calls)
                if args.phase == "factor" and "control" in calls:
                    results[band]["comparisons"] = {
                        name: comparison(results[band], "control", name)
                        for name in calls
                        if name != "control"
                    }
            print("DONE", band, flush=True)
            # Preserve each completed band even if a later phase is refused.
            partial = args.output.with_suffix(".partial.json")
            partial.write_text(
                json.dumps(
                    dict(
                        phase=args.phase,
                        source_sha256=before,
                        driver_sha256=driver_hashes,
                        results=results,
                    ),
                    indent=2,
                )
                + "\n"
            )
    if fingerprint(ROOT) != before:
        raise RuntimeError("runtime changed during capture")
    data = dict(
        phase=args.phase,
        split=args.split,
        source_sha256=before,
        control_sha256=digest(HERE / "p38_r2_baseline.json"),
        corpus_sha256=digest(corpus_path),
        driver_sha256=driver_hashes,
        implementation=platform.python_implementation(),
        python=platform.python_version(),
        results=results,
        budgets=dict(
            work=WORK,
            wall_cpu_small=0.2,
            wall_cpu_30d=2,
            collector_wall_cpu=30,
            owned_bytes=MEMORY,
        ),
        config={
            b: {
                v: asdict(
                    selected_config(
                        control["qs.siqs"]
                        if v == "control"
                        else current["qs.siqs"],
                        b,
                        v,
                    )
                )
                for v in variants
            }
            for b in args.bands.split(",")
        },
        max_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        * (1 if sys.platform == "darwin" else 1024),
        scope="Warmed complete calls include setup and validation. Censored "
        "times are cohort costs. Profiles are separate diagnostics. Owned "
        "workspace excludes interpreter/JIT RSS. CPU-heavy runs serialized.",
    )
    args.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
