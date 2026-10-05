"""Matched, validated P2/P3 cost controls against a frozen runtime checkout.

Invoke this file directly with --runtime-root to select an isolated runtime.
Cold pools, warmed execution and optional parent profiles are separate cases.
All timing gates require three seconds of validated warmup and nine samples;
unstable cases receive five seconds and fifteen samples, up to three attempts.
"""

import argparse
import ast
import cProfile
import hashlib
import importlib
import json
import platform
import random
import resource
import statistics
import sys
import tempfile
import time
from collections import Counter
from dataclasses import asdict, replace
from math import gcd, isqrt, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CORPUS = Path(__file__).parent / "inputs/performance_audit_corpus.json"
BASELINE = Path(__file__).parent / "inputs/performance_audit_baseline.json"


def checked_baseline(path):
    """Validate owned immutable source bytes before loading a control."""
    data = json.loads(path.read_text())
    for name, source in data["source"].items():
        parts = Path(name).parts
        if not (
            parts[0] == "v2"
            and len(parts) in (2, 3)
            and (len(parts) == 2 or parts[1] == "qs")
            and parts[-1].endswith(".py")
            and all(part not in (".", "..") for part in parts)
        ):
            raise ValueError("invalid baseline module path")
        if (
            hashlib.sha256(source.encode()).hexdigest()
            != data["source_sha256"][name]
        ):
            raise ValueError("corrupt baseline source")
    return data


def materialize_baseline(path, directory):
    for name, source in checked_baseline(path)["source"].items():
        destination = directory / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(source)


def apply_revert(variant):
    """Replace one diagnosed component with its owned exact control."""
    if variant == "none":
        return
    changes = {
        "preprocessing": (
            ("preprocessing", None, "integer_root"),
            ("portfolio", None, "_advance"),
        ),
        "preparation": (("qs.pipeline", "QSJob", "_solve"),),
        "serial_poll": (("qs.pipeline", "QSJob", "run"),),
        "sssf_trees": (("qs.sss", "SSSCollector", "__init__"),),
        "segments": (("stage_jobs", None, "peek_prime"),),
        "collisions": (("qs.sss", None, "collision_candidates"),),
        "division": (
            ("qs.sieve_collector", "SieveCollector", "_divide"),
            ("qs.sieve_collector", "SieveCollector", "_resieve"),
        ),
        "scoring": (
            ("qs.sieve_collector", "SieveCollector", "_sieve"),
            ("qs.sieve_collector", "SieveCollector", "_candidate_passes"),
        ),
        "matrix": (
            ("qs.linear_algebra", None, "filter_matrix"),
            ("qs.linear_algebra", None, "verify_dependency"),
            ("qs.extraction", None, "extract_dependency"),
        ),
    }
    data = checked_baseline(BASELINE)
    for name, owner, function in changes[variant]:
        source = data["source"]["v2/" + name.replace(".", "/") + ".py"]
        tree = ast.parse(source)
        nodes = (
            tree.body
            if owner is None
            else next(
                node.body
                for node in tree.body
                if isinstance(node, ast.ClassDef) and node.name == owner
            )
        )
        node = next(
            node
            for node in nodes
            if isinstance(node, ast.FunctionDef) and node.name == function
        )
        if owner is not None:
            # Standalone AST functions have no implicit __class__ cell.
            # Bind super explicitly to the same runtime owner and instance.
            for call in ast.walk(node):
                if (
                    isinstance(call, ast.Call)
                    and isinstance(call.func, ast.Name)
                    and call.func.id == "super"
                    and not call.args
                    and not call.keywords
                ):
                    call.args = [
                        ast.Name(id=owner, ctx=ast.Load()),
                        ast.Name(id="self", ctx=ast.Load()),
                    ]
            ast.fix_missing_locations(node)
        module = importlib.import_module("v2." + name)
        module.__dict__.setdefault("Counter", Counter)
        namespace = {}
        exec(
            compile(
                ast.Module(body=[node], type_ignores=[]),
                str(BASELINE) + ":" + name,
                "exec",
            ),
            module.__dict__,
            namespace,
        )
        setattr(
            module if owner is None else getattr(module, owner),
            function,
            namespace[function],
        )
    if variant == "preprocessing":
        sys.modules["v2.portfolio"].integer_root = sys.modules[
            "v2.preprocessing"
        ].integer_root
    elif variant == "preparation":
        from v2.qs import extraction, pipeline

        pipeline.prepare_relations = extraction.prepare_relations
    elif variant == "segments":
        sys.modules["v2.portfolio"].peek_prime = sys.modules[
            "v2.stage_jobs"
        ].peek_prime
    elif variant == "matrix":
        sys.modules["v2.qs.pipeline"].filter_matrix = sys.modules[
            "v2.qs.linear_algebra"
        ].filter_matrix


def verify_corpus(corpus):
    """Check Pocklington proofs independently of library primality labels."""
    verified = set()

    def verify(n):
        if n in verified:
            return
        proof = corpus["certificates"][str(n)]
        if proof["kind"] == "trial":
            if not 2 <= n < 2**16 or any(
                n % d == 0 for d in range(2, isqrt(n) + 1)
            ):
                raise AssertionError("invalid trial certificate")
        else:
            q, witness = proof["q"], proof["witness"]
            if not 2 <= q < n or (n - 1) % q:
                raise AssertionError("invalid certificate dependency")
            verify(q)
            if not (
                q * q > n
                and pow(witness, n - 1, n) == 1
                and gcd(pow(witness, (n - 1) // q, n) - 1, n) == 1
            ):
                raise AssertionError("invalid Pocklington certificate")
        verified.add(n)

    seen = set()
    for fixture in corpus["fixtures"]:
        for factor in fixture["factors"]:
            verify(factor)
        if prod(fixture["factors"]) != fixture["n"] or fixture["n"] in seen:
            raise AssertionError("invalid or repeated input")
        seen.add(fixture["n"])


def fingerprint(root):
    """Identify runtime bytes, including new accounting helpers."""
    paths = sorted((root / "v2").glob("*.py"))
    paths += sorted((root / "v2/qs").glob("*.py"))
    return {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in paths
    }


def make_config(case, band, modules):
    """Fixed controls; candidate policies have explicit case names."""
    if case.startswith("p2"):
        if case == "p2":
            return dict(rho_attempts=0, pm1_attempts=0, ecm_tiers=())
        choices = {
            "p2_full": {},
            "p2_pm1": dict(rho_attempts=0, ecm_tiers=()),
            "p2_rho": dict(pm1_attempts=0, ecm_tiers=()),
            "p2_ecm": dict(rho_attempts=0, pm1_attempts=0),
            "p2_full_chunk8": dict(chunk_size=8),
            "p2_full_chunk32": dict(chunk_size=32),
            "p2_full_gcd64": dict(gcd_batch=64),
            "p2_full_gcd256": dict(gcd_batch=256),
            "p2_full_rho128": dict(rho_batch=128),
            "p2_full_cache": dict(schedule_cache_bytes=2**20),
        }
        return choices[case]
    sieve = modules["sieve_collector"].SieveConfig
    options = sieve(
        score_policy="powers",
        division="bucket",
        residual_bound=10000,
        max_atoms=8192,
        max_relations=4096,
        max_partials=1024,
        memory_bytes=64 * 2**20,
    )
    if case.startswith("sss"):
        mode = "sssf" if case.startswith("sssf") else "sss"
        filtered = mode == "sssf" and "unfiltered" not in case
        return modules["sss"].SSSConfig(
            mode=mode,
            memory_bytes=(
                128 if "feasible" in case or "memory128" in case else 64
            )
            * 2**20,
            base_bound=(
                {"40d": 8000, "60d": 82000, "80d": 82000}[band]
                if "capacity" in case
                else (400 if band == "small" else 5000)
            ),
            search_rounds=128 if band == "small" else 4096,
            selection_size=7 if mode == "sssf" and "six" not in case else 6,
            filter_divisor=2 if mode == "sssf" else 10,
            filter_bound=(1000 if band == "small" else 10**14)
            if filtered
            else 0,
            collector=replace(
                options,
                score_policy="adaptive",
                division="roots",
                residual_bound=25_000_000 if "residual" in case else 10000,
            ),
        )
    if case.startswith(("parallel", "native")):
        params = dict(
            base_bound=(
                1000
                if "base1k" in case
                else (
                    10000
                    if "base10k" in case
                    else (200 if band == "small" else 400)
                )
            ),
            half_width=256 if band == "small" else 512,
            factor_count=1 if band == "small" else 3,
            family_count=4,
            pool_size=16,
            max_batch_atoms=512,
            assignment_work=50_000_000,
            collector=options,
        )
        if "chunk" in case:
            params["batch_width"] = 128 if "chunk128" in case else 256
        if "strict" in case:
            params["poll_interval"] = 1
        if "cap" in case:
            params["max_batch_atoms"] = 1024
        return modules["parallel"].ParallelConfig(**params)
    if case.startswith("siqs"):
        if "feasible" in case:
            return modules["siqs"].SIQSConfig(
                base_bound=5000 if "base5k" in case else 10000,
                multiplier=0 if "multiplier" in case else 1,
                half_width=2048 if "width2k" in case else 8192,
                max_half_width=8192,
                factor_count=5 if "five" in case else 4,
                family_count=64,
                pool_size=32,
                max_stalled=8192,
                max_trivial=4096,
                row_excess=32,
                batch_width=4096,
                memory_bytes=128 * 2**20,
                checkpoint_bytes=2**20,
                collector=replace(
                    options,
                    block_width=4096,
                    residual_bound=10000000,
                    max_atoms=65536,
                    max_relations=32768,
                    max_partials=32768,
                ),
            )
        return modules["siqs"].SIQSConfig(
            base_bound=200 if band == "small" else 5000,
            half_width=256 if band == "small" else 512,
            factor_count=1 if band == "small" else 4,
            family_count=64,
            pool_size=64,
            max_stalled=64,
            max_trivial=4096,
            memory_bytes=64 * 2**20,
            collector=options,
        )
    if case.startswith("collector"):
        return replace(
            options,
            division=case.split("_")[1],
            small_prime_cutoff=5 if "cutoff" in case else 0,
        )
    if case in ("collision", "matrix"):
        return None
    raise ValueError("unknown audit case: " + case)


def run_one(case, fixture, seed, config, modules, pool=None):
    """Include setup/classification and validate unresolved results."""
    utils = modules["utils"]
    n = fixture["n"]
    seconds = (
        30
        if "feasible" in case
        else (0.2 if fixture["band"] in ("40d", "60d", "80d") else 5)
    )
    budget = modules["budget"].Budget(
        work_limit=10**13 if "feasible" in case else 2_000_000_000,
        seconds=seconds,
        cpu_seconds=seconds,
    )
    began = time.perf_counter()
    divisor, stats, signature = None, {}, None
    factors, remaining, labels = [], [n], []
    try:
        if case.startswith("p2"):
            portfolio = modules["portfolio"]
            outcome = portfolio.factorize_bounded(
                n,
                seed=seed,
                budget=budget,
                config=portfolio.PortfolioConfig(**config),
            )
            if outcome.result.reconstruct() != n:
                raise AssertionError("portfolio reconstruction failed")
            reason = outcome.reason
            stats = dict(
                events=outcome.events,
                stage_seconds=outcome.checkpoint["payload"]["state"][
                    "stage_seconds"
                ],
            )
            factors = sorted(
                f.value
                for f in outcome.result.factors
                for _ in range(f.exponent)
            )
            remaining = list(outcome.result.remaining)
            labels = [
                f.certainty.value
                for f in outcome.result.factors
                for _ in range(f.exponent)
            ]
        elif case == "matrix":
            rows = (1 << 99999,) * 32
            matrix = modules["linear_algebra"].filter_matrix(
                rows, weight_two=True, budget=budget, memory_bytes=64 * 2**20
            )
            masks = (
                modules["linear_algebra"]
                .DependencySolver(matrix, budget=budget)
                .run()
            )
            for mask in masks:
                modules["linear_algebra"].verify_dependency(mask, rows)
            stats, reason, signature = (
                matrix.stats,
                "kernel_verified",
                str(masks),
            )
        elif case.startswith("collector") or case == "collision":
            base = modules["factor_base"].build_factor_base(
                n, bound=1000, budget=budget, memory_bytes=64 * 2**20
            )
            divisor = base.divisor
            reason = "factor_found" if divisor else "complete"
            if not divisor:
                polynomial = modules["polynomial"].qs_polynomial(
                    base.factor_base
                )
                if case == "collision":
                    collector = modules["sss"].SSSCollector(
                        polynomial,
                        base.factor_base,
                        search_config=modules["sss"].SSSConfig(),
                        budget=budget,
                        seed=seed,
                    )
                    candidates = collector.assignment(7)
                    signature = str(candidates)
                    stats = dict(candidates=len(candidates))
                else:
                    collector = modules["sieve_collector"].SieveCollector(
                        polynomial,
                        base.factor_base,
                        config=config,
                        budget=budget,
                    )
                    result = collector.collect(-256, 257)
                    divisor, reason, stats = (
                        result.divisor,
                        result.reason,
                        result.stats,
                    )
                    signature = str(
                        tuple(
                            sorted(
                                (a.position, a.sign, a.exponents, a.residual)
                                for a in result.atoms
                            )
                        )
                    )
        else:
            if case.startswith("parallel"):
                job = modules["parallel"].ParallelSIQSJob(
                    n, seed=seed, budget=budget, config=config
                )
                result = job.run(pool=pool, fixed_work="fixed" in case)
                if "fixed" in case and job.engine is not None:
                    signature = str(
                        tuple(
                            sorted(
                                (
                                    a.relation_id,
                                    a.sign,
                                    a.exponents,
                                    a.residual,
                                )
                                for a in job.engine.collector._atoms.values()
                            )
                        )
                    )
            elif case.startswith("native"):
                native = modules["siqs"].SIQSConfig(
                    base_bound=config.base_bound,
                    half_width=config.half_width,
                    factor_count=config.factor_count,
                    family_count=config.family_count,
                    pool_size=config.pool_size,
                    memory_bytes=config.parent_memory_bytes,
                    max_stalled=64,
                    max_trivial=4096,
                    collector=config.collector,
                )
                result = (
                    modules["siqs"]
                    .SIQSJob(n, seed=seed, budget=budget, config=native)
                    .run()
                )
            elif case.startswith("siqs"):
                result = (
                    modules["siqs"]
                    .SIQSJob(n, seed=seed, budget=budget, config=config)
                    .run()
                )
            else:
                result = (
                    modules["sss"]
                    .SSSJob(n, seed=seed, budget=budget, config=config)
                    .run()
                )
            divisor, reason, stats = (
                result.divisor,
                result.reason,
                result.stats,
            )
            if (divisor or 1) * result.cofactor != n:
                raise AssertionError("job reconstruction failed")
    except modules["budget"].BudgetExhaustedError:
        reason = budget.reason
    except MemoryError as error:
        reason = "memory_limit"
        stats["refusal_detail"] = str(error)
    if divisor is not None:
        if not utils.valid_divisor(divisor, n):
            raise AssertionError("invalid split")
        children = sorted((divisor, n // divisor))
        try:
            for child in children:
                budget.consume(child.bit_length() * 32)
                label = utils.classify_prime(child, rng=random.Random(seed))
                if label is utils.Primality.COMPOSITE:
                    raise AssertionError("nonterminal semiprime child")
                labels.append(label.value)
            factors, remaining = children, []
        except modules["budget"].BudgetExhaustedError:
            reason, labels = "classification_" + budget.reason, []
    if prod(factors) * prod(remaining) != n or (
        not remaining and factors != fixture["factors"]
    ):
        raise AssertionError("result differs from the independent corpus")
    return dict(
        id=fixture["id"],
        seed=seed,
        seconds=time.perf_counter() - began,
        cpu_seconds=budget.cpu_used,
        work_used=budget.used,
        completed=not remaining,
        factors=factors,
        remaining=remaining,
        certainty=labels,
        reason=reason,
        stats=stats,
        worker_rss_high_water_bytes=(
            list(pool.state["rss"])
            if pool is not None and pool.mode == "process"
            else []
        ),
        process_rss_high_water_bytes=int(
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        )
        * (1 if sys.platform == "darwin" else 1024),
        signature=None
        if signature is None
        else hashlib.sha256(signature.encode()).hexdigest(),
    )


def measure(call, args):
    attempts = []
    for attempt in range(3):
        began, warm_calls = time.perf_counter(), 0
        while time.perf_counter() - began < max(
            args.warmup_seconds, 3 if not attempt else 5
        ):
            call()
            warm_calls += 1
        warm_seconds = time.perf_counter() - began
        samples = []
        for _ in range(max(args.repetitions, 9 if not attempt else 15)):
            started = time.perf_counter()
            rows = call()
            samples.append(
                dict(seconds=time.perf_counter() - started, rows=rows)
            )
        times = [sample["seconds"] for sample in samples]
        median = statistics.median(times)
        q1, _, q3 = statistics.quantiles(times, n=4)
        drift = abs(
            statistics.median(times[:3]) / statistics.median(times[-3:]) - 1
        )
        spread = (q3 - q1) / median
        stable = drift <= 0.15 and spread <= 0.2
        attempts.append(
            dict(
                warmup_seconds=warm_seconds,
                warmup_calls=warm_calls,
                median_seconds=median,
                drift=drift,
                relative_iqr=spread,
                stable=stable,
                samples=samples,
            )
        )
        print(
            "sample gate",
            attempt + 1,
            len(samples),
            round(median, 6),
            "stable" if stable else "extending",
            flush=True,
        )
        if stable:
            break
    return dict(stable=stable, median_seconds=median, attempts=attempts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-root", type=Path, default=ROOT)
    parser.add_argument("--baseline-json", type=Path)
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument(
        "--cases",
        default=(
            "p2,collision,collector_bucket,collector_resieve,"
            "collector_bucket_cutoff,siqs,sss,sssf,native,"
            "parallel_serial,parallel_chunk_serial"
        ),
    )
    parser.add_argument("--bands", default="small,medium,30d")
    parser.add_argument(
        "--split", choices=("training", "held_out"), default="held_out"
    )
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument(
        "--polynomial-order",
        choices=("generated", "target"),
        default="generated",
    )
    parser.add_argument(
        "--revert",
        choices=(
            "none",
            "preprocessing",
            "preparation",
            "serial_poll",
            "sssf_trees",
            "segments",
            "collisions",
            "division",
            "scoring",
            "matrix",
        ),
        default="none",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 is required")
    if args.output.exists():
        parser.error("preserve existing captures; choose a new output")
    if args.baseline_json:
        with tempfile.TemporaryDirectory(prefix="factor-p361-") as directory:
            args.runtime_root = Path(directory)
            materialize_baseline(args.baseline_json, args.runtime_root)
            run(args)
    else:
        run(args)


def run(args):
    """Measure exactly one selected runtime without mutating its files."""
    args.runtime_root = args.runtime_root.resolve()
    if (
        "v2" in sys.modules
        and Path(sys.modules["v2"].__file__).resolve().parent.parent
        != args.runtime_root
    ):
        raise RuntimeError(
            "invoke this file directly to isolate a different runtime"
        )
    sys.path.insert(0, str(args.runtime_root))
    modules = {
        name: importlib.import_module(
            "v2."
            + (
                "qs."
                if name
                in (
                    "factor_base",
                    "polynomial",
                    "sieve_collector",
                    "sss",
                    "siqs",
                    "parallel",
                    "linear_algebra",
                )
                else ""
            )
            + name
        )
        for name in (
            "utils",
            "budget",
            "portfolio",
            "factor_base",
            "polynomial",
            "sieve_collector",
            "sss",
            "siqs",
            "parallel",
            "linear_algebra",
        )
    }
    corpus = json.loads(args.corpus.read_text())
    verify_corpus(corpus)
    apply_revert(args.revert)
    if args.polynomial_order == "target":
        from v2.qs import families

        original = families.family_assignments

        def ordered(base, half_width, **options):
            assignments = original(base, half_width, **options)
            target = families.a_target(base.n_prime, half_width)
            options["budget"].consume(
                len(assignments) * options["factor_count"]
            )
            return tuple(
                sorted(
                    assignments, key=lambda primes: abs(prod(primes) - target)
                )
            )

        modules["siqs"].family_assignments = ordered
        modules["parallel"].family_assignments = ordered
    before = fingerprint(args.runtime_root)
    data = dict(
        runtime=sys.version,
        machine=platform.machine(),
        source_sha256=before,
        corpus_sha256=hashlib.sha256(args.corpus.read_bytes()).hexdigest(),
        split=args.split,
        revert=args.revert,
        polynomial_order=args.polynomial_order,
        control_sha256=hashlib.sha256(BASELINE.read_bytes()).hexdigest()
        if args.revert != "none"
        else None,
        work_limit=2_000_000_000,
        seconds=5,
        larger_diagnostic_seconds=0.2,
        feasible_seconds=30,
        feasible_work_limit=10**13,
        results={},
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    for band in args.bands.split(","):
        fixtures = [
            f
            for f in corpus["fixtures"]
            if f["band"] == band and f["split"] == args.split
        ]
        if not fixtures:
            raise ValueError("band has no fixtures: " + band)
        for case in args.cases.split(","):
            print("starting", band, case, flush=True)
            config = make_config(case, band, modules)
            pool = None
            if case.startswith("parallel"):
                parts = case.split("_")
                mode = next(
                    (p for p in parts if p in ("serial", "thread", "process")),
                    "serial",
                )
                workers = (
                    int(parts[-1])
                    if parts[-1].isdigit()
                    else (1 if mode == "serial" else 2)
                )
                pool = modules["parallel"].CollectionPool(mode, workers)

            def cohort():
                if "cold" in case:
                    with modules["parallel"].CollectionPool(
                        mode, workers
                    ) as fresh:
                        return [
                            run_one(case, f, seed, config, modules, fresh)
                            for f in fixtures
                            for seed in corpus.get("seeds", [7, 29])
                        ]
                return [
                    run_one(case, f, seed, config, modules, pool)
                    for f in fixtures
                    for seed in corpus.get("seeds", [7, 29])
                ]

            try:
                result = measure(cohort, args)
                result["config"] = (
                    config
                    if config is None or isinstance(config, dict)
                    else asdict(config)
                )
                if args.profile:
                    profiler = cProfile.Profile()
                    profiler.runcall(cohort)
                    profiler.dump_stats(
                        str(args.output.with_suffix(f".{band}.{case}.prof"))
                    )
                data["results"][band + "/" + case] = result
                args.output.write_text(json.dumps(data, indent=2) + "\n")
                last = result["attempts"][-1]["samples"][0]["rows"]
                print(
                    band,
                    case,
                    round(result["median_seconds"], 6),
                    result["stable"],
                    sum(r["completed"] for r in last),
                    "/",
                    len(last),
                    flush=True,
                )
            finally:
                if pool:
                    pool.close()
    if before != fingerprint(args.runtime_root):
        raise AssertionError("runtime source changed during measurement")


if __name__ == "__main__":
    main()
