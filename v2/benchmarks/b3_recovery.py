"""Frozen lazy-recovery ablation; preserve the earlier B3 experiments."""

import argparse
import contextlib
import fcntl
import json
import os
import random
import resource
import signal
import subprocess
import sys
import time
import types
from math import lcm, prod
from pathlib import Path

from .. import arithmetic, ecm, ecm_chains, stage_jobs
from . import b3_coverage as coverage
from . import b3_production as study
from .a6_pm1 import assert_quiet, require_runtime
from .build_b3_inputs import write_new
from .build_c6_fast_inputs import LimitedRandom
from .build_phase_two_corpus import certified_prime, verify_certificates
from .prac_oracle import affine_multiply, twist_point

PROTOCOL = study.INPUTS / "controls/b3_recovery_protocol.json"
CORPUS = study.INPUTS / "corpora/b3_recovery_confirmation.json"
EAGER = study.INPUTS / "baselines/b3_eager_chains.json"
LAZY_STORE = ecm_chains.ChainPlans
SCOPES = (
    "preparation",
    "fresh1",
    "reuse8",
    "eviction",
    "portfolio",
    "resume_equal",
    "absolute",
    "ownership",
)

GROUPS = [
    ("training", "preparation", "python-int", False),
    ("training", "preparation", "gmpy2-mpz", False),
    ("training", "portfolio", "python-int", False),
    ("training", "absolute", "python-int", False),
    ("training", "ownership", "python-int", False),
    ("training", "portfolio", "python-int", True),
    ("training", "portfolio", "gmpy2-mpz", True),
    ("training", "portfolio", "gmpy2-mpz", False),
    ("confirmation", "portfolio", "python-int", False),
    ("confirmation", "portfolio", "gmpy2-mpz", False),
    ("training", "resume_equal", "python-int", False),
    ("training", "resume_equal", "gmpy2-mpz", False),
    ("confirmation", "resume_equal", "python-int", False),
    ("confirmation", "resume_equal", "gmpy2-mpz", False),
    ("training", "fresh1", "python-int", False),
    ("training", "reuse8", "python-int", False),
    ("training", "eviction", "python-int", False),
    ("confirmation", "preparation", "python-int", False),
    ("confirmation", "preparation", "gmpy2-mpz", False),
]


def freeze():
    """Freeze one candidate and fresh certified inputs before any timing."""
    if PROTOCOL.exists() or CORPUS.exists():
        raise ValueError("recovery protocol is immutable")
    generator, certificates, fixtures = LimitedRandom(116761), {}, []

    def prime(digits):
        for _ in range(64):
            value = certified_prime(
                (10**digits - 1).bit_length(), generator, certificates
            )
            if len(str(value)) == digits:
                return value
        raise RuntimeError("recovery prime draw cap")

    for digits in (40, 50, 60, 70, 80):
        for shape in ("balanced", "small10"):
            target = digits // 2 if shape == "balanced" else 10
            for _ in range(64):
                factors = [prime(target), prime(digits - target)]
                if (
                    len(str(prod(factors))) == digits
                    and len(set(factors)) == 2
                ):
                    break
            else:
                raise RuntimeError("recovery composite draw cap")
            fixtures.append(
                dict(id=f"{shape}_{digits}d", n=prod(factors), factors=factors)
            )
    for index, sizes in enumerate(((6, 6, 8), (6, 10, 20))):
        factors = [prime(size) for size in sizes]
        fixtures.append(
            dict(id=f"recursive_{index}", n=prod(factors), factors=factors)
        )
    verify_certificates(certificates)
    old = set()
    for path in CORPUS.parent.glob("*.json"):
        data = json.loads(path.read_text())
        if isinstance(data, dict):
            old.update(
                f["n"]
                for f in data.get("fixtures", [])
                if isinstance(f, dict) and "n" in f
            )
    values = {f["n"] for f in fixtures}
    if old & values or len(values) != 12:
        raise ValueError("recovery confirmation overlaps")
    data = dict(
        generation_seed=116761, certificates=certificates, fixtures=fixtures
    )
    if len(json.dumps(data).encode()) > 2**24:
        raise MemoryError("recovery corpus cap")
    write_new(CORPUS, data)
    paths = set(json.loads(coverage.PROTOCOL.read_text())["production_sha256"])
    paths.update(
        str(p.relative_to(study.ROOT))
        for p in (
            Path(__file__),
            CORPUS,
            EAGER,
            study.BASELINE,
            study.CORPUS,
            Path(study.__file__),
            Path(coverage.__file__),
        )
    )
    paths.update(
        str((Path(__file__).parent / name).relative_to(study.ROOT))
        for name in (
            "a6_pm1.py",
            "b4_common.py",
            "build_b3_inputs.py",
            "build_c6_fast_inputs.py",
            "build_phase_two_corpus.py",
            "prac_oracle.py",
        )
    )
    write_new(
        PROTOCOL,
        dict(
            control_commit="b956801",
            candidate=(
                "Only defer Executor construction until strict recovery; "
                "preserve all proofs, charges, memory and identities."
            ),
            sha256={
                name: study.digest(study.ROOT / name) for name in sorted(paths)
            },
            training_seeds=[31001, 38920],
            confirmation_seeds=[78515, 86434],
            scopes=list(SCOPES),
            groups=GROUPS,
            sampling=[[3, 9], [5, 18], [8, 27]],
            max_relative_iqr=0.15,
            bootstrap_seed=193001,
            bounds=[2000, 147396],
            curves=8,
            work=coverage.WORK_LIMIT,
            seconds=20,
            cpu_seconds=20,
            memory_bytes=33554432,
            phase_seconds=900,
            study_seconds=1800,
            worker_seconds=480,
            capture_bytes=67108864,
            generation=dict(
                random_draws=500000,
                decimal_attempts=64,
                seconds=60,
                output_bytes=16777216,
            ),
            preparation_repetitions=32,
            arms=(
                "Portfolio/resume: frozen mainline ladder, accepted eager B3 "
                "and lazy B3 in rotated/reversed order; stage scopes use "
                "the same three arms. Preparation: eager/lazy. Absolute: "
                "native/GMP on production routing and "
                "on identical Lucas records/readable kernels, stage reuse8, "
                "one worker. Ownership: lazy per-run versus one finite owner "
                "per timed cohort, plus a ladder control; its first miss "
                "remains inside the timer. "
                "Both extra scopes diagnose causes; they do not change APIs."
            ),
            validation=(
                "Reconstruct all results including unresolved. Check proper "
                "known divisors, deterministic work/outcomes and identical "
                "eager/lazy work. Unresolved cofactors require eight curves "
                "each. Resume cancels after the first certified chunk. "
                "Check independent affine stage outputs. Charge fresh setup, "
                "rebuilding, eviction and failed searches."
            ),
            acceptance=(
                "Require positive paired 95% intervals, CPU/halves and "
                "stability for setup on both backends; unchanged coverage "
                "and no >3% full-run regression interval against eager B3. "
                "Claim full gains only if training AND untouched confirmation "
                "beat each backend's own ladder under revised roadmap policy. "
                "Defaults stay off. No candidate tuning, chain search or new "
                "arithmetic formula."
            ),
            instrumentation=(
                "Separate counts and nested preparation/Executor times are "
                "causal diagnostics only, never accepted timing. Do not "
                "remeasure cold startup or conflate it with warmed execution."
            ),
        ),
    )


def inputs(split):
    settings = json.loads(PROTOCOL.read_text())
    for name, expected in settings["sha256"].items():
        if study.digest(study.ROOT / name) != expected:
            raise ValueError("recovery source/input changed: " + name)
    path = CORPUS if split == "confirmation" else study.CORPUS
    data = json.loads(path.read_text())
    verify_certificates(data["certificates"])
    fixtures = data["fixtures"]
    if split == "training":
        fixtures = [f for f in fixtures if f["split"] == "training"]
    if any(prod(f["factors"]) != f["n"] for f in fixtures):
        raise ValueError("recovery corpus reconstruction")
    return settings, fixtures


def eager_module(same_lucas=False):
    """Load accepted source with optional diagnostic substitutions."""
    data = json.loads(EAGER.read_text())
    source = data["source"]
    if study.hashlib.sha256(source.encode()).hexdigest() != data["sha256"]:
        raise ValueError("eager source snapshot corrupt")
    if same_lucas:
        replacements = (
            (
                'family = "prac" if backend == "python-int" else "lucas"',
                'family = "lucas"',
            ),
            ('native = backend.name == "python-int"', "native = False"),
            (
                'family = "prac-reduced" if backend == "python-int" '
                'else "lucas-tuple"',
                'family = "diagnostic-lucas-readable"',
            ),
        )
        for original, replacement in replacements:
            if source.count(original) != 1:
                raise ValueError("same-Lucas diagnostic substitution changed")
            source = source.replace(original, replacement)
    module = types.ModuleType("v2._b3_eager_control")
    module.__package__, module.__file__ = "v2", ecm_chains.__file__
    exec(compile(source, str(EAGER), "exec"), module.__dict__)
    return module


@contextlib.contextmanager
def selected(store):
    previous = ecm_chains.ChainPlans, study.ChainPlans
    ecm_chains.ChainPlans = study.ChainPlans = store
    try:
        yield
    finally:
        ecm_chains.ChainPlans, study.ChainPlans = previous


def stage_targets(fixtures, seeds, scope):
    curves = (
        3
        if scope == "eviction"
        else int(scope.replace("fresh", "").replace("reuse", ""))
    )
    targets = {}
    for fixture in fixtures:
        if not fixture["id"].startswith("balanced"):
            continue
        for seed in seeds:
            for curve in range(curves):
                bound = 1999 if scope == "eviction" and curve == 1 else 2000
                sigma = random.Random(seed + curve).randrange(6, 2**63)
                setup = ecm.setup_curve(fixture["n"], sigma)
                expected = []
                if setup.point is not None:
                    for prime in fixture["factors"]:
                        point, curve_a, curve_b = twist_point(
                            setup.point, setup.a24, prime
                        )
                        expected.append(
                            (
                                prime,
                                affine_multiply(
                                    lcm(*range(1, bound + 1)),
                                    point,
                                    prime,
                                    curve_a,
                                    curve_b,
                                ),
                            )
                        )
                targets[(fixture["id"], seed, curve)] = expected
    return targets


def preparation(store, backend, settings):
    rows = []
    for index in range(settings["preparation_repetitions"]):
        started = time.perf_counter()
        cache, budget = (
            store(ecm_chains.MIN_MEMORY_BYTES, backend, (2000,)),
            coverage.allowance(study.portfolio),
        )
        plan = cache.get(2000, budget)
        if len(plan.entries) != 303 or cache.misses != 1:
            raise AssertionError("incomplete prepared plan")
        rows.append(
            dict(
                fixture="preparation",
                complete=True,
                work=budget.used,
                owned_bytes=plan.owned_bytes,
                seconds=time.perf_counter() - started,
                index=index,
            )
        )
        del plan, cache
    return rows


def cohort(fixtures, seeds, backend, candidate, scope):
    """Validate full search on every unresolved recursive cofactor."""
    engine = study.portfolio if candidate else study.baseline()
    config = study.configuration(engine, backend, candidate)
    rows = []
    for fixture in fixtures:
        if scope == "resume_equal" and not fixture["id"].startswith(
            "balanced"
        ):
            continue
        for seed in seeds:
            started = time.perf_counter()
            options = dict(seed=seed)
            if scope == "resume_equal":
                options = dict(
                    checkpoint=coverage.equal_prefix_pause(
                        engine, fixture, config, seed
                    )
                )
            run = engine.factorize_bounded(
                fixture["n"],
                config=config,
                budget=coverage.allowance(engine),
                **options,
            )
            row = study.validate(run, fixture)
            counts = {
                str(value): sum(
                    event.get("stage") == "ecm" and event.get("n") == value
                    for event in run.events
                )
                for value in run.result.remaining
            }
            if not run.result.complete and (
                run.reason != "exhausted"
                or any(n != 8 for n in counts.values())
            ):
                raise AssertionError("incomplete unresolved curve schedule")
            rows.append(
                dict(
                    row,
                    seed=seed,
                    unresolved_curve_counts=counts,
                    root_curves=sum(
                        event.get("stage") == "ecm"
                        and event.get("n") == fixture["n"]
                        for event in run.events
                    ),
                    certified_pause_prefix=16
                    if scope == "resume_equal"
                    else None,
                    seconds=time.perf_counter() - started,
                )
            )
    return rows


def worker(args):
    settings, fixtures = inputs(args.split)
    seeds = settings[args.split + "_seeds"]
    eager = eager_module(args.scope == "absolute")
    scope = "reuse8" if args.scope == "absolute" else args.scope
    targets = (
        stage_targets(fixtures, seeds, scope)
        if scope in ("fresh1", "reuse8", "eviction")
        else {}
    )
    if args.scope in (
        "portfolio",
        "resume_equal",
        "fresh1",
        "reuse8",
        "eviction",
    ):
        arms = [
            ("ladder", None, args.backend),
            ("eager", eager.ChainPlans, args.backend),
            ("lazy", LAZY_STORE, args.backend),
        ]
    elif args.scope == "ownership":
        arms = [
            ("ladder", None, args.backend),
            ("lazy-per-run", LAZY_STORE, args.backend),
            ("lazy-cohort-owner", LAZY_STORE, args.backend),
        ]
    elif args.scope == "absolute":
        arms = [
            ("production-int", LAZY_STORE, "python-int"),
            ("production-gmp", LAZY_STORE, "gmpy2-mpz"),
            ("same-lucas-int", eager.ChainPlans, "python-int"),
            ("same-lucas-gmp", eager.ChainPlans, "gmpy2-mpz"),
        ]
    else:
        arms = [
            ("eager", eager.ChainPlans, args.backend),
            ("lazy", LAZY_STORE, args.backend),
        ]
    signatures, warmups, samples = {}, [], []

    def measure(index):
        name, store, backend = arms[index]
        if name == "lazy-cohort-owner":
            owned = []

            def factory(*values, **options):
                if not owned:
                    owned.append(LAZY_STORE(*values, **options))
                return owned[0]

            store = factory
        with selected(store or LAZY_STORE):
            started, cpu = time.perf_counter(), time.process_time()
            if scope == "preparation":
                rows = preparation(store, backend, settings)
            elif scope in ("portfolio", "resume_equal", "ownership"):
                rows = cohort(
                    fixtures,
                    seeds,
                    backend,
                    store is not None,
                    "portfolio" if scope == "ownership" else scope,
                )
            else:
                rows = study.stage_cohort(
                    fixtures, seeds, backend, store is not None, scope, targets
                )
            rows = arithmetic.canonical(rows)
            elapsed, used = (
                time.perf_counter() - started,
                time.process_time() - cpu,
            )
        signature = [
            {k: v for k, v in row.items() if k != "seconds"} for row in rows
        ]
        if name in signatures and signatures[name] != signature:
            raise AssertionError("non-deterministic outcome/work: " + name)
        signatures[name] = signature
        return dict(wall=elapsed, cpu=used, rows=rows)

    if args.scope == "absolute":
        comparisons = [(0, 1), (2, 3)]
    else:
        comparisons = [(0, 2), (1, 2)] if len(arms) == 3 else [(0, 1)]
    for seconds, count in settings["sampling"]:
        for arm in range(len(arms)):
            started, rounds = time.perf_counter(), 0
            while time.perf_counter() - started < seconds:
                measure(arm)
                rounds += 1
            warmups.append(
                dict(
                    arm=arms[arm][0],
                    seconds=time.perf_counter() - started,
                    cohorts=rounds,
                )
            )
        while len(samples) < count:
            order = list(range(len(arms)))
            shift = len(samples) % len(arms)
            order = order[shift:] + order[:shift]
            if len(samples) % 2:
                order.reverse()
            row = [None] * len(arms)
            for arm in order:
                row[arm] = measure(arm)
            samples.append(row)
            if args.output is not None:
                # Keep completed samples if the finite worker deadline wins.
                partial = args.output.with_suffix(".partial.json")
                temporary = partial.with_suffix(".tmp")
                temporary.write_text(
                    json.dumps(
                        dict(
                            incomplete=True,
                            scope=args.scope,
                            split=args.split,
                            backend=args.backend,
                            arms=[a[0] for a in arms],
                            samples=samples,
                            warmups=warmups,
                            protocol_sha256=study.digest(PROTOCOL),
                        )
                    )
                )
                temporary.replace(partial)
        summaries = {
            arms[a][0] + "_to_" + arms[b][0]: study.summarize(
                [[row[a], row[b]] for row in samples], settings
            )
            for a, b in comparisons
        }
        if all(summary["stable"] for summary in summaries.values()):
            break
    if (
        "eager" in signatures
        and "lazy" in signatures
        and signatures["eager"] != signatures["lazy"]
    ):
        raise AssertionError("lazy preparation changed outcomes/work")
    return dict(
        scope=args.scope,
        split=args.split,
        backend=args.backend,
        arms=[a[0] for a in arms],
        samples=samples,
        warmups=warmups,
        summaries=summaries,
        protocol_sha256=study.digest(PROTOCOL),
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    )


def diagnose(args):
    """Count and attribute preparation separately from accepted samples."""
    settings, fixtures = inputs(args.split)
    reports = []
    for name, module, store in (
        ("eager", eager_module(), None),
        ("lazy", ecm_chains, LAZY_STORE),
    ):
        store = store or module.ChainPlans
        counts = dict(
            plan_calls=0,
            plan_seconds=0.0,
            strict_calls=0,
            strict_seconds=0.0,
            chunk_calls=0,
            chunk_seconds=0.0,
            stage_one_calls=0,
            stage_one_seconds=0.0,
            stage_two_calls=0,
            stage_two_seconds=0.0,
        )
        original_init, original_executor = (
            module.ChainPlan.__init__,
            module.Executor,
        )

        def plan_init(self, *values, **options):
            started = time.perf_counter()
            try:
                return original_init(self, *values, **options)
            finally:
                counts["plan_calls"] += 1
                counts["plan_seconds"] += time.perf_counter() - started

        def executor(*values, **options):
            started = time.perf_counter()
            try:
                return original_executor(*values, **options)
            finally:
                counts["strict_calls"] += 1
                counts["strict_seconds"] += time.perf_counter() - started

        def counted(function, prefix):
            def call(*values, **options):
                started = time.perf_counter()
                try:
                    return function(*values, **options)
                finally:
                    counts[prefix + "_calls"] += 1
                    counts[prefix + "_seconds"] += (
                        time.perf_counter() - started
                    )

            return call

        original_execute = module.ChainPlan.execute
        original_one, original_two = (
            stage_jobs._stage_one,
            stage_jobs._stage_two,
        )
        module.ChainPlan.execute = counted(original_execute, "chunk")
        stage_jobs._stage_one = counted(original_one, "stage_one")
        stage_jobs._stage_two = counted(original_two, "stage_two")
        module.ChainPlan.__init__, module.Executor = plan_init, executor
        try:
            with selected(store):
                started = time.perf_counter()
                rows = cohort(
                    fixtures,
                    settings[args.split + "_seeds"],
                    args.backend,
                    True,
                    "portfolio",
                )
                reports.append(
                    dict(
                        arm=name,
                        instrumented=True,
                        wall=time.perf_counter() - started,
                        counts=counts,
                        rows=rows,
                    )
                )
        finally:
            module.ChainPlan.execute = original_execute
            stage_jobs._stage_one, stage_jobs._stage_two = (
                original_one,
                original_two,
            )
            module.ChainPlan.__init__, module.Executor = (
                original_init,
                original_executor,
            )
    return reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", action="store_true")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--diagnose", action="store_true")
    parser.add_argument("--scope", choices=SCOPES, default="portfolio")
    parser.add_argument(
        "--split", choices=("training", "confirmation"), default="training"
    )
    parser.add_argument(
        "--backend", choices=("python-int", "gmpy2-mpz"), default="python-int"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require_runtime()
    if args.freeze:
        signal.signal(
            signal.SIGALRM,
            lambda *unused: (_ for _ in ()).throw(
                TimeoutError("recovery generation cap")
            ),
        )
        signal.alarm(60)
        try:
            freeze()
        finally:
            signal.alarm(0)
        return
    if args.worker:
        print(json.dumps(diagnose(args) if args.diagnose else worker(args)))
        return
    if args.output is None or args.output.exists():
        parser.error("choose a new output path")
    settings, _ = inputs(args.split)
    owner = Path("/private/tmp/factor-performance-owner.json")
    with open("/private/tmp/factor-performance.lock", "a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert_quiet()
        owner.write_text(
            json.dumps(dict(owner="B3 recovery diagnosis", pid=os.getpid()))
        )
        try:
            command = [
                sys.executable,
                "-B",
                "-m",
                __spec__.name,
                "--worker",
                "--scope",
                args.scope,
                "--split",
                args.split,
                "--backend",
                args.backend,
            ]
            if args.diagnose:
                command.append("--diagnose")
            args.output.parent.mkdir(parents=True, exist_ok=True)
            command.extend(["--output", str(args.output)])
            result = subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
                timeout=settings["worker_seconds"],
            )
            data = json.loads(result.stdout)
            inputs(args.split)
            encoded = json.dumps(data, indent=2) + "\n"
            if len(encoded.encode()) > settings["capture_bytes"]:
                raise MemoryError("recovery capture cap")
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("x") as stream:
                stream.write(encoded)
            print(
                json.dumps(
                    data.get("summaries", data)
                    if isinstance(data, dict)
                    else data
                )
            )
        finally:
            owner.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
