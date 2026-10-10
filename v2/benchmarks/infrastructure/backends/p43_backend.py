"""Matched PyPy int/GMP study, including a frozen pre-P4.3 control.

Every warmup and sample checks its output. Import/startup is measured in
separate processes; setup/conversions remain inside the warmed workloads.
This bounded study does not calibrate a universal backend crossover.
"""

import argparse
import hashlib
import importlib
import importlib.abc
import importlib.util
import json
import math
import random
import resource
import statistics
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace

from .... import portfolio
from ....common import arithmetic, utils
from ....execution import stage_jobs
from ....execution.budget import Budget
from ....execution.schedules import SieveContext
from ....qs import SieveConfig, SIQSConfig, SIQSJob
from ....qs.linear_algebra import DependencySolver, filter_matrix
from ....qs.sss import SSSConfig, SSSJob
from ...suites.build_phase_two_corpus import verify_certificates
from ...suites.phase_one import environment
from ...support.paths import (
    BENCHMARK_ROOT,
    source_path,
)

INPUTS = BENCHMARK_ROOT / "inputs"
CORPUS = INPUTS / "corpora/p43_backend_corpus.json"
BASELINE = INPUTS / "baselines/p43_before_sources.json"


class _FrozenFinder(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    def __init__(self, data):
        self.sources = data["sources"]
        self.name = "_factor_p43_before"
        for path, source in self.sources.items():
            if (
                hashlib.sha256(source.encode()).hexdigest()
                != data["sha256"][path]
            ):
                raise ValueError("corrupt P4.3 baseline source")

    def _path(self, name):
        suffix = name.removeprefix(self.name).replace(".", "/")
        root = "v2" + suffix
        return (
            root + ".py"
            if root + ".py" in self.sources
            else root + "/__init__.py"
        )

    def find_spec(self, fullname, path=None, target=None):
        if fullname == self.name or fullname.startswith(self.name + "."):
            source_path = self._path(fullname)
            if source_path in self.sources:
                return importlib.util.spec_from_loader(
                    fullname,
                    self,
                    is_package=source_path.endswith("/__init__.py"),
                )
        return None

    def create_module(self, spec):
        return None

    def exec_module(self, module):
        path = self._path(module.__name__)
        module.__file__ = str(BASELINE) + ":" + path
        exec(
            compile(self.sources[path], module.__file__, "exec"),
            module.__dict__,
        )


def load_baseline():
    """Import hash-verified immutable sources, including their dependencies."""
    finder = _FrozenFinder(json.loads(source_path(BASELINE).read_text()))
    sys.meta_path.insert(0, finder)
    return SimpleNamespace(
        portfolio=importlib.import_module(finder.name + ".portfolio"),
        stage_jobs=importlib.import_module(finder.name + ".stage_jobs"),
        siqs=importlib.import_module(finder.name + ".qs.siqs"),
        sss=importlib.import_module(finder.name + ".qs.sss"),
        linear=importlib.import_module(finder.name + ".qs.linear_algebra"),
        ecm=importlib.import_module(finder.name + ".ecm"),
    )


def _budget():
    return Budget(work_limit=200_000, seconds=None, cpu_seconds=None)


def _result_signature(result):
    return (
        result.original,
        result.sign,
        tuple(
            (f.value, f.exponent, f.certainty.value) for f in result.factors
        ),
        result.remaining,
    )


def _validated_portfolio(run, fixture):
    """Check partial reconstruction and certainty against prime proofs."""
    result = run.result
    if result.reconstruct() != fixture["n"]:
        raise AssertionError("portfolio lost a cofactor")
    remaining = fixture["n"]
    certified = {prime for prime, _ in fixture["factors"]}
    for factor in result.factors:
        if factor.value not in certified:
            raise AssertionError("terminal factor is not independently prime")
        if factor.certainty.value == "proven_prime" and factor.value >= 2**64:
            raise AssertionError("backend upgraded probable-prime certainty")
        power = factor.value**factor.exponent
        if remaining % power:
            raise AssertionError("terminal factor multiplicity is incorrect")
        remaining //= power
    if math.prod(result.remaining) != remaining:
        raise AssertionError("unresolved cofactors do not reconstruct")
    return (_result_signature(result), run.reason, run.work_used)


def _portfolio_call(module, fixtures, seeds, backend=None):
    options = dict(
        trial_bound=100,
        rho_attempts=1,
        rho_evaluations=512,
        pm1_b1=100,
        pm1_b2=1000,
        ecm_tiers=((100, 1000, 2),),
        segment_size=128,
        trace_limit=32,
    )
    if backend is not None:
        options["backend"] = backend
    config = module.PortfolioConfig(**options)
    outputs = []
    for fixture in fixtures:
        for seed in seeds:
            run = module.factorize_bounded(
                fixture["n"], seed=seed, config=config, budget=_budget()
            )
            outputs.append(_validated_portfolio(run, fixture))
    return outputs


def _stage_call(kind, n, backend, module=stage_jobs):
    config = portfolio.PortfolioConfig(
        trial_bound=10,
        rho_evaluations=2048,
        chunk_size=16,
        gcd_batch=64,
        pm1_b1=200,
        pm1_b2=2000,
        ecm_tiers=((200, 2000, 1),),
    )
    integer = arithmetic.get_backend(backend).integer
    job = module.new_job(kind, integer(n), 104729, 200, 2000)
    context, budget = SieveContext(2001, segment_size=128), _budget()
    while not job["done"]:
        module.advance_job(job, budget, context, config)
    if job["factor"] is not None and not utils.valid_divisor(job["factor"], n):
        raise AssertionError("stage returned an invalid divisor")
    return arithmetic.canonical(job), budget.used


def _qs_call(mode, backend, before=None):
    module = (
        before.siqs
        if before is not None
        else SimpleNamespace(SIQSConfig=SIQSConfig, SIQSJob=SIQSJob)
    )
    # Configurations include setup, collection, filtering and extraction.
    options = dict(
        mode=mode,
        base_bound=200,
        half_width=256,
        collector=SieveConfig(
            residual_bound=500,
            max_atoms=4096,
            max_relations=4096,
            max_partials=4096,
        ),
    )
    if before is not None:
        options["collector"] = module.SieveConfig(
            **asdict(options["collector"])
        )
    else:
        options["backend"] = backend
    config = module.SIQSConfig(**options)
    budget = Budget(work_limit=20_000_000, seconds=None, cpu_seconds=None)
    job = module.SIQSJob(4001 * 5003, seed=73, config=config, budget=budget)
    result = job.run()
    if result.divisor is not None:
        if (
            not utils.valid_divisor(result.divisor, 4001 * 5003)
            or result.divisor * result.cofactor != 4001 * 5003
        ):
            raise AssertionError("QS output failed independent reconstruction")
    elif result.cofactor != 4001 * 5003:
        raise AssertionError("QS lost its unresolved cofactor")
    return result.reason, result.divisor, result.cofactor, budget.used


def _sss_call(mode, backend, before=None):
    module = (
        before.sss
        if before is not None
        else SimpleNamespace(SSSConfig=SSSConfig, SSSJob=SSSJob)
    )
    options = dict(mode=mode, base_bound=400, search_rounds=32)
    if before is None:
        options["backend"] = backend
    job = module.SSSJob(
        4001 * 4003,
        seed=73,
        config=module.SSSConfig(**options),
        budget=Budget(work_limit=20_000_000, seconds=None, cpu_seconds=None),
    )
    result = job.run()
    if result.divisor is not None:
        if (
            not utils.valid_divisor(result.divisor, 4001 * 4003)
            or result.divisor * result.cofactor != 4001 * 4003
        ):
            raise AssertionError("SSS output failed reconstruction")
    elif result.cofactor != 4001 * 4003:
        raise AssertionError("SSS lost its unresolved cofactor")
    return result.reason, result.divisor, result.cofactor, job.budget.used


def _matrix_call(rows, backend, before=None):
    integer = arithmetic.get_backend(backend).integer
    typed = tuple(integer(row) for row in rows)
    module = (
        before.linear
        if before
        else SimpleNamespace(
            filter_matrix=filter_matrix, DependencySolver=DependencySolver
        )
    )
    matrix = module.filter_matrix(
        typed,
        weight_two=True,
        budget=Budget(work_limit=10**9, seconds=None, cpu_seconds=None),
    )
    solver = module.DependencySolver(
        matrix, budget=Budget(work_limit=10**9, seconds=None, cpu_seconds=None)
    )
    masks = solver.run()
    return tuple(int(mask) for mask in masks), solver.xors


def _micro_call(operation, pairs, backend):
    integer = arithmetic.get_backend(backend).integer
    pairs = tuple((integer(a), integer(b)) for a, b in pairs)
    if operation == "gcd":
        return [arithmetic.gcd(a, b) for a, b in pairs]
    if operation == "powmod":
        return [arithmetic.pow(a, 65537, b) for a, b in pairs]
    if operation == "inverse":
        return [utils.modular_inverse(a, b) for a, b in pairs]
    if operation == "roots":
        return [arithmetic.integer_root(a, 7) for a, _ in pairs]
    if operation == "division":
        return [arithmetic.divexact(a * b, b) for a, b in pairs]
    raise ValueError(operation)


def _interval(left, right):
    """Conditional bootstrap time-change interval for isolated arm samples."""
    generator = random.Random(4305)
    changes = []
    for _ in range(2000):
        a = statistics.median(generator.choice(left) for _ in left)
        b = statistics.median(generator.choice(right) for _ in right)
        changes.append(100 * (b / a - 1))
    changes.sort()
    return changes[50], changes[1949]


def _wire(output):
    return json.dumps(
        output,
        sort_keys=True,
        separators=(",", ":"),
        default=arithmetic.json_integer,
    )


def _workload(spec):
    """Resolve one arm in its own interpreter so GMP cannot alter int JITs."""
    backend = spec["backend"]
    before = load_baseline() if backend == "before-int" else None
    name = "python-int" if before is not None else backend
    kind = spec["kind"]
    if kind == "micro":
        return lambda: _micro_call(spec["operation"], spec["pairs"], name)
    if kind == "matrix":
        return lambda: _matrix_call(tuple(spec["rows"]), name, before)
    if kind == "stage":
        module = before.stage_jobs if before is not None else stage_jobs
        return lambda: _stage_call(spec["stage"], spec["n"], name, module)
    if kind == "qs":
        return lambda: _qs_call(spec["mode"], name, before)
    if kind == "sss":
        return lambda: _sss_call(spec["mode"], name, before)
    if kind == "portfolio":
        module = before.portfolio if before is not None else portfolio
        return lambda: _portfolio_call(
            module, spec["fixtures"], spec["seeds"], None if before else name
        )
    raise ValueError("unknown workload")


def _measure_worker(request):
    call = _workload(request["spec"])
    expected = request["expected"]
    started, count = time.perf_counter(), 0
    while (
        count == 0 or time.perf_counter() - started < request["warmup_seconds"]
    ):
        if _wire(call()) != expected:
            raise AssertionError("invalid isolated warmup output")
        count += 1
    warmup_seconds = time.perf_counter() - started
    samples, target = [], request["repetitions"]
    while len(samples) < target:
        started = time.perf_counter()
        output = call()
        elapsed = time.perf_counter() - started
        if _wire(output) != expected:
            raise AssertionError("invalid isolated timed output")
        samples.append(elapsed)
        if len(samples) == target and target < 27:
            if statistics.stdev(samples) / statistics.mean(samples) > 0.15:
                target = min(27, target + 9)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return dict(
        warmup_seconds=warmup_seconds,
        warmup_calls=count,
        samples_seconds=samples,
        median_seconds=statistics.median(samples),
        p95_seconds=sorted(samples)[math.ceil(0.95 * len(samples)) - 1],
        relative_stdev=statistics.stdev(samples) / statistics.mean(samples),
        peak_rss_bytes=rss if sys.platform == "darwin" else rss * 1024,
    )


def measure(
    name, specs, expected, repetitions, warmup_seconds, *, reverse=False
):
    """Isolate JITs, include conversions and verify every timed output."""
    row = dict(name=name, arms={})
    labels = list(specs)
    if reverse:
        labels.reverse()
    for label in labels:
        request = dict(
            spec=specs[label],
            expected=_wire(expected),
            repetitions=repetitions,
            warmup_seconds=warmup_seconds,
        )
        child = subprocess.run(
            [
                sys.executable,
                "-m",
                "v2.benchmarks.infrastructure.backends.p43_backend",
                "--measure-worker",
            ],
            input=json.dumps(request),
            text=True,
            capture_output=True,
            check=True,
        )
        row["arms"][label] = json.loads(child.stdout)
    control = row["arms"]["python-int"]
    row["changes_vs_int"] = {
        label: dict(
            median_percent=100
            * (arm["median_seconds"] / control["median_seconds"] - 1),
            bootstrap_95_percent=_interval(
                control["samples_seconds"], arm["samples_seconds"]
            ),
        )
        for label, arm in row["arms"].items()
        if label != "python-int"
    }
    print(
        name,
        {
            label: round(arm["median_seconds"], 6)
            for label, arm in row["arms"].items()
        },
        flush=True,
    )
    return row


def cold_samples(backends, repetitions):
    """Time new interpreters and validate their full factorization output."""
    samples = {name: [] for name in backends}
    for index in range(repetitions):
        offset = index % len(backends)
        for name in backends[offset:] + backends[:offset]:
            started = time.perf_counter()
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "v2.benchmarks.infrastructure.backends.p43_backend",
                    "--cold-worker",
                    name,
                ],
                check=True,
                text=True,
                capture_output=True,
            )
            elapsed = time.perf_counter() - started
            value = json.loads(result.stdout)
            if value != [
                626100403,
                [[25013, 1, "proven_prime"], [25031, 1, "proven_prime"]],
                [],
            ]:
                raise AssertionError("invalid cold factorization")
            samples[name].append(elapsed)
    return {
        name: dict(
            samples_seconds=values, median_seconds=statistics.median(values)
        )
        for name, values in samples.items()
    }


def run(repetitions=9, warmup_seconds=3):
    if not hasattr(sys, "pypy_version_info") or sys.version_info[:2] != (
        3,
        11,
    ):
        raise RuntimeError("measure only PyPy Python 3.11")
    if repetitions < 9 or warmup_seconds < 3:
        raise ValueError(
            "require at least nine samples and three warmup seconds"
        )
    corpus = json.loads(source_path(CORPUS).read_text())
    verify_certificates(corpus["certificates"])
    for fixture in corpus["fixtures"]:
        if math.prod(p**e for p, e in fixture["factors"]) != fixture["n"]:
            raise ValueError("corrupt corpus reconstruction")
    backends = ["python-int"]
    unavailable = {}
    try:
        arithmetic.get_backend("gmpy2-mpz")
        backends.append("gmpy2-mpz")
    except RuntimeError as error:
        unavailable["gmpy2-mpz"] = str(error)
    report = dict(
        environment=environment(),
        unavailable=unavailable,
        identities={
            name: arithmetic.get_backend(name).identity for name in backends
        },
        corpus_sha256=hashlib.sha256(
            source_path(CORPUS).read_bytes()
        ).hexdigest(),
        baseline_sha256=hashlib.sha256(
            source_path(BASELINE).read_bytes()
        ).hexdigest(),
        cold=cold_samples(backends, repetitions),
        rows=[],
    )

    def compare(name, specification, *, frozen=True):
        specs = {name: dict(specification, backend=name) for name in backends}
        if frozen:
            specs["before-int"] = dict(specification, backend="before-int")
        expected = _workload(
            specs["before-int"] if frozen else specs["python-int"]
        )()
        if specification["kind"] == "matrix":
            # Oracle work stays outside the timed matrix kernel.
            for mask in expected[0]:
                parity = 0
                for index, row in enumerate(specification["rows"]):
                    if mask >> index & 1:
                        parity ^= row
                if parity or not mask:
                    raise AssertionError("invalid independent matrix oracle")
        report["rows"].append(
            measure(
                name,
                specs,
                expected,
                repetitions,
                warmup_seconds,
                reverse=len(report["rows"]) % 2 == 1,
            )
        )
        return expected

    generator = random.Random(corpus["seed"])
    modulus = 2**255 - 19
    pairs = tuple(
        (generator.getrandbits(256) or 1, modulus) for _ in range(128)
    )
    for operation in ("gcd", "powmod", "inverse", "roots", "division"):
        compare(
            "arithmetic_" + operation,
            dict(kind="micro", operation=operation, pairs=pairs),
            frozen=False,
        )
    rows = tuple(generator.getrandbits(256) for _ in range(320))
    compare("matrix_256_columns_320_rows", dict(kind="matrix", rows=rows))
    for band in ("balanced_50d", "balanced_80d"):
        n = next(f["n"] for f in corpus["fixtures"] if f["band"] == band)
        for kind in ("rho", "pm1", "ecm"):
            compare(
                f"complete_{kind}_stage_{n.bit_length()}bit",
                dict(kind="stage", stage=kind, n=n),
            )
    for mode in ("qs", "mpqs", "siqs"):
        compare(mode + "_setup_to_extraction", dict(kind="qs", mode=mode))
    for mode in ("sss", "sssf"):
        compare(mode + "_setup_to_extraction", dict(kind="sss", mode=mode))
    fixtures = corpus["fixtures"]
    seeds = corpus["seeds"]
    expected = compare(
        "portfolio_11_inputs_3_seeds",
        dict(kind="portfolio", fixtures=fixtures, seeds=seeds),
    )
    report["portfolio_outcomes"] = [
        dict(id=fixture["id"], seed=seed, result=output)
        for (fixture, seed), output in zip(
            ((f, s) for f in fixtures for s in seeds), expected
        )
    ]
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    report["process_peak_rss_bytes"] = (
        rss if sys.platform == "darwin" else rss * 1024
    )
    report["limits"] = dict(
        measurement=(
            "each arm has an isolated PyPy JIT; alternate arm order by case; "
            "validation outside timed samples"
        ),
        core_budget=1,
        portfolio_work_limit=200_000,
        portfolio_seconds=None,
        portfolio_cpu_seconds=None,
        memory_policy=(
            "same finite engine workspace caps; "
            "each isolated arm reports its process peak RSS"
        ),
        interpretation=(
            "fixed-work complete/partial outcomes; "
            "no universal default promotion or thread/GIL claim"
        ),
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--measure-worker", action="store_true")
    parser.add_argument("--cold-worker", choices=("python-int", "gmpy2-mpz"))
    args = parser.parse_args()
    if args.measure_worker:
        print(json.dumps(_measure_worker(json.load(sys.stdin))))
        return
    if args.cold_worker:
        from ....factor import factorize

        result = factorize(626100403, seed=7, backend=args.cold_worker)
        print(
            json.dumps(
                [
                    result.original,
                    [
                        [f.value, f.exponent, f.certainty.value]
                        for f in result.factors
                    ],
                    result.remaining,
                ]
            )
        )
        return
    if args.output is None:
        parser.error("--output is required for preserved benchmark evidence")
    report = run(args.repetitions, args.warmup_seconds)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, default=arithmetic.json_integer) + "\n"
    )
    print(args.output)


if __name__ == "__main__":
    main()
