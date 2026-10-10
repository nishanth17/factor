"""Matched A10 primality/factoring evidence on PyPy Python 3.11 only.

Required inputs are immutable sources and independently verified proofs.
Every warmup/sample validates outputs. Run under the shared machine lock.
Native reference timings use the same bases, not a weaker probable test.
"""

import argparse
import builtins
import fcntl
import hashlib
import importlib
import importlib.abc
import importlib.util
import json
import math
import os
import random
import statistics
import subprocess
import sys
import time
import types
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from .... import portfolio
from ....common import arithmetic, utils
from ...suites.phase_one import environment
from ...support.paths import (
    BENCHMARK_ROOT,
    REPOSITORY_ROOT,
    source_path,
)
from .a10_inputs import CORPUS, REPORTED_PRIME, load_corpus

BASELINE = BENCHMARK_ROOT / "inputs/baselines/a10_before_sources.json"
PROTOCOL = BENCHMARK_ROOT / "inputs/corpora/a10_protocol.json"


def load_protocol():
    """Reject changed required inputs before measuring the frozen protocol."""
    protocol = json.loads(source_path(PROTOCOL).read_text())
    for path, key in (
        (BASELINE, "control_sha256"),
        (CORPUS, "corpus_sha256"),
    ):
        if (
            hashlib.sha256(source_path(path).read_bytes()).hexdigest()
            != protocol[key]
        ):
            raise ValueError("A10 input differs from frozen protocol")
    return protocol


@contextmanager
def performance_window():
    """Fail immediately if another factor performance owner holds the lock."""
    path = Path("/private/tmp/factor-performance.lock")
    owner_path = Path("/private/tmp/factor-performance-owner.json")
    with path.open("a+") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        metadata = dict(
            owner="A10",
            pid=os.getpid(),
            cwd=str(Path.cwd()),
            started_utc=datetime.now(timezone.utc).isoformat(),
        )
        owner_path.write_text(json.dumps(metadata) + "\n")
        try:
            yield
        finally:
            if (
                json.loads(source_path(owner_path).read_text()).get("pid")
                == os.getpid()
            ):
                owner_path.unlink()
            fcntl.flock(handle, fcntl.LOCK_UN)


class FrozenSources(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    name = "_factor_a10_before"

    def __init__(self, data):
        self.data = data
        for path, source in data["sources"].items():
            if (
                hashlib.sha256(source.encode()).hexdigest()
                != (data["sha256"][path])
            ):
                raise ValueError("corrupt A10 control source")

    def path(self, name):
        root = "v2" + name.removeprefix(self.name).replace(".", "/")
        return (
            root + ".py"
            if root + ".py" in self.data["sources"]
            else root + "/__init__.py"
        )

    def find_spec(self, fullname, path=None, target=None):
        if fullname == self.name or fullname.startswith(self.name + "."):
            source_path = self.path(fullname)
            if source_path in self.data["sources"]:
                return importlib.util.spec_from_loader(
                    fullname,
                    self,
                    is_package=source_path.endswith("/__init__.py"),
                )
        return None

    def create_module(self, spec):
        return None

    def exec_module(self, module):
        path = self.path(module.__name__)
        module.__file__ = str(BASELINE) + ":" + path
        exec(
            compile(self.data["sources"][path], module.__file__, "exec"),
            module.__dict__,
        )


def load_control():
    finder = FrozenSources(json.loads(source_path(BASELINE).read_text()))
    sys.meta_path.insert(0, finder)
    return types.SimpleNamespace(
        utils=importlib.import_module(finder.name + ".utils"),
        portfolio=importlib.import_module(finder.name + ".portfolio"),
        budget=importlib.import_module(finder.name + ".budget"),
    )


def factoring_candidate(name, module):
    """Isolate a powering challenger across complete portfolio calls."""
    data = json.loads(source_path(BASELINE).read_text())
    root = REPOSITORY_ROOT
    data["sources"] = {
        path: (source_path(root / path)).read_text()
        for path in data["sources"]
    }
    data["sha256"] = {
        path: hashlib.sha256(source.encode()).hexdigest()
        for path, source in data["sources"].items()
    }
    finder = FrozenSources(data)
    finder.name = "_factor_a10_candidate_" + name
    sys.meta_path.insert(0, finder)
    engine = importlib.import_module(finder.name + ".portfolio")
    if name == "early_one":
        engine.utils._strong_probable_prime = module._strong_probable_prime
    else:
        engine.utils.pow = module.pow
    return engine


def load_v1_adapter():
    """Emulate only v1 utils' Python-2 dependencies, not native v1 timing."""
    data = json.loads(source_path(BASELINE).read_text())
    source = data["v1_utils"]
    if hashlib.sha256(source.encode()).hexdigest() != data["v1_sha256"]:
        raise ValueError("corrupt v1 utility snapshot")
    module = types.ModuleType("_factor_a10_v1_emulated")
    module.__dict__["xrange"] = range
    # utils uses fractions.gcd only in its unrelated GCD helper. Avoid any
    # global monkeypatch of fractions or Python-2 factoring emulation.
    source = source.replace(
        "import fractions", "from math import gcd as _a10_math_gcd"
    )
    source = source.replace("fractions.gcd(a, b)", "_a10_math_gcd(a, b)")
    exec(
        compile(source, str(BASELINE) + ":v1/utils.py", "exec"),
        module.__dict__,
    )
    return module


def binary_power(base, exponent, modulus):
    value = 1
    while exponent:
        if exponent & 1:
            value = value * base % modulus
        base = base * base % modulus
        exponent >>= 1
    return value


def early_one(n, base, odd_part, shifts):
    base %= n
    if base in (0, 1):
        return True
    value = arithmetic.pow(base, odd_part, n)
    if value in (1, n - 1):
        return True
    for _ in range(shifts - 1):
        value = value * value % n
        if value == n - 1:
            return True
        if value == 1:
            return False
    return False


def candidates():
    """Clone helper globals so experiments never patch production modules."""
    source = source_path(Path(utils.__file__)).read_text()
    small_filter = (
        "    for prime in SMALL_PRIMES:\n"
        "        if n == prime:\n"
        "            return Primality.PROVEN\n"
        "        if n % prime == 0:\n"
        "            return Primality.COMPOSITE\n"
    )
    arms = {"accepted": utils}
    for name in (
        "native_pow",
        "binary_pow",
        "early_one",
        "extra_filters",
        "square_filter",
        "primorial_filter",
        "thirteen_bases",
    ):
        module = types.ModuleType("_factor_a10_" + name)
        module.__package__ = "v2"
        candidate_source = source
        if name == "primorial_filter":
            if source.count(small_filter) != 1:
                raise ValueError("A10 small-filter experiment needs review")
            candidate_source = source.replace(
                small_filter,
                "    if n in SMALL_PRIMES:\n"
                "        return Primality.PROVEN\n"
                "    if _a10_gcd(n, _a10_primorial) > 1:\n"
                "        return Primality.COMPOSITE\n",
            )
        elif name == "square_filter":
            marker = "    if not use_probabilistic and n < 41 * 41:\n"
            if source.count(marker) != 1:
                raise ValueError("A10 square-filter experiment needs review")
            candidate_source = source.replace(
                marker,
                "    if isqrt(n) ** 2 == n:\n"
                "        return Primality.COMPOSITE\n" + marker,
            )
        exec(
            compile(candidate_source, utils.__file__, "exec"), module.__dict__
        )
        if name == "native_pow":
            module.pow = builtins.pow
        elif name == "binary_pow":
            module.pow = binary_power
        elif name == "early_one":
            module._strong_probable_prime = early_one
        elif name == "thirteen_bases":
            dispatch = module.deterministic_bases

            def thirteen(n, dispatch=dispatch):
                if 2**64 <= n < utils.DETERMINISTIC_LIMIT:
                    return utils.THIRTEEN_PRIME_BASES
                return dispatch(n)

            module.deterministic_bases = thirteen
        elif name == "extra_filters":
            module.SMALL_PRIMES += (
                41,
                43,
                47,
                53,
                59,
                61,
                67,
                71,
                73,
                79,
                83,
                89,
                97,
            )
            # Trial divisors change; the independently supported bases do not.
            module.deterministic_bases = utils.deterministic_bases
        elif name == "primorial_filter":
            module._a10_gcd = math.gcd
            module._a10_primorial = math.prod(utils.SMALL_PRIMES)
        arms[name] = module
    return arms


def classification_work(module, fixtures, seed, *, native=False):
    rng = random.Random(seed)
    signature = []
    for fixture in fixtures:
        n = fixture["n"]
        value = arithmetic.get_backend("gmpy2-mpz").integer(n) if native else n
        result = module.classify_prime(value, rng=rng).value
        expected = (
            (
                "proven_prime"
                if n < module.DETERMINISTIC_LIMIT
                else "probable_prime"
            )
            if fixture["prime"]
            else "composite"
        )
        if result != expected:
            raise AssertionError((n, result, expected))
        signature.append(result)
    return signature


def sympy_work(fixtures):
    """Compare SymPy's MR with identical proven ranges; skip BPSW routing."""
    from sympy.ntheory.primetest import mr

    for fixture in fixtures:
        n = fixture["n"]
        bases = utils.deterministic_bases(n)
        if bases is None:
            continue
        result = bool(mr(n, bases))
        if result != fixture["prime"]:
            raise AssertionError("SymPy MR differs from independent oracle")
    return dict(certainty="same fixed bases within strict supported range")


def factoring_work(module, budget_module, fixtures, seeds):
    config = module.PortfolioConfig(
        trial_bound=30000, rho_attempts=0, pm1_attempts=0, ecm_tiers=()
    )
    signature = []
    for fixture in fixtures:
        for seed in seeds:
            run = module.factorize_bounded(
                fixture["n"],
                seed=seed,
                config=config,
                budget=budget_module.Budget(
                    work_limit=1000000, seconds=None, cpu_seconds=None
                ),
            )
            result = run.result
            if (
                result.reconstruct() != fixture["n"]
                or result.remaining
                or run.reason != "complete"
                or [(p.value, p.exponent) for p in result.factors]
                != [tuple(pair) for pair in fixture["factors"]]
            ):
                raise AssertionError("factoring result differs from proof")
            for factor in result.factors:
                expected = (
                    "proven_prime"
                    if factor.value < module.utils.DETERMINISTIC_LIMIT
                    else "probable_prime"
                )
                if factor.certainty.value != expected:
                    raise AssertionError("incorrect terminal certainty")
            signature.append(
                dict(
                    work=run.work_used,
                    certainty=[p.certainty.value for p in result.factors],
                )
            )
    return signature


def measure_group(functions, warmup, repetitions, batch_seconds=0.1):
    """Warm each arm, then interleave samples in frozen random order."""
    results = {}
    for name, function in functions.items():
        started, count = time.perf_counter(), 0
        while time.perf_counter() - started < warmup:
            function()
            count += 1
        elapsed = time.perf_counter() - started
        results[name] = dict(
            samples=[],
            warmup_seconds=elapsed,
            warmup_iterations=count,
            sample_iterations=max(1, int(batch_seconds * count / elapsed)),
        )
    order_rng = random.Random(20261009)
    target, orders = repetitions, []
    while True:
        while len(orders) < target:
            order = list(functions)
            order_rng.shuffle(order)
            orders.append(order)
            for name in order:
                result = results[name]
                before = time.perf_counter()
                for _ in range(result["sample_iterations"]):
                    signature = functions[name]()
                result["samples"].append(
                    (time.perf_counter() - before)
                    / result["sample_iterations"]
                )
                result["signature"] = signature
        for result in results.values():
            samples = result["samples"]
            median = statistics.median(samples)
            spread = (max(samples) - min(samples)) / median
            result.update(
                median=median,
                spread=spread,
                stable=spread <= 0.15,
                sample_orders=orders,
            )
        if all(r["stable"] for r in results.values()) or target >= 45:
            return results
        target = min(45, target + 9)


def analyze_capture(path):
    """Conditional paired-round bootstrap; fixed inputs are not resampled."""
    capture = json.loads(source_path(path).read_text())
    results = capture["results"]
    output = dict(
        capture_sha256=hashlib.sha256(
            source_path(path).read_bytes()
        ).hexdigest(),
        bootstrap_seed=20261009,
        bootstrap_draws=10000,
        comparisons={},
    )
    for key, result in results.items():
        name, cohort = key.split("/")
        for baseline in ("control", "accepted"):
            baseline_key = baseline + "/" + cohort
            if name == baseline or baseline_key not in results:
                continue
            before = results[baseline_key]["samples"]
            after = result["samples"]
            if len(before) != len(after) or len(before) < 9:
                raise ValueError("unmatched A10 sample rounds")
            rng = random.Random(20261009)
            savings = []
            for _ in range(10000):
                indices = [rng.randrange(len(before)) for _ in before]
                savings.append(
                    100
                    * (
                        1
                        - statistics.median(after[i] for i in indices)
                        / statistics.median(before[i] for i in indices)
                    )
                )
            savings.sort()
            output["comparisons"][name + "_vs_" + baseline + "/" + cohort] = (
                dict(
                    median_saving_percent=100
                    * (1 - result["median"] / results[baseline_key]["median"]),
                    conditional_95_percent=[savings[250], savings[9749]],
                    stable=result["stable"]
                    and results[baseline_key]["stable"],
                    samples=len(before),
                )
            )
    return output


def cold_startups(repetitions):
    """Measure process startup/import/one classification separately."""
    root = REPOSITORY_ROOT
    before = root / "v2/audit/a10/control"
    data = json.loads(source_path(BASELINE).read_text())
    FrozenSources(data)
    for path, source in data["sources"].items():
        output = before / path
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(source)
    results = {name: [] for name in ("control", "accepted")}
    order_rng = random.Random(20261009)
    for _ in range(repetitions):
        order = list(results)
        order_rng.shuffle(order)
        for name in order:
            cwd = before if name == "control" else root
            expected = (
                "probable_prime" if name == "control" else "proven_prime"
            )
            code = (
                "import random\nfrom v2.common.utils import classify_prime\n"
                f"result = classify_prime({REPORTED_PRIME}, "
                "rng=random.Random(7)).value\n"
                f"if result != {expected!r}:\n"
                "    raise AssertionError('invalid cold classification')\n"
            )
            env = dict(
                os.environ, PYTHONPATH=str(cwd), PYTHONDONTWRITEBYTECODE="1"
            )
            started = time.perf_counter()
            subprocess.run(
                [sys.executable, "-B", "-c", code],
                cwd=cwd,
                env=env,
                check=True,
                capture_output=True,
                timeout=30,
            )
            results[name].append(time.perf_counter() - started)
    return {
        name: dict(samples=samples, median=statistics.median(samples))
        for name, samples in results.items()
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--batch-seconds", type=float, default=0.1)
    parser.add_argument(
        "--split", choices=("training", "confirmation"), default="training"
    )
    parser.add_argument("--arms", nargs="+", default=["control", "accepted"])
    parser.add_argument("--factoring", action="store_true")
    parser.add_argument("--cold", action="store_true")
    parser.add_argument("--within-range", action="store_true")
    parser.add_argument("--analyze", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("use a fresh output path for A10 evidence")
    if not hasattr(sys, "pypy_version_info") or sys.version_info[:2] != (
        3,
        11,
    ):
        raise RuntimeError("A10 requires PyPy implementing Python 3.11")
    if (
        not math.isfinite(args.warmup_seconds)
        or args.warmup_seconds < 3
        or args.repetitions < 9
        or not math.isfinite(args.batch_seconds)
        or args.batch_seconds < 0.1
    ):
        raise ValueError("at least three seconds warmup and nine samples")
    if args.analyze:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(analyze_capture(args.analyze), indent=2) + "\n"
        )
        return
    data = load_corpus()
    protocol = load_protocol()
    control = load_control()
    arms = candidates()
    arms["control"] = control.utils
    output = dict(
        environment=environment(),
        protocol=protocol,
        arguments=vars(args).copy(),
        control_commit=json.loads(source_path(BASELINE).read_text())["commit"],
        corpus_sha256=hashlib.sha256(
            json.dumps(data, sort_keys=True).encode()
        ).hexdigest(),
        results={},
    )
    output["arguments"]["output"] = str(args.output)
    output["references"] = {}
    if "gmp_same_bases" in args.arms:
        import gmpy2

        output["references"]["gmp"] = dict(
            gmpy2=gmpy2.version(),
            gmp=gmpy2.mp_version(),
            extension=gmpy2.__file__,
        )
    if "sympy_same_bases" in args.arms:
        import sympy
        from sympy.external import gmpy as ground

        output["references"]["sympy"] = dict(
            version=sympy.__version__,
            module=sympy.__file__,
            ground_types=ground.GROUND_TYPES,
            integer_type=ground.MPZ.__module__ + "." + ground.MPZ.__name__,
            ground_override=os.environ.get("SYMPY_GROUND_TYPES"),
        )
    if args.cold:
        output["cold_startups"] = cold_startups(args.repetitions)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(output, indent=2) + "\n")
        return
    if "sympy_same_bases" in args.arms and not args.within_range:
        raise ValueError("reference comparisons require --within-range")
    if args.factoring:
        workloads = {
            "complete": [
                f
                for f in data["factoring"]
                if f["split"] in (args.split, "regression")
            ]
        }
    else:
        workloads = {}
        for fixture in data["fixtures"]:
            if fixture["split"] == args.split and (
                not args.within_range
                or fixture["n"] < utils.DETERMINISTIC_LIMIT
            ):
                workloads.setdefault(fixture["cohort"], []).append(fixture)
    for cohort, fixtures in workloads.items():
        functions = {}
        for name in args.arms:
            if args.factoring:
                if name == "control":
                    mod = control.portfolio
                elif name == "accepted":
                    mod = portfolio
                elif name in ("native_pow", "binary_pow", "early_one"):
                    mod = factoring_candidate(name, arms[name])
                else:
                    raise ValueError(
                        "unsupported complete factoring candidate"
                    )
                budgets = types.SimpleNamespace(Budget=mod.Budget)
                functions[name] = (
                    lambda mod=mod, budgets=budgets: factoring_work(
                        mod, budgets, fixtures, data["seeds"]
                    )
                )
            elif name == "sympy_same_bases":
                functions[name] = lambda: sympy_work(fixtures)
            else:
                module = utils if name == "gmp_same_bases" else arms[name]
                functions[name] = (
                    lambda module=module, name=name: classification_work(
                        module,
                        fixtures,
                        104729,
                        native=name == "gmp_same_bases",
                    )
                )
        measured = measure_group(
            functions,
            args.warmup_seconds,
            args.repetitions,
            args.batch_seconds,
        )
        for name, result in measured.items():
            key = name + "/" + cohort
            output["results"][key] = result
            print(key, result["median"], flush=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    with performance_window():
        main()
