"""P3.1-P3.3 independent audit and matched M26/M27 correction costs."""

import argparse
import hashlib
import json
import random
import resource
import subprocess
import sys
import time
import types
from pathlib import Path

from .. import prime_sieve, utils, work_budget
from ..budget import Budget
from .phase_one import environment
from .phase_three_pipeline import _complete, _config, _corpus, _record
from .phase_three_reference import _rss_bytes

ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = ROOT / "v2/benchmarks/qs_m26_baseline.json"
MODULES = (
    "factor_base",
    "polynomial",
    "families",
    "relations",
    "reference_collector",
    "power_sieve",
    "sieve_collector",
    "linear_algebra",
    "extraction",
    "pipeline",
)


def load_arm(name, current=False):
    """Compile hash-checked owned sources in memory, preserving class types."""
    snapshot = json.loads(SNAPSHOT.read_text())
    package = types.ModuleType(name)
    package.__path__ = []
    package.__package__ = name
    package.utils = utils
    package.prime_sieve = prime_sieve
    package.budget = sys.modules[Budget.__module__]
    package.work_budget = work_budget
    sys.modules[name] = package
    for child in ("utils", "prime_sieve", "budget", "work_budget"):
        sys.modules[name + "." + child] = getattr(package, child)
    qs = types.ModuleType(name + ".qs")
    qs.__path__ = []
    qs.__package__ = qs.__name__
    package.qs = qs
    sys.modules[qs.__name__] = qs
    for child in MODULES:
        key = "v2/qs/" + child + ".py"
        if not current and key not in snapshot["source"]:
            continue
        source = (
            (ROOT / key).read_text() if current else snapshot["source"][key]
        )
        if (
            not current
            and hashlib.sha256(source.encode()).hexdigest()
            != (snapshot["source_sha256"][key])
        ):
            raise ValueError("corrupt M26 source snapshot")
        module = types.ModuleType(qs.__name__ + "." + child)
        module.__package__ = qs.__name__
        module.__file__ = str(ROOT / key) if current else str(SNAPSHOT)
        setattr(qs, child, module)
        sys.modules[module.__name__] = module
        exec(compile(source, module.__file__, "exec"), module.__dict__)
    return package


def oracle_audit(arm):
    """Check roots, collected payloads and lifted kernels independently."""
    from ..tests.test_qs import reference_positions, unlimited_budget
    from ..tests.test_qs_pipeline import dense_kernel, span

    qs = arm.qs
    counters = {"root_sets": 0, "collection_windows": 0, "kernels": 0}
    for n in (101 * 137, 211 * 307, 307 * 401):
        for h in (1, 3, 9):
            base = qs.factor_base.build_factor_base(
                n, multiplier=h, bound=100, budget=unlimited_budget()
            ).factor_base
            polynomials = (
                qs.polynomial.qs_polynomial(base),
                qs.polynomial.mpqs_polynomial(
                    base, 31, budget=unlimited_budget()
                ),
            )
            for polynomial in polynomials:
                for entry in base.entries:
                    roots = qs.polynomial.polynomial_roots(
                        polynomial, base, entry, budget=unlimited_budget()
                    )
                    actual = (
                        tuple(range(entry.prime))
                        if roots.all_positions
                        else roots.roots
                    )
                    expected = tuple(
                        x
                        for x in range(entry.prime)
                        if polynomial.value(x) % entry.prime == 0
                    )
                    assert actual == expected
                    counters["root_sets"] += 1
                expected = reference_positions(polynomial, base, -31, 34, 97)
                for backend in ("list", "bytearray", "array"):
                    for division in ("full", "roots", "bucket", "resieve"):
                        for policy in (
                            "adaptive",
                            "conservative",
                            "candidate",
                        ):
                            config = qs.sieve_collector.SieveConfig(
                                score_backend=backend,
                                division=division,
                                score_policy=policy,
                                block_width=17,
                                small_prime_cutoff=11,
                                residual_bound=97,
                                max_atoms=4096,
                                max_partials=4096,
                                max_relations=4096,
                                memory_bytes=33554432,
                            )
                            worker = qs.sieve_collector.SieveCollector(
                                polynomial,
                                base,
                                config=config,
                                budget=unlimited_budget(),
                            )
                            result = worker.collect(-31, 34)
                            assert result.reason == "complete"
                            actual = {
                                atom.position: (
                                    atom.sign,
                                    atom.exponents,
                                    atom.residual,
                                )
                                for atom in result.atoms
                            }
                            assert actual == expected
                            counters["collection_windows"] += 1
    generator = random.Random(2733)
    for _ in range(80):
        rows = tuple(
            generator.randrange(1 << 12)
            for _ in range(generator.randrange(1, 11))
        )
        expected = dense_kernel(rows)
        for weight_two in (False, True):
            matrix = qs.linear_algebra.filter_matrix(
                rows, weight_two=weight_two, budget=unlimited_budget()
            )
            for pivot in ("highest", "lowest"):
                solver = qs.linear_algebra.DependencySolver(
                    matrix, pivot=pivot, budget=unlimited_budget()
                )
                assert span(solver.run()) == expected
                counters["kernels"] += 1
    return {"passed": True, "seed": 2733, **counters}


def cohort(arm, division="bucket"):
    """Complete a fresh fixed-policy cohort including setup and consumption."""
    return _complete(
        arm,
        _corpus(2735, 16),
        _config(arm, division=division, residual_bound=1),
        weight_two=True,
        pivot="highest",
        row_excess=2,
        batch_width=256,
    )


def main():
    """Measure only validated results; keep failed old cap cases unranked."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--cold-arm", choices=("before", "after"))
    args = parser.parse_args()
    if args.cold_arm:
        arm = load_arm("_m27_cold", args.cold_arm == "after")
        answers = cohort(arm)
        usage = resource.getrusage(resource.RUSAGE_SELF)
        print(
            json.dumps(
                {
                    "answers": answers,
                    "peak_rss_bytes": _rss_bytes(),
                    "cpu_to_output_seconds": usage.ru_utime + usage.ru_stime,
                }
            )
        )
        return
    if args.output.exists() or args.warmup_seconds < 3 or args.repetitions < 9:
        parser.error("use a new capture, three-second warmup and nine samples")
    measured_environment = environment()
    arms = {
        name: load_arm("_m27_" + name, name == "after")
        for name in ("before", "after")
    }
    audit = oracle_audit(arms["after"])
    print("independent audit", audit, flush=True)
    expected = [case["factors"] for case in _corpus(2735, 16)]
    results = []
    for division in ("roots", "bucket"):
        for name, arm in arms.items():
            results.append(
                _record(
                    name + "_complete_" + division,
                    lambda arm=arm, division=division: cohort(arm, division),
                    expected,
                    args,
                )
            )
    rows = tuple((index * 37) % 256 for index in range(64))
    expected_kernel = None
    for name, arm in arms.items():

        def matrix_call(arm=arm):
            matrix = arm.qs.linear_algebra.filter_matrix(
                rows,
                weight_two=True,
                budget=Budget(work_limit=200000000),
            )
            solver = arm.qs.linear_algebra.DependencySolver(
                matrix, budget=Budget(work_limit=200000000)
            )
            return solver.run()

        value = matrix_call()
        if expected_kernel is None:
            expected_kernel = value
        assert value == expected_kernel
        results.append(_record(name + "_matrix", matrix_call, value, args))
    cold = []
    for name in arms:
        for _ in range(9):
            start = time.perf_counter()
            prior = resource.getrusage(resource.RUSAGE_CHILDREN)
            child = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "v2.benchmarks.phase_three_audit",
                    "--output",
                    str(args.output),
                    "--cold-arm",
                    name,
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=30,
            )
            sample = json.loads(child.stdout)
            assert sample["answers"] == expected
            usage = resource.getrusage(resource.RUSAGE_CHILDREN)
            sample.update(
                arm=name,
                lifecycle_seconds=time.perf_counter() - start,
                lifecycle_cpu_seconds=(
                    usage.ru_utime
                    + usage.ru_stime
                    - prior.ru_utime
                    - prior.ru_stime
                ),
            )
            cold.append(sample)
    assert (
        environment()["source_sha256"]
        == (measured_environment["source_sha256"])
    )
    args.output.write_text(
        json.dumps(
            {
                "milestone": "M27 / P3.1-P3.3 audit and corrections",
                "environment": measured_environment,
                "command": sys.orig_argv,
                "snapshot": str(SNAPSHOT),
                "snapshot_sha256": hashlib.sha256(
                    SNAPSHOT.read_bytes()
                ).hexdigest(),
                "independent_audit": audit,
                "results": results,
                "cold_samples": cold,
                "corpus": _corpus(2735, 16),
                "scope": "Matched complete QS and cold lifecycle",
                "ownership": "Reserve simultaneous private workspace",
                "rss_scope": "Warm RSS is shared; cold RSS is per child",
                "cap_bug": "Failed old output is unranked",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
