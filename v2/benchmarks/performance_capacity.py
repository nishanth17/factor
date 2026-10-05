"""Reproduce inspected worker refusals and isolate norm-bound reservations."""

import argparse
import ast
import json
from math import isqrt, prod
from pathlib import Path
from types import SimpleNamespace

from ..budget import Budget
from ..qs import parallel
from ..qs.relations import verify_atomic
from . import performance_audit as audit


def verify_trial(corpus):
    """Independently prove every small prime in the historical P3.6 set."""
    seen = set()

    for fixture in corpus["fixtures"]:
        for factor in fixture["factors"]:
            if not 2 <= factor < 2**24 or any(
                factor % d == 0 for d in range(2, isqrt(factor) + 1)
            ):
                raise AssertionError("invalid historical trial proof")

        if prod(fixture["factors"]) != fixture["n"] or fixture["n"] in seen:
            raise AssertionError("invalid historical reconstruction")
        seen.add(fixture["n"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new capture path")
    if audit.sys.implementation.name != "pypy":
        raise RuntimeError("PyPy implementing Python 3.11 is required")
    if audit.sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 is required")
    root = Path(__file__).resolve().parents[2]
    before = audit.fingerprint(root)
    corpus = (
        Path(__file__).parent / "inputs/corpora/phase_three_p36_corpus.json"
    )
    audit.verify_corpus = verify_trial
    audit.run(
        SimpleNamespace(
            runtime_root=root,
            baseline_json=None,
            corpus=corpus,
            split="held_out",
            bands="medium",
            repetitions=9,
            warmup_seconds=3,
            profile=False,
            polynomial_order="generated",
            revert="none",
            cases="native_base1k,parallel_base1k_serial,"
            "parallel_base1k_cap_serial,parallel_base1k_chunk_serial",
            output=args.output,
        )
    )
    data = json.loads(args.output.read_text())
    fixture = next(
        f
        for f in json.loads(corpus.read_text())["fixtures"]
        if f["band"] == "medium" and f["split"] == "held_out"
    )
    # This frozen post-chunk/pre-norm implementation differs only in the
    # reservation being compared; workers still execute the current code.
    snapshot = json.loads(
        (
            Path(__file__).parent
            / "inputs/baselines"
            / "p38_r1_measured_sources.json"
        ).read_text()
    )
    source = snapshot["source"]["v2/qs/parallel.py"]
    expected = snapshot["source_sha256"]["v2/qs/parallel.py"]
    if audit.hashlib.sha256(source.encode()).hexdigest() != expected:
        raise ValueError("corrupt dense-reservation control")
    owner = next(
        n
        for n in ast.parse(source).body
        if isinstance(n, ast.ClassDef) and n.name == "ParallelSIQSJob"
    )
    method = next(
        n
        for n in owner.body
        if isinstance(n, ast.FunctionDef) and n.name == "_memory"
    )
    namespace = {}
    exec(
        compile(
            ast.Module(body=[method], type_ignores=[]),
            "immutable-dense-reservation",
            "exec",
        ),
        parallel.__dict__,
        namespace,
    )
    current = parallel.ParallelSIQSJob._memory
    config = parallel.ParallelConfig(
        base_bound=10000,
        half_width=512,
        factor_count=3,
        batch_width=256,
        family_count=4,
        pool_size=16,
        max_batch_atoms=1024,
    )
    data["reservation_control_sha256"] = expected
    data["reservations"] = {}

    try:
        with parallel.CollectionPool("process", 4) as pool:
            for name, method in (
                ("dense", namespace["_memory"]),
                ("norm_bound", current),
            ):
                parallel.ParallelSIQSJob._memory = method

                def run():
                    budget = Budget(
                        work_limit=10**10, seconds=5, cpu_seconds=5
                    )
                    job = parallel.ParallelSIQSJob(
                        fixture["n"], config=config, budget=budget
                    )

                    result = job.run(
                        pool=pool, max_assignments=1, fixed_work=True
                    )

                    if (result.divisor or 1) * result.cofactor != job.n:
                        raise AssertionError("capacity reconstruction failed")
                    for atom in job.engine.collector._atoms.values():
                        verify_atomic(atom, job.base, residual_bound=10000)
                    return [dict(reason=result.reason, stats=result.stats)]

                data["reservations"][name] = audit.measure(
                    run, SimpleNamespace(warmup_seconds=3, repetitions=9)
                )
    finally:
        parallel.ParallelSIQSJob._memory = current

    if before != audit.fingerprint(root):
        raise AssertionError("source changed during measurement")
    args.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
