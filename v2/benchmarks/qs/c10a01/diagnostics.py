"""Four forced policies and two bounded capacity checks; no complete pilots."""

import argparse
import fcntl
import json
import random
import subprocess
import time
from pathlib import Path

from ....execution.budget import Budget, BudgetExhaustedError
from ....qs.extraction import prepare_relations
from ....qs.linear_algebra import (
    DependencySolver,
    filter_matrix,
    matrix_workspace,
    verify_dependency,
)
from ..c10a04.accounting import Ledger, digest, output_lock
from . import runner
from .configuration import OWNED, RSS, TIERS

PLAN = Path(__file__).with_name("diagnostics.json")


def matrix_rows(count, columns, kind):
    """Construct independent bounded fixtures with known duplicate kernels."""
    if kind == "cycle":
        return tuple(
            (1 << (i % columns)) | (1 << ((i + 1) % columns))
            for i in range(count)
        )
    generator = random.Random(2026110201)
    unique = tuple(
        sum(1 << j for j in generator.sample(range(columns), 35))
        for _ in range(count - 64)
    )
    return unique + unique[:64]


def matrix_check(spec, budget):
    """Check admission before allocation, then exact original-row kernels."""
    rows = matrix_rows(spec["rows"], spec["columns"], spec["kind"])
    reserve = matrix_workspace(len(rows), spec["columns"])
    refused = Budget(work_limit=10**12, seconds=None, cpu_seconds=None)
    try:
        filter_matrix(rows, budget=refused)
    except MemoryError:
        assert refused.used == 0
    else:
        raise AssertionError("matrix reserve did not reject the smaller grant")
    start, cpu, before = time.monotonic(), time.process_time(), budget.used
    matrix = filter_matrix(rows, weight_two=True, budget=budget)
    filtered = dict(
        wall=time.monotonic() - start,
        cpu=time.process_time() - cpu,
        work=budget.used - before,
        **matrix.stats,
    )
    solver = DependencySolver(matrix, budget=budget)
    start, cpu, before = time.monotonic(), time.process_time(), budget.used
    reason = "complete"
    try:
        solver.run()
    except BudgetExhaustedError:
        reason = budget.reason
    # Check the original rows directly, independently of elimination steps.
    for mask in solver.dependencies:
        assert verify_dependency(mask, rows)
    return dict(
        reason=reason,
        reserve=reserve,
        refused_below_reserve=True,
        filtering=filtered,
        elimination=dict(
            wall=time.monotonic() - start,
            cpu=time.process_time() - cpu,
            work=budget.used - before,
            next_row=solver.next_row,
            xors=solver.xors,
            dependencies=len(solver.dependencies),
        ),
    )


def memory_call(fixture, spec):
    """Retain collector/preparation while the synthetic matrix coexists."""
    digits = spec["digits"]
    start, cpu = time.monotonic(), time.process_time()
    row = runner.run_one(
        fixture["n"], digits, "control", 44, 370000000000, 264
    )
    runner.verify_outcome(row, fixture)
    job = runner.LAST_JOB
    budget = job.budget
    collector = job.engine.collector
    before, prep_start, prep_cpu = (
        budget.used,
        time.monotonic(),
        time.process_time(),
    )
    prepared = prepare_relations(
        collector.matrix_relations,
        job.base,
        collector._atoms,
        budget=budget,
        retained_workspace_bytes=collector._workspace,
    )
    preparation = dict(
        work=budget.used - before,
        wall=time.monotonic() - prep_start,
        cpu=time.process_time() - prep_cpu,
        rows=len(prepared.rows),
        shared=prepared.shared_workspace_bytes,
        private=prepared.workspace_bytes - prepared.shared_workspace_bytes,
    )
    simultaneous = (
        collector._workspace
        + preparation["private"]
        + matrix_workspace(spec["rows"], spec["columns"])
    )
    assert simultaneous < job.config.memory_bytes
    matrix = matrix_check(spec, budget)
    return dict(
        work=budget.used,
        wall=time.monotonic() - start,
        cpu=time.process_time() - cpu,
        setup_collection=row,
        preparation=preparation,
        matrix=matrix,
        simultaneous_native_reserve=simultaneous,
        tier_owned=OWNED[TIERS[digits]],
        rss_high_water_bytes=runner.rss(),
        rss_cap=RSS[TIERS[digits]],
        rss_scope="process including matrix and collection warmup",
    )


def memory_warm(fixture, spec):
    start, cpu = time.monotonic(), time.process_time()
    collection = runner.warm(fixture["n"], spec["digits"], "control")
    budget = Budget(
        work_limit=290000000000 - collection["work"],
        seconds=20,
        cpu_seconds=20,
    )
    matrix_start, matrix_cpu = time.monotonic(), time.process_time()
    small = dict(spec, rows=1024, columns=257)
    calls = 0
    while (
        min(time.monotonic() - matrix_start, time.process_time() - matrix_cpu)
        < 3
    ):
        result = matrix_check(small, budget)
        assert result["reason"] == "complete"
        calls += 1
    return dict(
        work=collection["work"] + budget.used,
        wall=time.monotonic() - start,
        cpu=time.process_time() - cpu,
        matrix_calls=calls,
        collection=collection,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("forced", "memory"))
    parser.add_argument("--id", required=True)
    args = parser.parse_args()
    runner.runtime()
    plan = json.loads(PLAN.read_text())
    spec = next(r for r in plan[args.action] if r["id"] == args.id)
    corpus = json.loads(runner.CORPUS.read_text())
    runner.validate(corpus)
    fixture = next(
        r
        for r in corpus["fixtures"]
        if r["digits"] == spec["digits"] and r["id"].endswith("_balanced")
    )
    source = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=runner.ROOT, text=True
    ).strip()
    subprocess.run(
        ["git", "diff", "--exit-code", "HEAD", "--", "v2"],
        cwd=runner.ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    identity = dict(
        source=source,
        plan=digest(PLAN),
        corpus=digest(runner.CORPUS),
        spec=spec,
    )
    print(json.dumps(runner.host_admission(spec["digits"])), flush=True)
    bootstrap = json.loads(runner.BOOTSTRAP.read_text())
    with open("/private/tmp/factor-performance.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with output_lock(runner.OUTPUT):
            ledger = Ledger(
                runner.OUTPUT, digest(runner.BOOTSTRAP), bootstrap["limits"]
            )
            tag = args.action + "-" + args.id
            if (runner.OUTPUT / ("warm-" + tag + ".json")).exists():
                raise RuntimeError("refuse replay of a process warmup")
            ledger.call(
                "warmup",
                "warm-" + tag,
                35,
                290000000000,
                identity,
                lambda: memory_warm(fixture, spec)
                if args.action == "memory"
                else runner.warm(
                    fixture["n"],
                    spec["digits"],
                    spec["arm"],
                    multiplier=spec["multiplier"],
                ),
            )
            result = ledger.call(
                args.action,
                tag,
                50,
                370000000000,
                identity,
                lambda: memory_call(fixture, spec)
                if args.action == "memory"
                else runner.run_one(
                    fixture["n"],
                    spec["digits"],
                    spec["arm"],
                    45,
                    370000000000,
                    33 * spec["polynomials"],
                    multiplier=spec["multiplier"],
                ),
            )
            if args.action == "forced":
                runner.verify_outcome(result, fixture)
                ledger.call(
                    "analysis",
                    tag + "-coordinates",
                    3,
                    10**10,
                    identity,
                    lambda: runner.describe(runner.LAST_JOB),
                )
            print(
                json.dumps(
                    {
                        k: v
                        for k, v in result.items()
                        if k
                        not in ("stats", "configuration", "setup_collection")
                    }
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
