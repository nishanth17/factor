"""Bounded fixed-prefix replay of C1's complete retained census."""

import argparse
import hashlib
import json
import time
from collections import Counter
from math import prod
from pathlib import Path

from ....execution.budget import Budget
from ....qs.large_primes import LargePrimeForest
from ..b1.b1_calibration import save
from ..c1.c1_followup import compact_report, validate_records
from . import c9_study as study


def replay(records, n, aligned, budget):
    """Compare finite FIFO retention and optimistic no-eviction rank."""
    forest, kept, parents, sizes = LargePrimeForest(), set(), {}, {}
    counts = Counter()

    def root(v):
        while parents[v] != v:
            budget.consume(1)
            v = parents[v]
        return v

    for index, record in enumerate(records):
        budget.consume(1)
        lp = record["lp"]
        if len(lp) == 1 and lp[0] > 3000000:
            counts["above_endpoint_slp"] += 1
            if aligned:
                continue
        if len(lp) == 2 and (
            max(lp) > 3000000 or record["residual"] > 115200000000
        ):
            raise ValueError(
                "census outside frozen C1 product/endpoint domain"
            )
        left, right = tuple([1] * (2 - len(lp)) + lp)
        for v in (left, right):
            if v not in parents:
                parents[v], sizes[v] = v, 1
        a, b = root(left), root(right)
        if a == b:
            counts["no_eviction_cycle_rank"] += 1
        else:
            if sizes[a] < sizes[b]:
                a, b = b, a
            parents[b] = a
            sizes[a] += sizes.pop(b)
        identity = str(index)
        plan = forest.plan(left, right, identity, 8192, budget)
        if plan.cycle is not None:
            kept.add(index)
            kept.update(int(i) for i in plan.cycle)
            counts["emitted_cycles"] += 1
            counts["cycle_atoms"] += len(plan.cycle) + 1
        counts["evictions"] += len(plan.evicted)
        forest.commit(plan)
    selected = kept | {int(i) for i in forest.edges}
    # Input records stay live while forest references and analysis coexist.
    # Reserve raw-record objects, two forest maps and replay union-find first.
    retained_reserve = sum(2048 + 128 * len(r["exponents"]) for r in records)
    replay_reserve = 32768 + 4096 * len(forest.edges) + 512 * len(parents)
    available = 512 * 2**20 - retained_reserve - replay_reserve
    report = compact_report(
        [records[i] for i in sorted(selected)],
        n,
        budget=budget,
        memory_bytes=available,
        first_factor=True,
    )
    return dict(
        counts=counts,
        retained=len(selected),
        raw_record_reserve=retained_reserve,
        replay_reserve=replay_reserve,
        pipeline=report,
    )


def validate_report(report, fixture):
    """Check every offline outcome against independently certified fixtures."""
    for arm in ("current", "aligned"):
        result = report[arm]["pipeline"]
        factors, remaining = result["factors"], result["remaining"]
        if prod(factors) * prod(remaining) != fixture["n"]:
            raise AssertionError("offline census did not reconstruct")
        if Counter(factors) - Counter(fixture["factors"]):
            raise AssertionError("offline factors disagree with corpus")
        if result["complete"] != (not remaining):
            raise AssertionError("offline completeness disagrees with output")
        for divisor in result.get("proper_divisors", ()):
            if not study.utils.valid_divisor(divisor, fixture["n"]):
                raise AssertionError(
                    "offline census returned improper divisor"
                )


def run(retained, output):
    start, cpu = time.monotonic(), time.process_time()
    control = json.loads(study.CONTROL.read_text())
    reports = []
    budget = Budget(work_limit=10**13, seconds=120, cpu_seconds=120)
    for index, fixture in enumerate(control["fixtures"]):
        path = retained / f"50-{index}-records.json"
        if path.stat().st_size > 32 * 2**20:
            raise ValueError("retained census byte cap")
        records = json.loads(path.read_text())
        validate_records(records, fixture["n"], budget)
        high = Counter(
            r["residual"]
            for r in records
            if len(r["lp"]) == 1 and r["residual"] > 3000000
        )
        reports.append(
            dict(
                id=fixture["id"],
                source=str(path),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                records=len(records),
                last_block=max(r["block"] for r in records),
                high_slp_records=sum(high.values()),
                high_slp_endpoints=len(high),
                high_slp_independent_pairs=sum(v - 1 for v in high.values()),
                current=replay(records, fixture["n"], False, budget),
                aligned=replay(records, fixture["n"], True, budget),
            )
        )
        validate_report(reports[-1], fixture)
        print(
            json.dumps(
                {
                    k: reports[-1][k]
                    for k in (
                        "id",
                        "records",
                        "high_slp_records",
                        "high_slp_independent_pairs",
                    )
                }
            ),
            flush=True,
        )
    save(
        output,
        dict(
            reports=reports,
            seconds=time.monotonic() - start,
            cpu_seconds=time.process_time() - cpu,
            work=budget.used,
            limitation="Retained C1 fixed-prefix counterfactual; acquisition, "
            "scoring and full-call cost are not measured here. "
            "No-eviction rank is optimistic, not a memory-feasible arm.",
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retained", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    study.runtime()
    with study.machine_window():
        run(args.retained, args.output)


if __name__ == "__main__":
    main()
