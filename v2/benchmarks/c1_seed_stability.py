"""Predeclared fixed-seed assessment; preserve pooled flags and all samples."""

import json
from pathlib import Path

from .c1_implementation import SEEDS, unstable


def assess(directory):
    groups = {}
    for path in sorted(Path(directory).glob("c1_fresh_*.json")):
        if "warmup" in path.name or "stability" in path.name:
            continue
        row = json.loads(path.read_text())
        if "sample" not in row:
            continue
        groups.setdefault(row["digits"], {}).setdefault(
            row["id"], {}
        ).setdefault(row["arm"], []).append(row)
    report = {}
    for digits, cases in groups.items():
        strata, costs = {}, {}
        for identity, arms in cases.items():
            strata[identity] = {}
            for arm, rows in arms.items():
                by_seed = {}
                for seed in SEEDS:
                    selected = sorted(
                        (r for r in rows if r["seed"] == seed),
                        key=lambda r: r["sample"],
                    )
                    by_seed[str(seed)] = dict(
                        count=len(selected),
                        stable=len(selected) >= 9 and not unstable(selected),
                    )
                    costs.setdefault(arm, {}).setdefault(seed, []).extend(
                        r["capped_seconds"] for r in selected
                    )
                pooled = unstable(rows) if len(rows) >= 9 else True
                conditional = all(item["stable"] for item in by_seed.values())
                strata[identity][arm] = dict(
                    pooled_unstable=pooled,
                    fixed_seed=by_seed,
                    conditional_stable=conditional,
                    stable_under_clarification=not pooled or conditional,
                )
        directions = {}
        for arm in costs:
            if arm == "slp":
                continue
            directions[arm] = {
                str(seed): 1 - sum(costs[arm][seed]) / sum(costs["slp"][seed])
                for seed in SEEDS
                if costs[arm][seed] and costs["slp"][seed]
            }
        report[str(digits)] = dict(
            strata=strata, per_seed_capped_cost_benefit=directions
        )
    return report
