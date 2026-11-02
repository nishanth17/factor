"""Historical C10a-04 schedule and preliminary stopping-rule checks."""

import random
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
CORPUS = ROOT / "v2/benchmarks/inputs/corpora/c10a04_short_fresh.json"
ARMS = ("control", "compact")
SEEDS = (7, 29, 47)


def schedule(corpus):
    """Cross every input with each seed before any outcome-dependent stop."""
    fixtures = corpus["fixtures"]
    if len(fixtures) != 3 or any(f["digits"] != 60 for f in fixtures):
        raise ValueError("expected three fresh 60-digit inputs")
    result = []
    for sample in range(9):
        order = list(ARMS)
        random.Random(104202613 + sample).shuffle(order)
        result.append(
            (sample, fixtures[sample % 3], SEEDS[sample // 3], order)
        )
    return result


def decision(rows):
    """Return descriptive trends, never a sequential significance claim."""
    groups = defaultdict(dict)
    for row in rows:
        if row["arm"] in groups[row["sample"]]:
            raise ValueError("duplicate cell")
        groups[row["sample"]][row["arm"]] = row
    if any(set(group) != set(ARMS) for group in groups.values()):
        raise ValueError("interim requires complete matched pairs")
    count = len(groups)
    if count not in (3, 6, 9):
        return "continue"

    pairs = [groups[index] for index in sorted(groups)]
    if not any(pair["compact"]["complete"] for pair in pairs):
        return "stop-no-compact-completions"
    if count == 6:
        by_input = defaultdict(list)
        for pair in pairs:
            by_input[pair["control"]["id"]].append(pair)
        if (
            len(by_input) == 3
            and all(len(group) == 2 for group in by_input.values())
            and all(pair["compact"]["complete"] for pair in pairs)
            and all(
                pair["compact"]["capped_seconds"]
                < pair["control"]["capped_seconds"]
                for pair in pairs
            )
            and all(
                sum(p["compact"]["capped_seconds"] for p in group)
                <= 0.8 * sum(p["control"]["capped_seconds"] for p in group)
                for group in by_input.values()
            )
        ):
            return "stop-preliminary-consistent-compact-trend"
    return "finished-nine-pairs" if count == 9 else "continue"
