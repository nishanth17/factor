"""Pure candidate ranking and explicit regime keys for F4 selection."""

import math
import random
import statistics

WORKER_COUNTS = (0, 1, 2, 4, 6, 8, 10, 12, 14, 16)


def regime_key(digits, tiers, backend="python-int"):
    """Use observable plans only; no factor labels or outcomes enter a key."""
    if type(digits) is not int or not 2 <= digits <= 1234:
        raise ValueError("invalid decimal-size key")
    plan = tuple(tuple(tier) for tier in tiers)
    if not plan or len(plan) > 32:
        raise ValueError("finite numerical schedule required")
    for tier in plan:
        if len(tier) != 3:
            raise ValueError("tier must contain bounds and quota")
        b1, b2, curves = tier
        if any(type(value) is not int for value in tier):
            raise ValueError("integer bounds and quotas required")
        if not 2 <= b1 <= b2 < 2**32 or not 0 <= curves <= 4096:
            raise ValueError("tier outside declared search space")
    if sum(tier[2] for tier in plan) > 4096:
        raise ValueError("curve quota exceeds executor contract")
    if backend != "python-int":
        raise ValueError("unmeasured backend")
    return digits, plan, backend


def matched_ratio(candidate, control, *, seed=20261108, draws=4000):
    """Conditional paired-block interval; this is not input independence."""
    if set(candidate) != set(control) or len(candidate) < 3:
        raise ValueError("matched timing blocks required")
    logs = []
    for key in sorted(candidate):
        pair = candidate[key], control[key]
        if any(
            type(value) not in (float, int)
            or not math.isfinite(value)
            or value <= 0
            for value in pair
        ):
            raise ValueError("finite positive timings required")
        logs.append(math.log(pair[0] / pair[1]))
    if type(draws) is not int or not 100 <= draws <= 10_000:
        raise ValueError("finite bootstrap allowance required")
    generator = random.Random(seed)
    samples = sorted(
        math.exp(statistics.mean(generator.choices(logs, k=len(logs))))
        for _ in range(draws)
    )
    return dict(
        ratio=math.exp(statistics.mean(logs)),
        interval=[
            samples[int(0.025 * draws)],
            samples[int(0.975 * draws) - 1],
        ],
        blocks=len(logs),
        scope="conditional timing blocks, not independent moduli",
    )


def shortlist(rows, *, target=4, relative_band=0.15):
    """Keep fast candidates and adjacent counts for interaction checking.

    Near ties expand beyond the target rather than being silently truncated.
    Three-sample screens cannot justify tight statistical elimination. Failed
    and censored candidates remain explicit, never fast successful medians.
    """
    if type(target) is not int or not 1 <= target <= 16:
        raise ValueError("finite shortlist size required")
    if not 0 <= relative_band <= 1:
        raise ValueError("invalid exploration band")
    valid = []
    for row in rows:
        if row.get("status") != "ok":
            continue
        value = row["seconds"]
        if (
            type(value) not in (float, int)
            or not math.isfinite(value)
            or value <= 0
            or type(row["count"]) is not int
            or row["count"] not in WORKER_COUNTS
        ):
            raise ValueError("invalid selection observation")
        valid.append(row)
    if not valid:
        return []
    ranked = sorted(valid, key=lambda row: (row["seconds"], row["count"]))
    best = ranked[0]
    winners = [
        row
        for row in ranked
        if row["seconds"] <= best["seconds"] * (1 + relative_band)
    ]
    for row in ranked[:target]:
        if row not in winners:
            winners.append(row)
    count_position = WORKER_COUNTS.index(best["count"])
    neighbors = WORKER_COUNTS[max(0, count_position - 1) : count_position + 2]
    # Neighbor coverage is a declared expansion, not an automatic rejection.
    for count in neighbors:
        choices = [row for row in ranked if row["count"] == count]
        if choices and not any(row["count"] == count for row in winners):
            winners.append(choices[0])
    serial = next((row for row in ranked if row["count"] == 0), None)
    if serial is not None and not any(row["count"] == 0 for row in winners):
        winners.append(serial)
    return winners


def censored_quantile(observations, probability):
    """Report a common-grant quantile only when completion identifies it."""
    if not observations or not 0 < probability <= 1:
        raise ValueError("observations and probability required")
    grants = {row["grant"] for row in observations}
    if len(grants) != 1:
        raise ValueError("this summary requires a common latency grant")
    grant = grants.pop()
    if not math.isfinite(grant) or grant <= 0:
        raise ValueError("invalid latency grant")
    completed = []
    for row in observations:
        if not math.isfinite(row["seconds"]) or row["seconds"] < 0:
            raise ValueError("invalid latency observation")
        if row["complete"]:
            if row["seconds"] > grant:
                raise ValueError("completed outcome exceeds grant")
            completed.append(row["seconds"])
    position = max(0, math.ceil(probability * len(observations)) - 1)
    complete = sorted(completed)
    identified = position < len(complete)
    return dict(
        seconds=complete[position] if identified else None,
        lower_bound=None if identified else grant,
        completion=len(complete) / len(observations),
        inputs=len(observations),
    )
