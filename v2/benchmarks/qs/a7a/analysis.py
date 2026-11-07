"""Paired modulus uncertainty and censored descriptive endpoints."""

import math
import random
import statistics
from collections import defaultdict


def quantile_bound(rows, probability):
    """Equal administrative caps: identify ranks only among completions."""
    if not rows or not 0 < probability <= 1:
        raise ValueError("nonempty finite quantile required")
    caps = {row["cap"] for row in rows}
    if len(caps) != 1:
        raise ValueError("stratify unequal censoring caps before quantiles")
    cap = caps.pop()
    rank = math.ceil(probability * len(rows))
    times = sorted(row["seconds"] for row in rows if row["complete"])
    if rank > len(times):
        return dict(value=None, lower_bound=cap, identified=False)
    return dict(value=times[rank - 1], identified=True)


def paired_interval(rows, *, nominee, resamples=4096):
    """Bootstrap independent moduli; repeated seeds stay in their cluster."""
    if not 100 <= resamples <= 10000:
        raise ValueError("finite bootstrap allowance required")
    groups = defaultdict(dict)
    for row in rows:
        key = row["fixture"], row["seed"], row["arm"]
        if key in groups[row["fixture"]]:
            raise ValueError("repeated observation identity")
        groups[row["fixture"]][key] = row
    differences, savings, losses = [], [], []
    for fixture, group in groups.items():
        seeds = {key[1] for key in group}
        completed, costs = [], []
        for seed in seeds:
            left = group.get((fixture, seed, "S0"))
            right = group.get((fixture, seed, nominee))
            if left is None or right is None:
                raise ValueError("missing paired attempt")
            losses.append(left["complete"] and not right["complete"])
            completed.append(int(right["complete"]) - int(left["complete"]))
            # Capped service cost retains failures; it is not latency among
            # successful calls and does not identify an uncensored quantile.
            costs.append(
                (left["seconds"] if left["complete"] else left["cap"])
                - (right["seconds"] if right["complete"] else right["cap"])
            )
        differences.append(statistics.mean(completed))
        savings.append(statistics.mean(costs))
    if not differences:
        raise ValueError("no independent pairs")
    rng = random.Random(202611040706)
    samples = []
    for _ in range(resamples):
        indices = [rng.randrange(len(savings)) for _ in savings]
        samples.append(statistics.mean(savings[i] for i in indices))
    samples.sort()
    count = len(groups)
    return dict(
        independent_moduli=count,
        capped_time_saving=statistics.mean(savings),
        saving_interval95=[
            samples[int(0.025 * (resamples - 1))],
            samples[int(0.975 * (resamples - 1))],
        ],
        completion_difference=statistics.mean(differences),
        # Zero observed losses cannot prove a zero population loss. With nine
        # independent pairs this one-sided upper bound is still about 28%.
        zero_loss_upper95=1 - 0.05 ** (1 / count) if not any(losses) else None,
        p90_scope=(
            "Nine moduli do not establish a precise p90 "
            "or zero-loss guardrail."
        ),
    )
