"""Matched modulus-level uncertainty with explicit censored outcomes."""

import math
import random
import statistics

from .selection import censored_quantile


def _validate_arms(candidate, control):
    if set(candidate) != set(control) or not candidate:
        raise ValueError("matched independent modulus identities required")
    for identity in candidate:
        left, right = candidate[identity], control[identity]
        if not left or len(left) != len(right) or len(left) > 100:
            raise ValueError("matched finite timing repetitions required")
        for observations in (left, right):
            for row in observations:
                if set(row) != {"seconds", "complete", "grant"}:
                    raise ValueError("invalid latency observation")
                if type(row["complete"]) is not bool:
                    raise ValueError("invalid completion flag")
                if (
                    type(row["seconds"]) not in (float, int)
                    or not math.isfinite(row["seconds"])
                    or row["seconds"] < 0
                ):
                    raise ValueError("invalid observed latency")
                if (
                    type(row["grant"]) not in (float, int)
                    or not math.isfinite(row["grant"])
                    or row["grant"] <= 0
                ):
                    raise ValueError("invalid common grant")
                if row["complete"] and row["seconds"] > row["grant"]:
                    raise ValueError("completed observation exceeds grant")
    grants = {
        row["grant"]
        for arms in (candidate, control)
        for observations in arms.values()
        for row in observations
    }
    if len(grants) != 1:
        raise ValueError("one latency comparison needs equal grants")


def _completion_and_restricted(arm, identities):
    # Each independent modulus has equal weight, regardless of sample count.
    completion, restricted = [], []
    for identity in identities:
        observations = arm[identity]
        completion.append(
            statistics.mean(row["complete"] for row in observations)
        )
        restricted.append(
            statistics.mean(
                row["seconds"] if row["complete"] else row["grant"]
                for row in observations
            )
        )
    return statistics.mean(completion), statistics.mean(restricted)


def compare_cohort(candidate, control, *, seed=20261110, draws=2000):
    """Resample moduli, keeping their seeds/repetitions together.

    Restricted latency charges failures/timeouts the common grant. It is a
    finite-horizon endpoint, not an estimate of unknown completion time. Raw
    quantiles stay censored when completion cannot identify them.
    """
    _validate_arms(candidate, control)
    if type(draws) is not int or not 100 <= draws <= 10_000:
        raise ValueError("finite bootstrap draw count required")
    identities = sorted(candidate)
    if len(identities) < 3:
        raise ValueError("at least three independent moduli required")
    # Keeping equal repetition counts also gives each modulus equal weight
    # in the pooled descriptive latency quantiles below.
    counts = {len(observations) for observations in candidate.values()}
    counts.update(len(observations) for observations in control.values())
    if len(counts) != 1:
        raise ValueError("cohort quantiles require equal repetition counts")
    c_completion, c_restricted = _completion_and_restricted(
        candidate, identities
    )
    s_completion, s_restricted = _completion_and_restricted(
        control, identities
    )
    if s_restricted <= 0:
        raise ValueError("zero control latency cannot define a ratio")
    # Cache per-modulus summaries: repetitions remain one cluster, and the
    # bootstrap never repeats their decoding or mean calculations.
    clusters = {
        identity: (
            _completion_and_restricted(candidate, [identity]),
            _completion_and_restricted(control, [identity]),
        )
        for identity in identities
    }
    generator = random.Random(seed)
    differences, ratios = [], []
    for _ in range(draws):
        chosen = generator.choices(identities, k=len(identities))
        complete_c = sum(clusters[key][0][0] for key in chosen)
        complete_s = sum(clusters[key][1][0] for key in chosen)
        restricted_c = sum(clusters[key][0][1] for key in chosen)
        restricted_s = sum(clusters[key][1][1] for key in chosen)
        differences.append((complete_c - complete_s) / len(identities))
        if restricted_s <= 0:
            raise ValueError("bootstrap control latency is zero")
        ratios.append(restricted_c / restricted_s)
    differences.sort()
    ratios.sort()
    lower, upper = int(0.025 * draws), int(0.975 * draws) - 1
    arms = {}
    for name, arm, completion, restricted in (
        ("candidate", candidate, c_completion, c_restricted),
        ("control", control, s_completion, s_restricted),
    ):
        observations = [
            row for identity in identities for row in arm[identity]
        ]
        arms[name] = dict(
            completion=completion,
            restricted_mean_seconds=restricted,
            median=censored_quantile(observations, 0.5),
            p90=censored_quantile(observations, 0.9),
        )
    return dict(
        independent_moduli=len(identities),
        repetitions_per_modulus=counts.pop(),
        scope="modulus-cluster bootstrap; censored finite-horizon latency",
        completion_difference=c_completion - s_completion,
        completion_interval=[differences[lower], differences[upper]],
        restricted_latency_ratio=c_restricted / s_restricted,
        restricted_ratio_interval=[ratios[lower], ratios[upper]],
        arms=arms,
        tail_warning=(
            "Small cohort: descriptive p90, substantial tail uncertainty."
            if len(identities) < 100
            else "Descriptive p90 requires tail/censoring uncertainty."
        ),
    )
