"""Offline empirical stopping sequences; no hidden-factor runtime inputs."""

import math

LADDERS = {
    "wide": ((2000, 147396, 64),),
    "compact": ((2000, 50000, 64),),
    "deep": ((2000, 147396, 16), (11000, 250000, 4)),
}
ENDPOINTS = {
    "wide": (8, 16, 32, 64),
    "compact": (8, 16, 32, 64),
    "deep": (8, 16, 17, 18, 20),
}


def prefix_tiers(ladder, count):
    """Compile an observed prefix into ordinary independent finite tiers."""
    if type(count) is not int or count < 0:
        raise ValueError("curve count must be a nonnegative integer")
    tiers = []
    for b1, b2, curves in LADDERS[ladder]:
        taken = min(count, curves)
        if taken:
            tiers.append((b1, b2, taken))
        count -= taken
    if count:
        raise ValueError("prefix exceeds observed ladder")
    return tuple(tiers)


def fit_stopping_sequence(records, endpoints, *, minimum_subjects=12):
    """Compare SIQS now with a measured block and its failure continuation.

    Records describe first useful root splits, not independent Bernoulli
    trials. Their weights come from the frozen offline workload. Probabilities
    are conditional on surviving every preceding curve. This empirical
    recurrence can justify a low-yield block whose later continuation pays;
    it never combines outcomes from unobserved mixed-bound sequences.
    """
    if not records or tuple(sorted(set(endpoints))) != tuple(endpoints):
        raise ValueError("nonempty records and increasing endpoints required")
    if not endpoints or endpoints[0] <= 0 or minimum_subjects < 1:
        raise ValueError("positive endpoints and subject floor required")
    for record in records:
        hit = record["hit"]
        required = endpoints[-1] if hit is None else hit
        numbers = (
            record["weight"],
            record["qs_cost"],
            record["recursive_cost"],
            record["setup"],
            *record["curve_costs"],
        )
        if (
            record["weight"] <= 0
            or any(not math.isfinite(x) or x < 0 for x in numbers)
            or len(record["curve_costs"]) != required
            or hit is not None
            and not 1 <= hit <= endpoints[-1]
        ):
            raise ValueError("invalid or censored sequence observation")

    states, next_cost = [], 0.0
    starts = (0, *endpoints[:-1])
    for start, end in reversed(tuple(zip(starts, endpoints))):
        risk = [r for r in records if r["hit"] is None or r["hit"] > start]
        if not risk:
            states.append(dict(start=start, end=end, choice="siqs", cost=0.0))
            next_cost = 0.0
            continue
        weight = sum(r["weight"] for r in risk)
        success = [r for r in risk if r["hit"] is not None and r["hit"] <= end]
        success_weight = sum(r["weight"] for r in success)
        probability = success_weight / weight
        stop_cost = sum(r["weight"] * r["qs_cost"] for r in risk) / weight
        block_cost = (
            sum(
                r["weight"]
                * (
                    sum(r["curve_costs"][start:end])
                    + (r["setup"] if start == 0 else 0)
                )
                for r in risk
            )
            / weight
        )
        recovery_cost = (
            sum(r["weight"] * r["recursive_cost"] for r in success) / weight
        )
        if end == endpoints[-1]:
            failed = [r for r in risk if r not in success]
            failed_weight = weight - success_weight
            next_cost = (
                sum(r["weight"] * r["qs_cost"] for r in failed) / failed_weight
                if failed_weight
                else 0.0
            )
        continuation = (
            block_cost + recovery_cost + (1 - probability) * next_cost
        )
        subjects = len({r["subject"] for r in risk})
        proceed = subjects >= minimum_subjects and continuation < stop_cost
        next_cost = continuation if proceed else stop_cost
        states.append(
            dict(
                start=start,
                end=end,
                subjects=subjects,
                conditional_success=probability,
                block_cost=block_cost,
                success_recovery_cost=recovery_cost,
                siqs_cost=stop_cost,
                continuation_cost=continuation,
                choice="continue" if proceed else "siqs",
                cost=next_cost,
            )
        )
    states.reverse()
    count = 0
    for state in states:
        if state["choice"] == "siqs":
            break
        count = state["end"]
    return dict(curves=count, expected_cpu=states[0]["cost"], states=states)


def fit_bands(records_by_ladder):
    """Fit two observable-size tables; factor labels only supplied weights."""
    result = {}
    for band in (30, 40):
        alternatives = {}
        for ladder, records in records_by_ladder.items():
            selected = [r for r in records if r["band"] == band]
            if selected:
                alternatives[ladder] = fit_stopping_sequence(
                    selected, ENDPOINTS[ladder]
                )
        if not alternatives:
            raise ValueError("missing calibration band")
        winner = min(
            alternatives,
            key=lambda name: (
                alternatives[name]["expected_cpu"],
                alternatives[name]["curves"],
                name,
            ),
        )
        count = alternatives[winner]["curves"]
        result[str(band)] = dict(
            ladder=winner,
            curves=count,
            tiers=prefix_tiers(winner, count),
            alternatives=alternatives,
        )
    return result
