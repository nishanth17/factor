"""Observable complete-cost estimates for the opt-in ECM engagement rule."""

import math


def calibrated_costs(
    tiers,
    *,
    serial_cold_body,
    serial_cold_cli,
    serial_warm,
    parallel_cold_cli,
    parallel_warm,
    parallel_warm_parent,
    pool_start,
    retained_calls=9,
):
    """Fit local fixed/variable costs; never inspect search outcomes.

    Cold and warm estimates use their measured serial full-manifest cost.
    Parent CPU plus measured serialization/planning wall time supplies a
    warm-overhead estimate. The remaining parallel cost is fitted as variable
    work. This is an engagement model, not a success-probability model or a
    claim that short prefixes scale linearly. Whole-call confirmation remains
    necessary, especially for early factors and recursive completion.
    """
    values = (
        serial_cold_body,
        serial_cold_cli,
        serial_warm,
        parallel_cold_cli,
        parallel_warm,
        parallel_warm_parent,
        pool_start,
    )
    if any(
        type(value) not in (int, float)
        or not math.isfinite(value)
        or value < 0
        for value in values
    ):
        raise ValueError("finite nonnegative measured costs required")
    if min(serial_cold_body, serial_warm, parallel_warm) <= 0:
        raise ValueError("positive validated work timings required")
    if serial_cold_cli < serial_cold_body:
        raise ValueError("serial CLI cost cannot exclude its factoring body")
    if type(retained_calls) is not int or not 1 <= retained_calls <= 64:
        raise ValueError("finite declared retained-call horizon required")
    if not isinstance(tiers, (list, tuple)) or not 1 <= len(tiers) <= 32:
        raise ValueError("bounded numerical tiers required")
    for tier in tiers:
        if (
            not isinstance(tier, (list, tuple))
            or len(tier) != 3
            or any(type(value) is not int for value in tier)
            or not 2 <= tier[0] <= tier[1] < 2**32
            or not 0 <= tier[2] <= 4096
        ):
            raise ValueError("invalid numerical tier")
    curves = sum(tier[2] for tier in tiers)
    if not 1 <= curves <= 4096:
        raise ValueError("positive finite planned curve quota required")

    # Remove common interpreter/import/output cost from the parallel intercept.
    # Cold calls still pay worker spawn, metadata, JIT and shutdown in full.
    shared_cli = serial_cold_cli - serial_cold_body
    warm_overhead = min(parallel_warm_parent, parallel_warm)
    cold_overhead = max(
        warm_overhead,
        parallel_cold_cli - parallel_warm + warm_overhead - shared_cli,
    )
    started_overhead = max(warm_overhead, cold_overhead - pool_start)
    variable = max(0.0, parallel_warm - warm_overhead)
    result = {}
    for mode, serial, horizon in (
        ("cold", serial_cold_body, 1),
        ("warm", serial_warm, retained_calls),
    ):
        slope = variable / serial
        denominator = 1 - slope

        def threshold(overhead):
            if denominator <= 0:
                return None
            value = overhead / denominator
            if not math.isfinite(value):
                return None
            return max(0.001, value)

        # Declared retained reuse amortizes only preparation above warm cost;
        # the actual first call remains cold and is recorded/charged as cold.
        cold = warm_overhead + (cold_overhead - warm_overhead) / horizon
        pool = warm_overhead + (started_overhead - warm_overhead) / horizon
        result[mode] = dict(
            seconds_per_curve=tuple(serial / curves for _ in tiers),
            cold_seconds=threshold(cold),
            warm_seconds=threshold(warm_overhead),
            pool_seconds=threshold(pool),
        )
    return result
