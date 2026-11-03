"""Disjoint A partitions with an exact core-envelope window objective."""

from math import prod

from .count_policy import CountPolicyRefusedError


def feasible_core_window(
    core_primes, lower, upper, count, size, target, budget
):
    """Choose an entire core pool before enumerating its distinct subsets.

    Every core product is between its sorted smallest/largest count-products.
    A pair from the two disjoint completion sets is between their endpoint
    products. The objective bounds continuous-envelope miss for every core;
    it does not assert an attainable product at every intermediate integer.
    """
    size = min(size, len(core_primes))
    if size < count or not lower or not upper:
        raise CountPolicyRefusedError(
            "windowed pair has insufficient disjoint prime sets"
        )
    best = None
    for start in range(len(core_primes) - size + 1):
        budget.consume((count + 3) * (target.bit_length() + 1))
        first_stop = start + count
        smallest = prod(core_primes[start:first_stop])
        stop = start + size
        last_start = stop - count
        largest = prod(core_primes[last_start:stop])
        miss = max(
            0,
            target - smallest * lower[-1] * upper[-1],
            largest * lower[0] * upper[0] - target,
        )
        center = core_primes[start + (size - 1) // 2]
        key = miss, abs(center ** (count + 2) - target), start
        if best is None or key < best:
            best = key
    start = best[2]
    stop = start + size
    return core_primes[start:stop], best[0]
