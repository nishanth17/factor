"""A finite joint starting design, not an exhaustive or optimality claim."""

from itertools import combinations, product

from ....qs.collector_kernels import combine

# Disjoint two-dimensional subspaces give balanced four-level factors.
# The remaining nonzero contrasts give balanced binary factors. Every pair
# of raw factors is covered; higher interactions remain confounded.
CATEGORICAL = {
    "workers": ((1, 2), (4, 6, 8, 10)),
    "packet": ((4, 8), (4, 8, 16, 64)),
    "verification": ((16, 5), ("parent", "worker", "batch64", "batch256")),
    "stop": ((6, 9), ("legacy", "shared", "private_event", "semaphore")),
}
BINARY = (
    "reuse_state",
    "packed_transport",
    "coefficient_cache",
    "cpu_snapshot",
    "arrival_order",
    "lookahead",
    "density",
    "presieve",
    "unroll4",
    "scan_fused",
    "parity_two",
    "packed_queue",
    "wide_hits",
    "wide_inverses",
    "wide_transfers",
    "pruned_transfers",
    "power_plan",
    "sparse_queue",
    "large_nursery",
)


def initial_design():
    """Return 32 explicit configurations for the retained 64-polynomial family.

    This first design uses two levels for binary/buffer options. Later
    refinements retain the broader declared parameter space. CPU snapshots
    are inactive without packet leases; other inactive mechanisms remain
    visible in the design rather than being counted as executed coverage.
    """
    occupied = {
        mask
        for (left, right), _ in CATEGORICAL.values()
        for mask in (left, right, left ^ right)
    }
    masks = [mask for mask in range(1, 32) if mask not in occupied]
    if len(masks) != len(BINARY):
        raise AssertionError("joint design contrast allocation changed")

    rows = []
    for index in range(32):
        levels = {}
        for name, ((left, right), values) in CATEGORICAL.items():
            level = (index & left).bit_count() % 2
            level += 2 * ((index & right).bit_count() % 2)
            levels[name] = values[level]
        levels.update(
            (name, bool((index & mask).bit_count() % 2))
            for name, mask in zip(BINARY, masks)
        )
        stop = levels["stop"]
        verification = levels["verification"]
        features = dict(
            lease_budget=stop != "legacy",
            reuse_state=levels["reuse_state"],
            packed_transport=levels["packed_transport"],
            verification="batch"
            if verification.startswith("batch")
            else verification,
            verification_batch_atoms=256 if verification == "batch256" else 64,
            lookahead=8 if levels["lookahead"] else 0,
            whole_family=levels["packet"] == 64,
            private_stop=stop in ("private_event", "semaphore"),
            stop_backend="semaphore" if stop == "semaphore" else "event",
            cpu_snapshot=levels["cpu_snapshot"] and stop != "legacy",
            arrival_order=levels["arrival_order"],
            coefficient_cache=levels["coefficient_cache"],
        )
        marking = combine(
            name
            for name in BINARY
            if name
            in (
                "density",
                "presieve",
                "unroll4",
                "scan_fused",
                "parity_two",
                "packed_queue",
                "wide_hits",
                "wide_inverses",
                "wide_transfers",
                "pruned_transfers",
            )
            and levels[name]
        )
        rows.append(
            dict(
                name=f"joint_{index:02}",
                arm=f"C{levels['workers']}",
                polynomials_per_task=levels["packet"],
                features=features,
                collector_overrides=dict(
                    marking=marking,
                    sparse_queue_bytes=32 * 2**20
                    if levels["sparse_queue"]
                    else 0,
                    power_plan_bytes=2 * 2**20 if levels["power_plan"] else 0,
                ),
                nursery_mb=64 if levels["large_nursery"] else 4,
                design_levels=levels,
            )
        )
    return rows


def pair_coverage(rows):
    """Report raw design coverage; runtime activity needs separate evidence."""
    names = tuple(CATEGORICAL) + BINARY
    result = {}
    for left, right in combinations(names, 2):
        levels = {
            name: values if name in CATEGORICAL else (False, True)
            for name, values in (
                (left, CATEGORICAL[left][1] if left in CATEGORICAL else ()),
                (right, CATEGORICAL[right][1] if right in CATEGORICAL else ()),
            )
        }
        observed = {
            (row["design_levels"][left], row["design_levels"][right])
            for row in rows
        }
        possible = set(product(levels[left], levels[right]))
        result[f"{left}:{right}"] = dict(
            observed=len(observed & possible), possible=len(possible)
        )
    return result
