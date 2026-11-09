"""Literal configurations: resource tiers never reroute collector policies."""

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.benchmarks.infrastructure.resource_control import (
    historical_replace as replace,
)

from ....qs.sieve_collector import GrowingDoubleLargeSieveConfig
from ....qs.siqs import CountedSIQSConfig

# Above 60 these are declared exploratory configurations, not tuned presets.
COUNTS = {60: 5000, 62: 6000, 65: 8500, 70: 10500}
FACTORS = {60: 8, 62: 9, 65: 9, 70: 10}
TIERS = {60: 1, 62: 2, 65: 2, 70: 3}
OWNED = {1: 2**30, 2: 2**31, 3: 2**32}
RSS = {1: 3 * 2**30, 2: 5 * 2**30, 3: 8 * 2**30}
ROWS = {60: 12000, 62: 16000, 65: 24000, 70: 32000}
ATOMS = {60: 65536, 62: 98304, 65: 131072, 70: 196608}
PARTIALS = {60: 32768, 62: 49152, 65: 65536, 70: 98304}


def outside_job_reserve(digits):
    """Keep assignment ownership and bounded harness state inside the tier."""
    return 16 * 2**20 + 65536 + 512 * COUNTS[digits]


def configuration(digits, arm, *, multiplier=0):
    """Only the selected representation and ranking differ within a pair."""
    if digits not in COUNTS or arm not in (
        "control",
        "parity2",
        "representation",
    ):
        raise ValueError("unfrozen size or arm")
    config = adapt_resource_configuration(
        CountedSIQSConfig,
        base_count=COUNTS[digits],
        base_bound=250000,
        multiplier=multiplier,
        half_width=65536,
        factor_count=FACTORS[digits],
        family_count=100000,
        pool_size=64,
        assignment_policy="flyer",
        polynomials_per_family=64,
        max_half_width=65536,
        max_stalled=100000,
        max_trivial=100000,
        row_excess=32,
        batch_width=4096,
        filter_row_growth=32,
        sign_policy="swap_high",
        solve_policy="collection",
        memory_bytes=OWNED[TIERS[digits]] - outside_job_reserve(digits),
        checkpoint_bytes=1048576,
        collector=adapt_resource_configuration(
            GrowingDoubleLargeSieveConfig,
            sparse_queue_bytes=32 * 2**20,
            marking="sparse",
            block_width=4096,
            division="bucket",
            score_policy="powers",
            power_plan_bytes=2 * 2**20,
            residual_bound=3000000,
            large_prime_bound=3000000,
            large_product_bound=115200000000,
            max_atoms=ATOMS[digits],
            max_partials=PARTIALS[digits],
            max_relations=ROWS[digits],
            split_call_limit=1048576,
        ),
    )
    if arm != "control":
        config = replace(
            config,
            polynomial_policy="parity2",
            multiplier_policy="standard"
            if arm == "representation"
            else "parity2",
        )
    return config
