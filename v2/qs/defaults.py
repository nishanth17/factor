"""Size-scoped SIQS defaults; explicit configurations win."""

from ..common import utils
from .sieve_collector import GrowingDoubleLargeSieveConfig
from .siqs import CountedSIQSConfig, SIQSConfig

C9_WORK_LIMIT = 10**13
C9_SECONDS = 1200
C9_MEMORY_BYTES = 1024 * 1024 * 1024
C9_PORTFOLIO_MEMORY_BYTES = C9_MEMORY_BYTES + 16 * 1024 * 1024
C10_POWER_PLAN_BYTES = 2 * 1024 * 1024


def uses_c9_defaults(n, *, backend="python-int"):
    """Return the observable size/backend scope of fresh C9 default routing."""
    utils.require_integer(n)
    return backend == "python-int" and 10**49 <= abs(n) < 10**51


def uses_siqs_defaults(n, *, backend="python-int"):
    """Return native bands selected for fresh automatic jobs."""
    utils.require_integer(n)
    return backend == "python-int" and (
        10**39 <= abs(n) < 10**40 or 10**49 <= abs(n) < 10**60
    )


def default_siqs_config(n, *, backend="python-int", memory_bytes=None):
    """Select 40 and 50–60-digit bundles without changing explicit quotas.

    The user-directed 52–60 extension uses compact count policies; only 55
    training and preliminary fresh 60-digit evidence cover that extension.
    Other sizes/backends retain the reference configuration. Explicit memory
    remains binding, including the temporary charged count-resolution sieve.
    """
    if not uses_siqs_defaults(n, backend=backend):
        return SIQSConfig(backend=backend, memory_bytes=memory_bytes)

    forty_digits = abs(n) < 10**40
    owned_memory = C9_MEMORY_BYTES if memory_bytes is None else memory_bytes
    # A1a passed at 50; the user extends routing through 55 despite drift at
    # 55 and unmeasured intermediate sizes. Keep the tested owned allowance.
    selected_marking = (
        "density+presieve+scan_fused"
        if 10**49 <= abs(n) < 10**55 and owned_memory == C9_MEMORY_BYTES
        else "sparse"
    )
    # C10 isolated this plan at 40/50 digits; the later 51-digit default study
    # confirms the same bundle. Other owned allowances were not tested.
    plan_bytes = C10_POWER_PLAN_BYTES if owned_memory == C9_MEMORY_BYTES else 0

    larger = abs(n) >= 10**51
    config_type = CountedSIQSConfig if larger else SIQSConfig
    count_policy = (
        {"base_count": 2500 if abs(n) < 10**55 else 5000} if larger else {}
    )
    factor_count = 5 if forty_digits else 6
    if larger:
        factor_count = 7 if abs(n) < 10**55 else 8
    return config_type(
        **count_policy,
        backend=backend,
        base_bound=250000 if larger else 5000 if forty_digits else 30000,
        multiplier=0 if larger else 1,
        half_width=65536,
        factor_count=factor_count,
        family_count=100000,
        pool_size=64 if larger else 32,
        assignment_policy="flyer",
        polynomials_per_family=64,
        max_half_width=65536,
        max_stalled=100000,
        max_trivial=100000,
        row_excess=32,
        batch_width=4096,
        filter_row_growth=32,
        # C10a-05 passed fresh complete-cost gates at 50 digits. The later
        # 51-digit default study confirms this bundle's completion gain; 40
        # digits retains window pacing after its inconclusive comparison.
        solve_policy=(
            "collection"
            if not forty_digits and owned_memory == C9_MEMORY_BYTES
            else "window"
        ),
        memory_bytes=owned_memory,
        checkpoint_bytes=1048576,
        collector=GrowingDoubleLargeSieveConfig(
            marking=selected_marking,
            block_width=4096,
            division="bucket",
            score_policy="powers",
            power_plan_bytes=plan_bytes,
            residual_bound=25000000 if forty_digits else 3000000,
            large_prime_bound=500000 if forty_digits else 3000000,
            large_product_bound=3200000000 if forty_digits else 115200000000,
        ),
    )


def default_siqs_budget(n, *, backend="python-int"):
    """Use finite size-scoped grants only when the caller omits a budget."""
    from ..execution.budget import Budget

    if uses_siqs_defaults(n, backend=backend):
        return Budget(
            work_limit=C9_WORK_LIMIT,
            seconds=1800 if abs(n) >= 10**55 else C9_SECONDS,
            cpu_seconds=1800 if abs(n) >= 10**55 else C9_SECONDS,
        )
    return Budget()
