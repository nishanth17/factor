"""Coupled SIQS policies with retained size/resource presets."""

from dataclasses import replace

from ..common import utils
from ..execution.resources import removed_arguments
from .collector_kernels import combine, components
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
    """Return the calibrated preset/portfolio bands, not the policy scope."""
    utils.require_integer(n)
    return backend == "python-int" and (
        10**39 <= abs(n) < 10**40 or 10**49 <= abs(n) < 10**60
    )


def _coupled_policies(config):
    """Apply the selected C10a-01 bundle only to a fresh automatic config.

    Explicit constructors and saved configurations retain their own fields.
    Parity arithmetic and the two-power initializer keep exact fallbacks for
    unsupported discriminants, backends and existing score/kernel choices.
    """
    marking = combine((*components(config.collector.marking), "parity_two"))
    return replace(
        config,
        multiplier=0,
        polynomial_policy="parity2",
        multiplier_policy="parity2_wide",
        assignment_policy="windowed_pair",
        root_update_policy="batch",
        collector=replace(config.collector, marking=marking),
    )


@removed_arguments(("memory_bytes",))
def legacy_siqs_config(n, *, backend="python-int", recipe="calibrated"):
    """Select coupled policies at every size without enlarging grants.

    The existing 40 and 50–60-digit parameter/resource presets remain intact.
    Other sizes/backends use the reference parameters with the same coupled
    policies; those parameters are not calibrated larger-size presets.
    The explicit recipe preserves historical numerical/cache choices
    independently of root resource limits.
    """
    if not uses_siqs_defaults(n, backend=backend):
        return _coupled_policies(
            SIQSConfig(
                backend=backend,
                sign_policy="swap_high",
            )
        )

    forty_digits = abs(n) < 10**40
    if recipe not in ("calibrated", "reference"):
        raise ValueError("unknown SIQS recipe")
    calibrated = recipe == "calibrated"
    # A1a passed at 50; the user extends routing through 55 despite drift at
    # 55 and unmeasured intermediate sizes. Keep the tested owned allowance.
    selected_marking = (
        "density+presieve+scan_fused"
        if 10**49 <= abs(n) < 10**55 and calibrated
        else "sparse"
    )
    # C10 isolated this plan at 40/50 digits; the later 51-digit default study
    # confirms the same bundle. Other owned allowances were not tested.
    plan_bytes = C10_POWER_PLAN_BYTES if calibrated else 0

    larger = abs(n) >= 10**51
    config_type = CountedSIQSConfig if larger else SIQSConfig
    count_policy = (
        {"base_count": 2500 if abs(n) < 10**55 else 5000} if larger else {}
    )
    factor_count = 5 if forty_digits else 6
    if larger:
        factor_count = 7 if abs(n) < 10**55 else 8
    config = config_type(
        **count_policy,
        sign_policy="swap_high",
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
            "collection" if not forty_digits and calibrated else "window"
        ),
        collector=GrowingDoubleLargeSieveConfig(
            sparse_queue_bytes=32 * 1024 * 1024,
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
    return _coupled_policies(config)


@removed_arguments(("memory_bytes",))
def default_siqs_config(n, *, backend="python-int", recipe="calibrated"):
    """Select Y0 for native 50+ digit QS; retain explicit reference recipes.

    Between defined numerical geometries, inherit the nearest row, preferring
    the smaller size on a tie. Above 70 is extrapolated; existing 75/80
    geometry stays intact. Smaller sizes and other backends retain their
    existing numerical policy.
    """
    utils.require_integer(n)
    if recipe not in ("calibrated", "reference"):
        raise ValueError("unknown SIQS recipe")
    if backend != "python-int" or abs(n) < 10**49 or recipe == "reference":
        return legacy_siqs_config(n, backend=backend, recipe=recipe)

    from .production_policy import y0_config

    return y0_config(n)


def default_siqs_budget(n, *, backend="python-int"):
    """Create the default root owner independently of numerical presets."""
    from ..execution.budget import Budget

    return Budget()
