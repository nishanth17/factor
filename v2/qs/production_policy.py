"""Selected A7a Y0 parameters and digit-based process recipes."""

from dataclasses import dataclass

from ..common import utils


@dataclass(frozen=True)
class QSRecipe:
    """Nearest measured worker row, with explicit extrapolation metadata."""

    digits: int
    anchor: int
    workers: int
    packet: int

    @property
    def extrapolated(self):
        return self.digits > 70


ROWS = (
    (40, 4, 4),
    (50, 4, 8),
    (55, 4, 16),
    (60, 8, 16),
    (65, 14, 64),
    (70, 8, 64),
)


def qs_recipe(n):
    """Route by decimal input size; ties select the smaller anchor."""
    utils.require_integer(n)
    digits = len(str(abs(n)))
    if digits < 40:
        return QSRecipe(digits, 0, 0, 0)
    anchor, workers, packet = min(ROWS, key=lambda row: abs(row[0] - digits))
    return QSRecipe(digits, anchor, workers, packet)


def default_qs_workers(n, config):
    """Resolve the frozen process bundle without starting a pool.

    Explicit unsupported engines keep their serial implementation. Pools are
    constructed only when the portfolio reaches an eligible SIQS stage.
    """
    recipe = qs_recipe(n)
    if (
        not recipe.workers
        or config is None
        or (
            config.mode != "siqs"
            or config.backend != "python-int"
            or config.growth_steps
            or not config.shared_relations
            or not config.diverse
            or getattr(config.collector, "recovery_workers", 0)
        )
    ):
        return None

    from .current_parallel import QSWorkerConfig

    return QSWorkerConfig(
        workers=recipe.workers,
        polynomials_per_task=recipe.packet,
        max_batch_atoms=2048,
        lease_budget=True,
        reuse_state=True,
        state_cache_bytes=64 * 2**20,
        packed_transport=True,
        verification="worker",
        lookahead=4,
        private_stop=True,
        stop_backend="semaphore",
        coefficient_cache=True,
    )


def y0_config(n):
    """Reproduce the H3-adapted frozen Y0 numerical anchors exactly."""
    from .extended_policy import BoundedYieldSieveConfig
    from .siqs import CountedSIQSConfig

    digits = len(str(abs(n)))
    anchor = min((50, 55, 60, 65, 70), key=lambda size: abs(size - digits))
    base_count, factors = {
        50: (1800, 6),
        55: (2500, 7),
        60: (5000, 8),
        65: (8500, 9),
        70: (10500, 10),
    }[anchor]
    large = anchor == 70
    collector = BoundedYieldSieveConfig(
        block_width=4096,
        marking="presieve+unroll4+parity_two",
        division="bucket",
        score_policy="powers",
        sparse_queue_bytes=32 * 2**20,
        power_plan_bytes=2 * 2**20,
        residual_bound=3000000,
        large_prime_bound=3819584 if large else 3000000,
        large_product_bound=108837419273339392 if large else 8950715676137152,
        candidate_bound=115200000000,
        small_prime_bound=954896 if large else 415256,
        target_rows=32000 if large else None,
        max_partials=98304 if large else None,
        split_call_limit=1048576 if large else None,
        split_policy="rho-first-v3",
        batch_policy="off",
        threshold_epochs=0,
        recompute_prime_limit=0,
        solver_panel_bits=0,
    )
    return CountedSIQSConfig(
        base_count=base_count,
        factor_count=factors,
        multiplier=0,
        half_width=65536,
        family_count=100000,
        pool_size=64,
        assignment_policy="windowed_pair",
        polynomials_per_family=64,
        max_half_width=65536,
        max_stalled=100000,
        max_trivial=100000,
        row_excess=32,
        batch_width=4096,
        filter_row_growth=32,
        polynomial_policy="parity2",
        multiplier_policy="parity2_wide",
        root_update_policy="batch",
        sign_policy="swap_high",
        solve_policy="collection",
        collector=collector,
    )
