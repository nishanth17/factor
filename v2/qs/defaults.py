"""Confirmed C9 size-scoped SIQS defaults; explicit configurations win."""

from ..common import utils
from .sieve_collector import GrowingDoubleLargeSieveConfig
from .siqs import SIQSConfig

C9_WORK_LIMIT = 10**13
C9_SECONDS = 1200
C9_MEMORY_BYTES = 1024 * 1024 * 1024
C9_PORTFOLIO_MEMORY_BYTES = C9_MEMORY_BYTES + 16 * 1024 * 1024


def uses_c9_defaults(n, *, backend="python-int"):
    """Return the observable size/backend scope of fresh C9 default routing."""
    utils.require_integer(n)
    return backend == "python-int" and 10**49 <= abs(n) < 10**50


def default_siqs_config(n, *, backend="python-int", memory_bytes=None):
    """Select the accepted 50-digit bundle without changing explicit quotas.

    Other sizes/backends retain the reference configuration. None uses the
    confirmed 1 GiB owned allowance at 50 digits and the usual legacy default
    elsewhere; pass an integer to keep a caller's finite memory allowance.
    """
    if not uses_c9_defaults(n, backend=backend):
        return SIQSConfig(backend=backend, memory_bytes=memory_bytes)

    return SIQSConfig(
        backend=backend,
        base_bound=30000,
        half_width=65536,
        factor_count=6,
        family_count=100000,
        pool_size=32,
        assignment_policy="flyer",
        polynomials_per_family=64,
        max_half_width=65536,
        max_stalled=100000,
        max_trivial=100000,
        row_excess=32,
        batch_width=4096,
        filter_row_growth=32,
        memory_bytes=C9_MEMORY_BYTES if memory_bytes is None else memory_bytes,
        checkpoint_bytes=1048576,
        collector=GrowingDoubleLargeSieveConfig(
            block_width=4096,
            division="bucket",
            score_policy="powers",
            residual_bound=3000000,
            large_prime_bound=3000000,
            large_product_bound=115200000000,
        ),
    )


def default_siqs_budget(n, *, backend="python-int"):
    """Use the confirmed finite grant only when the caller omits a budget."""
    from ..execution.budget import Budget

    if uses_c9_defaults(n, backend=backend):
        return Budget(
            work_limit=C9_WORK_LIMIT,
            seconds=C9_SECONDS,
            cpu_seconds=C9_SECONDS,
        )
    return Budget()
