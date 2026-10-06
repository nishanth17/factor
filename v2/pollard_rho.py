"""Bounded Brent rho with independent retries and local GCD recovery."""

from dataclasses import dataclass

from . import arithmetic, constants, utils
from .arithmetic import gcd


@dataclass
class RhoStats:
    """Work consumed across attempts, including saturation recovery."""

    attempts: int = 0
    evaluations: int = 0
    gcd_calls: int = 0
    saturated_batches: int = 0


def _brent_attempt(
    n, start, offset, batch_size, max_evaluations, recovery_limit, stats
):
    y = start % n
    cycle_length = 1
    used = 0

    while used < max_evaluations:
        # Brent doubles the cycle length; x anchors the next comparison run.
        x = y
        advance = min(cycle_length, max_evaluations - used)
        for _ in range(advance):
            y = (y * y + offset) % n
        used += advance
        stats.evaluations += advance
        position = 0

        while position < cycle_length and used < max_evaluations:
            # Save the batch start so saturation can replay the same walk.
            saved = y
            count = min(
                batch_size, cycle_length - position, max_evaluations - used
            )
            product = 1
            for _ in range(count):
                y = (y * y + offset) % n
                product = product * abs(x - y) % n
            # No early exit occurs inside a batch: account once per chunk.
            used += count
            stats.evaluations += count
            divisor = gcd(product, n)
            stats.gcd_calls += 1
            if 1 < divisor < n:
                return int(divisor)
            if divisor == n:
                stats.saturated_batches += 1
                # Recovery evaluations also consume this attempt's allowance.
                recover = min(count, recovery_limit, max_evaluations - used)

                for _ in range(recover):
                    saved = (saved * saved + offset) % n
                    used += 1
                    stats.evaluations += 1
                    divisor = gcd(abs(x - saved), n)
                    stats.gcd_calls += 1
                    if 1 < divisor < n:
                        return int(divisor)

                return None

            position += count

        cycle_length <<= 1

    return None


def factorize_rho(
    n,
    verbose=False,
    *,
    seed=None,
    rng=None,
    max_attempts=constants.RHO_ATTEMPTS,
    max_evaluations=constants.RHO_EVALUATIONS,
    batch_size=constants.RHO_BATCH_SIZE,
    recovery_limit=constants.RHO_RECOVERY_LIMIT,
    stats=None,
    _known_composite=False,
    backend=None,
):
    """Return a proper divisor or None within explicit evaluation limits.

    The dispatcher supplies _known_composite only after its own classification;
    ordinary public calls still check primality. All outputs remain divisors.
    """
    utils.require_integer(n, minimum=1)
    engine = (
        arithmetic.backend_for(n)
        if backend is None
        else arithmetic.get_backend(backend)
    )
    n = engine.integer(n)
    utils.require_integer(max_attempts, "max_attempts", 0)
    utils.require_integer(max_evaluations, "max_evaluations", 0)
    utils.require_integer(batch_size, "batch_size", 1)
    utils.require_integer(recovery_limit, "recovery_limit", 0)
    generator = utils.resolve_rng(seed, rng)
    if not max_attempts or not max_evaluations or n == 1:
        return None
    if n % 2 == 0:
        return 2 if n > 2 else None
    if not _known_composite and utils.is_prime(n, rng=generator):
        return None
    work = stats if stats is not None else RhoStats()

    for _ in range(max_attempts):
        start = engine.integer(generator.randint(1, int(n) - 1))
        offset = engine.integer(generator.randint(1, int(n) - 1))
        work.attempts += 1
        if verbose:
            print(f"Rho attempt {work.attempts}, offset={offset}")
        divisor = _brent_attempt(
            n, start, offset, batch_size, max_evaluations, recovery_limit, work
        )
        if utils.valid_divisor(divisor, n):
            return int(divisor)

    return None
