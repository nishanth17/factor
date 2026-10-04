"""Two-stage Pollard p-1 with complete relation batches and bounded retries."""

from math import gcd

from . import constants, prime_sieve, utils


def compute_bounds(n):
    """Return conservative explicit defaults, independent of float logs."""
    utils.require_integer(n, minimum=1)
    return constants.PM1_B1, constants.PM1_B2


def _stage_one(n, base, primes, b1):
    residue = base % n
    for prime in primes:
        previous = residue
        residue = pow(residue, utils.prime_power(prime, b1), n)
        divisor = gcd(residue - 1, n)
        if divisor == n:
            # Replay only this prime's exponent units to recover a lost split.
            power = prime
            while power <= b1:
                previous = pow(previous, prime, n)
                found = gcd(previous - 1, n)
                if utils.valid_divisor(found, n):
                    return residue, found
                if found == n:
                    break
                power *= prime
            return residue, n
        if utils.valid_divisor(divisor, n):
            return residue, divisor
    return residue, 1


def _stage_two_terms(residue, n, primes):
    """Yield every a**q - 1, including the first stage-2 prime."""
    previous_prime = 0
    value = 1
    gap_cache = {}
    for prime in primes:
        gap = prime - previous_prime
        if gap not in gap_cache:
            gap_cache[gap] = pow(residue, gap, n)
        # a**q = a**previous_prime * a**gap; cache only this residue/modulus.
        value = value * gap_cache[gap] % n
        yield (value - 1) % n
        previous_prime = prime


def factorize_pm1(
    n,
    verbose=False,
    *,
    b1=None,
    b2=None,
    base=2,
    max_attempts=constants.PM1_ATTEMPTS,
    batch_size=constants.GCD_BATCH_SIZE,
):
    """Return a proper divisor or None; bounds include B1 and B2.

    Stage 2 targets a divisor p whose p-1 has one eligible large prime,
    not a prime divisor of n between the bounds.
    """
    utils.require_integer(n, minimum=1)
    default_b1, default_b2 = compute_bounds(n)
    b1 = default_b1 if b1 is None else b1
    b2 = default_b2 if b2 is None else b2
    utils.require_integer(b1, "b1", 2)
    utils.require_integer(b2, "b2", b1)
    utils.require_integer(base, "base", 2)
    utils.require_integer(max_attempts, "max_attempts", 0)
    utils.require_integer(batch_size, "batch_size", 1)

    if n == 1 or not max_attempts:
        return None
    if n % 2 == 0:
        return 2 if n > 2 else None
    if utils.is_prime(n):
        return None

    stage_one_primes = prime_sieve.prime_sieve(b1 + 1)
    stage_two_primes = None

    for attempt in range(max_attempts):
        candidate = base + attempt
        divisor = gcd(candidate, n)
        if utils.valid_divisor(divisor, n):
            return divisor
        if divisor == n:
            continue
        if verbose:
            print(f"p-1 base={candidate}, bounds={b1}/{b2}")

        residue, divisor = _stage_one(n, candidate, stage_one_primes, b1)
        if utils.valid_divisor(divisor, n):
            return divisor
        if divisor == n:
            continue
        if stage_two_primes is None:
            stage_two_primes = prime_sieve.segmented_sieve(b1 + 1, b2 + 1)

        terms = []
        saturated = False
        for term in _stage_two_terms(residue, n, stage_two_primes):
            terms.append(term)
            if len(terms) == batch_size:
                divisor, saturated = utils.batch_factor(terms, n)
                if divisor is not None:
                    return divisor
                terms.clear()
                if saturated:
                    break

        if not saturated and terms:
            divisor, _ = utils.batch_factor(terms, n)
            if divisor is not None:
                return divisor

    return None
