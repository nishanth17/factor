"""Exact integer helpers and explicit primality classifications.

Fixed-base guarantees have strict upper bounds; see the source comparison
in benchmarks/README.md. The wider ranges are exhaustive computational
results of Sorenson and Webster, https://arxiv.org/abs/1509.00864.
Survivors outside supported ranges remain probable primes.
"""

import random
import sys
from bisect import bisect_left, bisect_right
from enum import Enum

from . import arithmetic, constants
from .arithmetic import gcd, isqrt, pow

SMALL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
WORD_DETERMINISTIC_LIMIT = 2**64
TWELVE_BASE_LIMIT = 318665857834031151167461
DETERMINISTIC_LIMIT = 3317044064679887385961981
DETERMINISTIC_BASES = (2, 325, 9375, 28178, 450775, 9780504, 1795265022)
THIRTEEN_PRIME_BASES = SMALL_PRIMES + (41,)
PRIMALITY_POLICY = "mr13-strict-v1"
LEGACY_PRIMALITY_POLICY = "mr64-strict-v1"
_USE_PYPY_INVERSE = hasattr(sys, "pypy_version_info")
_USE_PYPY_SEARCH = _USE_PYPY_INVERSE


class Primality(str, Enum):
    """Evidence supported by a primality test, not a truthy numeric residue."""

    COMPOSITE = "composite"
    PROBABLE = "probable_prime"
    PROVEN = "proven_prime"


def require_integer(value, name="n", minimum=None):
    """Reject floats/booleans; optionally enforce an integer lower bound."""
    if not arithmetic.is_integer(value):
        raise TypeError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return value


def _binary_search_bisect(value, array, include_equal=False):
    """Use the native bisect backend for generic sorted sequences."""
    return (bisect_left if include_equal else bisect_right)(array, value)


def _binary_search_jit(value, array, include_equal=False):
    """JIT-friendly boundaries with duplicate-safe equality shortcuts."""
    size = len(array)
    if size == 0 or value < array[0]:
        return 0
    if array[-1] < value:
        return size

    low, high = 0, size - 1

    while low <= high:
        middle = (low + high) >> 1
        item = array[middle]
        if item == value:
            # Equality alone cannot choose a boundary when duplicates exist.
            if include_equal:
                if middle == 0 or array[middle - 1] < value:
                    return middle
                high = middle - 1
            else:
                if middle + 1 == size or value < array[middle + 1]:
                    return middle + 1
                low = middle + 1
        elif item < value:
            low = middle + 1
        else:
            high = middle - 1

    return low


# Select once: a runtime branch inside every search impedes PyPy inlining.
binary_search = (
    _binary_search_jit if _USE_PYPY_SEARCH else _binary_search_bisect
)


def extended_gcd(a, b):
    """Return (g, x, y) with g >= 0 and a*x + b*y == gcd(a, b)."""
    require_integer(a, "a")
    require_integer(b, "b")
    if arithmetic.is_mpz(a) or arithmetic.is_mpz(b):
        return arithmetic.gcdext(a, b)

    old_r, remainder = a, b
    old_x, x = 1, 0
    old_y, y = 0, 1
    while remainder:
        quotient = old_r // remainder
        old_r, remainder = remainder, old_r - quotient * remainder
        old_x, x = x, old_x - quotient * x
        old_y, y = y, old_y - quotient * y

    if old_r < 0:
        return -old_r, -old_x, -old_y
    return old_r, old_x, old_y


def xgcd(a, b):
    """Return the coefficient of a in Bezout's identity (legacy name)."""
    return extended_gcd(a, b)[1]


def modular_inverse(value, modulus):
    """Return the canonical inverse; raise ValueError for a nonunit.

    Factoring callers must GCD-check a denominator before calling this.
    """
    require_integer(value, "value")
    require_integer(modulus, "modulus", 2)
    if arithmetic.is_mpz(value) or arithmetic.is_mpz(modulus):
        return arithmetic.invert(value, modulus)

    if not _USE_PYPY_INVERSE:
        return pow(value, -1, modulus)
    # Track only the coefficient of value. divmod shares quotient/remainder
    # work; computing both Bezout coefficients is unnecessary for an inverse.
    a, b = modulus, value % modulus
    coefficient, next_coefficient = 0, 1

    while b:
        quotient, remainder = divmod(a, b)
        a, b = b, remainder
        coefficient, next_coefficient = (
            next_coefficient,
            coefficient - quotient * next_coefficient,
        )

    if a != 1:
        raise arithmetic.NonInvertibleError(a)
    return coefficient % modulus


def prime_power(prime, bound):
    """Return the largest power of prime <= bound using exact integers."""
    require_integer(prime, "prime", 2)
    require_integer(bound, "bound", prime)
    power = prime
    while power <= bound // prime:
        power *= prime
    return power


def valid_divisor(divisor, n):
    """Reject failure sentinels, endpoints, booleans, and nondivisors."""
    return (
        arithmetic.is_integer(divisor)
        and not isinstance(divisor, bool)
        and 1 < divisor < n
        and n % divisor == 0
    )


def batch_factor(terms, n):
    """Return (proper_factor_or_None, saturated) for a finite term batch.

    A saturated product is replayed term by term to recover mixed factors.
    If no term splits n, the caller must retry rather than claim success.
    """
    product = 1
    for term in terms:
        product = product * term % n
    divisor = gcd(product, n)
    if valid_divisor(divisor, n):
        return divisor, False
    if divisor == n:
        for term in terms:
            divisor = gcd(term, n)
            if valid_divisor(divisor, n):
                return divisor, True
        return None, True
    return None, False


def resolve_rng(seed=None, rng=None):
    """Use an injected generator or create a private seeded generator."""
    if seed is not None and rng is not None:
        raise ValueError("specify seed or rng, not both")
    return rng if rng is not None else random.Random(seed)


def is_prime_bf(n):
    """Exact trial-division primality; suitable only for small inputs."""
    require_integer(n)
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    for divisor in range(5, isqrt(n) + 1, 6):
        if n % divisor == 0 or n % (divisor + 2) == 0:
            return False
    return True


def deterministic_bases(n):
    """Return the supported fixed witnesses, or None outside their domain.

    Callers handle small divisors before running these witnesses.
    Equality belongs to the next range: the two wider endpoints are known
    strong pseudoprimes to the preceding set, not certifiable primes.
    """
    require_integer(n)
    if n < 2:
        return None
    if n < 9_080_191:
        return (31, 73)
    if n < 4_759_123_141:
        return (2, 7, 61)
    if n < WORD_DETERMINISTIC_LIMIT:
        return DETERMINISTIC_BASES
    if n < TWELVE_BASE_LIMIT:
        return SMALL_PRIMES
    if n < DETERMINISTIC_LIMIT:
        return THIRTEEN_PRIME_BASES
    return None


def _strong_probable_prime(n, base, odd_part, shifts):
    base %= n
    if base in (0, 1):
        return True
    value = pow(base, odd_part, n)
    if value in (1, n - 1):
        return True
    for _ in range(shifts - 1):
        value = value * value % n
        if value == n - 1:
            return True
    return False


def classify_prime(
    n,
    *,
    use_probabilistic=False,
    tolerance=constants.PRIMALITY_ROUNDS,
    rng=None,
):
    """Classify n; random rounds apply to every survivor in requested mode.

    Small-prime membership/divisibility proves trivial cases without rounds.
    Otherwise probabilistic mode uses exactly tolerance independently drawn
    witnesses for a survivor, even below the deterministic limit. The default
    uses the bounded deterministic test, then random rounds above its domain.
    A compositeness witness may stop either test early.
    """
    require_integer(n)
    require_integer(tolerance, "tolerance", 1)
    if n < 2:
        return Primality.COMPOSITE
    for prime in SMALL_PRIMES:
        if n == prime:
            return Primality.PROVEN
        if n % prime == 0:
            return Primality.COMPOSITE

    if not use_probabilistic and n < 41 * 41:
        return Primality.PROVEN
    odd_part = n - 1
    shifts = (odd_part & -odd_part).bit_length() - 1
    odd_part >>= shifts
    bases = None if use_probabilistic else deterministic_bases(n)
    deterministic = bases is not None
    if not deterministic:
        generator = rng if rng is not None else random.SystemRandom()
        bases = (generator.randint(2, n - 2) for _ in range(tolerance))

    for base in bases:
        if not _strong_probable_prime(n, base, odd_part, shifts):
            return Primality.COMPOSITE
    return Primality.PROVEN if deterministic else Primality.PROBABLE


def is_prime_fast(
    n,
    use_probabilistic=False,
    tolerance=constants.PRIMALITY_ROUNDS,
    *,
    rng=None,
):
    """Boolean convenience test; use classify_prime for certainty."""
    return (
        classify_prime(
            n,
            use_probabilistic=use_probabilistic,
            tolerance=tolerance,
            rng=rng,
        )
        is not Primality.COMPOSITE
    )


def is_prime(
    n,
    use_probabilistic=False,
    tolerance=constants.PRIMALITY_ROUNDS,
    *,
    rng=None,
):
    """Return a bool; classify_prime distinguishes proof from probability."""
    return is_prime_fast(n, use_probabilistic, tolerance, rng=rng)
