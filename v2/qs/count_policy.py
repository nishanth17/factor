"""Charged exact eligible-prime count policies for SIQS setup."""

from math import gcd

from ..common import prime_sieve, utils
from ..execution.resources import removed_arguments
from .factor_base import checked_target

ENDPOINT_CEILING = 250000


class CountPolicyRefusedError(ValueError):
    """The requested base count cannot fit the finite endpoint ceiling."""


@removed_arguments(("memory_bytes",))
def endpoint_for_count(n, multiplier, count, budget, *, ceiling):
    """Resolve the smallest exclusive bound, or return a proper divisor.

    Include 2 and singular multiplier primes as build_factor_base does. The
    temporary sieve is released before root/base allocation. Resume repeats
    this deterministic setup under the remaining budget rather than refunding
    its earlier work or retaining an uncharged partial sieve.
    """
    target = checked_target(n, multiplier)
    utils.require_integer(count, "base count", 1)
    utils.require_integer(ceiling, "endpoint ceiling", 3)
    if ceiling > ENDPOINT_CEILING:
        raise ValueError("count policy endpoint exceeds its finite ceiling")
    8192 + 64 * ceiling + 8 * target.bit_length()
    budget.consume(n.bit_length())
    divisor = gcd(n, multiplier)
    if utils.valid_divisor(divisor, n):
        return None, divisor
    if divisor != 1:
        raise ValueError("multiplier is divisible by n")

    budget.consume(ceiling)
    found = 0
    for prime in prime_sieve.prime_sieve(ceiling):
        budget.consume(target.bit_length())
        divisor = gcd(prime, n)
        if utils.valid_divisor(divisor, n):
            return None, divisor
        residue = target % prime
        if (
            prime == 2
            or residue == 0
            or pow(residue, (prime - 1) // 2, prime) == 1
        ):
            found += 1
            if found == count:
                return prime + 1, None

    raise CountPolicyRefusedError("base count is unreachable under ceiling")
