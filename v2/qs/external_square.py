"""Bounded classical MPQS coefficients with separately checked squares."""

import random

from ..common import arithmetic, utils
from ..common.arithmetic import gcd, isqrt, pow
from .polynomial import MAX_COEFFICIENT_BITS, Polynomial, a_target


class CoefficientExhaustedError(Exception):
    """The finite coefficient-candidate allowance found no usable root."""


class CoefficientDivisorError(Exception):
    """Coefficient construction encountered a proper divisor of n."""

    def __init__(self, divisor):
        self.divisor = divisor


def external_square(base, half_width, *, budget, cursor=0, trials=4096):
    """Return a lifted polynomial, next cursor and coefficient certainty.

    Candidates q are 3 mod 4. A probable-prime filter is explicitly labelled;
    the construction depends on verified roots and an inverse, not that label.
    All accepted polynomials satisfy B²-N' = q²*C exactly, even if a composite
    were to pass the filter. No coefficient certainty promotes factors of n.
    """
    utils.require_integer(cursor, "coefficient cursor", 0)
    utils.require_integer(trials, "coefficient trials", 1)
    if trials > 65536 or cursor.bit_length() > MAX_COEFFICIENT_BITS // 2:
        raise ValueError(
            "external coefficient search exceeds its finite limit"
        )
    target = a_target(base.n_prime, half_width)
    start = max(3, isqrt(target), cursor)
    start += (3 - start) % 4
    for candidate in range(start, start + 4 * trials, 4):
        if 2 * candidate.bit_length() > MAX_COEFFICIENT_BITS:
            raise CoefficientExhaustedError("coefficient bit limit")
        budget.consume(
            64 * candidate.bit_length() ** 3 + base.n_prime.bit_length()
        )
        divisor = gcd(candidate, base.n)
        if utils.valid_divisor(divisor, base.n):
            raise CoefficientDivisorError(divisor)
        if gcd(candidate, base.n_prime) != 1:
            continue
        certainty = utils.classify_prime(
            candidate, rng=random.Random(candidate)
        )
        if certainty == utils.Primality.COMPOSITE:
            continue
        root = pow(base.n_prime, (candidate + 1) // 4, candidate)
        if (root * root - base.n_prime) % candidate:
            continue
        if gcd(2 * root, candidate) != 1:
            continue
        # Lift B=root+q*t: 2*root*t == (N'-root²)/q (mod q), making
        # B²-N' divisible by q² without relying on coefficient primality.
        lift = (
            arithmetic.divexact(base.n_prime - root * root, candidate)
            * utils.modular_inverse(2 * root, candidate)
        ) % candidate
        a = candidate * candidate
        b = root + candidate * lift
        polynomial = Polynomial(
            base.n, base.multiplier, a, min(b, a - b), candidate
        )
        return polynomial, candidate + 4, certainty.value

    raise CoefficientExhaustedError("coefficient candidate limit")
