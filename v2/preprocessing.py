"""Exact integer preprocessing; no arbitrary-size input enters a float."""

from math import isqrt

from . import utils


def integer_root(n, exponent):
    """Return floor(n**(1/exponent)) by integer Newton iteration."""
    utils.require_integer(n, minimum=0)
    utils.require_integer(exponent, "exponent", 1)
    if n < 2 or exponent == 1:
        return n
    if exponent >= n.bit_length():
        return 1
    estimate = 1 << ((n.bit_length() + exponent - 1) // exponent)
    while True:
        following = (
            (exponent - 1) * estimate + n // estimate ** (exponent - 1)
        ) // exponent
        if following >= estimate:
            return estimate
        estimate = following


def strip_twos(n):
    """Return (odd part, valuation at 2) for a positive integer."""
    utils.require_integer(n, minimum=1)
    exponent = (n & -n).bit_length() - 1
    return n >> exponent, exponent


def fermat_step(n, candidate):
    """Test one a in a²-n=b²; return a proper factor or None."""
    difference = candidate * candidate - n
    if difference < 0:
        raise ValueError("Fermat candidate is below ceil(sqrt(n))")
    root = isqrt(difference)
    divisor = candidate - root
    if root * root == difference and utils.valid_divisor(divisor, n):
        return divisor
    return None
