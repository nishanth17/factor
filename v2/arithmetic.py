"""Exact optional GMP arithmetic with coarse, instance-owned boundaries.

Algorithms select a representation at entry and retain it through operator
loops. Small indices and schedules stay native. No global backend switch or
per-multiply dispatch is needed; checkpoints and public results use integers.
"""

import builtins
import importlib
import math
import sys
from dataclasses import dataclass

_gmp = None


class BackendUnavailableError(RuntimeError):
    """The selected optional backend cannot run on this interpreter."""


def _load_gmp():
    global _gmp
    if _gmp is None:
        if not hasattr(sys, "pypy_version_info") or sys.version_info[:2] != (
            3,
            11,
        ):
            raise BackendUnavailableError(
                "GMP backend requires PyPy Python 3.11"
            )
        try:
            _gmp = importlib.import_module("gmpy2")
        except ImportError as error:
            raise BackendUnavailableError(
                "gmpy2 is unavailable on this PyPy; choose python-int "
                "or install gmpy2 for the same interpreter"
            ) from error
    return _gmp


def is_mpz(value):
    """Recognize immutable mpz only, including caller-imported instances."""
    if isinstance(value, int):
        return False
    if _gmp is not None:
        return isinstance(value, _gmp.mpz)
    if type(value).__module__ == "gmpy2" and type(value).__name__ == "mpz":
        return isinstance(value, _load_gmp().mpz)
    return False


def is_integer(value):
    return not isinstance(value, bool) and (
        isinstance(value, int) or is_mpz(value)
    )


def _require(value, name="value", minimum=None):
    if not is_integer(value):
        raise TypeError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")


@dataclass(frozen=True)
class ArithmeticBackend:
    """Stable identity and an explicit conversion at an engine boundary."""

    name: str

    def integer(self, value):
        _require(value)
        return (
            int(value) if self.name == "python-int" else _load_gmp().mpz(value)
        )

    @property
    def identity(self):
        if self.name == "python-int":
            return self.name
        module = _load_gmp()
        return f"gmpy2-mpz/{module.version()}/{module.mp_version()}"


def get_backend(name="python-int"):
    """Resolve explicit names; a missing GMP dependency never falls back."""
    if name not in ("python-int", "gmpy2-mpz"):
        raise ValueError("backend must be python-int or gmpy2-mpz")
    if name == "gmpy2-mpz":
        _load_gmp()
    return ArithmeticBackend(name)


def backend_for(value):
    return get_backend("gmpy2-mpz" if is_mpz(value) else "python-int")


def canonical(value):
    """Convert mpz recursively at JSON/public boundaries, preserving shape."""
    if is_mpz(value):
        return int(value)
    if isinstance(value, dict):
        return {canonical(key): canonical(item) for key, item in value.items()}
    if isinstance(value, list):
        return [canonical(item) for item in value]
    if isinstance(value, tuple):
        return tuple(canonical(item) for item in value)
    return value


def json_integer(value):
    """JSON default that encodes only exact mpz, rejecting other objects."""
    if is_mpz(value):
        return int(value)
    raise TypeError(f"{type(value).__name__} is not JSON serializable")


def gcd(*values):
    if any(is_mpz(value) for value in values):
        return _load_gmp().gcd(*values)
    return math.gcd(*values)


def pow(base, exponent, modulus=None):
    """Exact powering; negative modular exponents use checked inversion."""
    _require(base, "base")
    _require(exponent, "exponent")
    if modulus is None:
        if exponent < 0:
            raise ValueError("nonmodular exponent must be nonnegative")
        return builtins.pow(base, exponent)
    _require(modulus, "modulus", 1)
    if is_mpz(base) or is_mpz(modulus) or is_mpz(exponent):
        if exponent < 0:
            base = invert(base, modulus) if modulus > 1 else _load_gmp().mpz(0)
            exponent = -exponent
        return _load_gmp().powmod(base, exponent, modulus)
    if exponent < 0 and modulus > 1:
        base, exponent = invert(base, modulus), -exponent
    return builtins.pow(base, exponent, modulus)


def isqrt(value):
    _require(value, minimum=0)
    return _load_gmp().isqrt(value) if is_mpz(value) else math.isqrt(value)


class NonInvertibleError(ValueError):
    """A failed inverse retains its GCD, including full saturation."""

    def __init__(self, divisor):
        super().__init__("value is not invertible modulo modulus")
        self.divisor = divisor


def invert(value, modulus):
    _require(value)
    _require(modulus, "modulus", 2)
    if is_mpz(value) or is_mpz(modulus):
        try:
            return _load_gmp().invert(value, modulus)
        except ZeroDivisionError as error:
            raise NonInvertibleError(gcd(value, modulus)) from error
    # PyPy's integer Euclidean loop is the accepted inverse baseline.
    a, b = modulus, value % modulus
    coefficient, following = 0, 1
    while b:
        quotient, remainder = divmod(a, b)
        a, b = b, remainder
        coefficient, following = following, coefficient - quotient * following
    if a != 1:
        raise NonInvertibleError(a)
    return coefficient % modulus


def divexact(value, divisor):
    """Check divisibility before GMP divexact; never truncate a nondivisor."""
    _require(value)
    _require(divisor, "divisor")
    if divisor == 0:
        raise ZeroDivisionError("exact division by zero")
    if is_mpz(value) or is_mpz(divisor):
        module = _load_gmp()
        if not module.is_divisible(value, divisor):
            raise ValueError("exact division requires a divisor")
        return module.divexact(value, divisor)
    quotient, remainder = divmod(value, divisor)
    if remainder:
        raise ValueError("exact division requires a divisor")
    return quotient


def gcdext(a, b):
    """Return nonnegative GCD and exact Bezout coefficients."""
    _require(a, "a")
    _require(b, "b")
    if is_mpz(a) or is_mpz(b):
        return _load_gmp().gcdext(a, b)
    old_r, remainder = a, b
    old_x, x, old_y, y = 1, 0, 0, 1
    while remainder:
        quotient = old_r // remainder
        old_r, remainder = remainder, old_r - quotient * remainder
        old_x, x = x, old_x - quotient * x
        old_y, y = y, old_y - quotient * y
    return (-old_r, -old_x, -old_y) if old_r < 0 else (old_r, old_x, old_y)


def integer_root(value, exponent):
    """Floor root with integer Newton iteration or GMP iroot."""
    _require(value, minimum=0)
    _require(exponent, "exponent", 1)
    if value < 2 or exponent == 1:
        return value
    if exponent >= value.bit_length():
        return backend_for(value).integer(1)
    if is_mpz(value):
        return _load_gmp().iroot(value, int(exponent))[0]
    if exponent == 2:
        return math.isqrt(value)
    estimate = 1 << ((value.bit_length() + exponent - 1) // exponent)
    while True:
        following = (
            (exponent - 1) * estimate + value // estimate ** (exponent - 1)
        ) // exponent
        if following >= estimate:
            return estimate
        estimate = following
