"""Exact normalized QS/MPQS polynomials and exceptional modular roots."""

import hashlib
from dataclasses import dataclass, field
from math import isqrt

from .. import utils
from ..budget import Budget
from .factor_base import MAX_INPUT_BITS, checked_target

MAX_COEFFICIENT_BITS = MAX_INPUT_BITS
MAX_POSITION_BITS = MAX_INPUT_BITS


def checked_position(position):
    """Validate a signed, finite-size reference-collector position."""
    utils.require_integer(position, "position")
    if abs(position).bit_length() > MAX_POSITION_BITS:
        raise ValueError("position exceeds the reference bit limit")
    return position


@dataclass(frozen=True)
class Polynomial:
    """F(x)=A*x*x+2*B*x+C with U(x)^2-h*n=A*F(x), exactly.

    A is positive; B may be negative. C is derived only after checking
    divisibility. Immutable identities encode n,h,A,B in hexadecimal.
    square_coefficient records a known square factor of A, defaulting to 1.
    No floating-point approximation or primality assertion is made.
    """

    n: int
    multiplier: int
    a: int
    b: int
    square_coefficient: int = 1
    c: int = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self):
        """Validate finite coefficients and derive C and stable identity."""
        target = checked_target(self.n, self.multiplier)
        utils.require_integer(self.a, "a", 1)
        utils.require_integer(self.b, "b")
        if max(self.a.bit_length(), abs(self.b).bit_length()) > (
            MAX_COEFFICIENT_BITS
        ):
            raise ValueError("polynomial coefficient exceeds the bit limit")
        utils.require_integer(self.square_coefficient, "square_coefficient", 1)
        if self.square_coefficient.bit_length() > MAX_COEFFICIENT_BITS:
            raise ValueError("square coefficient exceeds the bit limit")
        if self.a % (self.square_coefficient * self.square_coefficient):
            raise ValueError("square coefficient squared must divide A")
        quotient, remainder = divmod(self.b * self.b - target, self.a)
        if remainder:
            raise ValueError("B*B-h*n must be divisible by A")
        object.__setattr__(self, "c", quotient)
        encoded = ":".join(
            hex(value) for value in (self.n, self.multiplier, self.a, self.b)
        )
        if self.square_coefficient != 1:
            encoded += ":square:" + hex(self.square_coefficient)
        object.__setattr__(
            self, "identity", hashlib.sha256(encoded.encode()).hexdigest()
        )

    @property
    def n_prime(self):
        """Return h*n as an exact integer."""
        return self.n * self.multiplier

    @property
    def supported_a(self):
        """A after removing the explicitly represented square coefficient."""
        return self.a // (self.square_coefficient * self.square_coefficient)

    def value(self, position):
        """Evaluate normalized F at a signed integer position."""
        checked_position(position)
        return (self.a * position + 2 * self.b) * position + self.c

    def u_value(self, position):
        """Return U=A*x+B; U squared is congruent to A*F modulo n."""
        checked_position(position)
        return self.a * position + self.b


def a_target(n_prime, half_width):
    """Return max(1,floor(sqrt(2*N_prime)/M)) using integer arithmetic.

    M is a positive half-width, not the number of positions in a block.
    This target is a design control; no tuned SIQS choice is implied.
    """
    utils.require_integer(n_prime, "n_prime", 1)
    utils.require_integer(half_width, "half_width", 1)
    if half_width.bit_length() > MAX_POSITION_BITS:
        raise ValueError("half_width exceeds the reference bit limit")
    if n_prime.bit_length() > MAX_INPUT_BITS + 20:
        raise ValueError("target exceeds the multiplied-input bit limit")
    return max(1, isqrt(2 * n_prime) // half_width)


def qs_polynomial(factor_base):
    """Return the A=1 QS polynomial centered at ceil(sqrt(h*n))."""
    center = isqrt(factor_base.n_prime)
    if center * center != factor_base.n_prime:
        center += 1
    return Polynomial(factor_base.n, factor_base.multiplier, 1, center)


def mpqs_polynomial(factor_base, half_width, *, budget=None, prime=None):
    """Choose A=q squared nearest the integer target and lift B modulo A.

    q is an odd factor-base prime not dividing h*n. Lift a cached root r
    with B=r+q*t, where 2*r*t=(h*n-r*r)/q modulo q; the coefficient is a
    unit. Choose the smaller of B and A-B reproducibly. This reference
    selects one polynomial, not SIQS families or a tuned production policy.
    Missing eligible primes or inversion failure raise ValueError.
    """
    target = a_target(factor_base.n_prime, half_width)
    budget = budget if budget is not None else Budget()
    budget.consume(len(factor_base.entries))
    eligible = [
        entry
        for entry in factor_base.entries
        if entry.prime != 2 and factor_base.n_prime % entry.prime != 0
    ]
    if not eligible:
        raise ValueError("MPQS requires an odd nonsingular factor-base prime")
    if prime is None:
        entry = min(
            eligible,
            key=lambda item: (abs(item.prime**2 - target), item.prime),
        )
    else:
        utils.require_integer(prime, "MPQS prime", 3)
        entry = next((item for item in eligible if item.prime == prime), None)
        if entry is None:
            raise ValueError("MPQS prime must be a nonsingular base member")

    prime, root = entry.prime, entry.square_roots[0]
    budget.consume(factor_base.n_prime.bit_length() + prime.bit_length() ** 2)
    inverse = utils.modular_inverse(2 * root, prime)
    correction = (
        ((factor_base.n_prime - root * root) // prime) * inverse % prime
    )
    a = prime * prime
    b = root + prime * correction
    b = min(b, a - b)
    return Polynomial(factor_base.n, factor_base.multiplier, a, b)


@dataclass(frozen=True)
class PolynomialRoots:
    """Canonical residues, or all positions for an identically zero value.

    all_positions=True uses an empty roots tuple instead of allocating p
    residues. An empty tuple with False means no positions are roots.
    """

    prime: int
    roots: tuple[int, ...]
    all_positions: bool = False


def polynomial_roots(polynomial, factor_base, entry, *, budget=None):
    """Derive roots of normalized F using checked cached target roots.

    Handle 2 by enumeration; p|A is linear or constant. Inversion is used
    only for nonzero coefficients modulo a proven prime. A failed inverse
    propagates ValueError; it never silently drops a root or yields a split.
    """
    if (polynomial.n, polynomial.multiplier) != (
        factor_base.n,
        factor_base.multiplier,
    ) or entry not in factor_base.entries:
        raise ValueError("polynomial and factor-base identity mismatch")

    budget = budget if budget is not None else Budget()
    prime = entry.prime
    budget.consume(prime.bit_length() ** 2)
    if prime == 2:
        roots = tuple(x for x in (0, 1) if polynomial.value(x) % 2 == 0)
        if len(roots) == 2:
            return PolynomialRoots(2, (), True)
        return PolynomialRoots(2, roots)

    if polynomial.a % prime == 0:
        coefficient = 2 * polynomial.b % prime
        constant = polynomial.c % prime
        if coefficient == 0:
            return PolynomialRoots(prime, (), constant == 0)
        inverse = utils.modular_inverse(coefficient, prime)
        return PolynomialRoots(prime, ((-constant * inverse) % prime,))

    inverse = utils.modular_inverse(polynomial.a, prime)
    roots = tuple(
        sorted(
            {
                (root - polynomial.b) * inverse % prime
                for root in entry.square_roots
            }
        )
    )
    return PolynomialRoots(prime, roots)
