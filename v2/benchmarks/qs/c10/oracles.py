"""Independent small arithmetic oracles for C10 collection hypotheses.

These bounded diagnostics do not call production root, score or factor-base
recovery routines. They are not an alternative factoring API.
"""

from collections import Counter
from dataclasses import dataclass
from math import gcd, prod


@dataclass(frozen=True)
class DiscriminantPolynomial:
    """P=a*x*x+b*x+c with (2*a*x+b)**2-D == 4*a*P exactly."""

    discriminant: int
    a: int
    b: int

    def __post_init__(self):
        values = (self.discriminant, self.a, self.b)
        if any(type(value) is not int for value in values):
            raise TypeError("oracle coefficients must be integers")
        if self.discriminant <= 0 or self.a <= 0 or self.a % 2 == 0:
            raise ValueError(
                "positive discriminant and odd positive A required"
            )
        if self.discriminant % 4 != 1 or self.b % 2 == 0:
            raise ValueError("normalized form requires D=1 mod 4 and odd B")
        if (self.b * self.b - self.discriminant) % (4 * self.a):
            raise ValueError("B is not a square root of D modulo 4A")

    @property
    def c(self):
        return (self.b * self.b - self.discriminant) // (4 * self.a)

    def value(self, x):
        return self.a * x * x + self.b * x + self.c

    def u_value(self, x):
        return 2 * self.a * x + self.b

    @classmethod
    def from_standard(cls, discriminant, a, b):
        """Choose the odd lift of a standard B, plus its position parity."""
        delta = 0 if b % 2 else 1
        return cls(discriminant, a, b + delta * a), delta

    def roots(self, prime):
        """Derive roots independently from enumerated small square roots."""
        if not 2 <= prime <= 257 or any(
            prime % d == 0 for d in range(2, prime)
        ):
            raise ValueError("oracle requires a prime at most 257")
        if prime == 2:
            return tuple(x for x in (0, 1) if self.value(x) % 2 == 0)
        if self.a % prime == 0:
            if self.b % prime == 0:
                return tuple(range(prime)) if self.c % prime == 0 else ()
            return ((-self.c * pow(self.b, -1, prime)) % prime,)

        inverse = pow(2 * self.a, -1, prime)
        square_roots = (
            r for r in range(prime) if (r * r - self.discriminant) % prime == 0
        )
        return tuple(
            sorted({(r - self.b) * inverse % prime for r in square_roots})
        )


def divide_over_base(value, primes):
    """Return sign, exact exponents and unresolved positive residual."""
    if not value or abs(value).bit_length() > 4096 or len(primes) > 128:
        raise ValueError("outside the finite nonzero oracle domain")
    if len(set(primes)) != len(primes) or any(p < 2 for p in primes):
        raise ValueError("invalid oracle base")
    residual, exponents = abs(value), []
    for prime in primes:
        exponent = 0
        while residual % prime == 0:
            residual //= prime
            exponent += 1
        exponents.append(exponent)
    sign = -1 if value < 0 else 1
    if (
        sign * residual * prod(p**e for p, e in zip(primes, exponents))
        != value
    ):
        raise AssertionError("oracle recovery did not reconstruct")
    return sign, tuple(exponents), residual


def product_census(assignments, target):
    """Distinguish unique families, A error and pairwise factor overlap."""
    if not 1 <= len(assignments) <= 256 or target < 1:
        raise ValueError("outside the finite assignment census")
    products = [prod(values) for values in assignments]
    errors = sorted(abs(a - target) for a in products)
    overlaps, adjacent = Counter(), Counter()
    sets = [set(values) for values in assignments]
    for i, left in enumerate(sets):
        for j in range(i):
            overlap = len(left & sets[j])
            overlaps[overlap] += 1
            if j == i - 1:
                adjacent[overlap] += 1
    return dict(
        count=len(assignments),
        unique_products=len(set(products)),
        target=target,
        minimum_a=min(products),
        maximum_a=max(products),
        error_sum=sum(errors),
        error_median_numerator=errors[(len(errors) - 1) // 2]
        + errors[len(errors) // 2],
        error_median_denominator=2,
        maximum_error=errors[-1],
        pair_overlap=dict(overlaps),
        adjacent_overlap=dict(adjacent),
    )


def threshold_losses(values, scores, thresholds, primes, bound, extras):
    """Audit deliberate report loss, not downstream usefulness or primality.

    Full division supplies an independent residual-magnitude oracle. A value
    within this bound may still fail C9 admission or never yield a dependency.
    This distinction prevents tuning on raw report count as a factor proxy.
    """
    if not 0 < len(values) <= 4096 or not (
        len(values) == len(scores) == len(thresholds)
    ):
        raise ValueError("invalid finite threshold census")
    if bound < 1 or not 1 <= len(extras) <= 16 or any(x < 0 for x in extras):
        raise ValueError("invalid threshold arms")
    eligible = [divide_over_base(v, primes)[2] <= bound for v in values]
    return {
        extra: dict(
            reports=sum(s >= t + extra for s, t in zip(scores, thresholds)),
            eligible=sum(eligible),
            lost=sum(
                ok and s < t + extra
                for ok, s, t in zip(eligible, scores, thresholds)
            ),
        )
        for extra in extras
    }


def multiplier_divisor(n, multiplier):
    """Keep an exceptional multiplier GCD separate from a usable choice."""
    divisor = gcd(n, multiplier)
    if 1 < divisor < n:
        return divisor
    return None
