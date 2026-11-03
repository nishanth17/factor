"""Explicit contracts for three proven endpoints and opaque square pairing."""

from .relations import MAX_RESIDUAL, AtomicRelation


class ResidualAtom(AtomicRelation):
    """A generalized residual; an opaque token is never labelled prime.

    Equal opaque tokens can cancel as R² without factoring R. Their atoms
    remain available for later exact endpoint conversion in the incidence
    collector. Verification must explicitly enable this contract.
    """

    def __init__(
        self, polynomial, position, sign, exponents, *, primes=(), opaque=None
    ):
        original = AtomicRelation(polynomial, position, sign, exponents)
        if not isinstance(primes, tuple) or len(primes) not in (1, 2, 3):
            if primes or opaque is None:
                raise ValueError(
                    "extended residual needs one to three endpoints"
                )
        if primes and opaque is not None:
            raise ValueError(
                "opaque and proven endpoint contracts are exclusive"
            )
        if any(
            type(p) is not int or not 2 <= p <= MAX_RESIDUAL for p in primes
        ):
            raise ValueError("extended endpoint exceeds its integer bound")
        if tuple(sorted(primes)) != primes:
            raise ValueError("extended endpoints must be sorted")
        if opaque is not None and (
            type(opaque) is not int or not 4 <= opaque <= MAX_RESIDUAL**3
        ):
            raise ValueError("opaque cofactor exceeds its integer bound")
        for key, value in original.__dict__.items():
            object.__setattr__(self, key, value)
        object.__setattr__(self, "large_primes", primes)
        object.__setattr__(self, "opaque", opaque)

    @property
    def residual_primes(self):
        # Generic arithmetic consumers use these as cancellation factors.
        # Only the explicit verifier decides which carry a prime certificate.
        return (self.opaque,) if self.opaque is not None else self.large_primes

    @property
    def residual_product(self):
        result = 1
        for value in self.residual_primes:
            result *= value
        return result

    @property
    def incidence_tokens(self):
        # Domain-separate opaque R from prime endpoint R in incidence space.
        return tuple(
            (self.opaque is not None, p) for p in self.residual_primes
        )
