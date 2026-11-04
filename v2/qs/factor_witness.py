"""Exact representation witnesses without decoding trusted compact factors."""

from .cold_factors import RecomputedFactors
from .indexed_relations import FactorRecord, IndexedFactors


def polynomial_snapshot(polynomial):
    """Bind every arithmetic/provenance field, including bypass mutations."""
    return (
        polynomial.n,
        polynomial.multiplier,
        polynomial.a,
        polynomial.b,
        polynomial.c,
        polynomial.square_coefficient,
        polynomial.identity,
    )


def factor_snapshot(factors):
    """Copy mutable storage; share only deeply immutable representations.

    Subclasses and foreign sequences retain the independent sequence witness.
    Recomputed payloads bind their own provenance, not just the enclosing atom.
    """
    if type(factors) is FactorRecord:
        if type(factors._primes) is not tuple:
            raise ValueError("packed factor base was modified")
        return "packed", factors._primes, len(factors), factors.encoded()
    if type(factors) is IndexedFactors:
        if (
            type(factors.primes) is not tuple
            or type(factors.indices) is not tuple
            or any(
                type(pair) is not tuple
                or len(pair) != 2
                or any(type(value) is not int for value in pair)
                for pair in factors.indices
            )
        ):
            raise ValueError("indexed factor representation was modified")
        return "indices", factors.primes, factors.indices
    if type(factors) is RecomputedFactors:
        if type(factors.small) is not tuple or any(
            type(prime) is not int for prime in factors.small
        ):
            raise ValueError("recomputed factor representation was modified")
        if type(factors.retained) is RecomputedFactors:
            raise ValueError("nested recomputation is not a packed record")
        return (
            "recomputed",
            polynomial_snapshot(factors.polynomial),
            factors.position,
            factors.small,
            factor_snapshot(factors.retained),
            len(factors),
        )
    return "sequence", tuple(factors)
