"""Tiny exhaustive QS/MPQS collection: exact division at every position."""

from dataclasses import dataclass
from math import gcd

from .. import utils
from ..budget import Budget, BudgetExhaustedError
from .factor_base import DEFAULT_MEMORY_BYTES
from .polynomial import checked_position
from .relations import AtomicRelation, checked_residual_bound, verify_atomic

MAX_BLOCK_WIDTH = 4096


@dataclass(frozen=True)
class CollectionResult:
    """Checked atoms and a finite stop boundary for the half-open block.

    next_position is the first unprocessed x. scanned counts committed
    positions, including rejected candidates and zeros. A proper divisor
    is distinct from relation yield. Storage is never evicted implicitly.
    """

    relations: tuple[AtomicRelation, ...]
    divisor: int | None
    scanned: int
    next_position: int
    reason: str
    zero_positions: tuple[int, ...]
    workspace_bytes: int


def collect_block(
    polynomial,
    factor_base,
    lo,
    hi,
    *,
    residual_bound=1,
    max_relations=MAX_BLOCK_WIDTH,
    memory_bytes=DEFAULT_MEMORY_BYTES,
    budget=None,
):
    """Collect every admissible relation in [lo,hi), including negative x.

    Full-only is the default. Optional small proven residuals support
    combination fixtures; there is no partial matching store or score sieve.
    Caps return a checked prefix and first unprocessed position. Temporary
    bigint and retained-atom workspace are reserved before work/publication.
    The reservation estimates owned bytes, not process RSS. Setup/root work
    shares the allowance only when the caller passes the same Budget.
    """
    checked_position(lo)
    checked_position(hi)
    checked_residual_bound(residual_bound)
    utils.require_integer(max_relations, "max_relations", 0)
    utils.require_integer(memory_bytes, "memory_bytes", 0)
    if not 0 <= hi - lo <= MAX_BLOCK_WIDTH:
        raise ValueError("block must be ordered and at most 4096 positions")
    if max_relations > MAX_BLOCK_WIDTH:
        raise ValueError("max_relations exceeds the reference limit")
    if (polynomial.n, polynomial.multiplier) != (
        factor_base.n,
        factor_base.multiplier,
    ):
        raise ValueError("polynomial and factor-base identity mismatch")
    budget = budget if budget is not None else Budget()
    relations, zero_positions = [], []
    position, reason, divisor = lo, "complete", None
    # A*x and its square are finite because coefficient/position bits are
    # capped. Reserve several simultaneous bigint temporaries at that bound.
    value_bits = 2 * (
        polynomial.a.bit_length()
        + max(abs(lo).bit_length(), abs(hi).bit_length())
        + abs(polynomial.b).bit_length()
        + factor_base.n_prime.bit_length()
        + 2
    )
    workspace = factor_base.workspace_bytes + 8192 + 8 * value_bits
    primes = factor_base.primes
    try:
        while position < hi:
            if len(relations) >= max_relations:
                reason = "relation_limit"
                break
            if workspace > memory_bytes:
                reason = "memory_limit"
                break

            budget.consume((len(primes) + 1) * value_bits)
            value = polynomial.a * polynomial.value(position)
            if value == 0:
                # Division of zero never terminates. Surface a square-target
                # GCD when proper, otherwise record the exceptional position.
                if workspace + 64 > memory_bytes:
                    reason = "memory_limit"
                    break
                zero_positions.append(position)
                workspace += 64
                candidate = gcd(polynomial.u_value(position), factor_base.n)
                position += 1
                if utils.valid_divisor(candidate, factor_base.n):
                    divisor, reason = candidate, "factor_found"
                    break
                continue

            remaining, exponents = abs(value), []
            for prime in primes:
                exponent = 0
                while remaining % prime == 0:
                    remaining //= prime
                    exponent += 1
                if exponent:
                    exponents.append((prime, exponent))

            atom = None
            if remaining <= residual_bound:
                budget.consume(remaining.bit_length() ** 2)
                if remaining == 1 or (
                    utils.classify_prime(remaining) == utils.Primality.PROVEN
                ):
                    atom = AtomicRelation(
                        polynomial,
                        position,
                        -1 if value < 0 else 1,
                        tuple(exponents),
                        remaining,
                    )
                    verify_atomic(
                        atom,
                        factor_base,
                        residual_bound=residual_bound,
                        budget=budget,
                    )

            if atom is not None:
                reserve = 1024 + 128 * len(exponents) + value_bits
                if workspace + reserve > memory_bytes:
                    reason = "memory_limit"
                    break
                relations.append(atom)
                workspace += reserve
            position += 1

    except BudgetExhaustedError:
        reason = budget.reason
    return CollectionResult(
        tuple(relations),
        divisor,
        position - lo,
        position,
        reason,
        tuple(zero_positions),
        workspace,
    )
