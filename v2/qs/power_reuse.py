"""Bounded exact higher-root transfers between polynomials sharing A."""

from ..common import utils
from ..common.arithmetic import pow
from .polynomial import Polynomial

MAX_TRANSFER_BYTES = 32 * 2**20
MAX_EXTENDED_TRANSFER_BYTES = 128 * 2**20


class RootTransfers:
    """Retain complete nonsingular root pairs, never window-pruned subsets.

    For fixed A and target, A*F(x)=(A*x+B)**2-target. Hence changing B by
    delta moves every root modulo m by -delta/A when A is invertible.
    Characteristic two and exceptional/singular primes keep ordinary lifting.
    """

    def __init__(
        self,
        byte_limit,
        *,
        capacity_limit=MAX_TRANSFER_BYTES,
        pruned_replay=False,
    ):
        utils.require_integer(byte_limit, "root transfer bytes", 0)
        utils.require_integer(capacity_limit, "root transfer capacity", 0)
        if capacity_limit > MAX_EXTENDED_TRANSFER_BYTES:
            raise ValueError("root transfer capacity exceeds 128 MiB")
        if byte_limit > capacity_limit:
            raise ValueError("root transfer bytes exceed the capacity ceiling")
        if type(pruned_replay) is not bool:
            raise ValueError("pruned replay must be Boolean")
        self.pruned_replay = pruned_replay
        self.pruned_replays = 0
        self.byte_limit = max(0, byte_limit - 65536)
        self.used = 0
        self.builds = self.replays = self.extra_work = 0
        self.key = self.polynomial = None
        self.b_words = 0
        self.values = {}

    def bind(self, polynomial):
        """Discard other-A data; bind only the checked immutable polynomial."""
        key = (
            (
                polynomial.n,
                polynomial.multiplier,
                polynomial.a,
                polynomial.square_coefficient,
            )
            if type(polynomial) is Polynomial
            else None
        )
        if key != self.key or key is None:
            self.values.clear()
            self.used = 0
        self.key = key
        self.polynomial = polynomial if key is not None else None
        self.b_words = (
            (abs(polynomial.b).bit_length() + 7) // 8
            if self.polynomial is not None
            else 0
        )

    def _charge(self, modulus):
        # Include map/tuple/limb capacity and old/new update coexistence.
        # The separate 64 KiB reserve bounds single-entry construction scratch.
        words = self.b_words + 4 * ((modulus.bit_length() + 7) // 8)
        return 2048 + 64 * words

    def shifted(self, polynomial, prime, modulus, initial, cost, budget):
        """Return a complete pair in the first-root traversal order."""
        if polynomial is not self.polynomial:
            return None
        previous = self.values.get(modulus)
        if previous is None:
            return None
        b, roots, inverse, charge, b_words = previous
        if b == polynomial.b:
            budget.consume(cost)
            budget.consume(cost)
            self.replays += 1
            return roots
        # B's width changes only at bind, not at every prime-power lookup.
        # Retaining a previous larger charge is safe and avoids recounting
        # invariant modulus storage or shrinking ownership on every shift.
        updated_words, updated_charge = b_words, charge
        if self.b_words > b_words:
            updated_words = self.b_words
            updated_charge += 64 * (updated_words - b_words)
            if self.used - charge + updated_charge > self.byte_limit:
                return None

        # Each root retains the original logical lift reservation. Reordering
        # by its first-level class preserves sparse event order after shifts.
        budget.consume(cost)
        delta = (polynomial.b - b) * inverse % modulus
        # Roots and delta are canonical, so each subtraction needs at most
        # one corrective addition rather than another division/remainder.
        first = roots[0] - delta
        if first < 0:
            first += modulus
        budget.consume(cost)
        second = roots[1] - delta
        if second < 0:
            second += modulus
        shifted = (
            (first, second) if first % prime == initial[0] else (second, first)
        )
        self.values[modulus] = (
            polynomial.b,
            shifted,
            inverse,
            updated_charge,
            updated_words,
        )
        self.used += updated_charge - charge
        self.replays += 1
        return shifted

    def shifted_class(
        self, polynomial, prime, modulus, initial, active, cost, budget
    ):
        """Select one active first-level class from a complete cached pair.

        Window pruning removes positions, never changes a surviving root's
        class modulo p. Replay the complete pair under the existing fixed-A
        identity, then select that class; retain both roots in the cache.
        Eligibility guarantees a unique nonsingular lift in each class.
        """
        if not self.pruned_replay:
            return None
        active %= prime
        if active not in initial:
            raise ValueError("active root is outside the initial classes")
        pair = self.shifted(polynomial, prime, modulus, initial, cost, budget)
        if pair is None:
            return None
        selected = tuple(root for root in pair if root % prime == active)
        if len(selected) != 1:
            raise ArithmeticError("cached root class is not unique")
        # shifted() funded both roots, although ordinary lifting needed only
        # the active one. Record that extra logical charge, alongside cache
        # construction, for exact work-accounting comparisons.
        self.extra_work += cost
        self.pruned_replays += 1
        return selected

    def remember(self, polynomial, parent_modulus, modulus, roots, budget):
        """Publish a complete pair only after its independent Hensel lift."""
        if polynomial is not self.polynomial or modulus in self.values:
            return
        charge = self._charge(modulus)
        if self.used + charge > self.byte_limit:
            return
        work = polynomial.a.bit_length() * modulus.bit_length() + 1
        if work > budget.work_limit - budget.used:
            # This inverse is optional; do not refuse an already completed
            # ordinary lift merely because its reusable copy is unaffordable.
            return
        budget.consume(work)
        self.extra_work += work
        previous = self.values.get(parent_modulus)
        if previous is None:
            inverse = pow(polynomial.a, -1, modulus)
        else:
            # Newton doubles the known p-adic inverse precision. For a
            # parent p**e, e>=1, that covers the next modulus p**(e+1).
            inverse = previous[2]
            inverse = inverse * (2 - polynomial.a * inverse) % modulus
        self.values[modulus] = (
            polynomial.b,
            roots,
            inverse,
            charge,
            self.b_words,
        )
        self.used += charge
        self.builds += 1
