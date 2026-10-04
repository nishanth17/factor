"""Verified relation preparation and bounded modular square extraction."""

from collections import Counter
from dataclasses import dataclass
from math import gcd

from .. import utils
from ..budget import Budget
from .factor_base import DEFAULT_MEMORY_BYTES, FactorBase
from .linear_algebra import MAX_MATRIX_ROWS, verify_dependency
from .relations import (
    AtomicRelation,
    CombinedRelation,
    parity_bits,
    verify_atomic,
    verify_combined,
)


@dataclass(frozen=True)
class PreparedRelations:
    """Unique checked full rows; duplicates retain their original indices."""

    relations: tuple
    rows: tuple
    source_indices: tuple
    duplicate_indices: tuple
    factor_base: FactorBase
    atoms: tuple
    workspace_bytes: int


def prepare_relations(
    relations,
    factor_base,
    atom_store,
    *,
    budget=None,
    memory_bytes=DEFAULT_MEMORY_BYTES,
):
    """Verify all inputs before removing exact provenance/payload duplicates.

    Distinct atoms with equal parity remain distinct rows. Partial atoms are
    not full matrix relations. Invalid provenance raises ValueError; work or
    memory refusal publishes nothing and never mutates the caller's store.
    """
    if not isinstance(relations, tuple):
        raise TypeError("relations must be an immutable tuple")
    if len(relations) > MAX_MATRIX_ROWS or len(atom_store) > MAX_MATRIX_ROWS:
        raise ValueError("relation/provenance row cap exceeded")
    utils.require_integer(memory_bytes, "memory_bytes", 0)
    reserve = factor_base.workspace_bytes + 32768
    for atom in atom_store.values():
        if not isinstance(atom, AtomicRelation):
            raise TypeError("atom store must contain atomic relations")
        reserve += 4096 + 256 * len(atom.exponents)
        reserve += 16 * (
            atom.polynomial.n_prime.bit_length()
            + abs(atom.position).bit_length()
        )
    for relation in relations:
        if not isinstance(relation, (AtomicRelation, CombinedRelation)):
            raise TypeError("matrix relation must be atomic or combined")
        reserve += 2048 + 256 * len(relation.exponents)
    if reserve > memory_bytes:
        raise MemoryError("relation provenance exceeds memory_bytes")
    budget = budget if budget is not None else Budget()
    unique, indices, duplicates, seen = [], [], [], set()
    for index, relation in enumerate(relations):
        if isinstance(relation, AtomicRelation):
            if relation.residual != 1:
                raise ValueError("partial atom is not a full matrix row")
            verify_atomic(relation, factor_base, budget=budget)
            identity = (
                "atomic",
                relation.relation_id,
                relation.sign,
                relation.exponents,
            )
        else:
            verify_combined(
                relation,
                factor_base,
                atom_store,
                budget=budget,
                memory_bytes=memory_bytes,
            )
            identity = (
                "combined",
                tuple(sorted(relation.atom_ids)),
                relation.u,
                relation.sign,
                relation.exponents,
                relation.square_correction,
            )
        if identity in seen:
            duplicates.append(index)
        else:
            seen.add(identity)
            unique.append(relation)
            indices.append(index)
    rows = tuple(parity_bits(relation, factor_base) for relation in unique)
    return PreparedRelations(
        tuple(unique),
        rows,
        tuple(indices),
        tuple(duplicates),
        factor_base,
        tuple(atom_store.values()),
        reserve,
    )


@dataclass(frozen=True)
class CongruenceResult:
    """A checked X²=Y² modulo n and both GCD outcomes; divisor may be None."""

    mask: int
    x: int
    y: int
    gcd_minus: int
    gcd_plus: int
    divisor: int | None


def extract_dependency(prepared, mask, *, budget=None):
    """Check original-row parity and even totals before modular extraction.

    Only bounded modular products are formed, including residual square
    corrections. Try both GCD signs and return only a checked proper divisor.
    Budget refusal raises without modifying prepared rows or extraction state.
    """
    verify_dependency(mask, prepared.rows)
    base = prepared.factor_base
    budget = budget if budget is not None else Budget()
    selected = [
        relation
        for index, relation in enumerate(prepared.relations)
        if mask & (1 << index)
    ]
    budget.consume(
        (len(base.entries) + 1) * len(selected) * (base.n.bit_length() + 1)
    )
    totals = Counter()
    x, correction, sign = 1, 1, 1
    for relation in selected:
        x = x * relation.u % base.n
        sign *= relation.sign
        if isinstance(relation, CombinedRelation):
            correction = correction * relation.square_correction % base.n
        for prime, exponent in relation.exponents:
            totals[prime] += exponent
    if sign != 1 or any(exponent % 2 for exponent in totals.values()):
        raise ValueError("dependency exponent totals are not all even")
    y = correction
    for prime, exponent in sorted(totals.items()):
        y = y * pow(prime, exponent // 2, base.n) % base.n
    if x * x % base.n != y * y % base.n:
        raise ValueError("dependency square congruence failed")
    minus, plus = gcd(x - y, base.n), gcd(x + y, base.n)
    divisor = next(
        (
            value
            for value in (minus, plus)
            if utils.valid_divisor(value, base.n)
        ),
        None,
    )
    return CongruenceResult(mask, x, y, minus, plus, divisor)


class DependencyExtractor:
    """Try finite dependencies, retaining position after trivial GCDs."""

    def __init__(self, prepared, dependencies, *, budget=None):
        """Accept at most 4096 masks; validate each during its trial."""
        if not isinstance(dependencies, tuple):
            raise TypeError("dependencies must be an immutable tuple")
        if len(dependencies) > MAX_MATRIX_ROWS:
            raise ValueError("dependency trial cap exceeded")
        self.prepared, self.dependencies = prepared, dependencies
        self.budget = budget if budget is not None else Budget()
        self.next_dependency = 0
        self.trials = []

    def run(self):
        """Return a proper divisor or None; refusal preserves the trial."""
        while self.next_dependency < len(self.dependencies):
            result = extract_dependency(
                self.prepared,
                self.dependencies[self.next_dependency],
                budget=self.budget,
            )
            self.trials.append(result)
            self.next_dependency += 1
            if result.divisor is not None:
                return result.divisor
        return None
