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
    """Unique checked rows and atom ownership, with shared reservations.

    shared_workspace_bytes covers the base and pinned atoms. A pipeline
    retaining their original collector counts that reservation once, then
    adds preparation's private containers and scratch to its live workspace.
    """

    relations: tuple
    rows: tuple
    source_indices: tuple
    duplicate_indices: tuple
    factor_base: FactorBase
    atoms: tuple
    workspace_bytes: int
    shared_workspace_bytes: int = 0


class _VerificationCache:
    """Bounded attestations for immutable rows and identical provenance.

    The owning collector reserves this storage and clears it before releasing
    its pinned relations. Public preparation never accepts caller attestations.
    A cache entry is published only after full exact verification succeeds.
    """

    def __init__(self, memory_bytes):
        self.memory_bytes = memory_bytes
        self.records = {}
        self.used = self.hits = self.misses = 0

    def clear(self):
        self.records.clear()
        self.used = 0

    def check(self, relation, base, store, budget, memory_bytes):
        record = self.records.get(id(relation))
        atoms = (
            tuple(store.get(identity) for identity in relation.atom_ids)
            if isinstance(relation, CombinedRelation)
            else ()
        )
        budget.consume(len(atoms) + 1)
        if record is not None and (
            record[0] is relation
            and record[1] is base
            and len(record[2]) == len(atoms)
            and all(a is b for a, b in zip(record[2], atoms))
        ):
            self.hits += 1
            return record[3]
        self.misses += 1
        _verify_row(relation, base, store, budget, memory_bytes)
        row = parity_bits(relation, base)
        reserve = 1024 + 256 * len(atoms) + 8 * ((row.bit_length() + 7) // 8)
        if record is None and self.used + reserve <= self.memory_bytes:
            self.records[id(relation)] = (relation, base, atoms, row)
            self.used += reserve
        return row


def _verify_row(relation, base, store, budget, memory_bytes):
    if isinstance(relation, AtomicRelation):
        if relation.residual != 1:
            raise ValueError("partial atom is not a full matrix row")
        verify_atomic(relation, base, budget=budget)
    else:
        verify_combined(
            relation, base, store, budget=budget, memory_bytes=memory_bytes
        )


def prepare_relations(
    relations,
    factor_base,
    atom_store,
    *,
    budget=None,
    memory_bytes=DEFAULT_MEMORY_BYTES,
    retained_workspace_bytes=0,
):
    """Verify and deduplicate immutable rows, with full public validation."""
    return _prepare_relations(
        relations,
        factor_base,
        atom_store,
        budget=budget,
        memory_bytes=memory_bytes,
        retained_workspace_bytes=retained_workspace_bytes,
    )


def _prepare_relations(
    relations,
    factor_base,
    atom_store,
    *,
    budget=None,
    memory_bytes=DEFAULT_MEMORY_BYTES,
    retained_workspace_bytes=0,
    verification_cache=None,
):
    """Verify all inputs before removing exact provenance/payload duplicates.

    Distinct atoms with equal parity remain distinct rows. Partial atoms are
    not full matrix relations. Invalid provenance raises ValueError; work or
    memory refusal publishes nothing and never mutates the caller's store.
    retained_workspace_bytes includes a collector retaining the same base
    and atoms; its private storage must coexist with preparation scratch.
    """
    if not isinstance(relations, tuple):
        raise TypeError("relations must be an immutable tuple")
    if len(relations) > MAX_MATRIX_ROWS or len(atom_store) > MAX_MATRIX_ROWS:
        raise ValueError("relation/provenance row cap exceeded")
    utils.require_integer(memory_bytes, "memory_bytes", 0)
    utils.require_integer(
        retained_workspace_bytes, "retained_workspace_bytes", 0
    )
    shared_reserve = factor_base.workspace_bytes
    for atom in atom_store.values():
        if not isinstance(atom, AtomicRelation):
            raise TypeError("atom store must contain atomic relations")
        shared_reserve += 4096 + 256 * len(atom.exponents)
        shared_reserve += 16 * (
            atom.polynomial.n_prime.bit_length()
            + abs(atom.position).bit_length()
        )
    reserve = shared_reserve + 32768
    for relation in relations:
        if not isinstance(relation, (AtomicRelation, CombinedRelation)):
            raise TypeError("matrix relation must be atomic or combined")
        reserve += 2048 + 256 * len(relation.exponents)
    simultaneous = retained_workspace_bytes + reserve - shared_reserve
    if max(reserve, simultaneous) > memory_bytes:
        raise MemoryError("relation provenance exceeds memory_bytes")
    budget = budget if budget is not None else Budget()
    unique, indices, duplicates, seen, rows = [], [], [], set(), []
    for index, relation in enumerate(relations):
        if verification_cache is None:
            _verify_row(
                relation, factor_base, atom_store, budget, memory_bytes
            )
            row = parity_bits(relation, factor_base)
        else:
            row = verification_cache.check(
                relation, factor_base, atom_store, budget, memory_bytes
            )
        if isinstance(relation, AtomicRelation):
            identity = (
                "atomic",
                relation.relation_id,
                relation.sign,
                relation.exponents,
            )
        else:
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
            rows.append(row)
    return PreparedRelations(
        tuple(unique),
        tuple(rows),
        tuple(indices),
        tuple(duplicates),
        factor_base,
        tuple(atom_store.values()),
        reserve,
        shared_reserve,
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
    selected, bits = [], mask
    while bits:
        bit = bits & -bits
        selected.append(prepared.relations[bit.bit_length() - 1])
        bits ^= bit
    budget.consume(
        (len(base.entries) + 1) * len(selected) * (base.n.bit_length() + 1)
    )
    totals = Counter()
    x, correction, sign = 1, 1, 1
    for relation in selected:
        x = x * relation.u % base.n
        sign *= relation.sign
        correction = (
            correction * getattr(relation, "square_correction", 1) % base.n
        )
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
        """Accept at most 65536 masks; validate each during its trial."""
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
