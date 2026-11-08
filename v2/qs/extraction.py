"""Verified relation preparation and bounded modular square extraction."""

from collections import Counter
from dataclasses import dataclass

from ..common import utils
from ..common.arithmetic import gcd, pow
from ..execution.budget import Budget
from ..execution.resources import removed_arguments
from .factor_base import FactorBase
from .linear_algebra import verify_dependency
from .relations import (
    AtomicRelation,
    CombinedRelation,
    parity_bits,
    verify_atomic,
    verify_combined,
)
from .storage import storage_limit


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
    row_identities: tuple = ()


class _VerificationCache:
    """Bounded attestations for immutable rows and identical provenance.

    The owning collector reserves this storage and clears it before releasing
    its pinned relations. Public preparation never accepts caller attestations.
    A cache entry is published only after full exact verification succeeds.
    """

    def __init__(self, memory_bytes, *, checked=None):
        self.memory_bytes = memory_bytes
        self.records = {}
        self.tested_dependencies = set()
        self.used = self.hits = self.misses = 0
        self.checked = checked

    def clear(self):
        self.records.clear()
        self.tested_dependencies.clear()
        self.used = 0

    @removed_arguments(("max_atoms", "max_rows", "memory_bytes"))
    def check(self, relation, base, store, budget):
        if self.checked is not None:
            if self.checked.base is not base:
                raise ValueError("checked preparation uses a foreign base")
            return self.checked.row(relation, store, budget)
        record = self.records.get(id(relation))
        atoms = (
            tuple(store.get(identity) for identity in relation.atom_ids)
            if isinstance(relation, CombinedRelation)
            else ()
        )
        budget.consume(len(atoms) + 1)
        # An object ID alone may be reused; pinned payloads must also match.
        if record is not None and (
            record[0] is relation
            and record[1] is base
            and len(record[2]) == len(atoms)
            and all(a is b for a, b in zip(record[2], atoms))
        ):
            self.hits += 1
            return record[3]

        self.misses += 1
        _verify_row(relation, base, store, budget)
        row = parity_bits(relation, base)
        reserve = 1024 + 256 * len(atoms) + 8 * ((row.bit_length() + 7) // 8)
        if record is None and self.used + reserve <= self.memory_bytes:
            self.records[id(relation)] = (relation, base, atoms, row)
            self.used += reserve
        return row

    def dependency_key(self, prepared, mask, budget):
        """Bind a tested selection to complete immutable row/base payloads."""
        verify_dependency(mask, prepared.rows)
        budget.consume(
            mask.bit_count() * (len(prepared.factor_base.entries) + 1) + 1
        )
        if len(prepared.row_identities) != len(prepared.rows):
            raise ValueError("prepared row identities are missing")
        identities, bits = [], mask

        while bits:
            bit = bits & -bits
            index = bit.bit_length() - 1
            identity = prepared.row_identities[index]
            if identity[0] is not prepared.factor_base or (
                identity[1] is not prepared.relations[index]
            ):
                raise ValueError("prepared row identity differs from payload")
            identities.append(identity)
            bits ^= bit

        return frozenset(identities)

    def remember_dependency(self, key):
        """Retain only checked trivial trials inside the shared cache cap."""
        reserve = 2048 + sum(512 + 128 * len(identity[2]) for identity in key)
        if key not in self.tested_dependencies and (
            self.used + reserve <= self.memory_bytes
        ):
            self.tested_dependencies.add(key)
            self.used += reserve


@removed_arguments(("max_atoms", "max_rows", "memory_bytes"))
def _verify_row(relation, base, store, budget):
    if isinstance(relation, AtomicRelation):
        if relation.residual != 1 or getattr(relation, "large_primes", ()):
            raise ValueError("partial atom is not a full matrix row")
        verify_atomic(relation, base, budget=budget)
    else:
        verify_combined(relation, base, store, budget=budget)


@removed_arguments(("max_atoms", "max_rows", "memory_bytes"))
def prepare_relations(
    relations,
    factor_base,
    atom_store,
    *,
    budget=None,
    retained_workspace_bytes=0,
):
    """Verify and deduplicate immutable rows, with full public validation."""
    return _prepare_relations(
        relations,
        factor_base,
        atom_store,
        budget=budget,
        retained_workspace_bytes=retained_workspace_bytes,
    )


@removed_arguments(("max_atoms", "max_rows", "memory_bytes"))
def _prepare_relations(
    relations,
    factor_base,
    atom_store,
    *,
    budget=None,
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
        if isinstance(relation, CombinedRelation):
            reserve += 128 * len(relation.atom_ids)

    # Shared atoms/base stay pinned once; private scratch must coexist.
    retained_workspace_bytes + reserve - shared_reserve
    budget = budget if budget is not None else Budget()
    unique, indices, duplicates, seen, rows = [], [], [], set(), []
    identities = []

    for index, relation in enumerate(relations):
        if verification_cache is None:
            _verify_row(relation, factor_base, atom_store, budget)
            row = parity_bits(relation, factor_base)
        else:
            row = verification_cache.check(
                relation, factor_base, atom_store, budget
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
            atoms = (
                tuple(atom_store[i] for i in relation.atom_ids)
                if isinstance(relation, CombinedRelation)
                else ()
            )
            # Native equality includes every immutable polynomial, exponent,
            # sign, residual, square correction and factor-base/root payload.
            identities.append((factor_base, relation, atoms))

    return PreparedRelations(
        tuple(unique),
        tuple(rows),
        tuple(indices),
        tuple(duplicates),
        factor_base,
        tuple(atom_store.values()),
        reserve,
        shared_reserve,
        tuple(identities),
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


def extract_dependency(
    prepared, mask, *, budget=None, product_policy="reference"
):
    """Check original-row parity and even totals before modular extraction.

    Only bounded modular products are formed, including residual square
    corrections. Try both GCD signs and return only a checked proper divisor.
    Budget refusal raises without modifying prepared rows or extraction state.
    """
    if product_policy not in ("reference", "joint-v2"):
        raise ValueError("unknown exponent-product policy")
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
        from .lazy_relations import LazyCombinedRelation

        if isinstance(relation, LazyCombinedRelation):
            relation.bind_budget(budget)
            root, square = relation.materialize()
        else:
            root, square = (
                relation.u,
                getattr(relation, "square_correction", 1),
            )
        x = x * root % base.n
        sign *= relation.sign
        correction = correction * square % base.n
        for prime, exponent in relation.exponents:
            totals[prime] += exponent

    if sign != 1 or any(exponent % 2 for exponent in totals.values()):
        raise ValueError("dependency exponent totals are not all even")
    if product_policy == "joint-v2":
        y = (
            correction
            * joint_exponent_product(totals, base.n, budget)
            % base.n
        )
    else:
        y = correction
        for prime, exponent in sorted(totals.items()):
            y = y * pow(prime, exponent // 2, base.n) % base.n
    if x * x % base.n != y * y % base.n:
        raise ValueError("dependency square congruence failed")
    # Either sign can expose a split when the other GCD is trivial.
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


def joint_exponent_product(totals, modulus, budget):
    """Bit-plane Horner product with an explicitly bounded accumulator."""
    powers = tuple(
        (p, exponent // 2) for p, exponent in sorted(totals.items())
    )
    width = max((exponent.bit_length() for _, exponent in powers), default=0)
    cap, result = max(256, 2 * modulus.bit_length()), 1
    for bit in range(width - 1, -1, -1):
        budget.consume(modulus.bit_length() ** 2 + len(powers))
        result = result * result % modulus
        product = 1
        for prime, exponent in powers:
            if (exponent >> bit) & 1:
                budget.consume(product.bit_length() * prime.bit_length())
                if product.bit_length() + prime.bit_length() > cap:
                    result = result * product % modulus
                    product = 1
                product *= prime
        result = result * product % modulus
    return result


class DependencyExtractor:
    """Try finite dependencies, retaining position after trivial GCDs."""

    def __init__(
        self,
        prepared,
        dependencies,
        *,
        budget=None,
        tested_cache=None,
        max_dependencies=None,
        product_policy="reference",
    ):
        """Honor an optional trial-count limit; validate each mask exactly."""
        if not isinstance(dependencies, tuple):
            raise TypeError("dependencies must be an immutable tuple")
        if max_dependencies is not None and len(dependencies) > storage_limit(
            max_dependencies
        ):
            raise ValueError("dependency trial cap exceeded")
        self.prepared, self.dependencies = prepared, dependencies
        self.budget = budget if budget is not None else Budget()
        self.next_dependency = 0
        self.trials = []
        self.tested_cache = tested_cache
        self.cache_skips = 0
        self.product_policy = product_policy

    def run(self):
        """Return a proper divisor or None; refusal preserves the trial."""
        while self.next_dependency < len(self.dependencies):
            key = None
            if self.tested_cache is not None:
                key = self.tested_cache.dependency_key(
                    self.prepared,
                    self.dependencies[self.next_dependency],
                    self.budget,
                )
                if key in self.tested_cache.tested_dependencies:
                    self.cache_skips += 1
                    self.next_dependency += 1
                    continue

            result = extract_dependency(
                self.prepared,
                self.dependencies[self.next_dependency],
                budget=self.budget,
                product_policy=self.product_policy,
            )
            if key is not None and result.divisor is None:
                self.tested_cache.remember_dependency(key)
            self.trials.append(result)
            self.next_dependency += 1
            if result.divisor is not None:
                return result.divisor

        return None
