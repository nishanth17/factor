"""Bounded flat provenance with deferred modular-root arithmetic.

The DAG has depth one: checked immutable original atoms are its only leaves.
Factor/parity metadata is eager; roots and square corrections are constructed
only for extraction or interchange. No recursively expanding merge histories.
"""

from collections import Counter
from dataclasses import FrozenInstanceError
from math import isqrt

from ..common.arithmetic import pow
from .factor_witness import factor_snapshot
from .relations import MAX_COMBINED_ATOMS, CombinedRelation


def structural_values(atoms, budget):
    factors, residuals = Counter(), Counter()
    sign = 1
    for atom in atoms:
        budget.consume(len(atom.exponents) + len(atom.residual_primes) + 1)
        sign *= atom.sign
        for prime, exponent in atom.exponents:
            factors[prime] += exponent
        residuals.update(atom.residual_primes)
    if any(getattr(atom, "opaque", None) is not None for atom in atoms):
        product = 1
        for atom in atoms:
            budget.consume(
                (product.bit_length() + 1) * atom.residual_product.bit_length()
            )
            product *= atom.residual_product
        budget.consume(product.bit_length() ** 2)
        root = isqrt(product)
        if root * root != product:
            raise ValueError("opaque cycle product is not square")
        residuals = {root: 2}
    if any(count % 2 for count in residuals.values()):
        raise ValueError("lazy cycle has an odd residual multiplicity")
    return (
        sign,
        tuple(sorted(factors.items())),
        tuple(sorted(residuals.items())),
    )


class LazyCombinedRelation(CombinedRelation):
    """Immutable leaves with separately charged arithmetic materialization."""

    def __init__(self, atoms, modulus, budget, *, factors=None):
        identities = tuple(atom.relation_id for atom in atoms)
        if not 1 <= len(atoms) <= MAX_COMBINED_ATOMS or (
            len(set(identities)) != len(atoms)
        ):
            raise ValueError("invalid lazy provenance")
        sign, exponents, residuals = structural_values(atoms, budget)
        for name, value in (
            ("atom_ids", identities),
            ("sign", sign),
            ("exponents", exponents if factors is None else factors),
            ("_leaves", tuple(atoms)),
            ("_modulus", modulus),
            ("_residuals", residuals),
            ("_budget", budget),
        ):
            object.__setattr__(self, name, value)

    def __setattr__(self, name, value):
        raise FrozenInstanceError("lazy relations are immutable")

    def __delattr__(self, name):
        raise FrozenInstanceError("lazy relations are immutable")

    def bind_budget(self, budget):
        # Continuation changes the allowance, never the mathematical payload.
        object.__setattr__(self, "_budget", budget)

    def _with_storage(self, atoms, factors, budget):
        """Rebind checked structural metadata to its private packed owners.

        Arena encoding changes representation only. Recomputing all factor
        sums and the opaque product square root here would repeat admission
        work. Public/restore verification still checks original arithmetic.
        """
        atoms = tuple(atoms)
        if tuple(atom.relation_id for atom in atoms) != self.atom_ids:
            raise ValueError("packed lazy provenance changed")
        budget.consume(len(atoms) + 1)
        result = object.__new__(type(self))
        for name, value in (
            ("atom_ids", self.atom_ids),
            ("sign", self.sign),
            ("exponents", factors),
            ("_leaves", atoms),
            ("_modulus", self._modulus),
            ("_residuals", self._residuals),
            ("_budget", budget),
        ):
            object.__setattr__(result, name, value)
        return result

    def materialize(self):
        bits = self._modulus.bit_length()
        self._budget.consume(
            (3 * len(self._leaves) + 2 * len(self._residuals) + 1) * bits**2
        )
        root, correction = 1, 1
        for atom in self._leaves:
            root = root * atom.u % self._modulus
            correction = correction * atom.square_correction % self._modulus
        for residual, count in self._residuals:
            correction = (
                correction
                * pow(residual, count // 2, self._modulus)
                % self._modulus
            )
        # The consumer uses both results once. No mutable arithmetic cache
        # can masquerade as a proof after an unsupported payload mutation.
        return root, correction

    @property
    def u(self):
        return self.materialize()[0]

    @property
    def square_correction(self):
        return self.materialize()[1]

    def structural_snapshot(self):
        return (
            self.atom_ids,
            self.sign,
            factor_snapshot(self.exponents),
            self._modulus,
            self._residuals,
            tuple(id(a) for a in self._leaves),
        )

    def verify_structure(self, atoms, budget):
        if len(atoms) != len(self._leaves) or any(
            atom is not leaf for atom, leaf in zip(atoms, self._leaves)
        ):
            raise ValueError("lazy provenance ownership changed")
        sign, factors, residuals = structural_values(atoms, budget)
        if (sign, factors, residuals) != (
            self.sign,
            tuple(self.exponents),
            self._residuals,
        ):
            raise ValueError("lazy structural payload changed")
