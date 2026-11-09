"""Base-bound compact factors and integer provenance in append-only slabs.

Published slab ranges are never overwritten. Records own their slab, so a
snapshot remains valid after collector eviction or replacement of the arena.
Only the current slab is additionally pinned by the allocator. Conservative
collector reservations cover slabs and construction copies separately.
"""

import base64
import binascii
import struct
from collections.abc import Sequence
from dataclasses import FrozenInstanceError

from .relations import AtomicRelation, CombinedRelation, _check_exponents
from .residual_relations import ResidualAtom

BLOCK_BYTES = 4096
BLOCK_OVERHEAD = 256
HANDLE = struct.Struct(">Q")


class IndexedFactors(Sequence):
    """Immutable recovery output retaining factor-base columns directly."""

    __slots__ = ("indices", "primes")

    def __init__(self, indices, primes):
        if not isinstance(primes, tuple):
            raise TypeError("indexed recovery requires immutable base primes")
        indices = tuple((index, exponent) for index, exponent in indices)
        previous = 0
        for index, exponent in indices:
            if type(index) is not int or not previous < index <= len(primes):
                raise ValueError("recovered indices must increase in the base")
            if type(exponent) is not int or not 1 <= exponent < 2**32:
                raise ValueError("invalid recovered exponent")
            previous = index
        object.__setattr__(self, "indices", indices)
        object.__setattr__(self, "primes", primes)

    def __setattr__(self, name, value):
        raise FrozenInstanceError("indexed recovery is immutable")

    def __delattr__(self, name):
        raise FrozenInstanceError("indexed recovery is immutable")

    def __len__(self):
        return len(self.indices)

    def __iter__(self):
        return ((self.primes[i - 1], e) for i, e in self.indices)

    def __getitem__(self, index):
        return tuple(self)[index]

    def __eq__(self, other):
        return tuple(self) == tuple(other)

    def __hash__(self):
        return hash(tuple(self))


def _uint(value):
    if not 0 <= value < 2**32:
        raise ValueError("compact integer exceeds u32")
    result = bytearray()
    while value >= 128:
        result.append((value & 127) | 128)
        value >>= 7
    result.append(value)
    return result


def _read_uint(data, offset, end):
    value, shift = 0, 0
    while offset < end and shift <= 28:
        part = data[offset]
        offset += 1
        value |= (part & 127) << shift
        if part < 128:
            if value >= 2**32 or (shift and part == 0):
                raise ValueError("noncanonical compact integer")
            return value, offset
        shift += 7
    raise ValueError("truncated or oversized compact integer")


def encode_factors(exponents, base):
    """Delta-code base indices; exponent one has no additional field."""
    if isinstance(exponents, IndexedFactors):
        if exponents.primes is not base.primes:
            raise ValueError("recovery indices belong to a different base")
        pairs = exponents.indices
    else:
        exponents = tuple(exponents)
        _check_exponents(exponents)
        pairs = tuple((base._columns.get(p), e) for p, e in exponents)
    data, previous = bytearray(), 0
    for index, exponent in pairs:
        if index is None:
            raise ValueError("compact prime is outside its factor base")
        data.extend(_uint(2 * (index - previous) + (exponent != 1)))
        if exponent != 1:
            data.extend(_uint(exponent))
        previous = index
    return bytes(data)


def decode_factors(data, primes, count):
    """Validate index, exponent, count and canonical termination exactly."""
    if type(count) is not int or not 0 <= count <= len(primes):
        raise ValueError("compact factor count exceeds its base")
    offset, index, end = 0, 0, len(data)
    for _ in range(count):
        token, offset = _read_uint(data, offset, end)
        delta, explicit = divmod(token, 2)
        index += delta
        if not delta or index > len(primes):
            raise ValueError("compact indices are not increasing in the base")
        exponent = 1
        if explicit:
            exponent, offset = _read_uint(data, offset, end)
            if exponent <= 1:
                raise ValueError("noncanonical explicit exponent")
        yield primes[index - 1], exponent
    if offset != end:
        raise ValueError("trailing compact factor data")


class FactorRecord(Sequence):
    """Immutable sequence interface without repeated tuple allocation."""

    __slots__ = ("_block", "_offset", "_length", "_count", "_primes", "_hash")

    def __init__(self, block, offset, length, count, primes):
        for name, value in zip(
            self.__slots__, (block, offset, length, count, primes)
        ):
            object.__setattr__(self, name, value)
        object.__setattr__(self, "_hash", None)

    def __setattr__(self, name, value):
        raise FrozenInstanceError("compact factors are immutable")

    def __delattr__(self, name):
        raise FrozenInstanceError("compact factors are immutable")

    def encoded(self):
        return bytes(
            memoryview(self._block)[self._offset : self._offset + self._length]
        )

    def __len__(self):
        return self._count

    def __iter__(self):
        # A readonly view avoids copying on ordinary iteration. No writable
        # view is exposed; bounds and canonical shape are checked by decoding.
        data = memoryview(self._block)[
            self._offset : self._offset + self._length
        ].toreadonly()
        return decode_factors(data, self._primes, self._count)

    def __getitem__(self, index):
        return tuple(self)[index]

    def __eq__(self, other):
        if self is other:
            return True
        if isinstance(other, FactorRecord) and self._primes is other._primes:
            return (
                self._count == other._count
                and self.encoded() == other.encoded()
            )
        if isinstance(other, (tuple, FactorRecord)):
            return tuple(self) == tuple(other)
        return NotImplemented

    def __hash__(self):
        # Retain one scalar, not a decoded duplicate of every sparse payload.
        # The private slab contract is immutable; proof reuse separately
        # compares bytes and catches unsupported mutation of that storage.
        if self._hash is None:
            object.__setattr__(self, "_hash", hash(tuple(self)))
        return self._hash


class RecordArena:
    """Slabs for factors and stable integer handles for owned provenance."""

    def __init__(self, base, *, config=None, budget=None):
        self.base = base
        self.block, self.used = bytearray(), 0
        self.allocated_bytes = 0
        self.shared_contexts = (
            getattr(config, "extension_policy", None) == "extended-v2"
        )
        self.polynomials = {}
        self.slots, self.atoms = {}, {}
        self.next_slot = 0
        self.recompute_limit = getattr(config, "recompute_prime_limit", 0)
        self.spill = None
        if getattr(config, "spill_directory", None) is not None:
            from .cold_factors import SpillStore

            self.spill = SpillStore(
                config.spill_directory,
                config.spill_cache_bytes,
                budget,
            )

    def factors(self, exponents, budget):
        budget.consume(4 * len(exponents) + 1)
        data = encode_factors(exponents, self.base)
        if self.spill is not None:
            from .cold_factors import ColdFactors

            self.spill.budget = budget
            reference = self.spill.append(data)
            return ColdFactors(
                self.spill, reference, self.base.primes, len(exponents)
            )
        if len(data) > len(self.block) - self.used:
            self.block = bytearray(max(BLOCK_BYTES, len(data)))
            self.allocated_bytes += BLOCK_OVERHEAD + len(self.block)
            self.used = 0
        offset = self.used
        self.block[offset : offset + len(data)] = data
        self.used += len(data)
        return FactorRecord(
            self.block, offset, len(data), len(exponents), self.base.primes
        )

    def plan(self, atom, row, budget):
        """Allocate unpublished records; no handle becomes visible yet."""
        if atom.relation_id in self.slots or self.next_slot >= 2**64:
            raise ValueError("duplicate or exhausted compact handle")
        from .residual_relations import ResidualAtom

        cls = (
            IndexedResidualAtom
            if isinstance(atom, ResidualAtom)
            else IndexedAtomicRelation
        )
        polynomial = self.polynomial(atom.polynomial)
        factors = atom.exponents
        small = (
            tuple(p for p, _ in factors if p <= self.recompute_limit)
            if self.recompute_limit
            else ()
        )
        if small:
            from .cold_factors import RecomputedFactors

            remaining = (
                IndexedFactors(
                    tuple(
                        (i, e)
                        for i, e in factors.indices
                        if self.base.primes[i - 1] > self.recompute_limit
                    ),
                    self.base.primes,
                )
                if isinstance(factors, IndexedFactors)
                else tuple(
                    (p, e) for p, e in factors if p > self.recompute_limit
                )
            )
            retained = self.factors(remaining, budget)
            stored = RecomputedFactors(
                polynomial, atom.position, small, retained, len(factors)
            )
        else:
            stored = self.factors(factors, budget)
        packed_atom = cls(atom, stored)
        object.__setattr__(packed_atom, "polynomial", polynomial)
        slot = self.next_slot
        packed_row = None
        if row is not None:
            from .lazy_relations import LazyCombinedRelation

            ids = tuple(
                slot if identity == atom.relation_id else self.slots[identity]
                for identity in row.atom_ids
            )
            budget.consume(len(ids))
            if isinstance(row, LazyCombinedRelation):
                leaves = tuple(
                    packed_atom if i == slot else self.atoms[i] for i in ids
                )
                packed_row = row._with_storage(
                    leaves,
                    self.factors(row.exponents, budget),
                    budget,
                )
            else:
                packed_row = IndexedCombinedRelation(
                    row,
                    self.factors(row.exponents, budget),
                    self,
                    b"".join(HANDLE.pack(i) for i in ids),
                )
        return packed_atom, packed_row, slot

    def commit(self, atom, slot):
        if slot != self.next_slot or atom.relation_id in self.slots:
            raise ValueError("stale compact publication")
        if self.shared_contexts:
            self.polynomials[atom.polynomial.identity] = atom.polynomial
        self.slots[atom.relation_id] = slot
        self.atoms[slot] = atom
        self.next_slot += 1

    def polynomial(self, polynomial):
        """Reuse a context only after checking all polynomial fields."""
        if not self.shared_contexts:
            return polynomial
        previous = self.polynomials.get(polynomial.identity)
        if previous is not None and previous != polynomial:
            raise ValueError(
                "polynomial identity has conflicting coefficients"
            )
        return polynomial if previous is None else previous

    def allocation_bound(self, counts):
        """Preflight worst-case new slabs before any encoding allocation."""
        if self.spill is not None:
            return 0
        # Each index/exponent occupies at most two canonical u32 varints.
        return sum(BLOCK_OVERHEAD + max(BLOCK_BYTES, 10 * n) for n in counts)

    def row(self, row, budget):
        from .lazy_relations import LazyCombinedRelation

        if isinstance(row, LazyCombinedRelation):
            return row._with_storage(
                tuple(self.atoms[self.slots[i]] for i in row.atom_ids),
                self.factors(row.exponents, budget),
                budget,
            )
        budget.consume(len(row.atom_ids))
        handles = b"".join(HANDLE.pack(self.slots[i]) for i in row.atom_ids)
        return IndexedCombinedRelation(
            row, self.factors(row.exponents, budget), self, handles
        )

    def forget(self, identity):
        slot = self.slots.pop(identity)
        del self.atoms[slot]


class _IndexedPayload:
    @property
    def exponents(self):
        return self._factors

    def __setattr__(self, name, value):
        raise FrozenInstanceError("compact relations are immutable")

    def __delattr__(self, name):
        raise FrozenInstanceError("compact relations are immutable")


class IndexedAtomicRelation(_IndexedPayload, AtomicRelation):
    def __init__(self, source, factors):
        for name in (
            "polynomial",
            "position",
            "sign",
            "residual",
            "large_primes",
            "relation_id",
        ):
            object.__setattr__(self, name, getattr(source, name))
        object.__setattr__(self, "_factors", factors)


class IndexedResidualAtom(_IndexedPayload, ResidualAtom):
    def __init__(self, source, factors):
        IndexedAtomicRelation.__init__(self, source, factors)
        object.__setattr__(self, "opaque", source.opaque)


class IndexedCombinedRelation(_IndexedPayload, CombinedRelation):
    def __init__(self, source, factors, arena, handles):
        for name in ("u", "sign", "square_correction"):
            object.__setattr__(self, name, getattr(source, name))
        object.__setattr__(self, "_factors", factors)
        object.__setattr__(self, "_arena", arena)
        object.__setattr__(self, "_handles", handles)

    @property
    def atom_ids(self):
        return tuple(
            self._arena.atoms[slot].relation_id
            for (slot,) in HANDLE.iter_unpack(self._handles)
        )


def export_factors(row, base):
    data = (
        row.exponents.encoded()
        if isinstance(row.exponents, FactorRecord)
        else encode_factors(row.exponents, base)
    )
    return {
        "index-v2": base64.b64encode(data).decode("ascii"),
        "count": len(row.exponents),
    }


def import_factors(record, base, budget):
    if not isinstance(record, dict) or set(record) != {"index-v2", "count"}:
        raise ValueError("invalid compact factor encoding")
    count, text = record["count"], record["index-v2"]
    if type(count) is not int or not 0 <= count <= len(base.entries):
        raise ValueError("compact factor count exceeds its base")
    if not isinstance(text, str) or len(text) > 4 * ((10 * count + 2) // 3):
        raise ValueError("compact factor encoding exceeds its bound")
    budget.consume(len(text) + 1)
    try:
        data = base64.b64decode(text, validate=True)
    except (ValueError, binascii.Error) as error:
        raise ValueError("invalid compact base64") from error
    if base64.b64encode(data).decode("ascii") != text:
        raise ValueError("noncanonical compact base64")
    return [list(pair) for pair in decode_factors(data, base.primes, count)]
