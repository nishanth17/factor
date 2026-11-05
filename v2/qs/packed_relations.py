"""Immutable packed exponent records, with lossless checked interchange."""

import base64
import binascii
import struct
from dataclasses import FrozenInstanceError

from .relations import AtomicRelation, CombinedRelation, _check_exponents

PAIR = struct.Struct(">II")


def encode_exponents(exponents):
    """Native factor-base primes and bounded norm valuations fit u32 pairs."""
    _check_exponents(exponents)
    if any(p >= 2**32 or e >= 2**32 for p, e in exponents):
        raise ValueError("packed exponent exceeds its u32 domain")
    return b"".join(PAIR.pack(p, e) for p, e in exponents)


def decode_exponents(data, maximum):
    """Check shape, ordering and positive entries before decoding."""
    if (
        type(data) is not bytes
        or len(data) % PAIR.size
        or (len(data) // PAIR.size > maximum)
    ):
        raise ValueError("packed exponent record exceeds its shape bound")
    exponents = tuple(PAIR.iter_unpack(data))
    _check_exponents(exponents)
    return exponents


class _PackedPayload:
    """No decoded tuple is retained; caller snapshots share immutable bytes."""

    @property
    def exponents(self):
        return tuple(PAIR.iter_unpack(self._packed_exponents))

    def __setattr__(self, name, value):
        raise FrozenInstanceError("packed relations are immutable")

    def __delattr__(self, name):
        raise FrozenInstanceError("packed relations are immutable")


class PackedAtomicRelation(_PackedPayload, AtomicRelation):
    def __init__(
        self,
        polynomial,
        position,
        sign,
        exponents,
        residual=1,
        *,
        large_primes=(),
    ):
        original = AtomicRelation(
            polynomial,
            position,
            sign,
            exponents,
            residual,
            large_primes=large_primes,
        )
        for name in (
            "polynomial",
            "position",
            "sign",
            "residual",
            "large_primes",
            "relation_id",
        ):
            object.__setattr__(self, name, getattr(original, name))
        object.__setattr__(
            self, "_packed_exponents", encode_exponents(exponents)
        )


class PackedCombinedRelation(_PackedPayload, CombinedRelation):
    def __init__(self, atom_ids, u, sign, exponents, square_correction):
        original = CombinedRelation(
            atom_ids, u, sign, exponents, square_correction
        )
        for name in ("atom_ids", "u", "sign", "square_correction"):
            object.__setattr__(self, name, getattr(original, name))
        object.__setattr__(
            self, "_packed_exponents", encode_exponents(exponents)
        )


def pack_relation(relation, budget):
    """Charge encoding before publishing the original identity."""
    if isinstance(relation, _PackedPayload):
        return relation
    exponents = relation.exponents
    budget.consume(2 * len(exponents) + 1)
    if isinstance(relation, AtomicRelation):
        return PackedAtomicRelation(
            relation.polynomial,
            relation.position,
            relation.sign,
            exponents,
            relation.residual,
            large_primes=relation.large_primes,
        )
    return PackedCombinedRelation(
        relation.atom_ids,
        relation.u,
        relation.sign,
        exponents,
        relation.square_correction,
    )


def export_exponents(relation):
    data = (
        relation._packed_exponents
        if isinstance(relation, _PackedPayload)
        else encode_exponents(relation.exponents)
    )
    return {"u32be64": base64.b64encode(data).decode("ascii")}


def import_exponents(record, maximum, budget):
    """Bound and check canonical base64 before decoding."""
    if not isinstance(record, dict) or set(record) != {"u32be64"}:
        raise ValueError("invalid packed exponent encoding")
    text = record["u32be64"]
    if not isinstance(text, str) or len(text) > 4 * ((8 * maximum + 2) // 3):
        raise ValueError("packed exponent text exceeds its bound")
    budget.consume(len(text) + 1)
    try:
        data = base64.b64decode(text, validate=True)
    except (ValueError, binascii.Error) as error:
        raise ValueError("invalid packed base64") from error
    if base64.b64encode(data).decode("ascii") != text:
        raise ValueError("noncanonical packed base64")
    return [list(pair) for pair in decode_exponents(data, maximum)]
