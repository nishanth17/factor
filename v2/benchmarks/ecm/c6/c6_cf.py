"""Bounded CF records and the paper's distinct three-point interpreter."""

import json
from dataclasses import dataclass
from math import gcd
from pathlib import Path

from ....common import prime_sieve, utils
from ....ecm import prac
from ...support.paths import (
    BENCHMARK_ROOT,
)
from . import c6_chains, c6_fast

DATA = BENCHMARK_ROOT / "inputs/controls/c6_cf_records.json"


def decode(values):
    """Decode upstream integer chains independently into differential SSA."""
    if type(values) is not list or not 2 <= len(values) <= 21:
        raise ValueError("CF chain exceeds finite depth")
    if any(type(v) is not int for v in values):
        raise ValueError("CF entries must be integers")
    if values == [1, 2]:
        return prac.Chain(2, ((0, 0, -1),), 1, "lucas"), b""
    if values[:3] != [1, 2, 3] or not 3 <= values[-1] <= 2000:
        raise ValueError("invalid bounded CF chain")
    operations, bits = [(0, 0, -1), (0, 1, 0)], []
    a, b, c = 0, 1, 2
    for index in range(3, len(values)):
        value = values[index]
        if value == values[b] + values[c]:
            operations.append((b, c, a))
            bits.append(0)
            a, b, c = b, c, index
        elif value == values[a] + values[c]:
            operations.append((a, c, b))
            bits.append(1)
            a, b, c = a, c, index
        else:
            raise ValueError("upstream chain is outside the CF family")
    chain = prac.Chain(values[-1], tuple(operations), c, "lucas")
    prac.verify_chain(chain)
    return chain, bytes(bits)


def verify_minimum(scalar, length):
    """Independent reverse-Euclid exhaustion, only for the bounded CF family.

    Each terminal tuple has a<b, a+b=scalar and gcd(a,b)=1. The predecessor
    is unique: (b-a,a) when b<2a, or (a,b-a) when b>2a. Exhaust every terminal
    pair, stopping a path once it cannot beat the supplied length. Reaching
    (1,2) establishes a CF path; no forward-search pruning is trusted here.
    """
    utils.require_integer(scalar, "scalar", 3)
    utils.require_integer(length, "length", 2)
    if scalar > 2000 or length > 20:
        raise ValueError("independent CF exhaustion exceeds bounds")
    best, nodes = length + 1, 0
    for first in range(1, (scalar + 1) // 2):
        a, b = first, scalar - first
        if gcd(a, b) != 1:
            continue
        count = 2
        while (a, b) != (1, 2) and count < best:
            nodes += 1
            if nodes > 40000:
                raise ValueError("CF verification node cap")
            if b < 2 * a:
                a, b = b - a, a
            else:
                a, b = a, b - a
            count += 1
        if (a, b) == (1, 2):
            best = min(best, count)
    if best != length:
        raise ValueError("claimed CF minimum is incorrect")
    return nodes


@dataclass(frozen=True)
class CFRecord:
    certified: c6_fast.CertifiedRecord
    prime: int
    exponent: int
    bits: bytes

    def __post_init__(self):
        if (
            type(self.prime) is not int
            or not 2 <= self.prime <= 2000
            or type(self.exponent) is not int
            or not 1 <= self.exponent <= 10
            or type(self.bits) is not bytes
            or len(self.bits) > 18
            or any(bit > 1 for bit in self.bits)
            or self.prime**self.exponent > 2000
        ):
            raise ValueError("invalid CF metadata bounds")
        values = [1, 2]
        if self.prime > 2:
            values.append(3)
            a, b, c = 1, 2, 3
            for bit in self.bits:
                a, b, c = (a, c, a + c) if bit else (b, c, b + c)
                values.append(c)
        elif self.bits:
            raise ValueError("doubling has no CF bits")
        chain, _ = decode(values)
        if chain.scalar != self.prime:
            raise ValueError("CF metadata scalar mismatch")
        power = self.prime**self.exponent
        expected = c6_chains.compact(c6_chains.compose(chain, power))
        if expected != self.certified.record:
            raise ValueError("three-point action differs from certified code")


class ThreePoint(c6_fast.FastRecord):
    """Compose CF prime actions using three persistent working points."""

    def __init__(self, record, backend):
        super().__init__(record.certified, backend)
        self.prime, self.exponent, self.bits = (
            record.prime,
            record.exponent,
            record.bits,
        )
        self.run = self.three

    def three(self, point, n, a24):
        add, double = self.backend.ecm.point_add, self.backend.ecm.point_double
        factor, masks = c6_fast._factor, self.masks
        product = factor(1, point, masks[0])
        index = 1
        for _ in range(self.exponent):
            a = point
            b = double(*a, n, a24)
            product = factor(product, b, masks[index])
            index += 1
            if self.prime == 2:
                point = b
                continue
            c = add(*a, *b, *a, n)
            product = factor(product, c, masks[index])
            index += 1
            for bit in self.bits:
                if bit:
                    a, b, c = a, c, add(*a, *c, *b, n)
                else:
                    a, b, c = b, c, add(*b, *c, *a, n)
                product = factor(product, c, masks[index])
                index += 1
            point = c
        return point, product


def load_catalog(path=DATA):
    with Path(path).open("rb") as stream:
        raw = stream.read(c6_chains.MAX_FILE_BYTES + 1)
    if len(raw) > c6_chains.MAX_FILE_BYTES:
        raise ValueError("CF catalog exceeds byte cap")
    data = json.loads(raw)
    if data["schema"] != 1 or data["bound"] != 2000:
        raise ValueError("unknown CF catalog")
    rows = data["primes"]
    if len(rows) > c6_chains.MAX_RECORDS:
        raise ValueError("too many CF prime records")
    result = {}
    for text, values in rows.items():
        prime = int(text)
        chain, bits = decode(values)
        if chain.scalar != prime:
            raise ValueError("CF scalar mismatch")
        power, exponent = prime, 1
        while power <= 2000:
            if power not in result and len(result) >= c6_chains.MAX_RECORDS:
                raise ValueError("CF retained-record cap")
            record = c6_chains.compact(c6_chains.compose(chain, power))
            certified = c6_fast.CertifiedRecord(
                record, c6_fast.frontier(record)
            )
            result[power] = CFRecord(certified, prime, exponent, bits)
            power *= prime
            exponent += 1
    if set(result) != {
        p**e
        for p in prime_sieve.prime_sieve(2001)
        for e in range(1, 11)
        if p**e <= 2000
    }:
        raise ValueError("incomplete CF catalog")
    if len(result) > c6_chains.MAX_RECORDS:
        raise ValueError("too many composed CF records")
    return result


def build_program(bound, mode, backend, batch=16, catalog=None):
    utils.require_integer(bound, "bound", 2)
    if bound > 2000 or mode not in ("tuple", "three"):
        raise ValueError("unsupported CF program")
    if catalog is None:
        catalog = load_catalog()
    entries = []
    for prime in prime_sieve.prime_sieve(bound + 1):
        power = utils.prime_power(prime, bound)
        record = catalog[power]
        action = (
            ThreePoint(record, backend)
            if mode == "three"
            else c6_fast.FastRecord(record.certified, backend)
        )
        entries.append((prime, power, action, catalog[prime].certified))
    return c6_fast.Program(entries, backend, batch)
