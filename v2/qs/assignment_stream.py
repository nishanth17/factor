"""Distinct SIQS assignments from a bounded pool and an integer cursor."""

import random
from bisect import bisect_left
from math import comb, gcd, prod

from .. import utils
from .families import MAX_A_FACTORS, _checksum, _identity
from .polynomial import a_target


def unrank_combination(size, count, rank):
    """Unrank a lexicographic subset without enumerating its predecessors."""
    if not 0 <= count <= size or not 0 <= rank < comb(size, count):
        raise ValueError("combination rank is out of range")
    selected, start = [], 0
    for left in range(count, 0, -1):
        for index in range(start, size - left + 1):
            block = comb(size - index - 1, left - 1)
            if rank < block:
                selected.append(index)
                start = index + 1
                break
            rank -= block
    return tuple(selected)


class AssignmentStream:
    """Indexable unique A sets; storage does not depend on the search quota.

    The flyer is chosen outside the core pool, so removing it recovers the
    unique core subset. This proves uniqueness without a growing seen set.
    Count is excluded from identity: increasing it preserves every prefix.
    """

    def __init__(
        self,
        base,
        half_width,
        *,
        factor_count,
        family_count,
        pool_size,
        seed,
        policy,
        budget,
        memory_bytes,
    ):
        for name, value, low, high in (
            ("factor_count", factor_count, 1, MAX_A_FACTORS),
            ("family_count", family_count, 1, 2**31 - 1),
            ("pool_size", pool_size, factor_count, 4096),
            ("seed", seed, 0, 2**64 - 1),
        ):
            utils.require_integer(value, name, low)
            if value > high:
                raise ValueError(f"{name} exceeds the stream limit")
        if policy not in ("nearest", "flyer"):
            raise ValueError("unknown streaming A policy")
        self.workspace_bytes = (
            base.workspace_bytes + 65536 + 512 * len(base.entries)
        )
        if self.workspace_bytes > memory_bytes:
            raise MemoryError("assignment stream exceeds memory_bytes")
        self.budget = budget
        self.target = a_target(base.n_prime, half_width)
        budget.consume(len(base.entries) * (factor_count + 1))
        eligible = tuple(p for p in base.primes if p != 2 and base.n_prime % p)
        self.pool = tuple(
            sorted(
                eligible, key=lambda p: (abs(p**factor_count - self.target), p)
            )[:pool_size]
        )
        self.core_count = factor_count - (
            policy == "flyer" and factor_count > 1
        )
        self.flyers = ()
        if self.core_count < factor_count:
            core = set(self.pool)
            self.flyers = tuple(p for p in eligible if p not in core)
            if not self.flyers:
                raise ValueError(
                    "flyer policy needs primes outside the core pool"
                )
        if len(self.pool) < self.core_count:
            raise ValueError("too few nonsingular primes for the requested A")
        self.total = comb(len(self.pool), self.core_count)
        self.count = min(family_count, self.total)
        generator = random.Random(seed)
        # Keep rank zero first, then permute all remaining ranks bijectively.
        modulus = max(1, self.total - 1)
        self.offset = generator.randrange(modulus)
        self.stride = 1
        for _ in range(128):
            budget.consume(1)
            candidate = generator.randrange(1, modulus + 1)
            if gcd(candidate, modulus) == 1:
                self.stride = candidate
                break
        self.identity = _checksum(
            [
                "assignment-stream-v1",
                _identity(base),
                self.target,
                factor_count,
                pool_size,
                seed,
                policy,
            ]
        )

    def __len__(self):
        return self.count

    def __getitem__(self, index):
        utils.require_integer(index, "assignment index", 0)
        if index >= self.count:
            raise IndexError("assignment cursor is exhausted")
        self.budget.consume(len(self.pool) * (self.core_count + 1))
        rank = (
            0
            if index == 0
            else 1
            + ((index - 1) * self.stride + self.offset) % (self.total - 1)
        )
        values = [
            self.pool[i]
            for i in unrank_combination(len(self.pool), self.core_count, rank)
        ]
        if self.flyers:
            product = prod(values)
            position = bisect_left(
                self.flyers, (self.target + product - 1) // product
            )
            lower, upper = max(0, position - 1), position + 1
            choices = self.flyers[lower:upper]
            values.append(
                min(choices, key=lambda p: (abs(product * p - self.target), p))
            )
        return tuple(sorted(values))


def assignment_identity(assignments):
    """Keep legacy checkpoint fingerprints and compact stream identities."""
    if isinstance(assignments, AssignmentStream):
        return assignments.identity
    if isinstance(assignments, range):
        return _checksum(["external-squares-v1"])
    return _checksum(assignments)
