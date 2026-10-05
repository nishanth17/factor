"""Capped product/remainder trees with complete smooth-part detection."""

from math import gcd

from .. import utils
from ..budget import Budget


def _tree_reservation(values, max_bits, max_nodes, memory_bytes):
    """Check every leaf before allocating the tree or multiplying nodes."""
    if not isinstance(values, tuple):
        raise TypeError("tree leaves must be a tuple")
    for name, value in (
        ("max_bits", max_bits),
        ("max_nodes", max_nodes),
        ("memory_bytes", memory_bytes),
    ):
        utils.require_integer(value, name, 0)

    if max_bits > 1048576 or max_nodes > 8192:
        raise ValueError("tree configuration exceeds the finite hard limits")
    if len(values) > max_nodes:
        raise MemoryError("product tree exceeds its node limit")
    bits = 0
    for value in values:
        utils.require_integer(value, "tree leaf", 1)
        bits += value.bit_length()
        if bits > max_bits:
            raise MemoryError("product tree exceeds its bit limit")

    nodes, width = 0, len(values)
    while width:
        nodes += width
        if width == 1:
            break
        width = (width + 1) // 2

    if bits > max_bits or nodes > max_nodes:
        raise MemoryError("product tree exceeds bit/node limits")
    reserve = 4096 + 192 * nodes
    reserve += 4 * bits * (len(values).bit_length() + 2)
    if reserve > memory_bytes:
        raise MemoryError("product/remainder tree exceeds memory_bytes")
    return reserve


def product_tree(
    values,
    *,
    budget=None,
    max_bits=262144,
    max_nodes=8192,
    memory_bytes=32 * 1024 * 1024,
):
    """Return immutable levels; empty input has no root, one is permitted."""
    _tree_reservation(values, max_bits, max_nodes, memory_bytes)
    budget = budget if budget is not None else Budget()
    budget.consume(len(values) + 1)
    levels, level = [], values

    while level:
        levels.append(level)
        if len(level) == 1:
            break
        next_level = []
        for index in range(0, len(level), 2):
            left = level[index]
            right = level[index + 1] if index + 1 < len(level) else 1
            budget.consume(left.bit_length() + right.bit_length())
            next_level.append(left * right)

        level = tuple(next_level)

    return tuple(levels)


class SmoothBatch:
    """Retain a squarefree base product; detect all valuations at its primes.

    At a positive leaf v, raise z mod v to 2**e >= bit_length(v).
    Each base-prime valuation in v is smaller than this exponent, so
    gcd(v, z**(2**e) mod v) is its complete smooth part. This detects
    smoothness, not its prime exponents: relation recovery remains separate.
    """

    def __init__(
        self,
        primes,
        *,
        budget=None,
        max_bits=262144,
        max_nodes=8192,
        memory_bytes=32 * 1024 * 1024,
    ):
        self.budget = budget if budget is not None else Budget()
        self.max_bits, self.max_nodes = max_bits, max_nodes
        self.memory_bytes = memory_bytes
        if not isinstance(primes, tuple):
            raise TypeError("primes must be a tuple")
        _tree_reservation(primes, max_bits, max_nodes, memory_bytes)
        previous = 1

        for prime in primes:
            self.budget.consume(prime.bit_length() ** 2)
            if prime <= previous or prime > 100000:
                raise ValueError("base primes must increase and be <= 100000")
            if utils.classify_prime(prime) is not utils.Primality.PROVEN:
                raise ValueError("smoothness base contains a composite")
            previous = prime

        tree = product_tree(
            primes,
            budget=self.budget,
            max_bits=max_bits,
            max_nodes=max_nodes,
            memory_bytes=memory_bytes,
        )
        self.radical = tree[-1][0] if tree else 1

    def residuals(self, values):
        """Return exact non-base parts in input order, including duplicates."""
        reserve = _tree_reservation(
            values, self.max_bits, self.max_nodes, self.memory_bytes
        )
        if reserve + 8 * self.radical.bit_length() > self.memory_bytes:
            raise MemoryError("radical and tree coexistence exceeds cap")
        tree = product_tree(
            values,
            budget=self.budget,
            max_bits=self.max_bits,
            max_nodes=self.max_nodes,
            memory_bytes=self.memory_bytes,
        )
        if not tree:
            return ()
        self.budget.consume(
            self.radical.bit_length() + tree[-1][0].bit_length()
        )
        remainders = (self.radical % tree[-1][0],)

        # Each child divides its parent: reduce the shared remainder down
        # the tree instead of dividing the full radical at every leaf.
        for level in reversed(tree[:-1]):
            children = []
            for index, modulus in enumerate(level):
                parent = remainders[index // 2]
                self.budget.consume(parent.bit_length() + modulus.bit_length())
                children.append(parent % modulus)
            remainders = tuple(children)

        output = []

        for value, remainder in zip(values, remainders):
            for _ in range((value.bit_length() - 1).bit_length()):
                self.budget.consume(2 * value.bit_length())
                remainder = remainder * remainder % value
            self.budget.consume(value.bit_length())
            output.append(value // gcd(value, remainder))

        return tuple(output)
