"""Capped product/remainder trees with complete smooth-part detection."""

from ..common import arithmetic, utils
from ..common.arithmetic import gcd
from ..execution.resources import removed_arguments, resolve_budget


@removed_arguments(("max_bits", "max_nodes", "memory_bytes"))
def _tree_reservation(values):
    """Check every leaf before allocating the tree or multiplying nodes."""
    if not isinstance(values, tuple):
        raise TypeError("tree leaves must be a tuple")

    bits = 0
    for value in values:
        utils.require_integer(value, "tree leaf", 1)
        bits += value.bit_length()

    nodes, width = 0, len(values)
    while width:
        nodes += width
        if width == 1:
            break
        width = (width + 1) // 2

    reserve = 4096 + 192 * nodes
    reserve += 4 * bits * (len(values).bit_length() + 2)
    return reserve


@removed_arguments(("max_bits", "max_nodes", "memory_bytes"))
def product_tree(values, *, budget=None):
    """Return immutable levels; empty input has no root, one is permitted."""
    _tree_reservation(values)
    budget = resolve_budget(budget)
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

    @removed_arguments(("max_bits", "max_nodes", "memory_bytes"))
    def __init__(self, primes, *, backend="python-int", budget=None):
        self.budget = resolve_budget(budget)
        if not isinstance(primes, tuple):
            raise TypeError("primes must be a tuple")
        _tree_reservation(primes)
        previous = 1

        for prime in primes:
            self.budget.consume(prime.bit_length() ** 2)
            if prime <= previous or prime > 100000:
                raise ValueError("base primes must increase and be <= 100000")
            if utils.classify_prime(prime) is not utils.Primality.PROVEN:
                raise ValueError("smoothness base contains a composite")
            previous = prime

        engine = arithmetic.get_backend(backend)
        tree = product_tree(
            tuple(engine.integer(prime) for prime in primes),
            budget=self.budget,
        )
        self.radical = tree[-1][0] if tree else 1

    def residuals(self, values):
        """Return exact non-base parts in input order, including duplicates."""
        _tree_reservation(values)
        tree = product_tree(
            values,
            budget=self.budget,
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
            output.append(arithmetic.divexact(value, gcd(value, remainder)))

        return tuple(output)
