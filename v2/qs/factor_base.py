"""Bounded factor bases and exact modular square roots for reference QS."""

from dataclasses import dataclass, field
from types import MappingProxyType

from .. import arithmetic, prime_sieve, utils
from ..arithmetic import gcd, pow
from ..budget import Budget

MAX_INPUT_BITS = 4096
MAX_MULTIPLIER = 1_000_000
MAX_FACTOR_BASE_BOUND = 1_000_000
DEFAULT_MEMORY_BYTES = 8 * 1024 * 1024


def checked_target(n, multiplier):
    """Return h*n for an odd n >= 3 and a bounded positive multiplier.

    Primes and powers are permitted for diagnostic fixtures; production
    callers should preprocess them. This does not classify n as composite.
    """
    utils.require_integer(n, "n", 3)
    utils.require_integer(multiplier, "multiplier", 1)
    if n % 2 == 0 or n.bit_length() > MAX_INPUT_BITS:
        raise ValueError("n must be odd and at most 4096 bits")
    if multiplier > MAX_MULTIPLIER:
        raise ValueError("multiplier exceeds the reference limit")
    return n * multiplier


def modular_square_roots(value, prime, *, budget=None):
    """Return sorted, distinct roots of value modulo a proven small prime.

    An empty tuple denotes a nonresidue. Tonelli-Shanks uses only integers;
    composite or oversized moduli raise ValueError. Work is cooperative and
    bounded by MAX_FACTOR_BASE_BOUND, including the nonresidue search.
    """
    utils.require_integer(value, "value")
    if abs(value).bit_length() > 4 * MAX_INPUT_BITS + 64:
        raise ValueError("square-root input exceeds the reference bit limit")
    utils.require_integer(prime, "prime", 2)
    if prime > MAX_FACTOR_BASE_BOUND:
        raise ValueError("prime exceeds the reference factor-base limit")
    budget = budget if budget is not None else Budget()
    budget.consume(prime.bit_length() ** 2)
    if utils.classify_prime(prime) != utils.Primality.PROVEN:
        raise ValueError("modulus must be a proven prime")
    value %= prime
    if prime == 2 or value == 0:
        return (value,)
    if pow(value, (prime - 1) // 2, prime) != 1:
        return ()
    if prime % 4 == 3:
        root = pow(value, (prime + 1) // 4, prime)
        return tuple(sorted((root, prime - root)))

    odd_part, shifts = prime - 1, 0
    while odd_part % 2 == 0:
        odd_part //= 2
        shifts += 1
    nonresidue = 2
    while pow(nonresidue, (prime - 1) // 2, prime) != prime - 1:
        budget.consume()
        nonresidue += 1
    root = pow(value, (odd_part + 1) // 2, prime)
    remainder = pow(value, odd_part, prime)
    correction = pow(nonresidue, odd_part, prime)

    # root² == value*remainder (mod prime); each correction preserves that
    # identity while lowering the power-of-two order of the remainder.
    while remainder != 1:
        budget.consume(shifts)
        index, squared = 0, remainder
        while squared != 1 and index < shifts:
            squared = squared * squared % prime
            index += 1
        if index == shifts:
            raise ArithmeticError("Tonelli-Shanks invariant failed")
        step = pow(correction, 1 << (shifts - index - 1), prime)
        root = root * step % prime
        correction = step * step % prime
        remainder = remainder * correction % prime
        shifts = index

    return tuple(sorted((root, prime - root)))


@dataclass(frozen=True)
class FactorBaseEntry:
    """A prime and the complete square-root set for the factor-base target."""

    prime: int
    square_roots: tuple[int, ...]


@dataclass(frozen=True)
class FactorBase:
    """Immutable primes below bound; sign occupies parity column zero.

    Entries include 2 and primes dividing h*n. workspace_bytes is a
    conservative owned-storage reservation, not a process RSS guarantee.
    """

    n: int
    multiplier: int
    bound: int
    entries: tuple[FactorBaseEntry, ...]
    _primes: tuple = field(init=False, repr=False, compare=False)
    _family_identity: str | None = field(
        default=None, init=False, repr=False, compare=False
    )
    _columns: object = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        """Reject mutable, unordered, composite, or incorrect root data."""
        target = checked_target(self.n, self.multiplier)
        if gcd(self.n, self.multiplier) != 1:
            raise ValueError("factor-base multiplier must be coprime to n")
        utils.require_integer(self.bound, "bound", 3)
        if self.bound > MAX_FACTOR_BASE_BOUND:
            raise ValueError("factor-base bound exceeds the reference limit")
        if not isinstance(self.entries, tuple):
            raise TypeError("entries must be a tuple")
        if not self.entries or self.entries[0].prime != 2:
            raise ValueError("factor base must include 2")
        previous = 1

        for entry in self.entries:
            prime = entry.prime
            utils.require_integer(prime, "prime", 2)
            if not previous < prime < self.bound:
                raise ValueError(
                    "factor-base primes must increase below bound"
                )
            if utils.classify_prime(prime) != utils.Primality.PROVEN:
                raise ValueError("factor-base entry is not a proven prime")
            residue = target % prime
            roots = entry.square_roots
            if not isinstance(roots, tuple):
                raise TypeError("square_roots must be a tuple")
            for root in roots:
                utils.require_integer(root, "root", 0)
            if roots != tuple(sorted(set(roots))) or any(
                root >= prime or root * root % prime != residue
                for root in roots
            ):
                raise ValueError("invalid factor-base square roots")

            if prime == 2 or residue == 0:
                expected_count = 1
            else:
                expected_count = 2
            if len(roots) != expected_count:
                raise ValueError("factor-base roots are incomplete")
            previous = prime

        object.__setattr__(
            self,
            "_columns",
            MappingProxyType(
                {entry.prime: i + 1 for i, entry in enumerate(self.entries)}
            ),
        )
        object.__setattr__(
            self, "_primes", tuple(entry.prime for entry in self.entries)
        )

    @property
    def n_prime(self):
        """Return the exact multiplied target."""
        return self.n * self.multiplier

    @property
    def primes(self):
        """Return primes in stable parity-column order, excluding sign."""
        return self._primes

    @property
    def workspace_bytes(self):
        """Reserve sieve/setup and immutable root storage conservatively."""
        return 8192 + 64 * self.bound + 8 * self.n_prime.bit_length()


@dataclass(frozen=True)
class FactorBaseBuild:
    """Either a usable base or a validated proper divisor found in setup."""

    factor_base: FactorBase | None
    divisor: int | None


def build_factor_base(
    n,
    *,
    multiplier=1,
    bound=100,
    budget=None,
    memory_bytes=DEFAULT_MEMORY_BYTES,
    backend=None,
):
    """Build roots for primes p < bound with gcd(h,n) checked first.

    Proper GCDs return a divisor instead of a base. A multiplier divisible
    by all of n raises ValueError. Memory failure occurs before allocation;
    budget exhaustion raises BudgetExhaustedError without publishing a base.
    """
    if backend is not None:
        n = arithmetic.get_backend(backend).integer(n)
    target = checked_target(n, multiplier)
    utils.require_integer(bound, "bound", 3)
    utils.require_integer(memory_bytes, "memory_bytes", 0)
    if bound > MAX_FACTOR_BASE_BOUND:
        raise ValueError("factor-base bound exceeds the reference limit")
    budget = budget if budget is not None else Budget()
    budget.consume(n.bit_length())
    divisor = gcd(multiplier, n)
    if utils.valid_divisor(divisor, n):
        return FactorBaseBuild(None, divisor)
    if divisor != 1:
        raise ValueError("multiplier is divisible by n")
    reserve = 8192 + 64 * bound + 8 * target.bit_length()
    if reserve > memory_bytes:
        raise MemoryError("factor-base workspace exceeds memory_bytes")
    budget.consume(bound)
    entries = []

    for prime in prime_sieve.prime_sieve(bound):
        budget.consume(n.bit_length())
        divisor = gcd(prime, n)
        if utils.valid_divisor(divisor, n):
            return FactorBaseBuild(None, divisor)
        roots = modular_square_roots(target, prime, budget=budget)
        if roots:
            entries.append(FactorBaseEntry(prime, roots))

    budget.consume(len(entries) * bound.bit_length() ** 2)
    base = FactorBase(n, multiplier, bound, tuple(entries))
    return FactorBaseBuild(base, None)
