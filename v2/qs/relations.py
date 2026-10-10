"""Immutable exact relations, checked provenance, and square corrections."""

import hashlib
from collections import Counter
from dataclasses import dataclass, field

from .. import arithmetic, utils
from ..arithmetic import gcd, pow
from ..budget import Budget
from .factor_base import DEFAULT_MEMORY_BYTES
from .polynomial import Polynomial, checked_position

MAX_RESIDUAL = 1_000_000_000_000
MAX_COMBINED_ATOMS = 256


def checked_residual_bound(bound):
    """Limit reference residuals to a small deterministic prime domain."""
    utils.require_integer(bound, "residual_bound", 1)
    if bound > MAX_RESIDUAL:
        raise ValueError("residual_bound exceeds the reference limit")
    return bound


def _check_exponents(exponents):
    """Validate immutable, sorted sparse positive prime/exponent pairs."""
    if not isinstance(exponents, tuple):
        raise TypeError("exponents must be a tuple")
    previous = 1

    for pair in exponents:
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise TypeError("exponent entries must be pairs in tuples")
        prime, exponent = pair
        utils.require_integer(prime, "prime", 2)
        utils.require_integer(exponent, "exponent", 1)
        if prime <= previous:
            raise ValueError("exponent primes must be strictly increasing")
        previous = prime


@dataclass(frozen=True)
class AtomicRelation:
    """Factorization of A*F(x), including a known square and proven residual.

    relation_id identifies the polynomial/position, not untrusted exponent
    content. Construction validates shape; verify_atomic validates the math
    before admission. The original U is derived from immutable provenance.
    """

    polynomial: Polynomial
    position: int
    sign: int
    exponents: tuple[tuple[int, int], ...]
    residual: int = 1
    large_primes: tuple[int, ...] = field(default=(), kw_only=True)
    relation_id: str = field(init=False)

    def __post_init__(self):
        """Reject mutable payloads and derive the stable atomic ID."""
        if not isinstance(self.polynomial, Polynomial):
            raise TypeError("polynomial must be a Polynomial")
        checked_position(self.position)
        utils.require_integer(self.sign, "sign")
        if self.sign not in (-1, 1):
            raise ValueError("sign must be -1 or 1")
        _check_exponents(self.exponents)
        utils.require_integer(self.residual, "residual", 1)
        if self.residual > MAX_RESIDUAL:
            raise ValueError("residual exceeds the reference limit")
        if not isinstance(self.large_primes, tuple) or len(
            self.large_primes
        ) not in (0, 2):
            raise TypeError("large_primes must be an empty tuple or a pair")
        if self.large_primes:
            if self.residual != 1:
                raise ValueError(
                    "a two-prime atom must retain unit scalar residual"
                )
            for prime in self.large_primes:
                utils.require_integer(prime, "large prime", 2)
                if prime > MAX_RESIDUAL:
                    raise ValueError("large prime exceeds the reference limit")
            if self.large_primes[0] > self.large_primes[1]:
                raise ValueError("large primes must be sorted")
        encoded = f"{self.polynomial.identity}:{hex(self.position)}"
        object.__setattr__(
            self, "relation_id", hashlib.sha256(encoded.encode()).hexdigest()
        )

    @property
    def residual_primes(self):
        """Return endpoints after verification; the scalar remains prime."""
        return self.large_primes or (
            (self.residual,) if self.residual != 1 else ()
        )

    @property
    def residual_product(self):
        return (
            self.large_primes[0] * self.large_primes[1]
            if self.large_primes
            else self.residual
        )

    @property
    def u(self):
        """Return the original A*x+B without storing duplicate state."""
        return self.polynomial.u_value(self.position)

    @property
    def square_correction(self):
        """Known square coefficient; exact verification checks its identity."""
        return self.polynomial.square_coefficient


def verify_atomic(
    relation,
    factor_base,
    *,
    residual_bound=1,
    large_prime_bound=0,
    large_product_bound=0,
    budget=None,
):
    """Return True for an exact relation; reject invalid math with ValueError.

    Verify the full integer identity, including A, sign and known square.
    Nonunit residuals must be proven primes <= residual_bound < 2**64.
    Exponent bounds are checked before powers to avoid unbounded allocation.
    BudgetExhaustedError propagates; no partially checked atom is admitted.
    """
    checked_residual_bound(residual_bound)
    polynomial = relation.polynomial
    if (polynomial.n, polynomial.multiplier) != (
        factor_base.n,
        factor_base.multiplier,
    ):
        raise ValueError("relation target differs from the factor base")

    if len(relation.exponents) > len(factor_base.entries):
        raise ValueError("relation has too many factor-base exponents")
    budget = budget if budget is not None else Budget()
    value = polynomial.a * polynomial.value(relation.position)
    budget.consume((len(relation.exponents) + 1) * abs(value).bit_length())
    if relation.u * relation.u - factor_base.n_prime != value:
        raise ValueError("polynomial identity failed")
    if value == 0:
        raise ValueError("zero is not a factorable relation")
    if relation.sign != (-1 if value < 0 else 1):
        raise ValueError("relation sign is incorrect")
    remaining = arithmetic.divexact(abs(value), relation.square_correction**2)
    primes = factor_base._columns

    for prime, exponent in relation.exponents:
        if prime not in primes:
            raise ValueError("exponent prime is outside the factor base")
        if exponent * (prime.bit_length() - 1) > remaining.bit_length():
            raise ValueError("exponent exceeds the relation bit bound")
        power = pow(prime, exponent)
        remaining, remainder = divmod(remaining, power)
        if remainder:
            raise ValueError("relation exponents do not divide its value")

    if remaining != relation.residual_product:
        raise ValueError("relation factorization is incomplete or incorrect")
    if relation.large_primes:
        checked_residual_bound(large_prime_bound)
        utils.require_integer(large_product_bound, "large_product_bound", 1)
        if large_product_bound > MAX_RESIDUAL**2:
            raise ValueError("large product bound exceeds the reference limit")
        if (
            remaining > large_product_bound
            or relation.large_primes[-1] > large_prime_bound
        ):
            raise ValueError("two-prime residual exceeds its bounds")
    if relation.residual > residual_bound:
        raise ValueError("residual exceeds residual_bound")
    for residual in relation.residual_primes:
        if residual in primes:
            raise ValueError("factor-base exponents must be fully recovered")
        budget.consume(residual.bit_length() ** 2)
        if utils.classify_prime(residual) != utils.Primality.PROVEN:
            raise ValueError("residual must be a proven prime")
    return True


def parity_bits(relation, factor_base):
    """Return sign/prime exponent parity with sign in bit zero.

    The caller must first verify the atomic or combined relation. Full
    exponents and square corrections remain stored for later extraction.
    """
    bits = arithmetic.backend_for(factor_base.n).integer(
        int(relation.sign < 0)
    )
    columns = factor_base._columns
    for prime, exponent in relation.exponents:
        if prime not in columns:
            raise ValueError("exponent prime is outside the factor base")
        if exponent % 2:
            bits |= 1 << columns[prime]

    return bits


@dataclass(frozen=True)
class CombinedRelation:
    """U^2 = sign*product(p^e)*square_correction^2 modulo n.

    atom_ids are immutable references into the caller's checked atom store.
    U and square_correction are canonical residues modulo n. This is a full
    relation: each nonunit residual occurs an even number of times.
    """

    atom_ids: tuple[str, ...]
    u: int
    sign: int
    exponents: tuple[tuple[int, int], ...]
    square_correction: int

    def __post_init__(self):
        """Validate finite, unique provenance and immutable field shapes."""
        if not isinstance(self.atom_ids, tuple):
            raise TypeError("atom_ids must be a tuple")
        if not 1 <= len(self.atom_ids) <= MAX_COMBINED_ATOMS:
            raise ValueError("combined provenance exceeds the atom limit")
        if len(set(self.atom_ids)) != len(self.atom_ids):
            raise ValueError("an atomic position cannot be used twice")
        if any(not isinstance(identity, str) for identity in self.atom_ids):
            raise TypeError("atomic IDs must be strings")
        utils.require_integer(self.u, "u", 0)
        utils.require_integer(self.square_correction, "square_correction", 0)
        utils.require_integer(self.sign, "sign")
        if self.sign not in (-1, 1):
            raise ValueError("sign must be -1 or 1")
        _check_exponents(self.exponents)


@dataclass(frozen=True)
class CombinationResult:
    """A verified full combination, or a proper residual GCD divisor."""

    relation: CombinedRelation | None
    divisor: int | None


def _combined_values(atoms, modulus):
    """Accumulate exact exponents and bounded modular square corrections."""
    exponents, residuals = Counter(), Counter()
    u, sign, correction = 1, 1, 1

    for atom in atoms:
        u = u * atom.u % modulus
        sign *= atom.sign
        correction = correction * atom.square_correction % modulus
        for prime, exponent in atom.exponents:
            exponents[prime] += exponent
        if atom.large_primes:
            residuals.update(atom.large_primes)
        elif atom.residual != 1:
            residuals[atom.residual] += 1

    for residual, count in sorted(residuals.items()):
        if count % 2:
            raise ValueError("combined residual multiplicities must be even")
        # Paired residuals enter the square root, not factor-base parity.
        correction = correction * pow(residual, count // 2, modulus) % modulus
    return u, sign, tuple(sorted(exponents.items())), correction


def combined_storage_reserve(exponent_count, modulus_bits):
    """Bound a retained sparse pair, including IDs and modular values."""
    return 4096 + 256 * exponent_count + 16 * modulus_bits


def _combination_workspace(atoms, factor_base, memory_bytes):
    """Reserve referenced atoms, exponent counters, and modular temporaries."""
    utils.require_integer(memory_bytes, "memory_bytes", 0)
    reserve = factor_base.workspace_bytes + 8192
    reserve += 256 * len(factor_base.entries)
    for atom in atoms:
        polynomial = atom.polynomial
        numeric_bits = (
            polynomial.a.bit_length()
            + abs(polynomial.b).bit_length()
            + abs(atom.position).bit_length()
            + factor_base.n_prime.bit_length()
        )
        reserve += 2048 + 128 * len(atom.exponents) + 4 * numeric_bits

    if reserve > memory_bytes:
        raise MemoryError(
            "combination/provenance workspace exceeds memory_bytes"
        )
    return reserve


def combine_relations(
    atoms, factor_base, *, budget=None, memory_bytes=DEFAULT_MEMORY_BYTES
):
    """Combine a bounded explicit tuple of checked atoms, without inversion.

    Retain residual square corrections rather than dividing U by residuals.
    GCD-check residuals before combination and validate every returned split.
    Matching/storage policy is left to P3.2; this helper does no searching.
    Referenced atoms count toward a conservative owned-memory reservation;
    refusal raises MemoryError before combination or relation publication.
    """
    if not isinstance(atoms, tuple):
        raise TypeError("atoms must be a tuple")
    if not 1 <= len(atoms) <= MAX_COMBINED_ATOMS:
        raise ValueError("combined provenance exceeds the atom limit")
    identities = tuple(atom.relation_id for atom in atoms)
    if len(set(identities)) != len(identities):
        raise ValueError("an atomic position cannot be used twice")
    _combination_workspace(atoms, factor_base, memory_bytes)
    budget = budget if budget is not None else Budget()
    for atom in atoms:
        verify_atomic(
            atom,
            factor_base,
            residual_bound=MAX_RESIDUAL,
            large_prime_bound=MAX_RESIDUAL,
            large_product_bound=MAX_RESIDUAL**2,
            budget=budget,
        )
    for atom in atoms:
        for residual in atom.residual_primes:
            budget.consume(factor_base.n.bit_length())
            divisor = gcd(residual, factor_base.n)
            if utils.valid_divisor(divisor, factor_base.n):
                return CombinationResult(None, divisor)
            if divisor != 1:
                raise ValueError("residual is a nonunit with no proper split")

    budget.consume(len(atoms) * (len(factor_base.entries) + 1))
    values = _combined_values(atoms, factor_base.n)
    combined = CombinedRelation(identities, *values)
    store = {atom.relation_id: atom for atom in atoms}
    verify_combined(
        combined, factor_base, store, budget=budget, memory_bytes=memory_bytes
    )
    return CombinationResult(combined, None)


def verify_combined(
    relation,
    factor_base,
    atom_store,
    *,
    budget=None,
    memory_bytes=DEFAULT_MEMORY_BYTES,
):
    """Recheck every atom, provenance fields, and the combined congruence.

    Missing or altered provenance and plausible but incorrect modular
    payloads raise ValueError. Verification never forms an integer product
    of all atomic values. The original atom store must remain available.
    """
    budget = budget if budget is not None else Budget()
    atoms = []
    for identity in relation.atom_ids:
        atom = atom_store.get(identity)
        if atom is None or atom.relation_id != identity:
            raise ValueError("missing or mismatched atomic provenance")
        atoms.append(atom)

    _combination_workspace(atoms, factor_base, memory_bytes)
    for atom in atoms:
        verify_atomic(
            atom,
            factor_base,
            residual_bound=MAX_RESIDUAL,
            large_prime_bound=MAX_RESIDUAL,
            large_product_bound=MAX_RESIDUAL**2,
            budget=budget,
        )
    budget.consume(len(atoms) * (len(factor_base.entries) + 1))
    values = _combined_values(atoms, factor_base.n)
    if values != (
        relation.u,
        relation.sign,
        relation.exponents,
        relation.square_correction,
    ):
        raise ValueError("combined fields differ from checked provenance")

    right = relation.sign % factor_base.n
    for prime, exponent in relation.exponents:
        right = right * pow(prime, exponent, factor_base.n) % factor_base.n
    right = right * pow(relation.square_correction, 2, factor_base.n)
    if relation.u * relation.u % factor_base.n != right % factor_base.n:
        raise ValueError("combined congruence failed")
    return True
