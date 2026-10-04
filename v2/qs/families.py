"""Exact bounded SIQS CRT families, Gray transitions and reusable roots."""

import hashlib
import json
import math
import random
from dataclasses import dataclass
from math import comb, prod

from .. import utils
from ..budget import Budget
from .factor_base import DEFAULT_MEMORY_BYTES
from .polynomial import Polynomial, PolynomialRoots, a_target, polynomial_roots

MAX_A_FACTORS = 32
MAX_FAMILIES = 64
MAX_FAMILY_POOL = 128
FAMILY_CHECKPOINT_VERSION = 1
MAX_CHECKPOINT_BYTES = 4096


def _identity(base):
    """Bind cached roots to the exact target, ordered base and endpoints."""
    payload = [
        base.n,
        base.multiplier,
        base.bound,
        [[entry.prime, entry.square_roots] for entry in base.entries],
    ]
    encoded = json.dumps(payload, separators=(",", ":"))
    return hashlib.sha256(encoded.encode()).hexdigest()


def _checksum(payload):
    """Canonical integrity marker; this is not an authentication signature."""
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode()).hexdigest()


def _checked_resources(resources):
    """Bound serialized resource fields before hashing untrusted state."""
    if not isinstance(resources, dict) or set(resources) != {
        "work_used",
        "wall_used",
        "cpu_used",
    }:
        raise ValueError("invalid family checkpoint resource record")
    work = utils.require_integer(resources["work_used"], "work_used", 0)
    if work.bit_length() > 64:
        raise ValueError("family checkpoint work exceeds 64 bits")
    for name in ("wall_used", "cpu_used"):
        value = resources[name]
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or (isinstance(value, int) and value.bit_length() > 64)
            or not math.isfinite(value)
            or value < 0
        ):
            raise ValueError("invalid family checkpoint time record")
    return resources


def family_assignments(
    base,
    half_width,
    *,
    factor_count=3,
    family_count=16,
    pool_size=16,
    seed=7,
    budget=None,
    memory_bytes=DEFAULT_MEMORY_BYTES,
):
    """Return reproducible distinct squarefree A sets near the integer target.

    Select a finite pool by distance from p**factor_count to the target;
    start with its nearest set, then sample diverse sets with a local seed.
    This is an untuned control, not a production parameter or digit table.
    Fewer sets may be returned if the pool or finite attempt cap is exhausted.
    """
    for name, value, minimum, maximum in (
        ("factor_count", factor_count, 1, 8),
        ("family_count", family_count, 1, MAX_FAMILIES),
        ("pool_size", pool_size, factor_count, MAX_FAMILY_POOL),
        ("seed", seed, 0, 2**64 - 1),
    ):
        utils.require_integer(value, name, minimum)
        if value > maximum:
            raise ValueError(f"{name} exceeds the family limit")
    utils.require_integer(memory_bytes, "memory_bytes", 0)
    reserve = base.workspace_bytes + 32768 + 256 * len(base.entries)
    reserve += family_count * factor_count * 256
    if reserve > memory_bytes:
        raise MemoryError("family assignment workspace exceeds memory_bytes")
    budget = budget if budget is not None else Budget()
    target = a_target(base.n_prime, half_width)
    budget.consume(len(base.entries) * (factor_count + 1))
    eligible = [
        entry.prime
        for entry in base.entries
        if entry.prime != 2 and base.n_prime % entry.prime
    ]
    pool = sorted(
        eligible, key=lambda prime: (abs(prime**factor_count - target), prime)
    )[:pool_size]
    if len(pool) < factor_count:
        raise ValueError("too few nonsingular primes for the requested A")
    wanted = min(family_count, comb(len(pool), factor_count))
    first = tuple(sorted(pool[:factor_count]))
    assignments, seen = [first], {first}
    generator = random.Random(seed)
    for _ in range(MAX_FAMILIES * 64):
        if len(assignments) == wanted:
            break
        budget.consume(factor_count + 1)
        values = tuple(sorted(generator.sample(pool, factor_count)))
        if values not in seen:
            assignments.append(values)
            seen.add(values)
    return tuple(assignments)


def verify_polynomial_roots(polynomial, base, roots, *, budget=None):
    """Certify a supplied complete root tuple without recomputing inverses.

    Root cardinality and exact modular evaluation certify nonsingular
    quadratic roots; linear/constant and characteristic-two cases retain
    explicit complete enumeration. Mutable or incomplete data is rejected.
    """
    if not isinstance(roots, tuple) or len(roots) != len(base.entries):
        raise ValueError("cached root tuple has the wrong base shape")
    if (polynomial.n, polynomial.multiplier) != (base.n, base.multiplier):
        raise ValueError("cached roots belong to a different target")
    budget = budget if budget is not None else Budget()
    for entry, item in zip(base.entries, roots):
        prime = entry.prime
        budget.consume(prime.bit_length() ** 2)
        if not isinstance(item, PolynomialRoots):
            raise ValueError("cached root has the wrong type")
        utils.require_integer(item.prime, "cached prime", 2)
        if item.prime != prime:
            raise ValueError("cached root prime order mismatch")
        if type(item.all_positions) is not bool or not isinstance(
            item.roots, tuple
        ):
            raise ValueError("cached roots must have immutable typed fields")
        for root in item.roots:
            utils.require_integer(root, "root", 0)
            if root >= prime:
                raise ValueError("cached root is outside its modulus")
        if item.roots != tuple(sorted(set(item.roots))):
            raise ValueError("cached roots must be distinct and increasing")
        a, b, c = (
            polynomial.a % prime,
            polynomial.b % prime,
            polynomial.c % prime,
        )
        if prime == 2:
            expected = tuple(
                x for x in (0, 1) if (a * x * x + 2 * b * x + c) % 2 == 0
            )
            all_positions = len(expected) == 2
            expected = () if all_positions else expected
            if (item.roots, item.all_positions) != (expected, all_positions):
                raise ValueError(
                    "cached characteristic-two roots are incomplete"
                )
            continue
        if a == 0 and 2 * b % prime == 0:
            if item.roots or item.all_positions != (c == 0):
                raise ValueError(
                    "cached constant-polynomial roots are incomplete"
                )
            continue
        count = 1 if a == 0 else len(entry.square_roots)
        if item.all_positions or len(item.roots) != count:
            raise ValueError("cached root cardinality is incomplete")
        if any((a * x * x + 2 * b * x + c) % prime for x in item.roots):
            raise ValueError("cached value is not a polynomial root")
    return True


@dataclass(frozen=True)
class FamilyStep:
    """One exact polynomial and complete roots, including actual B shift."""

    family_identity: str
    gray_index: int
    polynomial: Polynomial
    roots: tuple
    delta_b: int
    changed_prime: int | None


class PolynomialFamily:
    """Finite CRT family; one sign is fixed to avoid opposite-B duplicates.

    A is a product of distinct odd nonsingular factor-base primes. Its CRT
    terms choose every remaining sign in Gray order. B is recentered to the
    nearest residue to zero. Cached A inverses update nonsingular roots using
    the actual recentered B difference; 2 and p|A use their exact branches.
    """

    def __init__(
        self,
        base,
        a_primes,
        *,
        budget=None,
        memory_bytes=DEFAULT_MEMORY_BYTES,
    ):
        """Reserve root/CRT state, then cache exact modular inverses."""
        if not isinstance(a_primes, tuple):
            raise TypeError("A primes must be an immutable tuple")
        if not 1 <= len(a_primes) <= MAX_A_FACTORS:
            raise ValueError("A factor count exceeds the family limit")
        for prime in a_primes:
            utils.require_integer(prime, "A prime", 3)
        if a_primes != tuple(sorted(set(a_primes))):
            raise ValueError("A primes must be distinct and increasing")
        if any(
            p not in base.primes or base.n_prime % p == 0 for p in a_primes
        ):
            raise ValueError("A primes must be nonsingular base members")
        utils.require_integer(memory_bytes, "memory_bytes", 0)
        a = prod(a_primes)
        reserve = base.workspace_bytes + 32768 + 640 * len(base.entries)
        reserve += 512 * len(a_primes) + 16 * (
            a.bit_length() + base.n.bit_length()
        )
        if reserve > memory_bytes:
            raise MemoryError("family roots/CRT state exceeds memory_bytes")
        self.base, self.a_primes, self.a = base, a_primes, a
        self.budget = budget if budget is not None else Budget()
        self.memory_bytes, self.workspace_bytes = memory_bytes, reserve
        entries = {
            entry.prime: entry
            for entry in base.entries
            if entry.prime in a_primes
        }
        terms = []
        for prime in a_primes:
            self.budget.consume(a.bit_length() + prime.bit_length() ** 2)
            quotient = a // prime
            inverse = utils.modular_inverse(quotient, prime)
            terms.append(
                quotient * inverse * entries[prime].square_roots[0] % a
            )
        inverses = []
        for entry in base.entries:
            prime = entry.prime
            self.budget.consume(prime.bit_length() ** 2)
            inverses.append(
                None
                if prime == 2 or a % prime == 0
                else utils.modular_inverse(a, prime)
            )
        self.terms, self.inverses = tuple(terms), tuple(inverses)
        self.base_identity = _identity(base)
        self.identity = _checksum(
            [self.base_identity, a_primes, "center-zero-v1"]
        )
        self.count = 1 << (len(a_primes) - 1)
        self.next_index, self.current = 0, None

    def _make_step(self, index, previous):
        """Construct privately before publishing complete roots."""
        self.budget.consume(
            (len(self.base.entries) + len(self.terms) + 1)
            * (self.a.bit_length() + self.base.n.bit_length() + 1)
        )
        gray = index ^ (index >> 1)
        changed = None
        if previous is None:
            raw = self.terms[0]
            for bit, term in enumerate(self.terms[1:]):
                raw += -term if gray & (1 << bit) else term
        else:
            old_gray = previous.gray_index ^ (previous.gray_index >> 1)
            bit = (gray ^ old_gray).bit_length() - 1
            old_sign = -1 if old_gray & (1 << bit) else 1
            raw = previous.polynomial.b - 2 * old_sign * self.terms[bit + 1]
            changed = self.a_primes[bit + 1]
        b = (raw + self.a // 2) % self.a - self.a // 2
        polynomial = Polynomial(self.base.n, self.base.multiplier, self.a, b)
        delta = 0 if previous is None else b - previous.polynomial.b
        roots = []
        for offset, (entry, inverse) in enumerate(
            zip(self.base.entries, self.inverses)
        ):
            if inverse is None:
                roots.append(
                    polynomial_roots(
                        polynomial, self.base, entry, budget=self.budget
                    )
                )
                continue
            residues = (
                tuple(
                    (root - b) * inverse % entry.prime
                    for root in entry.square_roots
                )
                if previous is None
                else tuple(
                    (root - delta * inverse) % entry.prime
                    for root in previous.roots[offset].roots
                )
            )
            roots.append(PolynomialRoots(entry.prime, tuple(sorted(residues))))
        return FamilyStep(
            self.identity, index, polynomial, tuple(roots), delta, changed
        )

    def next(self):
        """Return the next checked-shape step, or None at finite exhaustion.

        Budget refusal retains the first unpublished index and prior root
        state. Retrying charges repeated private work under the same allowance.
        """
        if self.next_index == self.count:
            return None
        step = self._make_step(self.next_index, self.current)
        self.current = step
        self.next_index += 1
        return step

    def checkpoint(self):
        """Save compact family progress/resources; root caches are rebuilt.

        This is family state only, not a full SIQS relation-job checkpoint.
        Retain consumed work/time when creating a resume Budget.
        """
        payload = {
            "version": FAMILY_CHECKPOINT_VERSION,
            "base_identity": self.base_identity,
            "a_primes": list(self.a_primes),
            "next_index": self.next_index,
            "resources": {
                "work_used": self.budget.used,
                "wall_used": self.budget.wall_used,
                "cpu_used": self.budget.cpu_used,
            },
        }
        _checked_resources(payload["resources"])
        if len(json.dumps(payload).encode()) > MAX_CHECKPOINT_BYTES:
            raise ValueError("family checkpoint exceeds its byte limit")
        return {"payload": payload, "sha256": _checksum(payload)}

    @classmethod
    def from_checkpoint(
        cls,
        base,
        checkpoint,
        *,
        budget,
        memory_bytes=DEFAULT_MEMORY_BYTES,
    ):
        """Verify progress and reconstruct charged caches on resume."""
        if not isinstance(checkpoint, dict) or set(checkpoint) != {
            "payload",
            "sha256",
        }:
            raise ValueError("invalid family checkpoint envelope")
        payload = checkpoint["payload"]
        if not isinstance(payload, dict) or set(payload) != {
            "version",
            "base_identity",
            "a_primes",
            "next_index",
            "resources",
        }:
            raise ValueError("invalid family checkpoint payload")
        utils.require_integer(payload["version"], "checkpoint version", 1)
        if payload["version"] != FAMILY_CHECKPOINT_VERSION:
            raise ValueError("unsupported family checkpoint version")
        if (
            not isinstance(payload["base_identity"], str)
            or len(payload["base_identity"]) != 64
            or not isinstance(checkpoint["sha256"], str)
            or len(checkpoint["sha256"]) != 64
        ):
            raise ValueError("invalid family checkpoint identity shape")
        if payload["base_identity"] != _identity(base):
            raise ValueError("family checkpoint factor-base identity mismatch")
        primes = payload["a_primes"]
        if (
            not isinstance(primes, list)
            or not 1 <= len(primes) <= MAX_A_FACTORS
        ):
            raise ValueError("invalid family checkpoint A factors")
        for prime in primes:
            utils.require_integer(prime, "checkpoint A prime", 3)
            if prime not in base.primes:
                raise ValueError("checkpoint A prime is outside its base")
        index = utils.require_integer(payload["next_index"], "next_index", 0)
        if index > 1 << (len(primes) - 1):
            raise ValueError("family checkpoint index exceeds its limit")
        resources = _checked_resources(payload["resources"])
        if len(json.dumps(payload).encode()) > MAX_CHECKPOINT_BYTES:
            raise ValueError("family checkpoint exceeds its byte limit")
        if checkpoint["sha256"] != _checksum(payload):
            raise ValueError("family checkpoint integrity mismatch")
        prior = Budget(
            work_limit=budget.work_limit,
            used=resources["work_used"],
            prior_wall=resources["wall_used"],
            prior_cpu=resources["cpu_used"],
        )
        if (
            budget.used < prior.used
            or budget.wall_used < prior.prior_wall
            or budget.cpu_used < prior.prior_cpu
        ):
            raise ValueError("resume budget must retain consumed resources")
        family = cls(
            base, tuple(primes), budget=budget, memory_bytes=memory_bytes
        )
        if index:
            previous = (
                family._make_step(index - 2, None) if index > 1 else None
            )
            family.current = family._make_step(index - 1, previous)
        family.next_index = index
        return family
