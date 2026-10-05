"""Independent arithmetic, exhaustive roots, and corrupted QS provenance."""

import io
import unittest
from contextlib import redirect_stdout
from dataclasses import FrozenInstanceError, replace
from math import isqrt
from unittest.mock import patch

from v2.budget import Budget, BudgetExhaustedError
from v2.qs import (
    AtomicRelation,
    FactorBase,
    FactorBaseEntry,
    Polynomial,
    a_target,
    build_factor_base,
    collect_block,
    combine_relations,
    modular_square_roots,
    mpqs_polynomial,
    parity_bits,
    polynomial_roots,
    qs_polynomial,
    verify_atomic,
    verify_combined,
)


def reference_primes(bound):
    """Generate small primes by independent trial division."""
    return tuple(
        value
        for value in range(2, bound)
        if all(value % divisor for divisor in range(2, isqrt(value) + 1))
    )


def reference_factors(value):
    """Fully factor a small absolute value without the Factor utilities."""
    remaining, divisor, factors = abs(value), 2, []
    if remaining == 0:
        raise ValueError("zero cannot be trial factored")
    while divisor * divisor <= remaining:
        exponent = 0
        while remaining % divisor == 0:
            remaining //= divisor
            exponent += 1
        if exponent:
            factors.append((divisor, exponent))
        divisor += 1

    if remaining > 1:
        factors.append((remaining, 1))
    return tuple(factors)


def reference_positions(polynomial, base, lo, hi, residual_bound=1):
    """Enumerate admissible values from U squared and independent factors."""
    positions = {}

    for position in range(lo, hi):
        u = polynomial.a * position + polynomial.b
        value = u * u - base.n * base.multiplier
        if value == 0:
            continue
        factors = reference_factors(value)
        outside = [(p, e) for p, e in factors if p not in base.primes]
        if outside and not (
            len(outside) == 1
            and outside[0][1] == 1
            and outside[0][0] <= residual_bound
        ):
            continue

        residual = outside[0][0] if outside else 1
        exponents = tuple((p, e) for p, e in factors if p in base.primes)
        positions[position] = (-1 if value < 0 else 1, exponents, residual)

    return positions


def unlimited_budget():
    """Use finite work but no timer noise for deterministic oracle sweeps."""
    return Budget(work_limit=100_000_000, seconds=None, cpu_seconds=None)


class FactorBaseTests(unittest.TestCase):
    """Root completeness, half-open bounds, multipliers, and setup caps."""

    def test_all_small_modular_roots(self):
        """Exhaust residues, including the Tonelli-Shanks p=1 mod 4 path."""
        for prime in reference_primes(150):
            for value in range(-prime, prime):
                expected = tuple(
                    x for x in range(prime) if x * x % prime == value % prime
                )
                self.assertEqual(modular_square_roots(value, prime), expected)

    def test_factor_bases_against_enumeration(self):
        """Include 2, zero roots from h, and only quadratic residues."""
        for multiplier in (1, 2, 3, 5, 9, 15):
            for bound in (3, 7, 30, 31, 32, 100):
                base = build_factor_base(
                    101 * 103, multiplier=multiplier, bound=bound
                ).factor_base
                expected = []

                for prime in reference_primes(bound):
                    roots = tuple(
                        x
                        for x in range(prime)
                        if x * x % prime == base.n_prime % prime
                    )
                    if roots:
                        expected.append((prime, roots))

                self.assertEqual(
                    [
                        (entry.prime, entry.square_roots)
                        for entry in base.entries
                    ],
                    expected,
                )

    def test_setup_gcd_splits_are_proper(self):
        """Both multiplier and factor-base factors surface explicitly."""
        for kwargs, expected in (
            ({"multiplier": 7, "bound": 3}, 7),
            ({"bound": 20}, 7),
        ):
            result = build_factor_base(77, **kwargs)

            self.assertIsNone(result.factor_base)
            self.assertEqual(result.divisor, expected)
            self.assertEqual(77 % result.divisor, 0)

        with self.assertRaises(ValueError):
            build_factor_base(77, multiplier=77)

    def test_setup_limits_and_invalid_roots(self):
        """Fail before publishing oversized, mutable, or false base data."""
        for n in (2, 4, 2**4096 + 1):
            with self.assertRaises(ValueError):
                build_factor_base(n)
        with self.assertRaises(TypeError):
            build_factor_base(True)
        with self.assertRaises(MemoryError):
            build_factor_base(10403, memory_bytes=1)
        with self.assertRaises(BudgetExhaustedError):
            build_factor_base(10403, budget=Budget(work_limit=0))
        for prime in (1, 4, 9, 100001):
            with self.assertRaises(ValueError):
                modular_square_roots(1, prime)
        with self.assertRaises(ValueError):
            FactorBase(10403, 1, 7, (FactorBaseEntry(2, (0,)),))
        with self.assertRaises(ValueError):
            FactorBase(
                10403,
                1,
                7,
                (FactorBaseEntry(2, (1,)), FactorBaseEntry(3, (1,))),
            )

        with self.assertRaises(TypeError):
            FactorBase(10403, 1, 7, [FactorBaseEntry(2, (1,))])


class PolynomialTests(unittest.TestCase):
    """Exact targets, MPQS congruences, root degeneracy, and translations."""

    def test_exact_identity_and_integer_targets(self):
        """Large inputs beyond float range retain exact floor comparisons."""
        target = 10**1000 + 123
        for width in (1, 2, 31, 10**600):
            value = a_target(target, width)
            if value > 1:
                self.assertLessEqual((value * width) ** 2, 2 * target)

            self.assertGreater(((value + 1) * width) ** 2, 2 * target)

        polynomial = Polynomial(target, 1, 1, isqrt(target) + 1)
        for position in (-(2**200), -1, 0, 2**200):
            self.assertEqual(
                polynomial.u_value(position) ** 2 - target,
                polynomial.a * polynomial.value(position),
            )

        with self.assertRaises(ValueError):
            Polynomial(10403, 1, 7, 2)
        with self.assertRaises(ValueError):
            a_target(10, 0)
        with self.assertRaises(ValueError):
            Polynomial(10403, 1, 2**4096, 1)
        with self.assertRaises(ValueError):
            polynomial.value(2**4096)

    def test_exhaustive_normalized_polynomial_roots(self):
        """Enumerate x mod p independently, including p|A and repeated A."""
        for multiplier in (1, 2, 3, 9):
            base = build_factor_base(
                10403, multiplier=multiplier, bound=40
            ).factor_base

            for a in range(1, 33):
                for b in range(a):
                    if (b * b - base.n_prime) % a:
                        continue
                    for offset in (-2, 0, 2):
                        polynomial = Polynomial(
                            10403, multiplier, a, b + offset * a
                        )

                        for entry in base.entries:
                            answer = polynomial_roots(polynomial, base, entry)
                            expected = tuple(
                                x
                                for x in range(entry.prime)
                                if polynomial.value(x) % entry.prime == 0
                            )
                            actual = (
                                tuple(range(entry.prime))
                                if answer.all_positions
                                else answer.roots
                            )

                            self.assertEqual(actual, expected)

    def test_degenerate_and_linear_roots(self):
        """Distinguish all/none/one roots rather than skipping primes in A."""
        base = build_factor_base(10403, multiplier=9, bound=40).factor_base
        entry = next(entry for entry in base.entries if entry.prime == 3)

        self.assertTrue(
            polynomial_roots(
                Polynomial(10403, 9, 3, 0), base, entry
            ).all_positions
        )
        answer = polynomial_roots(Polynomial(10403, 9, 9, 0), base, entry)

        self.assertEqual(answer.roots, ())
        self.assertFalse(answer.all_positions)
        base = build_factor_base(10403, bound=40).factor_base
        entry = next(entry for entry in base.entries if entry.prime == 7)

        self.assertEqual(
            len(
                polynomial_roots(Polynomial(10403, 1, 7, 1), base, entry).roots
            ),
            1,
        )

    def test_inversion_failure_and_identity_mismatch(self):
        """A failed inverse must propagate instead of manufacturing roots."""
        base = build_factor_base(10403, bound=40).factor_base
        polynomial = qs_polynomial(base)
        entry = next(entry for entry in base.entries if entry.prime == 7)
        with patch(
            "v2.qs.polynomial.utils.modular_inverse",
            side_effect=ValueError("nonunit"),
        ):
            with self.assertRaisesRegex(ValueError, "nonunit"):
                polynomial_roots(polynomial, base, entry)

        other = build_factor_base(10403, multiplier=3, bound=40).factor_base
        with self.assertRaises(ValueError):
            polynomial_roots(polynomial, other, other.entries[0])

    def test_reference_mpqs_selection_and_lift(self):
        """Check nearest prime-square A and B against exhaustive residues."""
        for n, multiplier in ((10403, 1), (10403, 9), (1022117, 3)):
            base = build_factor_base(
                n, multiplier=multiplier, bound=40
            ).factor_base

            for width in (1, 7, 16, 128):
                polynomial = mpqs_polynomial(base, width)
                target = a_target(base.n_prime, width)
                eligible = [
                    p
                    for p in reference_primes(40)
                    if p != 2
                    and base.n_prime % p
                    and any(x * x % p == base.n_prime % p for x in range(p))
                ]
                prime = min(eligible, key=lambda p: (abs(p * p - target), p))

                self.assertEqual(polynomial.a, prime * prime)
                roots = [
                    b
                    for b in range(polynomial.a)
                    if (b * b - base.n_prime) % polynomial.a == 0
                ]

                self.assertEqual(polynomial.b, min(roots))
                result = collect_block(polynomial, base, -17, 20)

                self.assertEqual(
                    {
                        atom.position: (
                            atom.sign,
                            atom.exponents,
                            atom.residual,
                        )
                        for atom in result.relations
                    },
                    reference_positions(polynomial, base, -17, 20),
                )

        base = build_factor_base(10403, bound=3).factor_base
        with self.assertRaises(ValueError):
            mpqs_polynomial(base, 16)


class CollectorTests(unittest.TestCase):
    """Exhaustive relation admission, signed blocks, caps, and zeros."""

    def setUp(self):
        """Use a target with no factor below the diagnostic base bound."""
        self.base = build_factor_base(10403, bound=40).factor_base
        self.polynomial = qs_polynomial(self.base)

    def test_exhaustive_qs_and_mpqs_collectors(self):
        """Trial factor U squared directly and compare every admitted atom."""
        for multiplier in (1, 3, 9):
            base = build_factor_base(
                10403, multiplier=multiplier, bound=40
            ).factor_base
            polynomials = [qs_polynomial(base)]
            for a in (3, 7, 9, 11, 49):
                for b in range(a):
                    if (b * b - base.n_prime) % a == 0:
                        polynomials.append(Polynomial(10403, multiplier, a, b))
                        break

            for polynomial in polynomials:
                for bound in (1, 500):
                    result = collect_block(
                        polynomial,
                        base,
                        -43,
                        58,
                        residual_bound=bound,
                        budget=unlimited_budget(),
                    )

                    self.assertEqual(result.reason, "complete")
                    actual = {
                        atom.position: (
                            atom.sign,
                            atom.exponents,
                            atom.residual,
                        )
                        for atom in result.relations
                    }

                    self.assertEqual(
                        actual,
                        reference_positions(polynomial, base, -43, 58, bound),
                    )

    def test_half_open_empty_blocks_and_tails(self):
        """Split blocks preserve the complete sequence and stable atom IDs."""
        full = collect_block(self.polynomial, self.base, -13, 28)
        pieces = []

        for lo, hi in (
            (-13, -6),
            (-6, 1),
            (1, 8),
            (8, 15),
            (15, 22),
            (22, 28),
        ):
            pieces.extend(
                collect_block(self.polynomial, self.base, lo, hi).relations
            )

        self.assertEqual(tuple(pieces), full.relations)
        for position in (-13, 0, 28):
            result = collect_block(
                self.polynomial, self.base, position, position
            )

            self.assertEqual(result.relations, ())
            self.assertEqual(result.next_position, position)
            self.assertEqual(result.scanned, 0)
            self.assertEqual(result.reason, "complete")

    def test_zero_value_before_division(self):
        """A square target splits; an injected improper GCD can only skip."""
        base = build_factor_base(121, bound=7).factor_base
        result = collect_block(qs_polynomial(base), base, 0, 1)

        self.assertEqual(result.divisor, 11)
        self.assertEqual(result.reason, "factor_found")
        self.assertEqual(result.zero_positions, (0,))
        # Fault injection covers a full-modulus GCD without weakening the
        # coprime multiplier contract to manufacture an invalid target.
        with patch("v2.qs.reference_collector.gcd", return_value=121):
            result = collect_block(qs_polynomial(base), base, 0, 1)

        self.assertIsNone(result.divisor)
        self.assertEqual(result.reason, "complete")
        self.assertEqual(result.zero_positions, (0,))

    def test_work_cancel_time_relation_and_storage_caps(self):
        """Every refused position leaves a checked prefix and a finite stop."""
        for budget, reason in (
            (Budget(work_limit=0), "work_limit"),
            (Budget(seconds=0), "wall_limit"),
            (Budget(cpu_seconds=0), "cpu_limit"),
            (Budget(cancelled=lambda: True), "cancelled"),
        ):
            result = collect_block(
                self.polynomial, self.base, -10, 10, budget=budget
            )

            self.assertEqual(result.reason, reason)
            self.assertEqual(result.next_position, -10)
            self.assertEqual(result.relations, ())

        for changes, reason in (
            ({"max_relations": 0}, "relation_limit"),
            ({"memory_bytes": 0}, "memory_limit"),
        ):
            result = collect_block(
                self.polynomial, self.base, -10, 10, **changes
            )

            self.assertEqual(result.reason, reason)
            self.assertEqual(result.next_position, -10)

        prefix = collect_block(
            self.polynomial, self.base, -50, 50, max_relations=2
        )

        self.assertEqual(prefix.reason, "relation_limit")
        suffix = collect_block(
            self.polynomial, self.base, prefix.next_position, 50
        )
        full = collect_block(self.polynomial, self.base, -50, 50)

        self.assertEqual(prefix.relations + suffix.relations, full.relations)
        cap = prefix.workspace_bytes
        limited = collect_block(
            self.polynomial, self.base, -50, 50, memory_bytes=cap
        )

        self.assertEqual(limited.reason, "memory_limit")
        self.assertLessEqual(limited.workspace_bytes, cap)
        for atom in limited.relations:
            self.assertTrue(verify_atomic(atom, self.base))

    def test_invalid_blocks_and_quiet_calls(self):
        """Reject reversed, oversized, or noninteger blocks; no diagnostics."""
        for lo, hi in ((1, 0), (0, 4097)):
            with self.assertRaises(ValueError):
                collect_block(self.polynomial, self.base, lo, hi)
        with self.assertRaises(TypeError):
            collect_block(self.polynomial, self.base, 0.0, 1)
        with self.assertRaises(ValueError):
            collect_block(
                self.polynomial, self.base, 0, 1, residual_bound=2**64
            )
        stream = io.StringIO()
        with redirect_stdout(stream):
            build_factor_base(10403, bound=40)
            collect_block(self.polynomial, self.base, -10, 10)

        self.assertEqual(stream.getvalue(), "")

    def test_work_refusals_preserve_the_exact_prefix(self):
        """Interrupt division/verification and continue at the refused x."""
        full = collect_block(self.polynomial, self.base, -20, 21)

        for work in range(0, 20000, 997):
            budget = Budget(work_limit=work, seconds=None, cpu_seconds=None)
            prefix = collect_block(
                self.polynomial, self.base, -20, 21, budget=budget
            )
            suffix = collect_block(
                self.polynomial, self.base, prefix.next_position, 21
            )

            self.assertEqual(
                prefix.relations + suffix.relations, full.relations
            )
            self.assertLessEqual(budget.used, work)


class RelationTests(unittest.TestCase):
    """Full exponent recovery, bounded powers, and checked sums."""

    def setUp(self):
        """Retain small full/partial atoms for corruption and pairing tests."""
        self.base = build_factor_base(10403, bound=40).factor_base
        self.polynomial = qs_polynomial(self.base)
        self.atoms = collect_block(
            self.polynomial, self.base, -100, 101, residual_bound=500
        ).relations

    def test_repeated_factors_a_exponents_and_sign_parity(self):
        """Compare full A*F factorization, sign bit, and repeated powers."""
        polynomial = Polynomial(10403, 1, 49, 8)
        result = collect_block(
            polynomial, self.base, -20, 21, residual_bound=500
        )

        self.assertTrue(result.relations)
        for atom in result.relations:
            self.assertGreaterEqual(dict(atom.exponents)[7], 2)
            self.assertEqual(
                parity_bits(atom, self.base) & 1, int(atom.sign < 0)
            )
            self.assertTrue(verify_atomic(atom, self.base, residual_bound=500))

    def test_atomic_corruption_and_bounded_power(self):
        """Reject wrong signs, missing factors, and oversized exponents."""
        atom = next(atom for atom in self.atoms if atom.exponents)

        for corrupted in (
            replace(atom, sign=-atom.sign),
            replace(atom, exponents=()),
            replace(atom, exponents=((2, 10**100),)),
            replace(atom, residual=4),
            replace(atom, position=atom.position + 1),
        ):
            with self.assertRaises(ValueError):
                verify_atomic(corrupted, self.base, residual_bound=500)

        partial = next(atom for atom in self.atoms if atom.residual > 1)
        with self.assertRaises(ValueError):
            verify_atomic(partial, self.base)
        with self.assertRaises(ValueError):
            AtomicRelation(self.polynomial, 0, 1, ((3, 1), (2, 1)))
        with self.assertRaises(TypeError):
            AtomicRelation(self.polynomial, 0, 1, [(2, 1)])
        with self.assertRaises(FrozenInstanceError):
            atom.sign = -atom.sign

    def matching_pair(self):
        """Find two explicitly supplied same-residual atoms for a fixture."""
        for left in self.atoms:
            if left.residual <= 1 or 10403 % left.residual == 0:
                continue
            for right in self.atoms:
                if (
                    left.relation_id != right.relation_id
                    and left.residual == right.residual
                ):
                    return left, right

        self.fail("fixture has no matching residual pair")

    def test_combination_square_correction_and_full_provenance(self):
        """Check exact two-norm product independently for the small fixture."""
        pair = self.matching_pair()
        combined = combine_relations(pair, self.base).relation
        store = {atom.relation_id: atom for atom in pair}

        self.assertTrue(verify_combined(combined, self.base, store))
        left_value = pair[0].u ** 2 - self.base.n_prime
        right_value = pair[1].u ** 2 - self.base.n_prime
        product = left_value * right_value
        expected = combined.sign * combined.square_correction**2
        for prime, exponent in combined.exponents:
            expected *= prime**exponent

        self.assertEqual(expected, product)
        self.assertEqual(combined.square_correction, pair[0].residual)
        self.assertEqual(combined.u, pair[0].u * pair[1].u % 10403)
        self.assertEqual(
            parity_bits(combined, self.base),
            parity_bits(pair[0], self.base) ^ parity_bits(pair[1], self.base),
        )

    def test_corrupt_combination_and_atomic_provenance(self):
        """A valid-looking congruence cannot replace the checked atom data."""
        pair = self.matching_pair()
        combined = combine_relations(pair, self.base).relation
        store = {atom.relation_id: atom for atom in pair}

        for corrupted in (
            replace(combined, u=(combined.u + 1) % 10403),
            replace(combined, square_correction=1),
            replace(combined, sign=-combined.sign),
            replace(combined, exponents=()),
            replace(combined, atom_ids=("unknown",)),
        ):
            with self.assertRaises(ValueError):
                verify_combined(corrupted, self.base, store)

        store[pair[0].relation_id] = replace(pair[0], sign=-pair[0].sign)
        with self.assertRaises(ValueError):
            verify_combined(combined, self.base, store)
        with self.assertRaises(ValueError):
            combine_relations((pair[0], pair[0]), self.base)
        with self.assertRaises(ValueError):
            combine_relations((pair[0],), self.base)
        with self.assertRaises(ValueError):
            combine_relations(
                (replace(pair[0], residual=4), pair[1]), self.base
            )

    def test_residual_gcd_and_nonunit_branches(self):
        """A matching residual that splits n surfaces a proper divisor."""
        atom = AtomicRelation(self.polynomial, 1, 1, ((2, 1),), 103)
        result = combine_relations((atom,), self.base)

        self.assertEqual(result.divisor, 103)
        self.assertIsNone(result.relation)
        base = build_factor_base(101, bound=20).factor_base
        polynomial = Polynomial(101, 1, 1, 101)
        atom = AtomicRelation(polynomial, 0, 1, ((2, 2), (5, 2)), 101)
        with self.assertRaisesRegex(ValueError, "nonunit"):
            combine_relations((atom,), base)

    def test_full_atom_combination_and_budget_refusal(self):
        """A full atom keeps correction one; verification is chargeable."""
        atom = next(atom for atom in self.atoms if atom.residual == 1)
        result = combine_relations((atom,), self.base)

        self.assertEqual(result.relation.square_correction, 1)
        self.assertEqual(result.relation.exponents, atom.exponents)
        with self.assertRaises(BudgetExhaustedError):
            verify_atomic(atom, self.base, budget=Budget(work_limit=0))
        with self.assertRaises(MemoryError):
            combine_relations((atom,), self.base, memory_bytes=0)
        with self.assertRaises(MemoryError):
            verify_combined(
                result.relation,
                self.base,
                {atom.relation_id: atom},
                memory_bytes=0,
            )

    def test_exact_composite_residual_is_rejected(self):
        """A correct integer identity cannot certify a composite residual."""
        for position in range(-30, 31):
            value = self.polynomial.u_value(position) ** 2 - self.base.n_prime
            factors = reference_factors(value)
            exponents = tuple(
                (p, e) for p, e in factors if p in self.base.primes
            )
            residual = 1
            for prime, exponent in factors:
                if prime not in self.base.primes:
                    residual *= prime**exponent
            if residual > 1 and len(reference_factors(residual)) == 1:
                if reference_factors(residual)[0][1] == 1:
                    continue
            if 1 < residual <= 1_000_000:
                atom = AtomicRelation(
                    self.polynomial,
                    position,
                    -1 if value < 0 else 1,
                    exponents,
                    residual,
                )
                with self.assertRaisesRegex(ValueError, "proven prime"):
                    verify_atomic(atom, self.base, residual_bound=1_000_000)
                return

        self.fail("fixture has no exact composite residual")

    def test_multiple_partial_pairs_use_modular_corrections(self):
        """Four residual occurrences contribute r squared and four atom IDs."""
        pair = self.matching_pair()
        shifted = Polynomial(10403, 1, 1, self.polynomial.b + 1)
        extras = tuple(
            AtomicRelation(
                shifted,
                atom.position - 1,
                atom.sign,
                atom.exponents,
                atom.residual,
            )
            for atom in pair
        )
        atoms = pair + extras
        result = combine_relations(atoms, self.base).relation

        self.assertEqual(
            result.square_correction, pair[0].residual ** 2 % 10403
        )
        self.assertEqual(len(result.atom_ids), 4)
        self.assertTrue(
            verify_combined(
                result, self.base, {atom.relation_id: atom for atom in atoms}
            )
        )
        with self.assertRaises(ValueError):
            combine_relations((pair[0],) * 257, self.base)


if __name__ == "__main__":
    unittest.main()
