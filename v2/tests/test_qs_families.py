"""Independent CRT residues, Gray roots, resource refusal and resume."""

import copy
import json
import unittest
from dataclasses import replace
from functools import partial
from math import prod

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.qs import (
    PolynomialFamily,
    QSJob,
    SieveCollector,
    SieveConfig,
    build_factor_base,
    collect_block,
    family_assignments,
    polynomial_roots,
)
from v2.qs.families import _checksum
from v2.tests.test_qs import unlimited_budget


def resumed_budget(checkpoint):
    """Carry consumed resources and allow charged cache reconstruction."""
    resources = checkpoint["payload"]["resources"]
    return Budget(
        work_limit=200_000_000,
        used=resources["work_used"],
        prior_wall=resources["wall_used"],
        prior_cpu=resources["cpu_used"],
        seconds=None,
        cpu_seconds=None,
    )


class FamilyTests(unittest.TestCase):
    """Families preserve every identity through Gray signs and recentering."""

    def base(self, multiplier=1):
        return build_factor_base(
            10403,
            multiplier=multiplier,
            bound=100,
            budget=unlimited_budget(),
        ).factor_base

    def test_reproducible_distinct_bounded_assignments(self):
        base = self.base()
        options = dict(factor_count=3, family_count=16, pool_size=8, seed=29)
        first = family_assignments(base, 64, **options)

        self.assertEqual(first, family_assignments(base, 64, **options))
        self.assertEqual(len(first), len(set(first)))
        self.assertEqual(len(first), 16)
        for values in first:
            self.assertEqual(values, tuple(sorted(set(values))))
            self.assertTrue(all(p in base.primes for p in values))
        one = family_assignments(
            base, 64, factor_count=3, family_count=64, pool_size=3
        )

        self.assertEqual(len(one), 1)

    def test_crt_gray_roots_against_full_and_exhaustive_oracles(self):
        saw_recenter = False

        for h in (1, 2, 3, 9):
            base = self.base(h)

            for count in (1, 2, 3, 4):
                values = family_assignments(
                    base,
                    64,
                    factor_count=count,
                    family_count=1,
                    budget=unlimited_budget(),
                )[0]
                family = PolynomialFamily(
                    base, values, budget=unlimited_budget()
                )
                seen, prior = set(), None

                for index in range(family.count):
                    step = family.next()
                    polynomial = step.polynomial

                    self.assertEqual(step.gray_index, index)
                    self.assertEqual(polynomial.a, prod(values))
                    self.assertEqual(
                        polynomial.b**2 % polynomial.a,
                        base.n_prime % polynomial.a,
                    )
                    self.assertNotIn(polynomial.b, seen)
                    seen.add(polynomial.b)
                    if prior is not None:
                        self.assertEqual(
                            step.delta_b, polynomial.b - prior.polynomial.b
                        )
                        gray = index ^ (index >> 1)
                        old_gray = (index - 1) ^ ((index - 1) >> 1)

                        self.assertEqual((gray ^ old_gray).bit_count(), 1)
                        bit = (gray ^ old_gray).bit_length() - 1
                        old_sign = -1 if old_gray & (1 << bit) else 1
                        raw_delta = -2 * old_sign * family.terms[bit + 1]
                        saw_recenter |= raw_delta != step.delta_b

                    for entry, roots in zip(base.entries, step.roots):
                        expected = polynomial_roots(
                            polynomial, base, entry, budget=unlimited_budget()
                        )

                        self.assertEqual(roots, expected)
                        exhaustive = tuple(
                            x
                            for x in range(entry.prime)
                            if polynomial.value(x) % entry.prime == 0
                        )
                        actual = (
                            tuple(range(entry.prime))
                            if (roots.all_positions)
                            else roots.roots
                        )

                        self.assertEqual(actual, exhaustive)

                    for x in (-64, -1, 0, 1, 64):
                        self.assertEqual(
                            polynomial.u_value(x) ** 2 - base.n_prime,
                            polynomial.a * polynomial.value(x),
                        )

                    prior = step

                used = family.budget.used

                self.assertIsNone(family.next())
                self.assertIsNone(family.next())
                self.assertEqual(family.budget.used, used)

        self.assertTrue(saw_recenter)

    def test_refused_next_keeps_gray_and_root_state(self):
        base = self.base()
        primes = family_assignments(base, 64, family_count=1)[0]
        family = PolynomialFamily(base, primes, budget=unlimited_budget())
        first = family.next()
        used = family.budget.used
        family.budget = Budget(
            work_limit=used, used=used, seconds=None, cpu_seconds=None
        )
        with self.assertRaises(BudgetExhaustedError):
            family.next()

        self.assertEqual(family.next_index, 1)
        self.assertIs(family.current, first)
        family.budget = Budget(
            work_limit=200_000_000, used=used, seconds=None, cpu_seconds=None
        )
        second = family.next()

        self.assertEqual(second.gray_index, 1)

    def test_checkpoint_roundtrip_and_charged_root_reconstruction(self):
        base = self.base(3)
        primes = family_assignments(base, 64, factor_count=4, family_count=1)[
            0
        ]
        family = PolynomialFamily(base, primes, budget=unlimited_budget())
        for _ in range(3):
            family.next()
        checkpoint = json.loads(json.dumps(family.checkpoint()))
        budget = resumed_budget(checkpoint)

        restored = PolynomialFamily.from_checkpoint(
            base, checkpoint, budget=budget
        )

        self.assertGreater(
            budget.used, checkpoint["payload"]["resources"]["work_used"]
        )
        self.assertEqual(restored.current, family.current)
        while family.next_index < family.count:
            self.assertEqual(restored.next(), family.next())

        self.assertIsNone(restored.next())
        with self.assertRaises(ValueError):
            PolynomialFamily.from_checkpoint(
                base, checkpoint, budget=unlimited_budget()
            )

    def test_checkpoint_corruption_and_invalid_progress(self):
        base = self.base()
        primes = family_assignments(base, 64, family_count=1)[0]
        family = PolynomialFamily(base, primes, budget=unlimited_budget())
        family.next()
        original = family.checkpoint()
        corrupt = copy.deepcopy(original)
        corrupt["payload"]["next_index"] += 1
        with self.assertRaises(ValueError):
            PolynomialFamily.from_checkpoint(
                base, corrupt, budget=resumed_budget(original)
            )
        for key, value in (
            ("version", 3),
            ("next_index", 999),
            ("base_identity", "wrong"),
            ("a_primes", [7, 7]),
        ):
            corrupt = copy.deepcopy(original)
            corrupt["payload"][key] = value
            corrupt["sha256"] = _checksum(corrupt["payload"])
            with self.assertRaises(ValueError):
                PolynomialFamily.from_checkpoint(
                    base, corrupt, budget=resumed_budget(original)
                )

    def test_checkpoint_rejects_unbounded_fields_before_hashing(self):
        base = self.base()
        primes = family_assignments(base, 64, family_count=1)[0]
        family = PolynomialFamily(base, primes, budget=unlimited_budget())
        original = family.checkpoint()

        for key, value in (
            ("a_primes", [7] * 4097),
            ("base_identity", "x" * 4097),
            ("next_index", 1 << 200),
            (
                "resources",
                {"work_used": 1 << 200, "wall_used": 0, "cpu_used": 0},
            ),
            (
                "resources",
                {"work_used": 0, "wall_used": float("nan"), "cpu_used": 0},
            ),
        ):
            corrupt = copy.deepcopy(original)
            corrupt["payload"][key] = value
            with self.assertRaises(ValueError):
                PolynomialFamily.from_checkpoint(
                    base, corrupt, budget=unlimited_budget()
                )

        family.budget = Budget(work_limit=1 << 200, used=1 << 199)
        with self.assertRaises(ValueError):
            family.checkpoint()

    def test_family_input_memory_and_work_caps(self):
        base = self.base(3)
        for values in ((), (7, 7), (2,), (3,), (99991,), (17, 7)):
            with self.assertRaises(ValueError):
                PolynomialFamily(base, values)
        with self.assertRaises(TypeError):
            PolynomialFamily(base, [7, 13])
        values = family_assignments(base, 64, family_count=1)[0]
        with self.assertRaises(MemoryError):
            PolynomialFamily(base, values, memory_bytes=1)
        with self.assertRaises(MemoryError):
            family_assignments(base, 64, memory_bytes=1)
        with self.assertRaises(BudgetExhaustedError):
            PolynomialFamily(base, values, budget=Budget(work_limit=0))
        for kwargs in (
            {"factor_count": 9},
            {"family_count": 65},
            {"pool_size": 129},
        ):
            with self.assertRaises(ValueError):
                family_assignments(base, 64, **kwargs)

    def test_family_polynomials_use_existing_checked_extraction(self):
        n = 2003 * 8353
        budget = unlimited_budget()
        base = build_factor_base(n, bound=200, budget=budget).factor_base
        primes = family_assignments(base, 256, family_count=1, budget=budget)[
            0
        ]
        family = PolynomialFamily(base, primes, budget=budget)
        splits = []

        while (step := family.next()) is not None:
            job = QSJob(
                step.polynomial,
                base,
                -256,
                257,
                budget=budget,
                config=SieveConfig(residual_bound=1, memory_bytes=33554432),
                collector_class=partial(
                    SieveCollector, precomputed_roots=step.roots
                ),
            )

            result = job.run()
            if result.divisor is not None:
                self.assertEqual(result.divisor * result.cofactor, n)
                self.assertEqual(
                    {result.divisor, result.cofactor}, {2003, 8353}
                )
                splits.append(result.divisor)

        self.assertTrue(splits)

    def test_cached_collector_against_exhaustive_payloads(self):
        base = self.base(9)
        primes = family_assignments(base, 64, factor_count=3, family_count=1)[
            0
        ]
        family = PolynomialFamily(base, primes, budget=unlimited_budget())

        while (step := family.next()) is not None:
            reference = collect_block(
                step.polynomial,
                base,
                -31,
                34,
                residual_bound=97,
                budget=unlimited_budget(),
                memory_bytes=33554432,
            )
            worker = SieveCollector(
                step.polynomial,
                base,
                precomputed_roots=step.roots,
                budget=unlimited_budget(),
                config=SieveConfig(residual_bound=97, memory_bytes=33554432),
            )

            run = worker.collect(-31, 34)

            def payload(atoms):
                """Consume exact payloads independently of score internals."""
                return {
                    atom.position: (atom.sign, atom.exponents, atom.residual)
                    for atom in atoms
                }

            self.assertEqual(payload(run.atoms), payload(reference.relations))

    def test_cached_collector_rejects_wrong_or_missing_roots(self):
        base = self.base()
        primes = family_assignments(base, 64, family_count=1)[0]
        step = PolynomialFamily(base, primes, budget=unlimited_budget()).next()
        changed = list(step.roots)
        index = next(
            i
            for i, item in enumerate(changed)
            if item.prime > 2 and len(item.roots) == 2
        )
        changed[index] = replace(
            changed[index], roots=changed[index].roots[:1]
        )

        for roots in (tuple(changed), step.roots[:-1], list(step.roots)):
            with self.assertRaises(ValueError):
                SieveCollector(
                    step.polynomial,
                    base,
                    precomputed_roots=roots,
                    budget=unlimited_budget(),
                )

        wrong_type = (replace(step.roots[0], prime=2.0),) + step.roots[1:]
        with self.assertRaises(TypeError):
            SieveCollector(
                step.polynomial,
                base,
                precomputed_roots=wrong_type,
                budget=unlimited_budget(),
            )
