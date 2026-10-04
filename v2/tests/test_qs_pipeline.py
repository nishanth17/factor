"""Independent dense kernels, provenance corruption and complete QS splits."""

import gc
import random
import unittest
import weakref
from dataclasses import replace
from unittest.mock import patch

from v2.budget import Budget, BudgetExhaustedError
from v2.qs import (
    DependencyExtractor,
    DependencySolver,
    QSJob,
    SieveCollector,
    SieveConfig,
    build_factor_base,
    extract_dependency,
    filter_matrix,
    prepare_relations,
    qs_polynomial,
    verify_dependency,
)
from v2.tests.test_qs import unlimited_budget
from v2.tests.test_qs_sieve import collector


def dense_kernel(rows):
    """Enumerate selection vectors using independent boolean column sums."""
    columns = max((row.bit_length() for row in rows), default=0)
    return {
        mask
        for mask in range(1 << len(rows))
        if all(
            sum(
                bool(rows[index] & (1 << column))
                for index in range(len(rows))
                if mask & (1 << index)
            )
            % 2
            == 0
            for column in range(columns)
        )
    }


def span(masks):
    """Enumerate the finite span of a small computed basis."""
    values = {0}
    for mask in masks:
        values |= {value ^ mask for value in tuple(values)}
    return values


class MatrixTests(unittest.TestCase):
    """Filtering and pivot kernels agree with brute-force column equations."""

    def test_dense_kernel_random_and_cascades(self):
        generator = random.Random(331)
        fixtures = [(0, 0), (1, 3, 6, 4), (3, 5, 6), (1, 1, 2, 2), (3, 3), ()]
        fixtures += [
            tuple(generator.randrange(64) for _ in range(8)) for _ in range(60)
        ]
        for rows in fixtures:
            expected = dense_kernel(rows)
            for weight_two in (False, True):
                matrix = filter_matrix(
                    rows,
                    weight_two=weight_two,
                    budget=unlimited_budget(),
                )
                for pivot in ("highest", "lowest"):
                    solver = DependencySolver(
                        matrix,
                        pivot=pivot,
                        budget=unlimited_budget(),
                    )
                    self.assertEqual(span(solver.run()), expected)

    def test_weight_two_zero_rows_and_singleton_cascade(self):
        matrix = filter_matrix((3, 5, 6), weight_two=True)
        self.assertTrue(matrix.zero_dependencies)
        self.assertEqual(matrix.stats["weight_two_merges"], 2)
        cascade = filter_matrix((1, 3, 6, 12))
        self.assertEqual(cascade.rows, ())
        self.assertEqual(cascade.stats["singletons_removed"], 4)
        self.assertGreater(cascade.stats["rounds"], 1)

    def test_solver_refusal_and_resume(self):
        matrix = filter_matrix((3, 5, 6))
        solver = DependencySolver(matrix, budget=Budget(work_limit=0))
        with self.assertRaises(BudgetExhaustedError):
            solver.run()
        self.assertEqual(solver.next_row, 0)
        solver.budget = unlimited_budget()
        self.assertEqual(
            span(solver.run()), dense_kernel(matrix.original_rows)
        )
        self.assertEqual(solver.run(), tuple(solver.dependencies))

    def test_limits_and_bad_masks(self):
        for mask in (0, 8, 1):
            with self.assertRaises(ValueError):
                verify_dependency(mask, (3, 5, 6))
        with self.assertRaises(MemoryError):
            filter_matrix((3, 5, 6), memory_bytes=1)
        with self.assertRaises(ValueError):
            filter_matrix((1 << 100001,))
        with self.assertRaises(ValueError):
            filter_matrix((0,) * 4097)


class ExtractionTests(unittest.TestCase):
    """Checked congruences retain signs, exact exponents and corrections."""

    def setUp(self):
        self.base = build_factor_base(10403, bound=40).factor_base
        self.polynomial = qs_polynomial(self.base)
        self.run = collector(self.polynomial, self.base).collect(-128, 129)
        self.store = {atom.relation_id: atom for atom in self.run.atoms}
        self.relations = self.run.full_relations + self.run.combined_relations

    def test_distinct_equal_parity_and_exact_duplicates(self):
        prepared = prepare_relations(
            self.relations + self.relations[:1],
            self.base,
            self.store,
            budget=unlimited_budget(),
        )
        self.assertEqual(prepared.duplicate_indices, (len(self.relations),))
        self.assertEqual(len(prepared.relations), len(self.relations))
        groups = {}
        for index, row in enumerate(prepared.rows):
            groups.setdefault(row, []).append(index)
        equal = next(
            indices[:2] for indices in groups.values() if len(indices) >= 2
        )
        mask = sum(1 << index for index in equal)
        congruence = extract_dependency(
            prepared, mask, budget=unlimited_budget()
        )
        self.assertEqual(
            congruence.x**2 % self.base.n, congruence.y**2 % self.base.n
        )

    def test_combined_corrections_and_signed_dependencies(self):
        prepared = prepare_relations(
            self.relations, self.base, self.store, budget=unlimited_budget()
        )
        solver = DependencySolver(
            filter_matrix(prepared.rows), budget=unlimited_budget()
        )
        masks = solver.run()
        self.assertTrue(masks)
        divisors = []
        for mask in masks:
            result = extract_dependency(
                prepared, mask, budget=unlimited_budget()
            )
            selected = [
                row
                for index, row in enumerate(prepared.relations)
                if mask & (1 << index)
            ]
            # Independently reconstruct the square side from full signed
            # integer relation values, including every residual correction.
            from math import isqrt, prod

            product = prod(
                row.sign
                * prod(p**e for p, e in row.exponents)
                * getattr(row, "square_correction", 1) ** 2
                for row in selected
            )
            square = isqrt(product)
            self.assertEqual(square * square, product)
            self.assertEqual(square % self.base.n, result.y)
            if result.divisor:
                divisors.append(result.divisor)
        self.assertTrue(divisors)

    def test_corrupted_provenance_and_fields(self):
        combined = self.run.combined_relations[0]
        missing = dict(self.store)
        del missing[combined.atom_ids[0]]
        for relation, store in (
            (combined, missing),
            (
                replace(
                    combined, square_correction=combined.square_correction + 1
                ),
                self.store,
            ),
            (
                replace(
                    self.run.full_relations[0],
                    sign=-self.run.full_relations[0].sign,
                ),
                self.store,
            ),
        ):
            with self.assertRaises(ValueError):
                prepare_relations((relation,), self.base, store)
        partial = next(atom for atom in self.run.atoms if atom.residual != 1)
        with self.assertRaises(ValueError):
            prepare_relations((partial,), self.base, self.store)

    def test_prepared_rows_pin_atomic_provenance(self):
        """Releasing caller dictionaries retains checked source atoms."""
        store = dict(self.store)
        prepared = prepare_relations(
            self.relations, self.base, store, budget=unlimited_budget()
        )
        store.clear()
        identities = {atom.relation_id for atom in prepared.atoms}
        for relation in prepared.relations:
            if hasattr(relation, "atom_ids"):
                self.assertTrue(set(relation.atom_ids) <= identities)

    def test_trivial_trial_then_factor_and_refusal(self):
        prepared = prepare_relations(
            self.relations, self.base, self.store, budget=unlimited_budget()
        )
        masks = DependencySolver(
            filter_matrix(prepared.rows), budget=unlimited_budget()
        ).run()
        trivial = next(
            mask
            for mask in masks
            if extract_dependency(
                prepared, mask, budget=unlimited_budget()
            ).divisor
            is None
        )
        useful = next(
            mask
            for mask in masks
            if extract_dependency(
                prepared, mask, budget=unlimited_budget()
            ).divisor
            is not None
        )
        extractor = DependencyExtractor(
            prepared, (trivial, useful), budget=Budget(work_limit=0)
        )
        with self.assertRaises(BudgetExhaustedError):
            extractor.run()
        self.assertEqual(extractor.next_dependency, 0)
        extractor.budget = unlimited_budget()
        divisor = extractor.run()
        self.assertEqual(extractor.next_dependency, 2)
        self.assertEqual(self.base.n % divisor, 0)
        self.assertGreater(divisor, 1)
        self.assertLess(divisor, self.base.n)


class PipelineTests(unittest.TestCase):
    """Small balanced runs reconstruct input and preserve refused state."""

    def test_balanced_fixtures_and_stopping_controls(self):
        for p, q in ((101, 137), (211, 307), (1009, 1237), (4001, 5003)):
            n = p * q
            base = build_factor_base(n, bound=100).factor_base
            for weight_two in (False, True):
                for row_excess in (0, 2):
                    job = QSJob(
                        qs_polynomial(base),
                        base,
                        -256,
                        513,
                        weight_two=weight_two,
                        row_excess=row_excess,
                        batch_width=64,
                        budget=unlimited_budget(),
                        config=SieveConfig(memory_bytes=32 * 1024 * 1024),
                    )
                    result = job.run()
                    self.assertEqual(result.reason, "factor_found", result)
                    self.assertEqual(result.divisor * result.cofactor, n)
                    self.assertEqual({result.divisor, result.cofactor}, {p, q})
                    self.assertEqual(job.run(), result)

    def test_collection_and_solver_resume(self):
        base = build_factor_base(101 * 137, bound=100).factor_base
        job = QSJob(
            qs_polynomial(base),
            base,
            -128,
            129,
            budget=unlimited_budget(),
            config=SieveConfig(residual_bound=1),
        )
        job.budget = Budget(work_limit=0)
        stopped = job.run()
        self.assertEqual(stopped.reason, "work_limit")
        self.assertEqual(stopped.cofactor, base.n)
        self.assertEqual(stopped.next_position, -128)
        job.budget = unlimited_budget()
        with patch(
            "v2.qs.pipeline.DependencySolver.run",
            side_effect=BudgetExhaustedError,
        ):
            job.budget.reason = "work_limit"
            stopped = job.run()
        self.assertEqual(stopped.reason, "work_limit")
        self.assertIsNotNone(job.solver)
        consumed = job.budget.used
        job.budget = Budget(
            work_limit=100000000, used=consumed, seconds=None, cpu_seconds=None
        )
        result = job.run()
        self.assertEqual(result.reason, "factor_found")
        self.assertEqual(result.divisor * result.cofactor, base.n)

    def test_empty_and_finite_unsuccessful_window(self):
        base = build_factor_base(101 * 137, bound=40).factor_base
        for lo, hi in ((0, 0), (999, 1000)):
            job = QSJob(
                qs_polynomial(base), base, lo, hi, budget=unlimited_budget()
            )
            result = job.run()
            self.assertEqual(result.reason, "window_exhausted")
            self.assertEqual(result.cofactor, base.n)
            self.assertEqual(job.run(), result)

    def test_final_solve_of_unchanged_deferred_store(self):
        """Window exhaustion still solves a store deferred for low excess."""
        base = build_factor_base(10403, bound=40).factor_base
        polynomial = qs_polynomial(base)
        run = collector(polynomial, base).collect(-128, 129)
        prepared = prepare_relations(
            run.full_relations + run.combined_relations,
            base,
            {atom.relation_id: atom for atom in run.atoms},
            budget=unlimited_budget(),
        )
        groups = {}
        for index, row in enumerate(prepared.rows):
            if row:
                groups.setdefault(row, []).append(index)
        indices = next(
            values[:2] for values in groups.values() if len(values) >= 2
        )
        job = QSJob(
            polynomial, base, 0, 0, row_excess=4096, budget=unlimited_budget()
        )
        job.collector._atoms = {atom.relation_id: atom for atom in run.atoms}
        job.collector._full = [prepared.relations[index] for index in indices]
        self.assertIsNone(job._solve(False))
        self.assertEqual(job.stats["solve_calls"], 0)
        job._solve(True)
        self.assertEqual(job.stats["solve_calls"], 1)

    def test_batch_snapshots_released_before_more_collection(self):
        """Collection cannot pin previous snapshots outside its reservation."""

        class ObservedCollector(SieveCollector):
            """Check lifetime independently with weak references."""

            previous = None

            def collect(self, lo, hi):
                """Require release before allocating the next batch."""
                # PyPy uses tracing GC; weakrefs do not clear immediately
                # after the last strong reference is dropped.
                gc.collect()
                if self.previous is not None and self.previous() is not None:
                    raise AssertionError("prior collection snapshot is pinned")
                result = super().collect(lo, hi)
                self.previous = weakref.ref(result)
                return result

        base = build_factor_base(104729, bound=100).factor_base
        job = QSJob(
            qs_polynomial(base),
            base,
            -64,
            65,
            batch_width=16,
            config=SieveConfig(residual_bound=1),
            collector_class=ObservedCollector,
            budget=unlimited_budget(),
        )
        result = job.run()
        self.assertEqual(result.reason, "window_exhausted")
