"""Sparse filtering cascades, bounded work and exact lifted kernels."""

import random
import unittest

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.qs.linear_algebra import (
    DependencySolver,
    filter_matrix,
    verify_dependency,
)
from v2.tests.test_qs_pipeline import dense_kernel, span


def allowance(work=10**12):
    """Finite generous work without unrelated timer noise."""
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


class IncrementalFilterTests(unittest.TestCase):
    """Independent kernels cover incidence updates and stale queues."""

    def test_long_cycle_finishes_inside_sparse_work_allowance(self):
        count = 1200
        rows = tuple(
            (1 << index) | (1 << ((index + 1) % count))
            for index in range(count)
        )
        budget = allowance(20_000_000)
        matrix = filter_matrix(
            rows, weight_two=True, budget=budget, memory_bytes=8 * 1024**2
        )

        self.assertIs(matrix.original_rows, rows)
        self.assertEqual(matrix.rows, ())
        self.assertEqual(matrix.zero_dependencies, ((1 << count) - 1,))
        self.assertEqual(matrix.stats["weight_two_merges"], count - 1)
        self.assertTrue(verify_dependency(matrix.zero_dependencies[0], rows))

    def test_mixed_cascades_and_queue_transitions_match_dense_kernel(self):
        generator = random.Random(310431)
        fixtures = [(0, 3, 5, 6, 1, 9), (0, 0, 1, 3, 6, 12), (3, 3, 5, 5)]
        fixtures += [
            tuple(generator.randrange(256) for _ in range(9))
            for _ in range(40)
        ]
        for rows in fixtures:
            expected = dense_kernel(rows)

            for weight_two in (False, True):
                matrix = filter_matrix(
                    rows, weight_two=weight_two, budget=allowance()
                )
                for pivot in ("highest", "lowest"):
                    solver = DependencySolver(
                        matrix, pivot=pivot, budget=allowance()
                    )

                    self.assertEqual(span(solver.run()), expected)

    def test_refused_private_filter_preserves_input_and_can_retry(self):
        rows = (3, 5, 6, 24, 40, 48)
        completed = allowance()
        filter_matrix(rows, weight_two=True, budget=completed)
        for work in (0, completed.used // 2, completed.used - 1):
            with self.assertRaises(BudgetExhaustedError):
                filter_matrix(rows, weight_two=True, budget=allowance(work))

            self.assertEqual(rows, (3, 5, 6, 24, 40, 48))
        matrix = filter_matrix(rows, weight_two=True, budget=allowance())

        self.assertEqual(span(matrix.zero_dependencies), dense_kernel(rows))

    def test_cancellation_during_incidence_build_preserves_input(self):
        calls = 0

        def cancelled():
            nonlocal calls
            calls += 1
            return calls >= 3

        rows = (3,) * 1000
        budget = Budget(
            work_limit=10**12,
            seconds=None,
            cpu_seconds=None,
            cancelled=cancelled,
        )
        with self.assertRaises(BudgetExhaustedError):
            filter_matrix(rows, weight_two=True, budget=budget)

        self.assertEqual(budget.reason, "cancelled")
        self.assertEqual(calls, 3)
        self.assertEqual(rows, (3,) * 1000)

    def test_memory_refusal_precedes_incidence_and_work(self):
        rows = (3, 5, 6)
        budget = allowance()
        with self.assertRaises(MemoryError):
            filter_matrix(rows, weight_two=True, budget=budget, memory_bytes=1)
        self.assertEqual(budget.used, 0)


if __name__ == "__main__":
    unittest.main()
