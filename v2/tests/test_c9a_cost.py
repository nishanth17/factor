"""Amortized panel lookup and bounded lookahead refusal regressions."""

import unittest
from collections import Counter
from dataclasses import asdict
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.qs.c9a import extended
from v2.common import utils
from v2.execution.budget import BudgetExhaustedError
from v2.qs.cofactor_policy import CofactorBatch, split_cofactor, split_residual
from v2.qs.linear_algebra import (
    DependencySolver,
    filter_matrix,
    verify_dependency,
)
from v2.qs.structured_matrix import panel_reduce
from v2.tests.test_qs_dlp import allowance
from v2.tests.test_qs_pipeline import dense_kernel, span


class CountingPivots(dict):
    scans = 0

    def __iter__(self):
        self.scans += 1
        return super().__iter__()


class ChunkAllowance:
    def __init__(self, limit=10000):
        self.limit, self.used, self.polls = limit, 0, []

    def consume(self, amount=1):
        if self.used + amount > self.limit:
            raise BudgetExhaustedError()
        self.polls.append(amount)
        self.used += amount


def lookahead_worker(limit=10000):
    owner = SimpleNamespace(
        budget=ChunkAllowance(limit),
        config=SimpleNamespace(
            smooth_batch_size=0, division="bucket", large_product_bound=100
        ),
        _scores=[int(i in (10, 70)) for i in range(200)],
        _piecewise_passes=lambda offset: True,
    )
    owner._divide_double = (
        lambda worker, x, value, residual, factors, stats, hint=None: (
            x,
            None,
        )
    )
    owner._divide = lambda x, offset, stats: owner._divide_double(
        owner, x, x + 1, 1, (), stats
    )
    return owner


class CostRegressionTests(unittest.TestCase):
    def test_cost_profile_changes_only_the_measured_expensive_choices(self):
        old = extended.configuration(
            60, admission_base_bound=103814, profile="bounded-report"
        )
        new = extended.configuration(
            60, admission_base_bound=103814, profile="cost-report"
        )
        expected = asdict(old)
        expected["collector"].update(
            solver_panel_bits=0, split_policy="rho-first-v3"
        )
        self.assertEqual(asdict(new), expected)

    def test_rho_first_short_circuits_squfof_and_keeps_its_fallback(self):
        checked = SimpleNamespace(
            prime=lambda p, budget: (
                utils.classify_prime(p) is utils.Primality.PROVEN
            )
        )
        events = []

        def rho(*args, **kwargs):
            events.append("rho")
            return None

        def squfof(*args, **kwargs):
            events.append("squfof")
            return 1009

        with (
            patch("v2.qs.cofactor_policy.factorize_rho", side_effect=rho),
            patch("v2.qs.cofactor_policy.squfof", side_effect=squfof),
        ):
            pair, _, _ = split_cofactor(
                1009 * 1000003, 2000000, allowance(), checked, "rho-first-v3"
            )
        self.assertEqual(pair, (1009, 1000003))
        self.assertEqual(events, ["rho", "rho", "squfof"])

        with (
            patch("v2.qs.cofactor_policy.factorize_rho", return_value=17),
            patch(
                "v2.qs.cofactor_policy.squfof",
                side_effect=AssertionError("unneeded fallback"),
            ),
        ):
            pair, _, _ = split_cofactor(
                17 * 19, 100, allowance(), checked, "rho-first-v3"
            )
        self.assertEqual(pair, (17, 19))

    def test_three_prime_recovery_reuses_a_paid_proper_split(self):
        checked = SimpleNamespace(
            prime=lambda p, budget: (
                utils.classify_prime(p) is utils.Primality.PROVEN
            )
        )
        inputs = []

        def rho(number, **kwargs):
            inputs.append(number)
            return {1353: 3, 451: 11}[number]

        with patch("v2.qs.cofactor_policy.factorize_rho", side_effect=rho):
            primes, _, _ = split_residual(
                1353, 100, allowance(), checked, "rho-first-v3", dispatch=False
            )
        self.assertEqual(primes, (3, 11, 41))
        self.assertEqual(inputs, [1353, 451])

    def test_panel_presence_avoids_scanning_or_sorting_pivots(self):
        pivots = CountingPivots(
            {1: (1, 1), 2: (3, 2), 8: (10, 4), 32: (34, 8), 128: (130, 16)}
        )
        solver = SimpleNamespace(
            pivots=pivots,
            _panel_pivots=171,
            pivot="highest",
            budget=allowance(),
            matrix=SimpleNamespace(
                stats={"panel_bits": 4, "working_columns": 8},
                original_rows=tuple(range(8)),
            ),
        )
        expected = panel_reduce(solver, 171, 7)
        table = solver._panel_cache
        for _ in range(100):
            self.assertEqual(panel_reduce(solver, 171, 7), expected)
        self.assertEqual(pivots.scans, 0)

        pivots[512] = (512, 32)
        solver._panel_pivots |= 512
        self.assertEqual(panel_reduce(solver, 171, 7), expected)
        self.assertEqual(pivots.scans, 0)
        self.assertIs(solver._panel_cache, table)

    def test_noncontiguous_panels_match_independent_dense_kernels(self):
        rows = (3, 9, 10, 129, 136, 130, 4097, 4104, 4098, 4224)
        expected = span(dense_kernel(rows))
        for pivot in ("highest", "lowest"):
            matrix = filter_matrix(rows, panel_bits=3, budget=allowance())
            solver = DependencySolver(matrix, pivot=pivot, budget=allowance())
            dependencies = solver.run()
            self.assertEqual(solver._panel_pivots, sum(solver.pivots))
            self.assertEqual(span(dependencies), expected)
            for mask in dependencies:
                verify_dependency(mask, rows)

    def test_lookahead_prepays_a_chunk_without_reordering_candidates(self):
        owner, batch = lookahead_worker(), CofactorBatch(2)
        self.assertEqual(
            batch.recover(owner, 10, 10, 200, 1, Counter()), (10, None)
        )
        self.assertEqual(owner.budget.polls, [64])
        self.assertEqual(
            batch.recover(owner, 70, 70, 200, 1, Counter()), (70, None)
        )
        self.assertEqual(owner.budget.polls, [64])

    def test_refused_lookahead_publishes_no_records_or_divider_change(self):
        owner, batch = lookahead_worker(limit=63), CofactorBatch(2)
        divider = owner._divide_double
        with self.assertRaises(BudgetExhaustedError):
            batch.recover(owner, 10, 10, 200, 1, Counter())
        self.assertFalse(batch.records)
        self.assertIs(owner._divide_double, divider)
        self.assertEqual(owner.budget.used, 0)


if __name__ == "__main__":
    unittest.main()
