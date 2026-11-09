"""Parent CPU publication shares the reservation's one checked boundary."""

import unittest
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import Budget
from v2.execution.budget import BudgetExhaustedError
from v2.qs.parallel import CollectionPool, _ParentBudget


def pool_view():
    """Exercise the real accounting methods without launching any process."""
    pool = CollectionPool.__new__(CollectionPool)
    pool.mode = "process"
    pool.executor = SimpleNamespace(_processes={})
    pool.accounted_cpu = 0.0
    pool.state = dict(
        cpu=[3.0, 5.0],
        cpu_base=[1.0, 2.0],
        parent_cpu=SimpleNamespace(value=0.0),
    )
    return pool


class ParentBudgetBoundaryTests(unittest.TestCase):
    def test_public_poll_still_checks_limits_after_accounting_cpu(self):
        pool = pool_view()
        budget = Budget(work_limit=100, cancelled=lambda: True)
        with self.assertRaises(BudgetExhaustedError):
            pool.poll(budget)
        self.assertEqual(budget.reason, "cancelled")
        self.assertEqual(budget.used, 0)
        self.assertEqual(budget.prior_cpu, 0.0)

    def test_one_raw_check_per_boundary_retains_every_work_unit(self):
        pool = pool_view()
        budget = Budget(work_limit=100)
        view = _ParentBudget(budget, pool)
        with patch.object(budget, "consume", wraps=budget.consume) as checked:
            for _ in range(64):
                view.consume(1)
            self.assertEqual(checked.call_count, 1)
            self.assertEqual(budget.used, 64)
            view.consume(1)
            self.assertEqual(checked.call_count, 2)
        self.assertEqual(budget.used, 65)
        self.assertEqual(budget.prior_cpu, 0.0)

    def test_exhausted_limits_refuse_before_any_work_is_committed(self):
        for options, reason in (
            (dict(work_limit=0), "work_limit"),
            (dict(work_limit=100, cpu_seconds=0), "cpu_limit"),
            (dict(work_limit=100, seconds=0), "wall_limit"),
            (dict(work_limit=100, cancelled=lambda: True), "cancelled"),
        ):
            budget = Budget(**options)
            with self.assertRaises(BudgetExhaustedError):
                _ParentBudget(budget, pool_view()).consume(1)
            self.assertEqual(budget.reason, reason)
            self.assertEqual(budget.used, 0)

    def test_cancellation_between_boundaries_keeps_the_same_action_bound(self):
        stop = [False]
        budget = Budget(work_limit=100, cancelled=lambda: stop[0])
        view = _ParentBudget(budget, pool_view())
        view.consume(1)
        stop[0] = True
        for _ in range(63):
            view.consume(1)
        with self.assertRaises(BudgetExhaustedError):
            view.consume(1)
        self.assertEqual(budget.used, 64)
        self.assertEqual(budget.reason, "cancelled")
