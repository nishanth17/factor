"""Exact reuse and bounded polling across the complete SIQS execution path."""

import unittest
from dataclasses import replace
from math import prod
from types import SimpleNamespace
from unittest.mock import patch

from v2.execution.budget import BudgetExhaustedError
from v2.execution.work_budget import PollingBudget
from v2.qs import (
    Polynomial,
    SieveCollector,
    SieveConfig,
    SIQSJob,
    build_factor_base,
)
from v2.qs.residual_kernel import MAX_RESIDUAL, ResidualKernel
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_c9a_extended import policy
from v2.tests.test_siqs import allowance, configuration


class SwitchReuseTests(unittest.TestCase):
    def make_job(self):
        job = SIQSJob(
            4001 * 5003,
            seed=73,
            config=configuration(),
            budget=allowance(),
        )
        job._setup()
        job._activate()
        return job

    def test_same_a_reuses_checked_factors_with_identical_collection(self):
        jobs = [self.make_job() for _ in range(2)]
        outputs = []
        charges = []
        for index, job in enumerate(jobs):
            worker = job.engine.collector
            worker.collect(-64, 0)
            step = job.family.next()
            previous = worker._a_exponents
            self.assertEqual(
                step.polynomial.supported_a, worker._a_factor_value
            )
            if index:
                worker._a_factor_value = None

            before = job.budget.used
            worker.set_polynomial(step.polynomial, step.roots)
            charges.append(job.budget.used - before)
            factors = worker._a_exponents
            self.assertEqual(
                prod(p**e for p, e in zip(job.base.primes, factors)),
                step.polynomial.supported_a,
            )
            if not index:
                self.assertIs(factors, previous)
            outputs.append(worker.collect(-64, 65))

        self.assertLess(charges[0], charges[1])
        self.assertEqual(outputs[0], outputs[1])

    def test_cached_a_does_not_skip_root_validation_or_publish_on_refusal(
        self,
    ):
        job = self.make_job()
        worker = job.engine.collector
        step = job.family.next()
        old = worker.polynomial, worker._roots, worker._a_exponents
        roots = list(step.roots)
        item = next(i for i, r in enumerate(roots) if r.prime > 2 and r.roots)
        roots[item] = replace(roots[item], roots=())
        with self.assertRaises(ValueError):
            worker.set_polynomial(step.polynomial, tuple(roots))
        self.assertEqual(
            (worker.polynomial, worker._roots, worker._a_exponents), old
        )

        worker.budget = Budget(work_limit=0, seconds=None, cpu_seconds=None)
        with self.assertRaises(BudgetExhaustedError):
            worker.set_polynomial(step.polynomial, step.roots)
        self.assertEqual(
            (worker.polynomial, worker._roots, worker._a_exponents), old
        )

        worker.budget = allowance()
        worker.set_polynomial(step.polynomial, step.roots)
        self.assertIs(worker._a_exponents, old[2])


class WholeJobPollingTests(unittest.TestCase):
    def test_callback_updates_reach_the_original_ledger(self):
        ledger = allowance()
        view = PollingBudget(ledger)
        view.consume(1)
        view.cancelled = lambda: True
        self.assertIs(view.cancelled, ledger.cancelled)
        with self.assertRaises(BudgetExhaustedError):
            view.consume(0)
        self.assertEqual(ledger.reason, "cancelled")
        self.assertEqual(ledger.used, 1)

    def test_setup_and_switches_share_the_bounded_view_then_restore(self):
        original = allowance()
        job = SIQSJob(104729, config=configuration(), budget=original)
        setup, activate = job._setup, job._activate
        stages = []

        def checked_setup():
            self.assertIsInstance(job.budget, PollingBudget)
            stages.append("setup")
            return setup()

        def checked_activate():
            self.assertIsInstance(job.budget, PollingBudget)
            stages.append("switch")
            return activate()

        with (
            patch.object(job, "_setup", checked_setup),
            patch.object(job, "_activate", checked_activate),
        ):
            job.run(max_blocks=10)

        self.assertIn("setup", stages)
        self.assertIn("switch", stages)
        self.assertIs(job.budget, original)
        self.assertIs(job.engine.budget, original)
        self.assertIs(job.engine.collector.budget, original)
        self.assertIs(job.family.budget, original)

    def test_setup_cancellation_has_the_same_finite_action_bound(self):
        calls = 0

        def cancelled():
            nonlocal calls
            calls += 1
            return calls >= 2

        original = Budget(
            work_limit=10000,
            seconds=None,
            cpu_seconds=None,
            cancelled=cancelled,
        )
        job = SIQSJob(104729, config=configuration(), budget=original)

        def setup():
            for _ in range(1000):
                job.budget.consume(7)

        with patch.object(job, "_setup", setup):
            result = job.run(max_blocks=1)

        self.assertEqual(result.reason, "cancelled")
        self.assertLessEqual(original.used, 64 * 7)
        self.assertEqual(result.cofactor, job.n)
        self.assertIs(job.budget, original)

    def test_work_refusal_and_replacement_budget_keep_resume_checked(self):
        job = SIQSJob(
            4001 * 5003,
            config=configuration(),
            budget=Budget(work_limit=100, seconds=None, cpu_seconds=None),
        )
        result = job.run(max_blocks=1)
        self.assertEqual(result.reason, "work_limit")
        self.assertLessEqual(job.budget.used, 100)

        restored = SIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )
        result = restored.run()
        self.assertEqual(result.reason, "factor_found")
        self.assertEqual(result.divisor * result.cofactor, job.n)


class CompactDerivedStateTests(unittest.TestCase):
    def test_encoded_prime_and_opaque_coordinates_remain_distinct(self):
        kernel = ResidualKernel()
        config = policy(cycle_reduction_choices=4)
        cycles = []
        for handle, token in enumerate(
            (
                (False, MAX_RESIDUAL),
                (True, 4),
                (False, MAX_RESIDUAL),
                (True, 4),
            )
        ):
            atom = SimpleNamespace(
                incidence_tokens=(token,), residual_primes=(), exponents=()
            )
            plan = kernel.plan_atom(atom, handle, config, allowance())
            if plan.cycle is not None:
                cycles.append(plan.cycle)
            kernel.commit_atom(plan, atom, config.cycle_reduction_choices)
        self.assertEqual(cycles, [(0, 2), (1, 3)])
        state = kernel.state()
        self.assertEqual(
            [row[0] for row in state["pivots"]],
            [[False, MAX_RESIDUAL], [True, 4]],
        )

    def test_packed_queue_preserves_large_singular_power_weights(self):
        b = 2**298 - 1
        n = b * b - 2**300
        polynomial = Polynomial(n, 1, 1, b)
        base = build_factor_base(n, bound=3).factor_base
        config = SieveConfig(
            residual_bound=1,
            block_width=1,
            division="bucket",
            score_policy="powers",
        )
        direct, packed = (
            SieveCollector(
                polynomial,
                base,
                config=replace(config, sparse_queue_bytes=size),
                budget=allowance(),
            )
            for size in (0, 2**20)
        )
        expected = direct.collect(2, 3)
        observed = packed.collect(2, 3)
        self.assertIsNotNone(packed._sparse_queue)
        self.assertEqual(observed.atoms, expected.atoms)
        self.assertEqual(observed.atoms[0].exponents, ((2, 301),))
