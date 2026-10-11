"""Cumulative reservation, handoff, campaign and checkpoint acceptance."""

import copy
import hashlib
import json
import subprocess
import sys
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.execution.allocation import (
    ECMAllocation,
    HandoffRequiredError,
    PretestBudget,
    fallback_refusal,
)
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import SIQSConfig


def configuration(**changes):
    config = PortfolioConfig(
        trial_bound=5,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=((50, 1000, 3),),
        max_input_bits=128,
        segment_size=16,
        memory_bytes=64 * 2**20,
        siqs=SIQSConfig(
            base_bound=100, half_width=64, memory_bytes=32 * 2**20
        ),
        allocation=ECMAllocation(
            "pretest", pretest_work=0, fallback_work=1000
        ),
    )
    return replace(config, **changes)


def ledger(work=2_000_000, **changes):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None, **changes)


def reseal(checkpoint):
    encoded = json.dumps(
        checkpoint["payload"], sort_keys=True, separators=(",", ":")
    )
    checkpoint["sha256"] = hashlib.sha256(encoded.encode()).hexdigest()
    return checkpoint


class AllocationBudgetTests(unittest.TestCase):
    def test_atomic_work_reservation_and_prior_credit(self):
        budget = ledger(100, used=30)
        policy = ECMAllocation("pretest", pretest_work=80, fallback_work=40)
        view = PretestBudget(budget, policy, fallback=True)

        view.consume(30)
        with self.assertRaisesRegex(HandoffRequiredError, "reserved_work"):
            view.consume(1)
        self.assertEqual(budget.used, 60)
        with self.assertRaisesRegex(HandoffRequiredError, "pretest_work"):
            view.consume(21)
        self.assertEqual(budget.used, 60)

    def test_clock_reservations_and_insufficient_admission(self):
        policy = ECMAllocation(
            "pretest",
            pretest_work=100,
            fallback_seconds=2,
            fallback_cpu_seconds=2,
        )
        for changes, expected in (
            (dict(seconds=3, prior_wall=1.5, cpu_seconds=None), "wall"),
            (dict(cpu_seconds=3, prior_cpu=1.5, seconds=None), "cpu"),
        ):
            budget = Budget(**changes)
            with self.assertRaisesRegex(
                HandoffRequiredError, "reserved_" + expected
            ):
                PretestBudget(budget, policy, fallback=True).consume(0)
            self.assertEqual(
                fallback_refusal(budget, policy),
                "insufficient_fallback_" + expected,
            )

    def test_cancellation_precedes_handoff(self):
        budget = ledger(cancelled=lambda: True)
        view = PretestBudget(
            budget, ECMAllocation("pretest", pretest_work=0), fallback=False
        )
        with self.assertRaisesRegex(BudgetExhaustedError, "cancelled"):
            view.consume(1)
        self.assertEqual(budget.used, 0)

    def test_policy_validation_and_memory_reserve(self):
        for changes in (
            dict(mode="automatic"),
            dict(pretest_work=-1),
            dict(pretest_seconds=float("inf")),
            dict(fallback_seconds=True),
        ):
            values = dict(mode="pretest", pretest_work=0)
            values.update(changes)
            with self.assertRaises((ValueError, TypeError)):
                ECMAllocation(**values)
        with self.assertRaises(ValueError):
            ECMAllocation("campaign", pretest_work=1)
        with self.assertRaises(ValueError):
            configuration(siqs=None)
        with self.assertRaises(MemoryError):
            configuration(memory_bytes=32 * 2**20)


class AllocationPortfolioTests(unittest.TestCase):
    def test_real_handoff_and_quiet_reconstruction(self):
        import contextlib
        import io

        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            run = factorize_bounded(
                1009 * 1013, config=configuration(), budget=ledger()
            )

        self.assertTrue(run.result.complete)
        self.assertEqual(run.result.reconstruct(), 1009 * 1013)
        self.assertEqual(output.getvalue(), "")
        self.assertEqual(run.events[0]["outcome"], "pretest_work")
        self.assertTrue(any(event["stage"] == "siqs" for event in run.events))
        self.assertEqual(run.checkpoint["payload"]["version"], 12)

    def test_insufficient_fallback_is_resumable_without_pretest_restart(self):
        config = configuration(
            allocation=ECMAllocation(
                "pretest", pretest_work=0, fallback_work=500_000
            )
        )
        first = factorize_bounded(
            1009 * 1013, config=config, budget=ledger(1000)
        )
        self.assertEqual(first.reason, "insufficient_fallback_work")
        self.assertEqual(
            first.checkpoint["payload"]["state"]["current"]["stage"], "siqs"
        )

        resumed = factorize_bounded(
            1009 * 1013, checkpoint=first.checkpoint, budget=ledger()
        )
        self.assertTrue(resumed.result.complete)
        self.assertGreater(resumed.work_used, first.work_used)
        self.assertEqual(
            sum(e["stage"] == "handoff" for e in resumed.events), 1
        )

    def test_live_fallback_resume_keeps_admission_and_spent_resources(self):
        config = configuration(
            allocation=ECMAllocation(
                "pretest", pretest_work=0, fallback_work=10_000
            )
        )
        first = factorize_bounded(
            1009 * 1013, config=config, budget=ledger(30_000)
        )
        self.assertEqual(first.reason, "work_limit")
        current = first.checkpoint["payload"]["state"]["current"]
        self.assertTrue(current["fallback_admitted"])
        self.assertIn("siqs_checkpoint", current)

        restored = factorize_bounded(
            1009 * 1013, checkpoint=first.checkpoint, budget=ledger()
        )
        self.assertTrue(restored.result.complete)
        self.assertGreater(restored.work_used, first.work_used)
        self.assertGreaterEqual(restored.wall_seconds, first.wall_seconds)
        self.assertGreaterEqual(restored.cpu_seconds, first.cpu_seconds)

    def test_recursion_cannot_reset_pretest_ceiling(self):
        config = configuration()
        seen = []

        class SplitJob:
            def __init__(self, n, *, seed, config, budget):
                self.n, self.seed = n, seed
                self.finished_reason = None
                seen.append((n, budget.used))
                self.budget = budget

            def run(self, **unused):
                from types import SimpleNamespace

                self.budget.consume(100)
                divisor = next(p for p in (1009, 1013) if self.n % p == 0)
                return SimpleNamespace(
                    divisor=divisor, stats={}, reason="factor"
                )

        number = 1009 * 1013 * 1019
        with patch("v2.portfolio.SIQSJob", SplitJob):
            run = factorize_bounded(number, config=config, budget=ledger())
        self.assertTrue(run.result.complete)
        self.assertEqual(len(seen), 2)
        self.assertGreater(seen[1][1], seen[0][1])
        self.assertEqual(sum(e["stage"] == "handoff" for e in run.events), 2)
        self.assertFalse(any(e["stage"] == "ecm" for e in run.events))

    def test_partial_curve_handoff_retains_coverage_and_charge(self):
        config = configuration(
            siqs=None,
            allocation=ECMAllocation("pretest", pretest_work=1000),
            ecm_tiers=((2000, 147396, 32),),
        )
        run = factorize_bounded(
            1_000_003 * 1_000_033, config=config, budget=ledger()
        )
        self.assertFalse(run.result.complete)
        self.assertLessEqual(run.work_used, 1000)
        self.assertEqual(run.events[-1]["stage"], "handoff")
        self.assertEqual(run.events[-1]["outcome"], "pretest_work")
        self.assertEqual(run.result.reconstruct(), 1_000_003 * 1_000_033)

    def test_campaign_exhaustion_and_legacy_schema(self):
        config = configuration(
            siqs=None,
            allocation=ECMAllocation("campaign"),
            ecm_tiers=((2, 2, 1),),
        )
        run = factorize_bounded(
            1_000_003 * 1_000_033, config=config, budget=ledger()
        )
        self.assertEqual(run.reason, "exhausted")
        self.assertEqual(run.events[-1]["outcome"], "schedule_exhausted")
        legacy = factorize_bounded(
            1009 * 1013,
            config=replace(config, allocation=None),
            budget=ledger(),
        )
        self.assertLess(legacy.checkpoint["payload"]["version"], 12)
        self.assertNotIn("allocation", legacy.checkpoint["payload"]["config"])
        restored = factorize_bounded(
            1009 * 1013,
            config=replace(config, allocation=None),
            checkpoint=legacy.checkpoint,
            budget=ledger(),
        )
        self.assertEqual(restored.result, legacy.result)

    def test_policy_corruption_and_override_are_rejected(self):
        config = configuration()
        run = factorize_bounded(
            1009 * 1013, config=config, budget=ledger(1000)
        )
        corrupt = copy.deepcopy(run.checkpoint)
        corrupt["payload"]["allocation"] = "another-policy"
        with self.assertRaises(ValueError):
            factorize_bounded(
                1009 * 1013,
                config=config,
                checkpoint=reseal(corrupt),
                budget=ledger(),
            )
        with self.assertRaises(ValueError):
            factorize_bounded(
                1009 * 1013,
                config=replace(
                    config, allocation=ECMAllocation("pretest", pretest_work=1)
                ),
                checkpoint=run.checkpoint,
                budget=ledger(),
            )

    def test_prime_and_power_do_not_use_relation_engine(self):
        config = configuration(
            allocation=ECMAllocation("pretest", pretest_work=10000)
        )
        for number in (1009, 1009**3, -(1009**2), 2**30):
            with patch(
                "v2.portfolio.SIQSJob", side_effect=AssertionError("fallback")
            ):
                run = factorize_bounded(number, config=config, budget=ledger())
            self.assertTrue(run.result.complete)
            self.assertEqual(run.result.reconstruct(), number)

    def test_cli_explicit_campaign_and_invalid_override(self):
        command = [sys.executable, "-B", "-m", "v2.factor", "1022117"]
        result = subprocess.run(
            command
            + [
                "--ecm-policy",
                "campaign",
                "--ecm-tier",
                "50,1000,2",
                "--verbose",
            ],
            capture_output=True,
            text=True,
        )
        self.assertIn("Portfolio:", result.stdout)
        invalid = subprocess.run(
            command + ["--ecm-policy", "pretest"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(invalid.returncode, 2)
        self.assertIn("requires --pretest-work", invalid.stderr)
