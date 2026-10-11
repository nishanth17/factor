"""Cumulative reservation, handoff, campaign and checkpoint acceptance."""

import copy
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
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
            max_input_bits=256,
            ecm_chain_mode="off",
            allocation=ECMAllocation("pretest", pretest_work=10_000),
            ecm_tiers=((200, 7700, 32),),
        )
        number = (2**61 - 1) * (2**89 - 1)
        run = factorize_bounded(number, config=config, budget=ledger())
        self.assertFalse(run.result.complete)
        self.assertEqual(run.reason, "pretest_exhausted")
        self.assertLessEqual(run.work_used, 10_000)
        self.assertEqual(run.events[-1]["stage"], "handoff")
        self.assertEqual(run.events[-1]["outcome"], "pretest_work")
        partial = next(
            e
            for e in run.events
            if e["stage"] == "ecm" and e["outcome"] == "handoff"
        )
        self.assertGreater(partial["work"], 0)
        self.assertIn("phase", partial)
        self.assertIn("cursor", partial)
        self.assertEqual(run.result.reconstruct(), number)

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
            allocation=ECMAllocation("pretest", pretest_work=0)
        )
        for number in (1009, 1009**3, -(1009**2), 2**30):
            with patch(
                "v2.portfolio.SIQSJob", side_effect=AssertionError("fallback")
            ):
                run = factorize_bounded(number, config=config, budget=ledger())
            self.assertTrue(run.result.complete)
            self.assertEqual(run.result.reconstruct(), number)

    def test_direct_handoff_avoids_unused_schedule_setup(self):
        with patch(
            "v2.portfolio._execution_context",
            side_effect=AssertionError("unused ECM context"),
        ):
            run = factorize_bounded(
                1009 * 1013, config=configuration(), budget=ledger()
            )
        self.assertTrue(run.result.complete)

    def test_reserved_work_reaches_actual_fallback(self):
        # A campaign can end before an unfinished curve consumes the floor.
        config = configuration(
            allocation=ECMAllocation("campaign", fallback_work=400_000),
            ecm_tiers=((2000, 147396, 32),),
        )

        def expensive_action(job, budget, context, config):
            budget.consume(100_000)

        with patch("v2.portfolio.advance_job", expensive_action):
            run = factorize_bounded(
                1009 * 1013,
                config=config,
                seed=7,
                budget=ledger(450_000),
            )
        handoff = next(e for e in run.events if e["stage"] == "handoff")
        self.assertEqual(handoff["outcome"], "reserved_work")
        self.assertGreaterEqual(handoff["remaining_work"], 400_000)
        self.assertTrue(
            any(e["stage"] == "siqs" for e in run.events)
            or "siqs_seed"
            in (run.checkpoint["payload"]["state"]["current"] or {})
        )
        self.assertEqual(run.result.reconstruct(), 1009 * 1013)

    def test_cancelled_curve_resumes_deterministic_assignments(self):
        config = configuration(
            siqs=None,
            allocation=ECMAllocation("campaign"),
            ecm_tiers=((50, 1000, 3),),
            max_input_bits=256,
            ecm_chain_mode="off",
        )
        number = (2**61 - 1) * (2**89 - 1)
        uninterrupted = factorize_bounded(
            number, config=config, seed=7, budget=ledger()
        )
        budget = ledger()
        budget.cancelled = lambda: budget.used >= 8000
        first = factorize_bounded(number, config=config, seed=7, budget=budget)
        self.assertEqual(first.reason, "cancelled")
        restored = factorize_bounded(
            number, checkpoint=first.checkpoint, budget=ledger()
        )

        def assignments(run):
            return [
                (e.get("seed"), e.get("b1"), e.get("b2"), e["outcome"])
                for e in run.events
                if e["stage"] == "ecm"
            ]

        self.assertEqual(assignments(restored), assignments(uninterrupted))
        self.assertEqual(restored.result, uninterrupted.result)
        self.assertGreater(restored.work_used, uninterrupted.work_used)
        self.assertGreaterEqual(restored.wall_seconds, first.wall_seconds)

    def test_new_progress_metadata_rejects_resealed_corruption(self):
        first = factorize_bounded(
            1009 * 1013, config=configuration(), budget=ledger(1000)
        )
        for field, value in (
            ("fallback_admitted", "yes"),
            ("pending_handoff", "invented"),
        ):
            corrupt = copy.deepcopy(first.checkpoint)
            corrupt["payload"]["state"]["current"][field] = value
            with self.assertRaises(ValueError):
                factorize_bounded(
                    1009 * 1013, checkpoint=reseal(corrupt), budget=ledger()
                )

    def test_immutable_mainline_checkpoint_loads_with_identical_policy(self):
        from v2.benchmarks.ecm.c3.c3_study import baseline

        old = baseline()
        native_config = configuration(
            siqs=None,
            allocation=None,
            max_input_bits=256,
            ecm_chain_mode="off",
        )
        number = (2**61 - 1) * (2**89 - 1)
        values = dict(vars(native_config))
        values.pop("allocation")
        old_config = old.PortfolioConfig(**values)
        old_budget = sys.modules[old.__package__ + ".execution.budget"].Budget
        original = old.factorize_bounded(
            number,
            config=old_config,
            seed=7,
            budget=old_budget(work_limit=8000, seconds=None, cpu_seconds=None),
        )
        restored = factorize_bounded(
            number,
            config=native_config,
            checkpoint=original.checkpoint,
            budget=ledger(),
        )
        self.assertEqual(restored.result.reconstruct(), number)
        self.assertGreater(restored.work_used, original.work_used)

    def test_cli_policy_resume_restores_fallback_and_rejects_override(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "state.json"
            command = [sys.executable, "-B", "-m", "v2.factor"]
            first = subprocess.run(
                command
                + [
                    "1022117",
                    "--siqs",
                    "--ecm-policy",
                    "pretest",
                    "--pretest-work",
                    "0",
                    "--fallback-work",
                    "500000",
                    "--work-limit",
                    "1000",
                    "--checkpoint",
                    str(path),
                ],
                capture_output=True,
                text=True,
            )
            self.assertTrue(path.exists(), first.stderr)
            resumed = subprocess.run(
                command
                + [
                    "--resume",
                    str(path),
                    "--work-limit",
                    "200000000",
                    "--verbose",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(resumed.returncode, 0, resumed.stderr)
            self.assertIn("complete", resumed.stdout.lower())
            invalid = subprocess.run(
                command
                + [
                    "--resume",
                    str(path),
                    "--ecm-policy",
                    "pretest",
                    "--pretest-work",
                    "1",
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(invalid.returncode, 0)
            self.assertIn("incompatible", invalid.stderr)

    def test_both_sss_modes_restore_nested_collector_after_handoff(self):
        from v2.qs.sss import SSSConfig

        for mode in ("sss", "sssf"):
            config = configuration(
                siqs=None,
                sss=SSSConfig(
                    mode=mode, base_bound=100, memory_bytes=32 * 2**20
                ),
                allocation=ECMAllocation(
                    "pretest", pretest_work=0, fallback_work=10_000
                ),
            )
            first = factorize_bounded(
                41 * 43, config=config, budget=ledger(1000)
            )
            self.assertEqual(first.reason, "insufficient_fallback_work")
            restored = factorize_bounded(
                41 * 43, checkpoint=first.checkpoint, budget=ledger()
            )
            self.assertTrue(restored.result.complete)
            self.assertEqual(restored.result.reconstruct(), 41 * 43)
            self.assertGreater(restored.work_used, first.work_used)
            self.assertEqual(
                sum(e["stage"] == "handoff" for e in restored.events), 1
            )

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


class StudyContractTests(unittest.TestCase):
    def test_censored_fast_refusal_does_not_count_as_speed_gain(self):
        from v2.benchmarks.ecm.c3.analyze import compare, summarize

        rows = []
        for identity in ("first", "second"):
            for arm, wall in (("control", 5), ("quick8", 0.001)):
                for repetition in range(9):
                    rows.append(
                        dict(
                            id=identity,
                            seed=7,
                            arm=arm,
                            kind="balanced",
                            band=40,
                            complete=False,
                            proper_factor=False,
                            wall=wall,
                            cpu=wall,
                            work=10,
                            curves=0,
                            fallback_started=False,
                            memory_cap=100,
                            workspace_reserve=20,
                            fallback_owned_peak=0,
                            rss=200,
                            checkpoint_bytes=30,
                            cap_seconds=30,
                        )
                    )
        medians, unstable = summarize(rows)
        result = compare(medians, ["first", "second"], "quick8")
        self.assertEqual(unstable, [])
        self.assertEqual(result["saving"], 0)
        self.assertEqual(result["interval95"], [0, 0])
        self.assertEqual(result["inputs"], 2)
        self.assertEqual(result["candidate_completion"], 0)

    def test_historical_corpora_validate_and_do_not_enter_policy(self):
        from v2 import portfolio
        from v2.benchmarks.ecm.c3.c3_study import config_for, training

        fixtures = training()
        self.assertEqual(len(fixtures), 10)
        for fixture in fixtures:
            self.assertEqual(
                fixture["n"],
                __import__("math").prod(p**e for p, e in fixture["factors"]),
            )
        first = config_for(10**29 + 1, "quick8", portfolio)
        second = config_for(10**29 + 3, "quick8", portfolio)
        self.assertEqual(first, second)
        self.assertFalse(hasattr(first, "factors"))
