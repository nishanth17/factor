"""Round-two input/measurement controls without factoring timing claims."""

import copy
import unittest
from unittest.mock import patch

from v2 import portfolio
from v2.benchmarks.ecm.c3 import round2_corpus, round2_pilot


class RoundTwoControlsTests(unittest.TestCase):
    def test_certified_corpus_checks_actual_factor_stratum(self):
        corpus = dict(
            certificates={
                "1009": {"kind": "trial"},
                "1013": {"kind": "trial"},
            },
            fixtures=[
                dict(
                    id="small",
                    n=1009 * 1013,
                    factors=[(1009, 1), (1013, 1)],
                    digits=7,
                    small_digits=4,
                )
            ],
        )
        round2_corpus.validate_corpus(corpus)
        invalid = copy.deepcopy(corpus)
        invalid["fixtures"][0]["small_digits"] = 3
        with self.assertRaisesRegex(ValueError, "smaller-factor"):
            round2_corpus.validate_corpus(invalid)
        invalid = copy.deepcopy(corpus)
        invalid["fixtures"].append(invalid["fixtures"][0])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            round2_corpus.validate_corpus(invalid)

    def test_all_policies_construct_using_only_observable_size(self):
        for digits in (30, 40):
            for arm in round2_pilot.POLICIES:
                n = 10 ** (digits - 1) + 1
                config = round2_pilot.configuration(n, arm, portfolio)
                self.assertEqual(config.memory_bytes, 288 * 2**20)
                self.assertIsNotNone(config.siqs)
                self.assertEqual(
                    config.ecm_tiers, round2_pilot.POLICIES[arm]["tiers"]
                )
                if arm == "control":
                    self.assertIsNone(config.allocation)
                else:
                    self.assertEqual(
                        config.allocation.mode,
                        round2_pilot.POLICIES[arm]["mode"],
                    )
                if arm.startswith("economic"):
                    estimate = 0.2 if digits == 30 else 4.5
                    expected = (
                        estimate * round2_pilot.POLICIES[arm]["fraction"]
                    )
                    self.assertEqual(
                        config.allocation.pretest_seconds, expected
                    )
                    self.assertEqual(
                        config.allocation.pretest_cpu_seconds, expected
                    )

    def test_generation_has_a_real_draw_limit(self):
        generator = round2_corpus.GenerationRandom(7, draws=1)
        generator.getrandbits(2)
        with self.assertRaisesRegex(RuntimeError, "generation allowance"):
            generator.getrandbits(2)

    def test_instrumentation_preserves_executor_and_assignment(self):
        calls = []

        def execute(job, budget, context, config):
            calls.append((job["seed"], context, config))
            budget.used += 3
            job["done"] = True

        class Engine:
            advance_job = staticmethod(execute)

        class Ledger:
            used = 9

        job = dict(
            kind="ecm",
            n=1009 * 1013,
            seed=17,
            b1=50,
            b2=1000,
            done=False,
            factor=None,
            start_work=9,
        )
        with round2_pilot.attribution(Engine, True) as costs:
            Engine.advance_job(job, Ledger, "context", "config")
        self.assertEqual(calls, [(17, "context", "config")])
        self.assertEqual(costs[0]["work"], 3)
        self.assertEqual(costs[0]["outcome"], "exhausted")
        self.assertIs(Engine.advance_job, execute)

    def test_instrumentation_is_absent_for_performance_calls(self):
        with patch.object(portfolio, "advance_job") as execute:
            with round2_pilot.attribution(portfolio, False) as costs:
                self.assertIs(portfolio.advance_job, execute)
            self.assertEqual(costs, [])

    def test_historical_control_without_allocation_field(self):
        from v2.benchmarks.ecm.c3 import c3_study

        fixture = dict(
            id="small",
            kind="regression",
            n=1009 * 1013,
            factors=[(1009, 1), (1013, 1)],
        )
        row = round2_pilot.run_one(
            fixture, 17, "control", c3_study.baseline(), instrument=True
        )
        self.assertTrue(row["complete"])
        self.assertIsNone(row["allocation"])
        self.assertTrue(row["instrumented"])

    def test_reserved_charges_match_frozen_policy_boundaries(self):
        import hashlib
        import importlib
        import itertools
        import json
        import sys

        from v2.benchmarks.ecm.c3 import c3_study
        from v2.execution.allocation import ECMAllocation, PretestBudget
        from v2.execution.budget import Budget

        path = c3_study.INPUTS / "baselines/c3_round2_before_optimization.json"
        snapshot = json.loads(path.read_text())
        for name, source in snapshot["source"].items():
            self.assertEqual(
                hashlib.sha256(source.encode()).hexdigest(),
                snapshot["sha256"][name],
            )
        package = "_c3_before_reservation_optimization"
        if package not in sys.modules:
            sys.meta_path.insert(
                0, c3_study.SourceFinder(package, snapshot["source"])
            )
        old_view = importlib.import_module(
            package + ".execution.allocation"
        ).PretestBudget
        old_budget = importlib.import_module(
            package + ".execution.budget"
        ).Budget

        policies = [
            ECMAllocation("campaign"),
            ECMAllocation(
                "campaign",
                fallback_work=4,
                fallback_seconds=2,
                fallback_cpu_seconds=2,
            ),
            ECMAllocation(
                "pretest",
                pretest_work=8,
                pretest_seconds=3,
                pretest_cpu_seconds=3,
            ),
        ]

        def outcome(view, budget):
            outcomes = []
            for amount in (0, 5, 2, 9, 0):
                try:
                    view.consume(amount)
                    error = None
                except Exception as refusal:
                    error = (type(refusal).__name__, str(refusal))
                outcomes.append((error, budget.used, budget.reason))
            return outcomes

        with (
            patch("time.monotonic", return_value=100),
            patch("time.process_time", return_value=100),
        ):
            for policy, prior, clocks, fallback, cancel in itertools.product(
                policies,
                (0, 2, 4, 11),
                (None, 10),
                (False, True),
                (False, True),
            ):
                options = dict(
                    work_limit=20,
                    used=7,
                    prior_wall=prior,
                    prior_cpu=prior,
                    _wall_start=100,
                    _cpu_start=100,
                    seconds=clocks,
                    cpu_seconds=clocks,
                    cancelled=lambda: cancel,
                )
                new, old = Budget(**options), old_budget(**options)
                self.assertEqual(
                    outcome(
                        PretestBudget(new, policy, fallback=fallback), new
                    ),
                    outcome(old_view(old, policy, fallback=fallback), old),
                )
