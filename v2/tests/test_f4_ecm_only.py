"""ECM-only measurements reject fallback and invalid split publication."""

import unittest
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.ecm.b12_f4.selection_fixed_child import (
    MatchedSerialECM,
    SelectedECMOnly,
    checked_divisor,
    endpoint,
)


class ECMOnlyTests(unittest.TestCase):
    def test_no_fallback_is_explicit_and_fails_closed(self):
        self.assertEqual(
            endpoint(dict(endpoint="first", config=dict(siqs=None, sss=None))),
            "first",
        )
        for config in ({}, dict(siqs={}, sss=None), dict(siqs=None, sss={})):
            with self.assertRaises(ValueError):
                endpoint(dict(endpoint="first", config=config))
        self.assertEqual(endpoint({}), "fixed")
        with self.assertRaises(ValueError):
            endpoint(dict(endpoint="complete"))

    def test_proper_divisor_oracle(self):
        self.assertEqual(checked_divisor(221, 13), 13)
        self.assertEqual(checked_divisor(221, 17), 17)
        self.assertIsNone(checked_divisor(221, None))
        for value in (0, 1, 221, 442, 14, True, 13.0):
            with self.assertRaises(ValueError):
                checked_divisor(221, value)


class MatchedMasterTests(unittest.TestCase):
    def request(self):
        return dict(
            manifest="declared-search",
            assignments=[dict(n=221, seed=7, sigma=31, b1=20, b2=100)],
        )

    def test_master_arithmetic_owns_setup_and_preserves_full_config(self):
        from v2.tests.resource_stops import InterruptBudget as Budget

        owner = MatchedSerialECM(self.request())
        config = SimpleNamespace(ecm_tiers=((20, 100, 1),))
        programs = object()

        def advance(job, budget, context, actual_config):
            self.assertIs(context, programs)
            self.assertIs(actual_config, config)
            self.assertEqual(job["seed"], 7)
            job.update(done=True, sigma=31, factor=13)

        with (
            patch(
                "v2.portfolio._execution_context",
                return_value=(None, programs),
            ),
            patch("v2.execution.stage_jobs.advance_job", side_effect=advance),
        ):
            result = owner.run(221, config=config, budget=Budget())
        self.assertEqual(result.divisor, 13)
        self.assertFalse(result.stats["engaged"])
        self.assertEqual(result.stats["completed_curves"], 1)

    def test_wrong_sigma_or_invalid_factor_cannot_be_published(self):
        from v2.tests.resource_stops import InterruptBudget as Budget

        config = SimpleNamespace(ecm_tiers=((20, 100, 1),))
        for sigma, factor in ((32, 13), (31, 14), (31, 221)):

            def advance(job, *args):
                job.update(done=True, sigma=sigma, factor=factor)

            with (
                patch(
                    "v2.portfolio._execution_context",
                    return_value=(object(), None),
                ),
                patch(
                    "v2.execution.stage_jobs.advance_job", side_effect=advance
                ),
                self.assertRaises(ValueError),
            ):
                MatchedSerialECM(self.request()).run(
                    221, config=config, budget=Budget()
                )

    def test_declined_policy_includes_decision_and_serial_cost(self):
        from unittest.mock import Mock

        owner = SelectedECMOnly.__new__(SelectedECMOnly)
        owner.owner = Mock()
        owner.owner.run_portfolio.return_value = None
        owner.owner.last_stats = dict(engaged=False, selector_work=512)
        owner.serial = Mock(manifest="declared-search")
        owner.serial.run.return_value = SimpleNamespace(
            divisor=17, stats=dict(completed_curves=1)
        )
        result = owner.run(221, budget=object())
        owner.serial.run.assert_called_once()
        self.assertEqual(result.divisor, 17)
        self.assertEqual(result.stats["selector_work"], 512)
        self.assertEqual(result.stats["manifest"], "declared-search")
