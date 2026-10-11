"""Independent expected-cost examples for the offline sequence fitter."""

import copy
import unittest
from unittest.mock import patch

from v2 import portfolio
from v2.benchmarks.ecm.c3 import round2_train
from v2.benchmarks.ecm.c3.round2_sequence import (
    fit_stopping_sequence,
    prefix_tiers,
)


def observation(subject, costs, hit, qs=10, recovery=0):
    return dict(
        subject=subject,
        weight=1,
        curve_costs=costs,
        hit=hit,
        qs_cost=qs,
        recursive_cost=recovery,
        setup=0,
    )


class SequenceTests(unittest.TestCase):
    def test_failure_changes_conditional_yield(self):
        records = [
            observation("a", [1], 1),
            observation("b", [1, 1], 2),
            observation("c", [1, 1], None),
            observation("d", [1, 1], None),
        ]
        result = fit_stopping_sequence(records, (1, 2), minimum_subjects=1)
        self.assertEqual(result["curves"], 2)
        self.assertAlmostEqual(result["expected_cpu"], 6.75)
        self.assertEqual(result["states"][0]["conditional_success"], 0.25)
        self.assertAlmostEqual(
            result["states"][1]["conditional_success"], 1 / 3
        )

    def test_future_yield_can_pay_for_a_zero_yield_first_block(self):
        records = [observation("a", [1, 1], 2)]
        result = fit_stopping_sequence(records, (1, 2), minimum_subjects=1)
        self.assertEqual(result["curves"], 2)
        self.assertEqual(result["expected_cpu"], 2)
        self.assertEqual(result["states"][0]["conditional_success"], 0)

    def test_expensive_recovery_or_insufficient_subjects_stops(self):
        record = observation("a", [1], 1, recovery=20)
        self.assertEqual(
            fit_stopping_sequence([record], (1,), minimum_subjects=1)[
                "curves"
            ],
            0,
        )
        record["recursive_cost"] = 0
        self.assertEqual(fit_stopping_sequence([record], (1,))["curves"], 0)
        records = [copy.deepcopy(record) for _ in range(20)]
        self.assertEqual(fit_stopping_sequence(records, (1,))["curves"], 0)

    def test_weights_and_compilation_preserve_observed_paths(self):
        records = [observation("hit", [2], 1), observation("miss", [2], None)]
        records[0]["weight"] = 9
        fitted = fit_stopping_sequence(records, (1,), minimum_subjects=1)
        self.assertAlmostEqual(fitted["expected_cpu"], 3)
        self.assertEqual(
            prefix_tiers("deep", 18), ((2000, 147396, 16), (11000, 250000, 2))
        )
        self.assertEqual(prefix_tiers("wide", 0), ())
        with self.assertRaises(ValueError):
            prefix_tiers("deep", 21)

    def test_reject_censored_and_nonfinite_costs(self):
        with self.assertRaisesRegex(ValueError, "censored"):
            fit_stopping_sequence([observation("a", [1], None)], (1, 2))
        with self.assertRaises(ValueError):
            fit_stopping_sequence([observation("a", [float("nan")], 1)], (1,))

    def test_acceptance_call_has_no_root_instrumentation(self):
        with patch.object(portfolio, "new_job") as create:
            with round2_train.root_trace(123, False):
                self.assertIs(portfolio.new_job, create)

    def test_root_trace_records_only_the_selected_modulus(self):
        moments = iter((0, 1, 3, 4, 6, 8))

        def create(*args):
            return dict(kind=args[0], n=args[1])

        with (
            patch.object(portfolio, "new_job", side_effect=create),
            patch.object(portfolio, "_event"),
            patch.object(portfolio, "SIQSJob"),
            patch("time.process_time", side_effect=lambda: next(moments)),
        ):
            with round2_train.root_trace(123, True) as trace:
                portfolio.new_job("ecm", 321, 7, 20, 100)
                portfolio.new_job("ecm", 123, 7, 20, 100)
                portfolio._event(
                    {},
                    None,
                    stage="ecm",
                    n=123,
                    seed=7,
                    b1=20,
                    b2=100,
                    outcome="factor",
                )
        self.assertEqual(trace["curves"][0]["cpu"], 2)
        self.assertEqual(trace["recursive_cpu"], 1)

    def test_partial_handoff_does_not_become_a_completed_curve(self):
        with (
            patch.object(portfolio, "new_job", return_value={}),
            patch.object(portfolio, "_event"),
            patch.object(portfolio, "SIQSJob"),
        ):
            with round2_train.root_trace(123, True) as trace:
                portfolio.new_job("ecm", 123, 7, 20, 100)
                portfolio._event(
                    {},
                    None,
                    stage="ecm",
                    n=123,
                    outcome="handoff",
                    reason="reserved_cpu",
                )
        self.assertEqual(trace["curves"], [])
        self.assertTrue(trace["censored"])
