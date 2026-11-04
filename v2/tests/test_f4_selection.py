"""Independent contracts for finite selection and raw complete-call outputs."""

import copy
import math
import unittest

from v2.benchmarks.ecm.b12_f4.selection import (
    censored_quantile,
    matched_ratio,
    regime_key,
    shortlist,
)
from v2.benchmarks.ecm.b12_f4.selection_child import (
    summarize_stages,
    validate_reconstruction,
)
from v2.benchmarks.ecm.b12_f4.selection_inputs import (
    smaller_sizes,
    uniform_panel,
)


class SelectionTests(unittest.TestCase):
    def test_censoring_cannot_create_a_fast_tail(self):
        observations = [
            dict(seconds=1, complete=True, grant=10),
            dict(seconds=2, complete=True, grant=10),
            dict(seconds=0.1, complete=False, grant=10),
            dict(seconds=0.2, complete=False, grant=10),
        ]
        self.assertEqual(censored_quantile(observations, 0.5)["seconds"], 2)
        tail = censored_quantile(observations, 0.9)
        self.assertIsNone(tail["seconds"])
        self.assertEqual(tail["lower_bound"], 10)
        self.assertEqual(tail["completion"], 0.5)

    def test_all_censored_and_mixed_grants(self):
        row = dict(seconds=0.1, complete=False, grant=5)
        self.assertIsNone(censored_quantile([row], 0.5)["seconds"])
        with self.assertRaises(ValueError):
            censored_quantile([row, dict(row, grant=10)], 0.5)
        with self.assertRaises(ValueError):
            censored_quantile([dict(row, complete=True, seconds=6)], 0.5)

    def test_conditional_ratio_has_known_scale_and_finite_draws(self):
        control = {i: i + 1 for i in range(9)}
        candidate = {i: 2 * value for i, value in control.items()}
        result = matched_ratio(candidate, control, draws=100)
        self.assertAlmostEqual(result["ratio"], 2)
        for bound in result["interval"]:
            self.assertAlmostEqual(bound, 2)
        with self.assertRaises(ValueError):
            matched_ratio(candidate, {1: 1}, draws=100)
        with self.assertRaises(ValueError):
            matched_ratio(candidate, control, draws=10**9)

    def test_failure_and_timeout_are_not_selection_winners(self):
        rows = [
            dict(count=0, seconds=5, status="ok"),
            dict(count=4, seconds=2, status="ok"),
            dict(count=6, seconds=2.1, status="ok"),
            dict(count=8, seconds=0.01, status="censored"),
            dict(count=16, seconds=0.001, status="recorder_failure"),
        ]
        result = shortlist(rows, target=2)
        self.assertEqual(result[0]["count"], 4)
        self.assertIn(0, [row["count"] for row in result])
        self.assertNotIn(8, [row["count"] for row in result])
        self.assertNotIn(16, [row["count"] for row in result])

    def test_adjacent_count_is_retained_despite_one_slow_screen(self):
        rows = [
            dict(count=count, seconds=seconds, status="ok")
            for count, seconds in [(0, 10), (4, 2), (6, 1), (8, 2)]
        ]
        result = shortlist(rows, target=1)
        self.assertEqual({row["count"] for row in result}, {0, 4, 6, 8})

    def test_near_ties_are_not_truncated_at_the_screen_target(self):
        rows = [
            dict(count=count, seconds=1 + index / 1000, status="ok")
            for index, count in enumerate((0, 1, 2, 4, 6, 8, 10, 12, 14, 16))
        ]
        self.assertEqual(len(shortlist(rows, target=2)), 10)

    def test_invalid_timings_and_boolean_counts(self):
        for value in (0, -1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                shortlist([dict(count=4, seconds=value, status="ok")])
        with self.assertRaises(ValueError):
            shortlist([dict(count=True, seconds=1, status="ok")])

    def test_regime_keys_use_only_observable_plan_fields(self):
        key = regime_key(70, [(2000, 200000, 128)])
        self.assertEqual(key, (70, ((2000, 200000, 128),), "python-int"))
        self.assertNotEqual(key, regime_key(70, [(2000, 200000, 512)]))
        for tiers in ([], [(2000, 1000, 1)], [(2, 10, 4097)]):
            with self.assertRaises(ValueError):
                regime_key(70, tiers)


class InputPlanTests(unittest.TestCase):
    def test_factor_strata_are_mathematically_possible(self):
        for digits in (30, 40, 50, 60, 70, 80, 100):
            self.assertEqual(len(smaller_sizes(digits)), 4)
            for small in smaller_sizes(digits):
                self.assertLessEqual(small, digits // 2)

    def test_uniform_draws_are_reproducible_and_include_both_parities(self):
        rows = uniform_panel(2026110902, 12, 4)
        self.assertEqual(rows, uniform_panel(2026110902, 12, 4))
        self.assertEqual({row["n"] % 2 for row in rows}, {0, 1})
        self.assertEqual(
            {row["digits"] for row in rows}, {30, 40, 50, 60, 70, 80, 100}
        )
        self.assertTrue(
            all(len(str(row["n"])) == row["digits"] for row in rows)
        )
        self.assertEqual(len(rows), 52)


class StageAttributionTests(unittest.TestCase):
    def test_dropped_trace_does_not_become_no_ecm_success(self):
        result = summarize_stages([], 10, {"ecm": 0.1})
        self.assertTrue(result["ecm_reached"])
        self.assertEqual(result["ecm_success"], "unknown")
        self.assertIsNone(result["engaged"])

    def test_success_and_fallback_are_independent_of_trace_tail(self):
        events = [
            dict(
                stage="ecm_workers",
                outcome="factor_found",
                stats=dict(engaged=True),
            )
        ]
        result = summarize_stages(events, 10, {"siqs": 1})
        self.assertEqual(result["ecm_success"], "success")
        self.assertTrue(result["engaged"])
        self.assertTrue(result["fallback_reached"])


class RawResultTests(unittest.TestCase):
    def row(self):
        return dict(
            n=45,
            reconstruct=45,
            factors=[[3, 2, "proven_prime"]],
            remaining=[5],
            factorization_complete=False,
            endpoint="first",
            complete=True,
        )

    def test_raw_product_and_certainty_are_checked(self):
        validate_reconstruction(self.row())
        row = self.row()
        row["factors"][0][0] = 7
        with self.assertRaises(ValueError):
            validate_reconstruction(row)
        row = self.row()
        row["factors"][0][2] = "prime"
        with self.assertRaises(ValueError):
            validate_reconstruction(row)

    def test_prime_power_can_expose_a_proper_divisor(self):
        row = dict(
            self.row(),
            n=9,
            reconstruct=9,
            remaining=[],
            factorization_complete=True,
        )
        validate_reconstruction(row)

    def test_prime_exit_is_not_a_proper_divisor(self):
        row = dict(
            self.row(),
            n=3,
            reconstruct=3,
            factors=[[3, 1, "proven_prime"]],
            remaining=[],
            factorization_complete=True,
        )
        with self.assertRaises(ValueError):
            validate_reconstruction(row)

    def test_completion_cannot_hide_unresolved_cofactors(self):
        row = dict(self.row(), factorization_complete=True)
        with self.assertRaises(ValueError):
            validate_reconstruction(row)
        row = copy.deepcopy(self.row())
        row["remaining"] = [1, 45]
        with self.assertRaises(ValueError):
            validate_reconstruction(row)


if __name__ == "__main__":
    unittest.main()
