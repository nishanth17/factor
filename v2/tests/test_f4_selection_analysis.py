"""Independent censoring and modulus-clustering analysis oracles."""

import unittest

from v2.benchmarks.ecm.b12_f4.selection_analysis import (
    cohort_median_interval,
    compare_cohort,
    compare_pool_timings,
)


def row(seconds, complete=True):
    return dict(seconds=seconds, complete=complete, grant=10)


class CohortAnalysisTests(unittest.TestCase):
    def test_known_scale_and_repeated_seeds_do_not_create_moduli(self):
        control = {i: [row(2), row(2)] for i in range(6)}
        candidate = {i: [row(1), row(1)] for i in range(6)}
        result = compare_cohort(candidate, control, draws=100)
        self.assertEqual(result["independent_moduli"], 6)
        self.assertEqual(result["repetitions_per_modulus"], 2)
        self.assertEqual(result["restricted_ratio_interval"], [0.5, 0.5])
        self.assertEqual(result["completion_interval"], [0, 0])

    def test_fast_refusals_cannot_improve_a_latency_endpoint(self):
        control = {i: [row(1)] for i in range(5)}
        candidate = {i: [row(0.01, False)] for i in range(5)}
        result = compare_cohort(candidate, control, draws=100)
        self.assertEqual(result["restricted_latency_ratio"], 10)
        self.assertEqual(result["completion_difference"], -1)
        self.assertIsNone(result["arms"]["candidate"]["median"]["seconds"])
        self.assertEqual(result["arms"]["candidate"]["p90"]["lower_bound"], 10)

    def test_mixed_samples_remain_in_the_declared_population(self):
        control = {i: [row(2), row(2)] for i in range(4)}
        candidate = {i: [row(1), row(0.01, False)] for i in range(4)}
        result = compare_cohort(candidate, control, draws=100)
        self.assertEqual(result["arms"]["candidate"]["completion"], 0.5)
        self.assertEqual(result["restricted_latency_ratio"], 2.75)
        self.assertIsNone(result["arms"]["candidate"]["p90"]["seconds"])

    def test_unmatched_moduli_grants_repetitions_and_invalid_draws_reject(
        self,
    ):
        control = {i: [row(2)] for i in range(3)}
        with self.assertRaises(ValueError):
            compare_cohort({0: [row(1)]}, control)
        candidate = {i: [row(1)] for i in range(3)}
        candidate[0][0]["grant"] = 20
        with self.assertRaises(ValueError):
            compare_cohort(candidate, control)
        with self.assertRaises(ValueError):
            compare_cohort(control, control, draws=10**8)
        candidate = {i: [row(1)] for i in range(3)}
        candidate[0].append(row(1))
        with self.assertRaises(ValueError):
            compare_cohort(candidate, control)


class PoolTimingTests(unittest.TestCase):
    def test_known_scale_for_cold_and_warm_blocks(self):
        for control in ([[2]] * 9, [[2, 2, 2]] * 3):
            candidate = [[value / 2 for value in block] for block in control]
            result = compare_pool_timings(candidate, control, draws=100)
            self.assertEqual(result["ratio"], 0.5)
            self.assertEqual(result["interval"], [0.5, 0.5])
            self.assertEqual(result["candidate_pools"], len(control))
            self.assertEqual(result["independent_moduli"], 1)

    def test_nine_calls_in_one_pool_do_not_create_three_lifetimes(self):
        with self.assertRaises(ValueError):
            compare_pool_timings([[1] * 9], [[2] * 9], draws=100)
        with self.assertRaises(ValueError):
            compare_pool_timings([[1]] * 3, [[2]] * 3, draws=100)

    def test_invalid_or_unequal_blocks_are_rejected(self):
        control = [[2, 2, 2]] * 3
        for candidate in (
            [[0, 1, 1]] * 3,
            [[float("nan"), 1, 1]] * 3,
            [[True, 1, 1]] * 3,
            [[1, 1, 1], [1, 1], [1, 1, 1]],
            [[1]] * 9,
        ):
            with self.assertRaises(ValueError):
                compare_pool_timings(candidate, control, draws=100)
        with self.assertRaises(ValueError):
            compare_pool_timings(control, control, draws=100_000)


class CensoredMedianTests(unittest.TestCase):
    def test_exact_scale_and_modulus_cluster_count(self):
        control = {i: [row(2)] * 9 for i in range(6)}
        candidate = {i: [row(1)] * 9 for i in range(6)}
        result = cohort_median_interval(candidate, control, draws=100)
        self.assertEqual(result["ratio"], 0.5)
        self.assertEqual(result["interval"], [0.5, 0.5])
        self.assertEqual(result["independent_moduli"], 6)
        self.assertEqual(result["repetitions_per_modulus"], 9)

    def test_fast_failed_calls_keep_unbounded_median(self):
        control = {i: [row(2)] for i in range(4)}
        candidate = {i: [row(0.01, False)] for i in range(4)}
        result = cohort_median_interval(candidate, control, draws=100)
        self.assertIsNone(result["ratio"])
        self.assertEqual(result["ratio_bounds"], [5, None])
        self.assertEqual(result["interval"], [5, None])
        self.assertTrue(result["upper_unbounded"])

    def test_failed_control_gives_a_bound_without_inventing_a_ratio(self):
        control = {i: [row(0.01, False)] for i in range(4)}
        candidate = {i: [row(1)] for i in range(4)}
        result = cohort_median_interval(candidate, control, draws=100)
        self.assertIsNone(result["ratio"])
        self.assertEqual(result["ratio_bounds"], [0, 0.1])
        self.assertEqual(result["interval"], [0, 0.1])


if __name__ == "__main__":
    unittest.main()
