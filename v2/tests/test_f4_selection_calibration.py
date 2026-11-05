"""Independent complete-cost crossover and refusal oracles."""

import unittest

from v2.benchmarks.ecm.b12_f4.selection_calibration import calibrated_costs


class CalibrationTests(unittest.TestCase):
    def costs(self, **changes):
        values = dict(
            serial_cold_body=4,
            serial_cold_cli=4.2,
            serial_warm=4,
            parallel_cold_cli=2.4,
            parallel_warm=1.2,
            parallel_warm_parent=0.2,
            pool_start=0.6,
            retained_calls=9,
        )
        values.update(changes)
        return calibrated_costs(((2000, 200000, 20),), **values)

    def test_known_linear_cost_crossover_and_retained_amortization(self):
        # Serial = W; parallel = 1.2 + W/4, excluding common CLI cost 0.2.
        # Their cold crossover is W=1.6, while warm 0.2 + W/4 crosses 4/15.
        result = self.costs()
        self.assertAlmostEqual(result["cold"]["cold_seconds"], 1.6)
        self.assertAlmostEqual(result["cold"]["warm_seconds"], 4 / 15)
        self.assertAlmostEqual(result["cold"]["pool_seconds"], 0.8)
        self.assertEqual(result["cold"]["seconds_per_curve"], (0.2,))
        self.assertAlmostEqual(
            result["warm"]["cold_seconds"], (0.2 + 1 / 9) / 0.75
        )

    def test_no_scaling_or_invalid_observations_cannot_engage(self):
        result = self.costs(parallel_warm=4.5)
        for mode in result.values():
            self.assertIsNone(mode["cold_seconds"])
            self.assertIsNone(mode["warm_seconds"])
            self.assertIsNone(mode["pool_seconds"])
        for invalid in (float("nan"), float("inf"), -1, True):
            with self.assertRaises(ValueError):
                self.costs(pool_start=invalid)
        with self.assertRaises(ValueError):
            self.costs(serial_cold_cli=3)
        with self.assertRaises(ValueError):
            self.costs(retained_calls=0)

    def test_horizon_one_matches_one_shot_when_serial_rates_match(self):
        result = self.costs(retained_calls=1)
        self.assertEqual(result["cold"], result["warm"])

    def test_impossible_or_exhausted_plans_reject(self):
        for tiers in ((), ((2, 3, 0),), ((5, 2, 1),), ((2, 3, 4097),)):
            with self.assertRaises(ValueError):
                calibrated_costs(
                    tiers,
                    serial_cold_body=1,
                    serial_cold_cli=1,
                    serial_warm=1,
                    parallel_cold_cli=1,
                    parallel_warm=1,
                    parallel_warm_parent=0,
                    pool_start=0,
                )


if __name__ == "__main__":
    unittest.main()
