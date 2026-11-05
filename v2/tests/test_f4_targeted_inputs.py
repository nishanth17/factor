"""Prespecified cofactor-coverage strata without outcome-driven selection."""

import unittest

from v2.benchmarks.ecm.b12_f4.targeted_inputs import (
    CALIBRATION_DIGITS,
    build,
    enriched_plan,
    hard_plan,
)


class TargetedPlanTests(unittest.TestCase):
    def test_independent_moduli_and_observable_cofactor_regimes(self):
        rows = enriched_plan()
        self.assertEqual(len(rows), 12)
        self.assertEqual(len({r["id"] for r in rows}), 12)
        for digits in (60, 80, 100):
            band = [r for r in rows if r["digits"] == digits]
            self.assertEqual(len(band), 4)
            self.assertEqual(sum(r["two_exponent"] == 0 for r in band), 2)
            self.assertEqual(
                {r["algorithm_seed"] for r in band}, {7919, 104729}
            )
            for row in band:
                self.assertLessEqual(
                    row["smaller_cofactor_digits"], row["cofactor_digits"] // 2
                )
                self.assertLessEqual(row["cofactor_digits"], digits)
        self.assertEqual(
            {r["cofactor_digits"] for r in rows if r["two_exponent"]},
            {55, 75, 90},
        )

    def test_frontier_is_separate_and_mathematically_possible(self):
        rows = hard_plan()
        self.assertEqual(
            [(r["digits"], r["smaller_cofactor_digits"]) for r in rows],
            [(80, 22), (100, 25)],
        )
        self.assertTrue(all(r["population"] == "hard_frontier" for r in rows))
        self.assertTrue(all(r["two_exponent"] == 0 for r in rows))

    def test_calibration_and_confirmation_are_distinct_lanes(self):
        self.assertEqual(CALIBRATION_DIGITS, (55, 75, 81, 90))
        for kwargs in (
            dict(phase="unknown", seed=1),
            dict(phase="calibration", seed=-1),
            dict(phase="calibration", seed=1, seconds=1000),
            dict(phase="calibration", seed=1, draws=0),
        ):
            with self.assertRaises(ValueError):
                build(**kwargs)


if __name__ == "__main__":
    unittest.main()
