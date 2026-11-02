"""Frozen nominee and unrun pilot resource invariants."""

import unittest

from v2.benchmarks.qs.c10a01 import nominee_v4, plan_v4


class PilotV4FreezeTests(unittest.TestCase):
    def test_nominee_configs_and_receipt_bindings(self):
        frozen = nominee_v4.verify()
        self.assertEqual(len(frozen["calls"]), 8)
        self.assertEqual(frozen["pilot_status"], "unrun")
        self.assertEqual(frozen["control_commit"], plan_v4.START)
        self.assertEqual(frozen["policies"]["root_update_policy"], "batch")
        self.assertEqual(
            frozen["policies"]["multiplier_policy"], "parity2_wide"
        )
        for call in frozen["calls"]:
            config = call["configuration"]
            if call["arm"] == "control":
                self.assertNotIn("polynomial_policy", config)
                self.assertNotIn("multiplier_policy", config)
            else:
                self.assertEqual(config["polynomial_policy"], "parity2")
                self.assertEqual(config["multiplier_policy"], "parity2_wide")
                self.assertEqual(config["assignment_policy"], "windowed_pair")

    def test_exact_pilot_reservations(self):
        frozen = nominee_v4.verify()
        self.assertEqual(
            [frozen["stages"][stage]["total"]["wall"] for stage in "ABC"],
            [9200, 7600, 14800],
        )
        self.assertEqual(
            sum(frozen["stages"][stage]["total"]["wall"] for stage in "ABC"),
            31600,
        )
        for stage in "ABC":
            limits = frozen["stages"][stage]
            self.assertEqual(limits["total"]["wall"], limits["total"]["cpu"])
            self.assertEqual(
                limits["total"]["work"],
                sum(p["work"] for p in limits["phases"].values()),
            )


if __name__ == "__main__":
    unittest.main()
