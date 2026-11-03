"""Independent host/native drift rejection with unchanged SIQS timings."""

import copy
import unittest

from v2.benchmarks.qs.c10a01.regression65_analysis import timing_check


class TimingValidityTests(unittest.TestCase):
    def test_host_or_native_change_cannot_hide_behind_siqs_pass(self):
        power = dict(source="AC Power", low_power=0)
        sample = dict(wall=1, cpu=1)
        before = dict(
            host=dict(power=power, runtime={"qos": 33}),
            wall=12,
            cpu=12,
            samples=[sample.copy() for _ in range(9)],
            native_canary=dict(
                digest="fixed", samples=[sample.copy() for _ in range(9)]
            ),
        )
        self.assertTrue(
            timing_check(before, copy.deepcopy(before), power)["valid"]
        )
        after = copy.deepcopy(before)
        after["host"]["power"]["low_power"] = 1
        check = timing_check(before, after, power)
        self.assertTrue(check["checks"]["siqs_drift"])
        self.assertFalse(check["valid"])
        after = copy.deepcopy(before)
        for sample in after["native_canary"]["samples"]:
            sample["cpu"] = sample["wall"] = 2
        check = timing_check(before, after, power)
        self.assertTrue(check["checks"]["siqs_drift"])
        self.assertFalse(check["valid"])


if __name__ == "__main__":
    unittest.main()
