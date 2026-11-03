"""Power-state parsing and timing-validity tests without system mutation."""

import copy
import unittest

from v2.benchmarks.qs.c10a01 import host_state


class HostStateTests(unittest.TestCase):
    def test_active_profile_wins_over_other_profile(self):
        profiles = (
            "Battery Power:\n lowpowermode 0\nAC Power:\n lowpowermode 1\n"
        )
        self.assertEqual(
            host_state.parse_power("Now drawing from 'AC Power'\n", profiles),
            dict(source="AC Power", low_power=1),
        )
        self.assertEqual(
            host_state.parse_power(
                "Now drawing from 'Battery Power'\n", profiles
            ),
            dict(source="Battery Power", low_power=0),
        )
        with self.assertRaises(RuntimeError):
            host_state.parse_power("unavailable", profiles)

    def test_mode_or_runtime_change_invalidates_comparison(self):
        before = dict(
            power=dict(source="AC Power", low_power=0),
            runtime=dict(qos=33, nursery="4.0MB", environment={}),
        )
        self.assertTrue(host_state.matches(before, copy.deepcopy(before)))
        for key, field, value in (
            ("power", "low_power", 1),
            ("power", "source", "Battery Power"),
            ("runtime", "qos", 9),
            ("runtime", "nursery", "8.0MB"),
        ):
            after = copy.deepcopy(before)
            after[key][field] = value
            self.assertFalse(host_state.matches(before, after))


if __name__ == "__main__":
    unittest.main()
