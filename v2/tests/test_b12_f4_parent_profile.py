"""Independent nesting and refusal oracles for the diagnostic recorder."""

import unittest

from v2.benchmarks.ecm.b12_f4.parent_profile import (
    CoarseCounter,
    RecorderLimitError,
)


class ParentProfileTests(unittest.TestCase):
    def test_nested_inclusive_times_do_not_double_count_exclusive_cpu(self):
        cpu, wall = iter(range(4)), iter(range(4))
        recorder = CoarseCounter(
            cpu_clock=lambda: next(cpu), wall_clock=lambda: next(wall)
        )
        child = recorder.wrap("child", lambda: 17)
        parent = recorder.wrap("parent", lambda: child())
        recorder.enabled = True
        self.assertEqual(parent(), 17)
        self.assertEqual(recorder.rows["child"], [1, 1, 1, 1, 1])
        self.assertEqual(recorder.rows["parent"], [1, 3, 2, 3, 2])
        self.assertFalse(recorder.stack)

    def test_disabled_counter_preserves_behavior_without_clock_observation(
        self,
    ):
        recorder = CoarseCounter(cpu_clock=lambda: self.fail("unused clock"))
        self.assertEqual(recorder.wrap("f", lambda x: x + 1)(3), 4)
        self.assertFalse(recorder.rows)
        self.assertEqual(recorder.events, 0)

    def test_failure_unwinds_and_observer_exhaustion_has_a_distinct_type(self):
        recorder = CoarseCounter(maximum=1)
        recorder.enabled = True

        def fail():
            raise ValueError("production exception")

        wrapped = recorder.wrap("f", fail)
        with self.assertRaisesRegex(ValueError, "production exception"):
            wrapped()
        self.assertFalse(recorder.stack)
        with self.assertRaisesRegex(RecorderLimitError, "recorder_call_limit"):
            wrapped()
