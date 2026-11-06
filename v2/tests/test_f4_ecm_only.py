"""ECM-only measurements reject fallback and invalid split publication."""

import unittest

from v2.benchmarks.ecm.b12_f4.selection_fixed_child import (
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
