"""Regression for terminal-only materialization of frozen input ratios."""

import unittest

from v2.benchmarks.qs.a1b_terminal import geometric_mean


class A1bTerminalTests(unittest.TestCase):
    def test_geometric_mean_accepts_the_frozen_generator_call(self):
        self.assertAlmostEqual(geometric_mean(x for x in (1, 4, 16)), 4)
        self.assertAlmostEqual(geometric_mean((0.5, 2)), 1)


if __name__ == "__main__":
    unittest.main()
