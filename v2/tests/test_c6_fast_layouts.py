"""Identical-arithmetic coverage for the specialized three-point layout."""

import unittest

from v2 import prac
from v2.benchmarks import c6_chains, c6_fast
from v2.benchmarks.c6_fast_costs import ThreePoint
from v2.benchmarks.p41_campaign import PYTHON_BACKEND
from v2.benchmarks.prac_oracle import (
    affine_multiply,
    historical_points,
    matches,
)


class FastLayoutTests(unittest.TestCase):
    def test_identical_arithmetic_and_exceptional_recovery(self):
        subset = [
            c
            for c in c6_chains.load_lucas().values()
            if c6_chains.continued_fraction_bits(c) is not None
        ]
        self.assertEqual(len(subset), 45)
        controls = tuple(historical_points())
        for chain in subset:
            compact = c6_fast.FastRecord(
                c6_chains.compact(chain), PYTHON_BACKEND
            )
            three = ThreePoint(chain, PYTHON_BACKEND)
            ring = c6_fast.FastRecord(
                c6_chains.compact(chain, rolling=True), PYTHON_BACKEND
            )
            for modulus, curve_a, point in controls:
                a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
                expected = affine_multiply(
                    chain.scalar, point, modulus, curve_a
                )
                for action in (compact, three, ring):
                    self.assertTrue(
                        matches(
                            action((point[0], 1), modulus, a24),
                            expected,
                            modulus,
                        )
                    )
            for action in (compact, three, ring):
                with self.assertRaises(prac.NonunitPointError) as found:
                    action((5, 7), 35, 2)
                self.assertEqual(found.exception.factor, 7)
                self.assertEqual(action((0, 1), 101, 2), (0, 1))


if __name__ == "__main__":
    unittest.main()
