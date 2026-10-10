"""Residue equivalence, certified fusion and B4 control isolation."""

import unittest
from math import gcd

from v2.benchmarks.ecm.c6 import c6_b4, c6_cf, c6_fast
from v2.benchmarks.ecm.p41.p41_campaign import PYTHON_BACKEND
from v2.benchmarks.support.prac_oracle import (
    affine_multiply,
    matches,
    twist_point,
)
from v2.ecm import core as ecm
from v2.ecm import prac


class CombinedChainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.backend = c6_b4.load_backend(PYTHON_BACKEND)
        cls.catalog = c6_fast.load_catalog()
        cls.catalog["cf"] = {
            scalar: record.certified
            for scalar, record in c6_cf.load_catalog().items()
        }

    def test_pinned_ladder_exact_and_isolated(self):
        self.assertIs(PYTHON_BACKEND.ecm, ecm)
        self.assertIsNot(self.backend.ecm, ecm)
        for n in (101, 101 * 103, 1009**2):
            for x, z in ((0, 1), (1, 0), (3, 1), (101, 103)):
                for scalar in (0, 1, 2, 3, 31, 128, 2000):
                    self.assertEqual(
                        ecm.scalar_multiply(scalar, x, z, n, 7),
                        self.backend.ecm.scalar_multiply(scalar, x, z, n, 7),
                    )

    def test_kernels_are_exact_residues(self):
        for n in (35, 101, 1009**2):
            for x in range(13):
                for z in range(7):
                    p, q, r = (x, z), (z + 1, x + 1), (3, 2)
                    d = ecm.point_double(*p, n, 7)
                    a = ecm.point_add(*p, *q, *r, n)
                    self.assertEqual(c6_b4.point_double(*p, n, 7), d)
                    self.assertEqual(c6_b4.point_add(*p, *q, *r, n), a)
                    self.assertEqual(c6_b4.double_add(p, q, r, n, 7), (d, a))

    def test_every_record_exact_fusion_and_affine(self):
        fused_count = 0
        n, curve_a = 1009, 19 * pow(9, -1, 1009) % 1009
        a24 = (curve_a + 2) * pow(4, -1, n) % n
        for records in self.catalog.values():
            for scalar, certified in records.items():
                original = c6_fast.FastRecord(certified, PYTHON_BACKEND)
                fused = c6_b4.FusedRecord(certified, PYTHON_BACKEND)
                fused_count += sum(row[0] == 2 for row in fused.plan)
                for modulus, point in (
                    (n, (3, 1)),
                    (35, (5, 7)),
                    (1009**2, (3, 1)),
                ):
                    self.assertEqual(
                        fused.run(point, modulus, a24),
                        original.run(point, modulus, a24),
                    )
                expected = affine_multiply(scalar, (3, 7), n, curve_a)
                self.assertTrue(matches(fused((3, 1), n, a24), expected, n))
        self.assertGreater(fused_count, 0)

    def test_full_program_recovery_and_false_infinity(self):
        n, sigma = 33554520197234177, 2046841451
        setup = ecm.setup_curve(n, sigma)
        affine, a, b = twist_point(setup.point, setup.a24, n)
        expected = affine_multiply(ecm.stage_one_scalar(2000), affine, n, a, b)
        for family in ("prac", "lucas", "cf"):
            for mode in ("late", "reduced", "fused"):
                program = c6_b4.build_program(family, mode, self.backend, 16)
                point, factor = program(setup.point, n, setup.a24)
                if factor is None and point is not None:
                    self.assertTrue(matches(point, expected, n))
                elif factor is not None:
                    self.assertTrue(1 < factor < n and n % factor == 0)
                extra = {}
                _, factor = program((5, 7), 35, 2, extra)
                self.assertEqual(factor, 7)
                self.assertEqual(extra["block_replays"], 1)
                self.assertLessEqual(extra["record_replay_allowance"], 16)
                for _, _, action, _ in program.entries[:20]:
                    raw, guard = action.run((5, 7), 35, 2)
                    self.assertNotEqual(gcd(guard * raw[0] * raw[1], 35), 1)
                    with self.assertRaises(prac.NonunitPointError) as found:
                        action((5, 7), 35, 2)
                    self.assertEqual(found.exception.factor, 7)
