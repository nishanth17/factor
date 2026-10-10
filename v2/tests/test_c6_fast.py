"""Independent coverage, projective-action and bounded replay checks."""

import json
import unittest
from dataclasses import FrozenInstanceError, replace
from math import gcd

from v2.benchmarks.ecm.c6 import c6_chains as strict
from v2.benchmarks.ecm.c6 import c6_fast as fast
from v2.benchmarks.ecm.c6 import c6_study as study
from v2.benchmarks.ecm.p41.p41_campaign import PYTHON_BACKEND, gmp_backend
from v2.benchmarks.support.prac_oracle import (
    affine_multiply,
    historical_points,
    matches,
    twist_point,
)
from v2.ecm import core as ecm
from v2.ecm import prac


class FastChainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = fast.load_catalog()

    def test_every_integer_and_coverage_certificate(self):
        for family, records in self.catalog.items():
            self.assertLessEqual(len(records), 512)
            for scalar, certified in records.items():
                self.assertEqual(certified.record.scalar, scalar)
                self.assertTrue(strict.verify(certified.record))
                self.assertTrue(
                    fast.verify_frontier(certified.record, certified.masks)
                )
                self.assertEqual(
                    fast.frontier(certified.record), certified.masks
                )
                for index, mask in enumerate(certified.masks):
                    for bit in (1, 2):
                        if not mask & bit:
                            continue
                        broken = bytearray(certified.masks)
                        broken[index] ^= bit
                        with self.assertRaises(ValueError):
                            fast.CertifiedRecord(
                                certified.record, bytes(broken)
                            )
        certified = self.catalog["lucas"][13]
        with self.assertRaises(FrozenInstanceError):
            certified.masks = b""
        with self.assertRaises(ValueError):
            fast.CertifiedRecord(certified.record, b"\xff" * 512)
        with self.assertRaises(ValueError):
            fast.CertifiedRecord(
                replace(certified.record, scalar=17), certified.masks
            )

    def test_all_records_against_independent_affine(self):
        controls = tuple(historical_points())
        for records in self.catalog.values():
            for scalar, certified in records.items():
                actions = [
                    fast.FastRecord(certified, PYTHON_BACKEND, mode)
                    for mode in fast.MODES
                ]
                for modulus, curve_a, point in controls:
                    a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
                    expected = affine_multiply(scalar, point, modulus, curve_a)
                    for action in actions:
                        actual = action((point[0], 1), modulus, a24)
                        self.assertTrue(matches(actual, expected, modulus))

    def test_exhaustive_small_curves(self):
        for modulus in (5, 7, 11, 13):
            for curve_a in range(modulus):
                if (curve_a * curve_a - 4) % modulus == 0:
                    continue
                a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
                roots = {y * y % modulus: y for y in range(modulus)}
                for x in range(modulus):
                    rhs = (x**3 + curve_a * x * x + x) % modulus
                    if rhs not in roots:
                        continue
                    for scalar in (2, 3, 5, 7, 11, 13):
                        expected = affine_multiply(
                            scalar, (x, roots[rhs]), modulus, curve_a
                        )
                        for family in self.catalog:
                            for mode in fast.MODES:
                                action = fast.FastRecord(
                                    self.catalog[family][scalar],
                                    PYTHON_BACKEND,
                                    mode,
                                )
                                self.assertTrue(
                                    matches(
                                        action((x, 1), modulus, a24),
                                        expected,
                                        modulus,
                                    )
                                )

    def test_saturated_product_keeps_separate_factors(self):
        for records in self.catalog.values():
            for mode in fast.MODES:
                action = fast.FastRecord(records[3], PYTHON_BACKEND, mode)
                result, product = action.run((5, 7), 35, 2)
                self.assertEqual(gcd(product * result[0] * result[1], 35), 35)
                with self.assertRaises(prac.NonunitPointError) as found:
                    action((5, 7), 35, 2)
                self.assertEqual(found.exception.factor, 7)
                with self.assertRaises(ValueError):
                    action((0, 0), 35, 2)
                for point in ((5, 1), (1, 5), (5, 0)):
                    with self.assertRaises(prac.NonunitPointError) as found:
                        action(point, 35, 2)
                    self.assertEqual(found.exception.factor, 5)

    def test_composite_coverage_and_exact_strict_recovery(self):
        recovered = 0
        for n in (101 * 103, 1009 * 1013, 1009**2):
            for sigma in range(6, 18):
                setup = ecm.setup_curve(n, sigma)
                if setup.factor is not None or setup.retry:
                    continue
                for records in self.catalog.values():
                    for scalar in list(records)[:60]:
                        certified = records[scalar]
                        reference = strict.Executor(
                            certified.record, PYTHON_BACKEND
                        )
                        try:
                            expected = reference(setup.point, n, setup.a24)
                            factor = None
                        except prac.NonunitPointError as found:
                            expected, factor = None, found.factor
                        for mode in fast.MODES:
                            action = fast.FastRecord(
                                certified, PYTHON_BACKEND, mode
                            )
                            point, product = action.run(
                                setup.point, n, setup.a24
                            )
                            if expected is None:
                                self.assertNotEqual(
                                    gcd(product * point[0] * point[1], n), 1
                                )
                                with self.assertRaises(
                                    prac.NonunitPointError
                                ) as found:
                                    action(setup.point, n, setup.a24)
                                self.assertEqual(
                                    found.exception.factor, factor
                                )
                                recovered += 1
                            else:
                                self.assertEqual(
                                    action(setup.point, n, setup.a24), expected
                                )
        self.assertGreater(recovered, 0)

    def test_prime_squares_against_affine(self):
        for n in (1009**2, 10007**2):
            curve_a = 19 * pow(9, -1, n) % n
            a24 = (curve_a + 2) * pow(4, -1, n) % n
            for scalar, certified in list(self.catalog["lucas"].items())[:60]:
                try:
                    expected = affine_multiply(scalar, (3, 7), n, curve_a)
                except ValueError:
                    continue
                for mode in fast.MODES:
                    action = fast.FastRecord(certified, PYTHON_BACKEND, mode)
                    try:
                        result = action((3, 1), n, a24)
                    except prac.NonunitPointError as found:
                        self.assertTrue(1 < found.factor < n)
                        self.assertEqual(n % found.factor, 0)
                    else:
                        self.assertTrue(matches(result, expected, n))

    def test_false_infinity_and_batched_scalar_action(self):
        n, sigma = 33554520197234177, 2046841451
        setup = ecm.setup_curve(n, sigma)
        expected, curve_a, curve_b = twist_point(setup.point, setup.a24, n)
        expected = affine_multiply(
            ecm.stage_one_scalar(373), expected, n, curve_a, curve_b
        )
        for family in self.catalog:
            for mode in fast.MODES:
                for batch in (1, 16, 64):
                    program = fast.build_program(
                        373, family, mode, PYTHON_BACKEND, batch, self.catalog
                    )
                    point, factor = program(setup.point, n, setup.a24)
                    self.assertIsNone(factor)
                    self.assertNotEqual(point[1], 0)
                    self.assertTrue(matches(point, expected, n))

    def test_finite_block_replay_and_infinity(self):
        for mode in fast.MODES:
            program = fast.build_program(
                31, "lucas", mode, PYTHON_BACKEND, 16, self.catalog
            )
            extra = {}
            point, factor = program((0, 1), 101, 2, extra)
            self.assertIsNone(point)
            self.assertIsNone(factor)
            self.assertEqual(extra["block_replays"], 1)
            self.assertLessEqual(extra["record_replay_allowance"], 16)
            self.assertEqual(extra["prime_power_replays"], 1)
            self.assertEqual(extra["prime_units_replayed"], 1)
            point, factor = program((5, 7), 35, 2)
            self.assertIsNone(point)
            self.assertEqual(factor, 7)
        for kwargs in (dict(bound=2001), dict(batch=65), dict(batch=True)):
            values = dict(
                bound=31,
                family="lucas",
                mode="tuple",
                backend=PYTHON_BACKEND,
                batch=16,
                catalog=self.catalog,
            )
            values.update(kwargs)
            with self.assertRaises((ValueError, TypeError)):
                fast.build_program(**values)

    def test_full_stage_certificate_and_gcd_reduction(self):
        corpus = json.loads(study.control.CORPUS.read_text())
        case = next(c for c in corpus["fixtures"] if c["id"] == "balanced_40d")
        n, seed = case["n"], 41001
        from v2 import constants
        from v2.common import utils

        rng = utils.resolve_rng(seed, None)
        setup = ecm.setup_curve(n, rng.randint(6, constants.MAX_RANDOM_ECM))
        controls = study.affine_controls(n, tuple(case["factors"]), seed, 2000)
        for family in self.catalog:
            for mode in fast.MODES:
                program = fast.build_program(
                    2000, family, mode, PYTHON_BACKEND, 16, self.catalog
                )
                extra = {}
                point, factor = program(setup.point, n, setup.a24, extra)
                self.assertIsNone(factor)
                self.assertEqual(extra["block_replays"], 0)
                self.assertEqual(extra["guard_batches"], 19)
                for prime, expected in controls:
                    self.assertTrue(matches(point, expected, prime))

    def test_gmp_and_int_exact_action(self):
        try:
            backend = gmp_backend()
        except ImportError:
            self.skipTest("optional GMP unavailable")
        for family in self.catalog:
            for scalar, certified in list(self.catalog[family].items())[:60]:
                for mode in fast.MODES:
                    native = fast.FastRecord(certified, PYTHON_BACKEND, mode)
                    other = fast.FastRecord(certified, backend, mode)
                    expected = native((3, 1), 1009, 2)
                    actual = other(
                        (backend.integer(3), backend.integer(1)),
                        backend.integer(1009),
                        backend.integer(2),
                    )
                    self.assertEqual(tuple(map(int, actual)), expected)


if __name__ == "__main__":
    unittest.main()
