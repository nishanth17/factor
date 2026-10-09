"""C6 independent scalar, projective and finite-recovery contracts."""

import unittest
from dataclasses import replace
from math import gcd
from unittest.mock import patch

from v2 import ecm, prac, prime_sieve, utils
from v2.benchmarks import c6_chains as chains
from v2.benchmarks import c6_study as study
from v2.benchmarks.p41_campaign import PYTHON_BACKEND
from v2.benchmarks.prac_oracle import (
    affine_multiply,
    historical_points,
    matches,
    twist_point,
)


class CompactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = chains.load_lucas()

    def actions(self, prime, scalar=None):
        scalar = scalar or prime
        lucas = chains.compose(self.catalog[prime], scalar)
        for record in (
            chains.compact(prac.get_chain(scalar)),
            chains.compact(lucas),
            chains.compact(lucas, rolling=True),
        ):
            yield chains.Executor(record, PYTHON_BACKEND)

    def test_every_upstream_record_and_prime_power(self):
        for prime, chain in self.catalog.items():
            self.assertTrue(prac.verify_chain(chain))
            for bound in (31, 373, 2000):
                if prime > bound:
                    continue
                power = utils.prime_power(prime, bound)
                for action in self.actions(prime, power):
                    self.assertTrue(chains.verify(action.record))
                    self.assertLessEqual(action.record.slots, 16)
                    self.assertLessEqual(len(action.record.code), 2048)
        program = chains.build_program(2000, "lucas", PYTHON_BACKEND)
        records = {id(action) for row in program for action in row[2:]}
        self.assertLessEqual(len(records), 512)
        with self.assertRaises(ValueError):
            chains.build_program(2001, "lucas", PYTHON_BACKEND)

    def test_corruption_and_finite_workspace(self):
        record = chains.compact(self.catalog[13])
        for broken in (
            replace(record, scalar=17),
            replace(record, slots=17),
            replace(record, code=b"\0"),
            replace(record, scalar=True),
            replace(record, output=255),
            replace(record, code=b"\xff" * 4),
            replace(record, code=record.code * 513),
            replace(record, code=b"\0\0\0\0"),
        ):
            with self.assertRaises(ValueError):
                chains.Executor(broken, PYTHON_BACKEND)
        with self.assertRaises(ValueError):
            chains.decode_rows([{}] * 513)

    def test_three_point_interpreter(self):
        count = 0
        for chain in self.catalog.values():
            if chains.continued_fraction_bits(chain) is None:
                continue
            count += 1
            action = chains.ThreePointExecutor(chain, PYTHON_BACKEND)
            self.assertLessEqual(action.record.slots, 3)
            expected = affine_multiply(chain.scalar, (3, 293), 1009, 6)
            self.assertTrue(matches(action((3, 1), 1009, 2), expected, 1009))
            for point in ((0, 1), (7, 0)):
                result = action(point, 101, 2)
                self.assertNotEqual(result, (0, 0))
        self.assertGreater(count, 0)

    def test_all_records_independent_affine(self):
        for modulus, curve_a, point in historical_points():
            a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
            for prime in self.catalog:
                expected = affine_multiply(prime, point, modulus, curve_a)
                for action in self.actions(prime):
                    actual = action((point[0], 1), modulus, a24)
                    self.assertTrue(matches(actual, expected, modulus))

    def test_all_x_small_curves(self):
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
                    point = (x, roots[rhs])
                    for prime in (2, 3, 5, 7, 11, 13):
                        expected = affine_multiply(
                            prime, point, modulus, curve_a
                        )
                        for action in self.actions(prime):
                            actual = action((x, 1), modulus, a24)
                            self.assertTrue(matches(actual, expected, modulus))

    def test_composites_prime_powers_and_saturation(self):
        factors_seen = saturation_seen = checked = 0
        for primes in ((101, 103), (1009, 1013)):
            n = primes[0] * primes[1]
            for sigma in range(6, 25):
                setup = ecm.setup_curve(n, sigma)
                if setup.factor is not None or setup.retry:
                    continue
                for method in ("compact", "lucas", "rolling"):
                    point = setup.point
                    controls = [
                        twist_point(point, setup.a24, p) for p in primes
                    ]
                    program = chains.build_program(373, method, PYTHON_BACKEND)
                    for _, power, action, _ in program:
                        try:
                            point = action(point, n, setup.a24)
                        except prac.NonunitPointError as result:
                            self.assertTrue(
                                utils.valid_divisor(result.factor, n)
                            )
                            factors_seen += 1
                            break
                        controls = [
                            (affine_multiply(power, q, p, a, b), a, b)
                            for p, (q, a, b) in zip(primes, controls)
                        ]
                        for p, (expected, _, _) in zip(primes, controls):
                            self.assertTrue(matches(point, expected, p))
                        self.assertEqual(gcd(gcd(*point), n), 1)
                        checked += 1
                        if point[1] == 0:
                            saturation_seen += 1
                            break
        self.assertGreater(checked, 0)
        self.assertGreater(factors_seen, 0)
        # A prime-field infinity is valid but must trigger finite stage replay.
        program = chains.build_program(31, "lucas", PYTHON_BACKEND)
        extra = dict(prime_power_replays=0, prime_units_replayed=0)
        point, factor = study.stage_one(
            (0, 1), 101, 2, program, PYTHON_BACKEND, extra
        )
        self.assertIsNone(point)
        self.assertIsNone(factor)
        self.assertEqual(extra["prime_power_replays"], 1)
        self.assertEqual(extra["prime_units_replayed"], 1)

    def test_prime_squares(self):
        for n in (1009**2, 10007**2):
            curve_a = 19 * pow(9, -1, n) % n
            a24 = (curve_a + 2) * pow(4, -1, n) % n
            for prime in list(self.catalog)[:60]:
                try:
                    expected = affine_multiply(prime, (3, 7), n, curve_a)
                except ValueError:
                    continue
                for action in self.actions(prime):
                    try:
                        actual = action((3, 1), n, a24)
                    except prac.NonunitPointError as result:
                        self.assertTrue(utils.valid_divisor(result.factor, n))
                    else:
                        self.assertTrue(matches(actual, expected, n))

    def test_published_false_infinity(self):
        n = 33554520197234177
        setup = ecm.setup_curve(n, 2046841451)
        for method in ("compact", "lucas", "rolling"):
            point = setup.point
            expected, curve_a, curve_b = twist_point(point, setup.a24, n)
            for _, power, action, _ in chains.build_program(
                373, method, PYTHON_BACKEND
            ):
                expected = affine_multiply(
                    power, expected, n, curve_a, curve_b
                )
                point = action(point, n, setup.a24)
                self.assertTrue(matches(point, expected, n))
            self.assertNotEqual(point[1], 0)

    def test_nonunits_single_retry_and_discarded_factor(self):
        action = next(self.actions(3))
        for point in ((5, 1), (1, 5), (5, 0)):
            with self.assertRaises(prac.NonunitPointError) as result:
                action(point, 35, 2)
            self.assertEqual(result.exception.factor, 5)
        with self.assertRaises(ValueError):
            action((0, 0), 35, 2)
        with patch.object(ecm, "point_double", return_value=(0, 0)) as double:
            with self.assertRaises(prac.NonunitPointError) as result:
                action((3, 1), 101, 2)
            self.assertIsNone(result.exception.factor)
            self.assertEqual(double.call_count, 2)
        with patch.object(ecm, "point_double", return_value=(5, 1)):
            with patch.object(ecm, "point_add", return_value=(0, 0)):
                with self.assertRaises(prac.NonunitPointError) as result:
                    action((3, 1), 35, 2)
                self.assertEqual(result.exception.factor, 5)

    def test_backends_match_prime_power_action(self):
        try:
            backend = study.control.gmp_backend()
        except ImportError:
            self.skipTest("optional GMP unavailable")
        for p in prime_sieve.prime_sieve(100):
            for original in self.actions(p, utils.prime_power(p, 2000)):
                other = chains.Executor(original.record, backend)
                left = original((3, 1), 1009, 2)
                right = other(
                    (backend.integer(3), backend.integer(1)),
                    backend.integer(1009),
                    backend.integer(2),
                )
                self.assertEqual(tuple(map(int, right)), left)


if __name__ == "__main__":
    unittest.main()
