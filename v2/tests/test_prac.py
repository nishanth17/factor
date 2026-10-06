"""Integer proofs, independent group oracles and exceptional PRAC recovery."""

import random
import unittest
from dataclasses import replace
from math import gcd
from unittest.mock import patch

from v2 import ecm, prac, prime_sieve, utils
from v2.benchmarks.prac_oracle import (
    affine_add,
    affine_multiply,
    historical_points,
    matches,
    twist_point,
)


class ChainTests(unittest.TestCase):
    def test_integer_records_and_prime_power_schedules(self):
        scalars = set(range(10001))
        for bound in (31, 1000, 10000, 100000):
            scalars.update(
                utils.prime_power(p, bound)
                for p in prime_sieve.prime_sieve(bound + 1)
            )
        generator = random.Random(4101)
        scalars.update(generator.randrange(2**31, 2**32) for _ in range(200))
        for scalar in sorted(scalars):
            chain = prac.get_chain(scalar)
            self.assertTrue(prac.verify_chain(chain))
            self.assertLessEqual(len(chain.instructions), prac.MAX_CHAIN_STEPS)
            self.assertLessEqual(
                chain.cost(), prac._binary_chain(scalar).cost()
            )
            if chain.split:
                odd = scalar >> ((scalar & -scalar).bit_length() - 1)
                self.assertLess(odd // 2, chain.split)
                self.assertLess(chain.split, odd)
                self.assertEqual(gcd(odd, chain.split), 1)

    def test_boundary_scalars_and_bounded_cache(self):
        prac.clear_cache()
        for exponent in range(32):
            chain = prac.get_chain(1 << exponent)
            self.assertEqual(len(chain.instructions), exponent)
            self.assertTrue(all(step[2] == -1 for step in chain.instructions))
        for scalar in range(1000):
            prac.get_chain(scalar)
        info = prac.cache_info()
        self.assertEqual(info.currsize, prac.CACHE_SIZE)
        self.assertIs(prac.get_chain(997), prac.get_chain(997))
        self.assertIsNone(prac.get_chain(1 << 10000))
        self.assertEqual(prac.cache_info().currsize, info.currsize)
        for scalar in (-1, 1.5, True):
            with self.assertRaises((ValueError, TypeError)):
                prac.get_chain(scalar)
        for weight in (0, -1, 1.5, True, 2**32):
            with self.assertRaises((ValueError, TypeError)):
                prac.get_chain(97, add_cost=weight)
        self.assertTrue(prac.verify_chain(prac.get_chain(97, add_cost=17)))

    def test_corrupt_records_rejected_before_point_execution(self):
        chain = prac.get_chain(97)
        corrupt = (
            replace(chain, scalar=99),
            replace(chain, output=10000),
            replace(chain, output=True),
            replace(chain, instructions=list(chain.instructions)),
            replace(chain, instructions=((0, 5, 0),)),
            replace(chain, instructions=((0, 0, -2),)),
            replace(chain, instructions=((0, 0, True),)),
            replace(chain, instructions=((0, 0, -1), (0, 1, -1))),
            replace(chain, instructions=((0, 0, -1), (1, 1, 0))),
            replace(chain, instructions=((0, 0, -1),) * 513),
            replace(chain, split=97),
            replace(prac.get_chain(9), split=6),
        )
        for record in corrupt:
            with self.assertRaises(ValueError):
                prac.multiply(97, (3, 1), 1009, 2, None, None, chain=record)
        with self.assertRaises(ValueError):
            prac.multiply(99, (3, 1), 1009, 2, None, None, chain=chain)

    def test_independent_external_record(self):
        # A Lucas record for 13: 1,2,3,5,8,13. No PRAC metadata is trusted.
        chain = prac.Chain(
            13,
            ((0, 0, -1), (0, 1, 0), (1, 2, 0), (2, 3, 1), (3, 4, 2)),
            5,
            "lucas",
        )
        self.assertTrue(prac.verify_chain(chain))
        actual = prac.multiply(
            13, (3, 1), 1009, 2, ecm.point_add, ecm.point_double, chain=chain
        )
        self.assertTrue(
            matches(actual, affine_multiply(13, (3, 293), 1009, 6), 1009)
        )


class PointTests(unittest.TestCase):
    def test_all_x_points_on_small_nonsingular_curves(self):
        for modulus in (5, 7, 11, 13, 17):
            roots = {y * y % modulus: y for y in range(modulus)}
            for curve_a in range(modulus):
                if (curve_a * curve_a - 4) % modulus == 0:
                    continue
                a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
                for x in range(modulus):
                    square = (x**3 + curve_a * x * x + x) % modulus
                    if square not in roots:
                        continue
                    point = x, roots[square]
                    expected = None
                    for scalar in range(64):
                        if scalar:
                            expected = affine_add(
                                expected, point, modulus, curve_a
                            )
                        actual = ecm.multiply_prac(scalar, x, 1, modulus, a24)
                        self.assertTrue(
                            matches(actual, expected, modulus),
                            (modulus, curve_a, point, scalar),
                        )

    def test_original_16016_comparisons_including_exceptional_pairs(self):
        fixtures = list(historical_points())
        expected = [None] * len(fixtures)
        checks = 0
        fallback = prac._checked_ladder
        with patch.object(prac, "_checked_ladder", wraps=fallback) as recovery:
            for scalar in range(1001):
                for index, (modulus, curve_a, point) in enumerate(fixtures):
                    if scalar:
                        expected[index] = affine_add(
                            expected[index], point, modulus, curve_a
                        )
                    a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
                    actual = ecm.multiply_prac(
                        scalar, point[0], 1, modulus, a24
                    )
                    self.assertTrue(
                        matches(actual, expected[index], modulus),
                        (scalar, modulus, point, actual),
                    )
                    checks += 1
            self.assertGreater(recovery.call_count, 0)
        self.assertEqual(checks, 16016)

    def test_large_fields_scalars_and_projective_rescaling(self):
        generator = random.Random(4102)
        for modulus in (2**61 - 1, 2**127 - 1, 2**255 - 19, 2**521 - 1):
            point = (3, 7)
            curve_a = (49 - 27 - 3) * pow(9, -1, modulus) % modulus
            a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
            scalars = [0, 1, 2, 9, 2**31, 2**32 - 1, 2**32, 2**127 + 7]
            scalars += [generator.randrange(2**32) for _ in range(64)]
            for scalar in scalars:
                scale = generator.randrange(1, modulus)
                actual = ecm.multiply_prac(
                    scalar, 3 * scale, scale, modulus, a24
                )
                expected = affine_multiply(scalar, point, modulus, curve_a)
                self.assertTrue(matches(actual, expected, modulus))
                ladder = ecm.scalar_multiply(
                    scalar, 3 * scale, scale, modulus, a24
                )
                self.assertTrue(matches(ladder, expected, modulus))

    def test_complete_suyama_prime_power_schedules_over_composites(self):
        checked = factors = saturated = 0
        for primes in (
            (1009,),
            (1009, 1013),
            (10007, 10009),
            (2**61 - 1, 2**127 - 1),
        ):
            modulus = 1
            for prime in primes:
                modulus *= prime
            for sigma in range(6, 16):
                setup = ecm.setup_curve(modulus, sigma)
                if setup.factor or setup.retry:
                    continue
                for bound in (31, 128, 1000):
                    point = setup.point
                    controls = [
                        twist_point(point, setup.a24, p) for p in primes
                    ]
                    for prime in prime_sieve.prime_sieve(bound + 1):
                        power = utils.prime_power(prime, bound)
                        try:
                            point = ecm.multiply_prac(
                                power, *point, modulus, setup.a24
                            )
                        except prac.NonunitPointError as result:
                            self.assertTrue(
                                utils.valid_divisor(result.factor, modulus)
                            )
                            factors += 1
                            break
                        controls = [
                            (affine_multiply(power, q, p, a, b), a, b)
                            for p, (q, a, b) in zip(primes, controls)
                        ]
                        for p, (expected, _, _) in zip(primes, controls):
                            self.assertTrue(matches(point, expected, p))
                        self.assertEqual(gcd(gcd(*point), modulus), 1)
                        checked += 1
                        if point[1] == 0:
                            saturated += 1
                            break
        self.assertGreater(checked, 1000)
        self.assertGreater(factors, 0)
        self.assertGreater(saturated, 0)

    def test_prime_square_moduli(self):
        for modulus in (1009**2, 10007**2):
            curve_a = 19 * pow(9, -1, modulus) % modulus
            a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
            for scalar in range(80):
                try:
                    expected = affine_multiply(
                        scalar, (3, 7), modulus, curve_a
                    )
                except ValueError:
                    continue  # Affine division is not defined at a nonunit.
                try:
                    actual = ecm.multiply_prac(scalar, 3, 1, modulus, a24)
                except prac.NonunitPointError as result:
                    self.assertTrue(
                        utils.valid_divisor(result.factor, modulus)
                    )
                else:
                    self.assertTrue(matches(actual, expected, modulus))

    def test_gmp_ecm_published_false_infinity_regression(self):
        # GMP-ECM ecm.c warns that unguarded PRAC can falsely reach infinity
        # for this Suyama curve at B1=373: the order contains 23**2 > B1.
        modulus = 33554520197234177
        setup = ecm.setup_curve(modulus, 2046841451)
        self.assertIsNone(setup.factor)
        self.assertFalse(setup.retry)
        expected, curve_a, curve_b = twist_point(
            setup.point, setup.a24, modulus
        )
        actual = setup.point
        for prime in prime_sieve.prime_sieve(374):
            power = utils.prime_power(prime, 373)
            expected = affine_multiply(
                power, expected, modulus, curve_a, curve_b
            )
            actual = ecm.multiply_prac(power, *actual, modulus, setup.a24)
            self.assertTrue(matches(actual, expected, modulus))
        self.assertIsNotNone(expected)
        self.assertNotEqual(actual[1], 0)

    def test_nonunits_degenerate_states_and_single_retry(self):
        for x, z in ((5, 5), (1, 5), (5, 0)):
            with self.assertRaises(prac.NonunitPointError) as result:
                ecm.multiply_prac(9, x, z, 35, 2)
            self.assertEqual(result.exception.factor, 5)
        with self.assertRaises(ValueError):
            ecm.multiply_prac(9, 35, 70, 35, 2)
        for scalar in range(12):
            self.assertEqual(ecm.multiply_prac(scalar, 7, 0, 101, 2), (1, 0))
            self.assertEqual(
                ecm.multiply_prac(scalar, 0, 9, 101, 2),
                (0, 9) if scalar % 2 else (1, 0),
            )
        with patch.object(ecm, "point_double", return_value=(0, 0)) as double:
            with self.assertRaises(prac.NonunitPointError) as result:
                ecm.multiply_prac(9, 3, 1, 101, 2)
            self.assertIsNone(result.exception.factor)
            self.assertEqual(double.call_count, 2)
        # A common-zero output must not erase a retained nonunit X. Inject
        # the collapse to isolate this recovery contract from chain choice.
        with patch.object(ecm, "point_double", return_value=(5, 1)):
            with patch.object(ecm, "point_add", return_value=(0, 0)):
                with self.assertRaises(prac.NonunitPointError) as result:
                    ecm.multiply_prac(3, 3, 1, 35, 2)
                self.assertEqual(result.exception.factor, 5)
        self.assertFalse(matches((0, 0), None, 101))
        self.assertFalse(matches((0, 0), (3, 7), 101))
        self.assertFalse(matches((5, 10), (3, 7), 35))


if __name__ == "__main__":
    unittest.main()
