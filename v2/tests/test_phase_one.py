"""Independent oracles and adversarial checks for every Phase 1 gate."""

import hashlib
import io
import json
import random
import unittest
from bisect import bisect_left
from contextlib import redirect_stdout
from math import gcd, isqrt
from pathlib import Path
from unittest.mock import patch

from v2 import ecm, pollard_pm1, pollard_rho, prime_sieve, utils
from v2.factor import factorize, factorize_bf, print_factorization


def reference_primes(hi):
    """Full integer-index sieve, independent of odd/residue production code."""
    flags = bytearray(b"\x01") * max(hi, 2)
    flags[:2] = b"\x00\x00"
    for prime in range(2, isqrt(max(hi - 1, 0)) + 1):
        if flags[prime]:
            for composite in range(prime * prime, hi, prime):
                flags[composite] = 0
    return [value for value in range(2, hi) if flags[value]]


def affine_add(point, other, modulus, curve_a):
    """Independent full-coordinate group law on a small prime field."""
    if point is None:
        return other
    if other is None:
        return point
    x, y = point
    u, v = other
    if x == u and (y + v) % modulus == 0:
        return None
    if point == other:
        numerator = 3 * x * x + 2 * curve_a * x + 1
        denominator = 2 * y
    else:
        numerator, denominator = v - y, u - x
    slope = numerator * pow(denominator, -1, modulus) % modulus
    result_x = (slope * slope - curve_a - x - u) % modulus
    result_y = (slope * (x - result_x) - y) % modulus
    return result_x, result_y


class SequenceRandom:
    """Reproducible generator with draw counts for round/retry tests."""

    def __init__(self, values):
        self.values = iter(values)
        self.calls = 0

    def randint(self, low, high):
        value = next(self.values)
        if not low <= value <= high:
            raise AssertionError("test RNG returned an out-of-range value")
        self.calls += 1
        return value


class UtilityTests(unittest.TestCase):
    def test_exact_prime_powers(self):
        for prime in (2, 3, 5, 7, 11, 13, 17, 19):
            for exponent in range(1, 15):
                bound = prime**exponent
                self.assertEqual(utils.prime_power(prime, bound), bound)
                if exponent > 1:
                    self.assertEqual(
                        utils.prime_power(prime, bound - 1),
                        prime ** (exponent - 1),
                    )

    def test_bisect_empty_duplicates(self):
        self.assertEqual(utils.binary_search(3, []), 0)
        self.assertEqual(utils.binary_search(2, [1, 2, 2, 2, 3]), 4)
        self.assertEqual(utils.binary_search(2, [1, 2, 2, 2, 3], True), 1)

    def test_bezout_and_inverse_contracts(self):
        for a in range(-12, 13):
            for b in range(-12, 13):
                divisor, x, y = utils.extended_gcd(a, b)
                self.assertEqual(divisor, gcd(a, b))
                self.assertEqual(a * x + b * y, divisor)
                self.assertEqual(utils.xgcd(a, b), x)
        for modulus in (7, 60, 101, 1009):
            for value in range(1, 100):
                if gcd(value, modulus) == 1:
                    self.assertEqual(
                        value
                        * utils.modular_inverse(value, modulus)
                        % modulus,
                        1,
                    )
                else:
                    with self.assertRaises(ValueError):
                        utils.modular_inverse(value, modulus)

    def test_batch_recovery(self):
        self.assertEqual(utils.batch_factor([5, 7], 35), (5, True))
        self.assertEqual(utils.batch_factor([0, 1], 35), (None, True))
        self.assertEqual(utils.batch_factor([2, 3], 35), (None, False))
        self.assertEqual(utils.batch_factor([], 35), (None, False))

    def test_invalid_divisors(self):
        for value in (None, -1, 0, 1, True, 35, 4, 5.0):
            self.assertFalse(utils.valid_divisor(value, 35))
        self.assertTrue(utils.valid_divisor(7, 35))

    def test_primality_small_sweep(self):
        primes = set(reference_primes(20_001))
        for n in range(-3, 20_001):
            self.assertEqual(utils.is_prime(n), n in primes, n)
            self.assertIsInstance(utils.is_prime_fast(n), bool)

    def test_pseudoprimes_and_deterministic_boundary(self):
        for n in (
            341,
            561,
            1105,
            1729,
            3215031751,
            341550071728321,
            3825123056546413051,
            2**64 - 1,
            2**64 + 1,
            10**36 + 1,
        ):
            self.assertEqual(
                utils.classify_prime(n), utils.Primality.COMPOSITE
            )
        self.assertEqual(
            utils.classify_prime(2**64 - 59), utils.Primality.PROVEN
        )
        self.assertEqual(
            utils.classify_prime(2**127 - 1, rng=random.Random(1)),
            utils.Primality.PROBABLE,
        )

    def test_requested_probabilistic_rounds(self):
        for n in (1009, 2**127 - 1):
            for rounds in (1, 3, 9):
                rng = SequenceRandom([2] * rounds)
                result = utils.classify_prime(
                    n, use_probabilistic=True, tolerance=rounds, rng=rng
                )
                self.assertEqual(result, utils.Primality.PROBABLE)
                self.assertEqual(rng.calls, rounds)
        self.assertEqual(utils.classify_prime(29), utils.Primality.PROVEN)

    def test_primality_witness_stops_early(self):
        rng = SequenceRandom([2])
        self.assertEqual(
            utils.classify_prime(
                1009 * 1013, use_probabilistic=True, tolerance=8, rng=rng
            ),
            utils.Primality.COMPOSITE,
        )
        self.assertEqual(rng.calls, 1)

    def test_integer_argument_validation(self):
        for n in (1.5, "3", True):
            with self.assertRaises(TypeError):
                utils.classify_prime(n)
        for rounds in (0, -1):
            with self.assertRaises(ValueError):
                utils.classify_prime(101, tolerance=rounds)


class SieveTests(unittest.TestCase):
    def test_small_endpoints_and_tiny_atkin(self):
        for hi in range(-3, 501):
            expected = reference_primes(hi)
            for sieve in (
                prime_sieve.small_sieve,
                prime_sieve.prime_sieve,
                prime_sieve.sieve_of_atkin,
            ):
                self.assertEqual(sieve(hi), expected, (sieve.__name__, hi))

    def test_all_small_intervals(self):
        primes = reference_primes(130)
        for lo in range(-2, 101):
            for hi in range(-2, 103):
                expected = [p for p in primes if lo <= p < hi]
                self.assertEqual(prime_sieve.segmented_sieve(lo, hi), expected)

    def test_prime_square_is_marked_not_excluded(self):
        self.assertEqual(
            prime_sieve.segmented_sieve(3700, 3722), [3701, 3709, 3719]
        )
        for square in (25, 49, 121, 169, 289, 3721):
            self.assertNotIn(
                square, prime_sieve.segmented_sieve(square, square + 1)
            )

    def test_random_intervals_and_tiny_segments(self):
        primes = reference_primes(2_000_001)
        generator = random.Random(7301)
        for index in range(300):
            lo = generator.randrange(0, 1_990_000)
            hi = lo + generator.randrange(0, 5000)
            first = bisect_left(primes, lo)
            last = bisect_left(primes, hi)
            expected = primes[first:last]
            size = (1, 3, 16, 128, 65536)[index % 5]
            self.assertEqual(
                prime_sieve.segmented_sieve(lo, hi, segment_size=size),
                expected,
            )

    def test_dispatch_boundary_and_complete_atkin_values(self):
        reference = reference_primes(3_500_002)
        for hi in (3_499_999, 3_500_000, 3_500_001, 3_500_002):
            expected = reference[: bisect_left(reference, hi)]
            self.assertEqual(prime_sieve.prime_sieve(hi), expected)
        actual = prime_sieve.sieve_of_atkin(3_500_001)
        self.assertEqual(
            actual, reference[: bisect_left(reference, 3_500_001)]
        )
        self.assertEqual(len(actual), 250150)
        self.assertEqual(actual[17], 61)
        self.assertTrue(all(p < 3_500_001 for p in actual))

    def test_independent_atkin_state(self):
        first = prime_sieve.sieve_of_atkin(10000)
        prime_sieve.sieve_of_atkin(100)
        self.assertEqual(first, prime_sieve.sieve_of_atkin(10000))
        self.assertEqual(first, reference_primes(10000))


class RhoTests(unittest.TestCase):
    def test_legacy_first_attempt_failure_retries(self):
        # Legacy offset 499 modulo 35 is 9; y=4 saturates immediately.
        rng = SequenceRandom([4, 9, 2, 1])
        stats = pollard_rho.RhoStats()
        divisor = pollard_rho.factorize_rho(
            35,
            rng=rng,
            max_attempts=2,
            batch_size=34,
            max_evaluations=200,
            stats=stats,
        )
        self.assertTrue(utils.valid_divisor(divisor, 35))
        self.assertEqual(stats.attempts, 2)
        self.assertEqual(rng.calls, 4)

    def test_exhaustion_and_exact_work_caps(self):
        stats = pollard_rho.RhoStats()
        result = pollard_rho.factorize_rho(
            35,
            rng=SequenceRandom([4, 9] * 3),
            max_attempts=3,
            max_evaluations=5,
            batch_size=34,
            stats=stats,
        )
        self.assertIsNone(result)
        self.assertEqual(stats.attempts, 3)
        self.assertLessEqual(stats.evaluations, 15)
        stats = pollard_rho.RhoStats()
        self.assertIsNone(
            pollard_rho.factorize_rho(
                1009 * 1013,
                seed=2,
                max_attempts=3,
                max_evaluations=1,
                stats=stats,
            )
        )
        self.assertEqual(stats.evaluations, 3)

    def test_zero_budget_and_prime_inputs(self):
        for n in (1, 2, 3, 101):
            self.assertIsNone(pollard_rho.factorize_rho(n, seed=1))
        self.assertIsNone(pollard_rho.factorize_rho(35, max_attempts=0))
        self.assertIsNone(pollard_rho.factorize_rho(35, max_evaluations=0))

    def test_saturated_batch_recovery_is_bounded(self):
        stats = pollard_rho.RhoStats()
        result = pollard_rho.factorize_rho(
            35,
            rng=SequenceRandom([4, 9]),
            max_attempts=1,
            max_evaluations=4,
            batch_size=34,
            recovery_limit=0,
            stats=stats,
        )
        self.assertIsNone(result)
        self.assertEqual(stats.saturated_batches, 1)
        self.assertEqual(stats.evaluations, 2)

    def test_seeded_semiprimes(self):
        for n in (35, 77, 143, 10403, 25013 * 25031, 1000003 * 1000033):
            for seed in range(5):
                result = pollard_rho.factorize_rho(n, seed=seed)
                self.assertTrue(utils.valid_divisor(result, n), (n, seed))


class Pm1Tests(unittest.TestCase):
    def test_report_stage_two_counterexample(self):
        self.assertEqual(
            pollard_pm1.factorize_pm1(
                607 * 1019, b1=10, b2=200, max_attempts=1
            ),
            607,
        )

    def test_first_stage_two_prime_and_tail_batch(self):
        # 23-1=2*11: q=11 is the first eligible prime and the final term.
        self.assertEqual(
            pollard_pm1.factorize_pm1(
                23 * 47, b1=5, b2=11, max_attempts=1, batch_size=128
            ),
            23,
        )

    def test_every_prime_gap_term(self):
        n = 607 * 1019
        primes = reference_primes(201)
        primes = [p for p in primes if p > 10]
        residue = pow(2, 2520, n)
        expected = [(pow(residue, p, n) - 1) % n for p in primes]
        self.assertEqual(
            list(pollard_pm1._stage_two_terms(residue, n, primes)), expected
        )

    def test_stage_one_and_small_odd_input(self):
        self.assertTrue(utils.valid_divisor(pollard_pm1.factorize_pm1(9), 9))
        result = pollard_pm1.factorize_pm1(257 * 1019, b1=256, b2=256)
        self.assertTrue(utils.valid_divisor(result, 257 * 1019))

    def test_stage_one_saturation_replay(self):
        # After 2^8 mod 15 both factors saturate; the earlier 2^2 splits 3.
        residue, divisor = pollard_pm1._stage_one(15, 2, [2], 8)
        self.assertEqual(divisor, 3)
        self.assertEqual(residue, pow(2, 8, 15))

    def test_saturation_skips_continuation(self):
        with patch.object(pollard_pm1, "_stage_one", return_value=(1, 35)):
            with patch.object(pollard_pm1, "_stage_two_terms") as continuation:
                self.assertIsNone(
                    pollard_pm1.factorize_pm1(35, b1=5, b2=11, max_attempts=1)
                )
                continuation.assert_not_called()

    def test_stage_two_mixed_factor_replay(self):
        with patch.object(pollard_pm1, "_stage_one", return_value=(2, 1)):
            with patch.object(
                pollard_pm1, "_stage_two_terms", return_value=iter([5, 7])
            ):
                self.assertEqual(
                    pollard_pm1.factorize_pm1(
                        35, b1=5, b2=11, batch_size=2, max_attempts=1
                    ),
                    5,
                )

    def test_bounds_and_explicit_failure(self):
        for n in (1, 2, 3, 101):
            self.assertIsNone(pollard_pm1.factorize_pm1(n))
        self.assertIsNone(pollard_pm1.factorize_pm1(35, max_attempts=0))
        self.assertIsNone(
            pollard_pm1.factorize_pm1(1019 * 1237, b1=2, b2=2, max_attempts=1)
        )
        with self.assertRaises(ValueError):
            pollard_pm1.factorize_pm1(35, b1=10, b2=5)


class EcmTests(unittest.TestCase):
    def test_suyama_fixture(self):
        curve = ecm.setup_curve(101, 6)
        self.assertEqual(curve.a24, 93)
        self.assertEqual(
            curve.point[0] * pow(curve.point[1], -1, 101) % 101, 78
        )

    def test_setup_nonunit_and_singularity_branches(self):
        seen = set()
        for n in (9, 15, 21, 35, 77, 143, 1009 * 1013):
            for sigma in range(6, 60):
                curve = ecm.setup_curve(n, sigma)
                if curve.factor is not None:
                    self.assertTrue(utils.valid_divisor(curve.factor, n))
                    seen.add("factor")
                elif curve.retry:
                    seen.add("retry")
                else:
                    x, z = curve.point
                    self.assertNotEqual((x, z), (0, 0))
                    self.assertEqual(gcd(z, n), 1)
                    curve_a = (4 * curve.a24 - 2) % n
                    self.assertEqual(gcd(curve_a * curve_a - 4, n), 1)
                    seen.add("curve")
        self.assertEqual(seen, {"factor", "retry", "curve"})

    def test_ladder_and_prac_fallback_against_affine_oracle(self):
        modulus, curve_a = 1009, 6
        point = (3, 293)
        a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
        expected = None
        for scalar in range(200):
            if scalar:
                expected = affine_add(expected, point, modulus, curve_a)
            actual = ecm.scalar_multiply(scalar, point[0], 1, modulus, a24)
            self.assertNotEqual(actual, (0, 0), scalar)
            if expected is None:
                self.assertEqual(actual[1], 0)
            else:
                self.assertEqual(
                    (actual[0] - expected[0] * actual[1]) % modulus, 0, scalar
                )
            self.assertEqual(
                ecm.multiply_prac(scalar, point[0], 1, modulus, a24), actual
            )

    def test_invalid_scalar_and_projective_pair(self):
        with self.assertRaises(ValueError):
            ecm.scalar_multiply(-1, 3, 1, 1009, 2)
        with self.assertRaises(ValueError):
            ecm.scalar_multiply(5, 0, 0, 1009, 2)
        self.assertEqual(ecm.scalar_multiply(6, 1, 0, 1009, 2), (1, 0))

    def test_exact_stage_one_scalar(self):
        scalar = ecm.stage_one_scalar(243)
        self.assertEqual(scalar % 243, 0)
        self.assertEqual(scalar % 256, 128)
        for value in range(1, 244):
            self.assertEqual(scalar % value, 0)

    def test_zero_and_exact_curve_attempts(self):
        for curves in (0, 1, 2):
            stats = ecm.EcmStats()
            with patch.object(
                ecm, "setup_curve", return_value=ecm.CurveSetup(retry=True)
            ) as setup:
                self.assertIsNone(
                    ecm.factorize_ecm(
                        10403, b1=10, b2=100, max_curves=curves, stats=stats
                    )
                )
                self.assertEqual(setup.call_count, curves)
                self.assertEqual(stats.curves, curves)

    def test_saturated_stage_one_skips_stage_two(self):
        with patch.object(
            ecm, "setup_curve", return_value=ecm.CurveSetup((3, 1), 2)
        ):
            with patch.object(ecm, "scalar_multiply", return_value=(1, 10403)):
                with patch.object(ecm, "stage_two") as continuation:
                    stats = ecm.EcmStats()
                    self.assertIsNone(
                        ecm.factorize_ecm(
                            10403, b1=10, b2=100, max_curves=2, stats=stats
                        )
                    )
                    continuation.assert_not_called()
                    self.assertEqual(stats.stage_one_saturations, 2)

    def test_final_allowed_stage_two_success(self):
        with patch.object(
            ecm, "setup_curve", return_value=ecm.CurveSetup((3, 1), 2)
        ):
            with patch.object(ecm, "scalar_multiply", return_value=(3, 1)):
                with patch.object(
                    ecm, "stage_two", side_effect=[(None, True), (101, False)]
                ):
                    stats = ecm.EcmStats()
                    self.assertEqual(
                        ecm.factorize_ecm(
                            10403, b1=10, b2=100, max_curves=2, stats=stats
                        ),
                        101,
                    )
                    self.assertEqual(stats.stage_two_calls, 2)
                    self.assertEqual(stats.stage_two_saturations, 1)

    def test_stage_two_batch_recovery_and_tail(self):
        for terms, size in (([5, 7], 2), ([1, 5], 128)):
            with patch.object(
                ecm, "_stage_two_terms", return_value=iter(terms)
            ):
                divisor, _ = ecm.stage_two((3, 1), 35, 2, 10, [11, 13], size)
                self.assertEqual(divisor, 5)

    def test_tiny_b1_continuation(self):
        primes = [3, 5, 7]
        terms = list(ecm._stage_two_terms((3, 1), 1009, 2, 2, primes))
        self.assertEqual(
            terms,
            [ecm.scalar_multiply(prime, 3, 1, 1009, 2)[1] for prime in primes],
        )

    def test_seeded_real_curves(self):
        for n in (1009 * 1013, 10007 * 10009, 104729 * 1009):
            for seed in range(5):
                divisor = ecm.factorize_ecm(
                    n, b1=50, b2=1000, seed=seed, max_curves=32
                )
                self.assertTrue(utils.valid_divisor(divisor, n), (n, seed))

    def test_real_stage_two_only_success(self):
        # Fixed curve survives stage 1; q=17 reveals 1013, q=29 reveals 1009.
        n = 1009 * 1013
        setup = ecm.setup_curve(n, 6)
        point = ecm.scalar_multiply(
            ecm.stage_one_scalar(10), *setup.point, n, setup.a24
        )
        self.assertEqual(gcd(point[1], n), 1)
        primes = reference_primes(501)
        primes = [prime for prime in primes if prime > 10]
        divisor, _ = ecm.stage_two(point, n, setup.a24, 10, primes)
        self.assertEqual(divisor, 1013)
        direct = [
            gcd(ecm.scalar_multiply(prime, *point, n, setup.a24)[1], n)
            for prime in primes
        ]
        self.assertIn(1013, direct)

    def test_stage_two_terms_against_independent_field_oracle(self):
        modulus, curve_a, point = 1009, 6, (3, 293)
        a24 = 2
        multiples = [None]
        for _ in range(300):
            multiples.append(
                affine_add(multiples[-1], point, modulus, curve_a)
            )
        for b1 in (5, 10, 11, 20, 50):
            primes = [p for p in reference_primes(200) if p > b1]
            terms = list(
                ecm._stage_two_terms((3, 1), modulus, a24, b1, primes)
            )
            center = b1 if b1 % 2 else b1 - 1
            distance = min(isqrt(primes[-1]), (center - 1) // 2)
            step = 2 * distance
            for prime, term in zip(primes, terms):
                while prime > center + step:
                    center += step
                difference = prime - center
                expected = multiples[center][0] == multiples[difference][0]
                self.assertEqual(term == 0, expected, (b1, prime))


class FactorizationTests(unittest.TestCase):
    def test_reconstruct_all_small_positive_inputs(self):
        for n in range(1, 20001):
            result = factorize(n, seed=n)
            self.assertTrue(result.complete, n)
            self.assertTrue(result.proven, n)
            self.assertEqual(result.reconstruct(), n)

    def test_input_contract_and_formatting(self):
        with self.assertRaises(ValueError):
            factorize(0)
        for n in (2.5, "15", True):
            with self.assertRaises(TypeError):
                factorize(n)
        for n in (1, -1, -15):
            result = factorize(n)
            self.assertEqual(result.reconstruct(), n)
            self.assertTrue(print_factorization(n, result))
        self.assertEqual(print_factorization(1, factorize(1)), "1 = 1")

    def test_huge_trial_division_no_float(self):
        factors, remainder = factorize_bf(10**400)
        self.assertEqual(factors, [(2, 400), (5, 400)])
        self.assertEqual(remainder, 1)

    def test_forced_failure_preserves_all_factors(self):
        for cofactor in (25013 * 25031, 1000000000039 * 1000000000061):
            n = 2**3 * 3 * cofactor
            with patch("v2.factor.pollard_rho.factorize_rho", return_value=-1):
                with patch("v2.factor.ecm.factorize_ecm", return_value=None):
                    result = factorize(n)
                    self.assertFalse(result.complete)
                    self.assertEqual(result.remaining, (cofactor,))
                    self.assertEqual(result.reconstruct(), n)
                    self.assertIn("unresolved", print_factorization(n, result))

    def test_rho_failure_reaches_ecm(self):
        n = 25013 * 25031
        with patch("v2.factor.pollard_rho.factorize_rho", return_value=None):
            with patch(
                "v2.factor.ecm.factorize_ecm", return_value=25013
            ) as call:
                result = factorize(n)
                call.assert_called_once()
                self.assertTrue(result.complete)
                self.assertEqual(result.reconstruct(), n)

    def test_invalid_splitters_cannot_lose_cofactors(self):
        for value in (None, -1, 0, 1, True, 10403, 6, 101.0):
            with patch(
                "v2.factor.pollard_rho.factorize_rho", return_value=value
            ):
                with patch("v2.factor.ecm.factorize_ecm", return_value=value):
                    result = factorize(10403, level=2)
                    self.assertEqual(result.remaining, (10403,))
                    self.assertEqual(result.reconstruct(), 10403)

    def test_failed_child_preserves_partial_result(self):
        n = 1009 * 1013 * 1019
        with patch(
            "v2.factor.pollard_rho.factorize_rho", side_effect=[1009, None]
        ):
            with patch("v2.factor.ecm.factorize_ecm", return_value=None):
                result = factorize(n, level=2)
                self.assertEqual(result.remaining, (1013 * 1019,))
                self.assertEqual(result.factors[0].value, 1009)
                self.assertEqual(result.reconstruct(), n)

    def test_repeated_factor_multiplicities(self):
        n = 1009**3 * 1013**2
        result = factorize(n, seed=3)
        self.assertEqual(
            [(f.value, f.exponent) for f in result.factors],
            [(1009, 3), (1013, 2)],
        )
        self.assertEqual(result.reconstruct(), n)

    def test_probable_is_not_proven(self):
        n = 2**127 - 1
        result = factorize(n, seed=1)
        self.assertTrue(result.complete)
        self.assertFalse(result.proven)
        self.assertEqual(result.factors[0].certainty, utils.Primality.PROBABLE)
        self.assertIn("probable prime", print_factorization(n, result))

    def test_zero_work_preserves_cofactor(self):
        result = factorize(1009 * 1013, level=2, rho_attempts=0, ecm_curves=0)
        self.assertFalse(result.complete)
        self.assertEqual(result.remaining, (1009 * 1013,))

    def test_quiet_api(self):
        output = io.StringIO()
        with redirect_stdout(output):
            factorize(25013 * 25031, seed=1)
            pollard_pm1.factorize_pm1(607 * 1019, b1=10, b2=200)
            ecm.factorize_ecm(1009 * 1013, seed=1, b1=50, b2=1000)
        self.assertEqual(output.getvalue(), "")

    def test_preserved_baseline_hashes(self):
        root = Path(__file__).resolve().parents[2]
        manifest = json.loads(
            (root / "v2/audit/inputs/provenance.json").read_text()
        )
        for name, expected in manifest["original_source_sha256"].items():
            self.assertEqual(
                hashlib.sha256((root / "v1" / name).read_bytes()).hexdigest(),
                expected,
                name,
            )


if __name__ == "__main__":
    unittest.main()
