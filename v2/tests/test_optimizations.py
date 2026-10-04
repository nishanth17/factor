"""Regression tests for measured M9 optimizations and proof boundaries."""

import random
import unittest
from bisect import bisect_left, bisect_right
from unittest.mock import patch

from v2 import ecm, pollard_rho, utils
from v2.factor import factorize, factorize_bf


class OptimizationTests(unittest.TestCase):
    """Check avoided work without weakening arithmetic or result contracts."""

    def test_prime_skips_trial_and_splitters(self):
        """Exact and probable primes must avoid the expensive trial scan."""
        for n in (2147483647, 2**127 - 1):
            with patch("v2.factor.factorize_bf") as trial:
                with patch.object(pollard_rho, "factorize_rho") as rho:
                    with patch.object(ecm, "factorize_ecm") as curves:
                        result = factorize(n, seed=4)
            self.assertTrue(result.complete)
            self.assertEqual(result.reconstruct(), n)
            self.assertEqual(result.proven, n < 2**64)
            trial.assert_not_called()
            rho.assert_not_called()
            curves.assert_not_called()

    def test_trial_root_is_reused_and_squares_are_exact(self):
        """The loop computes one bound unless successful division shrinks n."""
        with patch.object(utils, "isqrt", wraps=utils.isqrt) as root:
            factors, remainder = factorize_bf(25013 * 25031)
        self.assertEqual(factors, [])
        self.assertEqual(remainder, 25013 * 25031)
        self.assertEqual(root.call_count, 1)
        for prime in (41, 257, 1009):
            self.assertEqual(factorize_bf(prime**2), ([(prime, 2)], 1))

    def test_repeated_child_reuses_local_classification(self):
        """Classify repeated children once and reuse composite evidence."""
        with patch.object(
            utils, "classify_prime", wraps=utils.classify_prime
        ) as classify:
            with patch.object(pollard_rho, "factorize_rho", return_value=1009):
                result = factorize(1009**2, level=2, seed=3)
        values = [call.args[0] for call in classify.call_args_list]
        self.assertEqual(values.count(1009**2), 1)
        self.assertEqual(values.count(1009), 1)
        self.assertEqual(result.factors[0].exponent, 2)
        self.assertTrue(result.proven)

    def test_probable_child_cache_preserves_certainty(self):
        """Repeated probable children never acquire a proof from the cache."""
        prime = 2**127 - 1
        with patch.object(
            utils, "classify_prime", wraps=utils.classify_prime
        ) as classify:
            with patch.object(ecm, "factorize_ecm", return_value=prime):
                result = factorize(prime**2, level=1, seed=8)
        values = [call.args[0] for call in classify.call_args_list]
        self.assertEqual(values.count(prime), 1)
        self.assertEqual(result.factors[0].exponent, 2)
        self.assertEqual(result.factors[0].certainty, utils.Primality.PROBABLE)
        self.assertTrue(result.complete)
        self.assertFalse(result.proven)
        self.assertEqual(result.reconstruct(), prime**2)

    def test_classification_cache_is_per_factorization(self):
        """No evidence or random-round result leaks into a different call."""
        with patch.object(
            utils, "classify_prime", wraps=utils.classify_prime
        ) as classify:
            factorize(2147483647, seed=1)
            factorize(2147483647, seed=2)
        self.assertEqual(classify.call_count, 2)

    def test_exact_paths_do_not_allocate_rng(self):
        """Proven primes and fully stripped composites need no random state."""
        with patch.object(utils, "resolve_rng") as generator:
            for n in (1, -1, 2147483647, 2**16 * 3**6 * 101):
                self.assertTrue(factorize(n, seed=7).proven)
        generator.assert_not_called()
        for n in (1, 2147483647, 35):
            with self.assertRaises(TypeError):
                factorize(n, seed=[])

    def test_bounded_mr_witness_counts_and_thresholds(self):
        """Strict witness-domain boundaries are checked independently."""
        for n, rounds in ((1000003, 2), (2147483647, 3), (2**61 - 1, 7)):
            with patch.object(
                utils,
                "_strong_probable_prime",
                wraps=utils._strong_probable_prime,
            ) as witness:
                self.assertEqual(
                    utils.classify_prime(n), utils.Primality.PROVEN
                )
            self.assertEqual(witness.call_count, rounds)
        for boundary in (9080191, 4759123141):
            for n in range(boundary - 3, boundary + 4):
                self.assertEqual(utils.is_prime(n), utils.is_prime_bf(n), n)
            self.assertFalse(utils.is_prime(boundary))
        generator = random.Random(991)
        for _ in range(500):
            n = generator.randrange(1681, 4759123141)
            self.assertEqual(utils.is_prime(n), utils.is_prime_bf(n), n)

    def test_dispatch_does_not_repeat_splitter_prime_checks(self):
        """Public splitters still check primes; the dispatcher reuses proof."""
        n = 25013 * 25031
        with patch.object(
            utils, "classify_prime", wraps=utils.classify_prime
        ) as classify:
            result = factorize(n, seed=7)
        self.assertTrue(result.proven)
        self.assertEqual(
            [call.args[0] for call in classify.call_args_list].count(n), 1
        )
        for splitter in (pollard_rho.factorize_rho, ecm.factorize_ecm):
            with patch.object(utils, "is_prime", wraps=utils.is_prime) as test:
                self.assertIsNone(splitter(1009, seed=1))
            test.assert_called_once()

    def test_inverse_backends_match_and_reject_nonunits(self):
        """Exercise both implementations on signed values and large moduli."""
        generator = random.Random(106)
        for backend in (False, True):
            with patch.object(utils, "_USE_PYPY_INVERSE", backend):
                for modulus in (7, 60, 101, 2**192 - 237, 2**512 - 159):
                    for _ in range(30):
                        value = generator.randrange(-2 * modulus, 2 * modulus)
                        try:
                            expected = pow(value, -1, modulus)
                        except ValueError:
                            with self.assertRaises(ValueError):
                                utils.modular_inverse(value, modulus)
                        else:
                            self.assertEqual(
                                utils.modular_inverse(value, modulus), expected
                            )

    def test_exact_square_split_and_neighbor(self):
        """Squares split without rho; a neighboring nonsquare remains whole."""
        prime = 1000003
        with patch.object(pollard_rho, "factorize_rho") as rho:
            with patch.object(ecm, "factorize_ecm") as curves:
                result = factorize(prime**2, level=2, seed=9)
        self.assertTrue(result.proven)
        self.assertEqual(result.factors[0].exponent, 2)
        rho.assert_not_called()
        curves.assert_not_called()
        with patch.object(pollard_rho, "factorize_rho", return_value=None):
            with patch.object(ecm, "factorize_ecm", return_value=None):
                result = factorize(prime**2 + 1, level=2, seed=9)
        self.assertEqual(result.remaining, (prime**2 + 1,))
        self.assertEqual(result.reconstruct(), prime**2 + 1)

    def test_search_backends_match_duplicate_boundaries(self):
        """Both backends match stdlib searches on random duplicate arrays."""
        generator = random.Random(108)
        for backend in (False, True):
            implementation = (
                utils._binary_search_jit
                if backend
                else utils._binary_search_bisect
            )
            with patch.object(utils, "binary_search", implementation):
                for _ in range(200):
                    array = sorted(generator.choices(range(-10, 11), k=50))
                    for value in range(-12, 13):
                        self.assertEqual(
                            utils.binary_search(value, array),
                            bisect_right(array, value),
                        )
                        self.assertEqual(
                            utils.binary_search(value, array, True),
                            bisect_left(array, value),
                        )
                for value in (0, float("inf"), float("nan")):
                    self.assertEqual(utils.binary_search(value, []), 0)

    def test_wheel_sieve_all_small_endpoints(self):
        """All six endpoint residues and nearby squares match trial proof."""
        from v2 import prime_sieve

        primes = [n for n in range(2, 1001) if utils.is_prime_bf(n)]
        for hi in range(-3, 1002):
            self.assertEqual(
                prime_sieve.small_sieve(hi), [p for p in primes if p < hi]
            )
