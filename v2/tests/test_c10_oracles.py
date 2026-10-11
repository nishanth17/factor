"""Independent normalized-polynomial and declared-loss counterexamples."""

import unittest

from v2.benchmarks.qs.c10.oracles import (
    DiscriminantPolynomial,
    divide_over_base,
    multiplier_divisor,
    product_census,
    threshold_losses,
)


class NormalizedPolynomialTests(unittest.TestCase):
    def test_explicit_normalization_and_nonsingular_two_lifts(self):
        polynomial, delta = DiscriminantPolynomial.from_standard(1009, 15, 2)
        self.assertEqual((polynomial.b, polynomial.c, delta), (17, -12, 1))
        for level in range(1, 8):
            modulus = 2**level
            roots = [
                x for x in range(modulus) if polynomial.value(x) % modulus == 0
            ]
            self.assertEqual(len(roots), 2)
            for root in roots:
                lifts = [
                    root + step * modulus
                    for step in (0, 1)
                    if polynomial.value(root + step * modulus) % (2 * modulus)
                    == 0
                ]
                self.assertEqual(len(lifts), 1)

    def test_identities_parity_and_exceptional_roots(self):
        cases = 0
        for discriminant in range(5, 130, 4):
            for a in (1, 3, 5, 7, 15, 21):
                for b in range(a):
                    if (b * b - discriminant) % a:
                        continue
                    normalized, delta = DiscriminantPolynomial.from_standard(
                        discriminant, a, b
                    )
                    standard_c = (b * b - discriminant) // a

                    for x in range(-9, 10):
                        y = 2 * x + delta
                        value = normalized.value(x)
                        self.assertEqual(
                            a * y * y + 2 * b * y + standard_c, 4 * value
                        )
                        self.assertEqual(
                            normalized.u_value(x) ** 2 - discriminant,
                            4 * a * value,
                        )
                        self.assertEqual(value % 2, (discriminant % 8 == 5))
                        if value:
                            divide_over_base(4 * a * value, (2, 3, 5, 7))

                    for prime in (2, 3, 5, 7, 11, 13):
                        expected = tuple(
                            x
                            for x in range(prime)
                            if normalized.value(x) % prime == 0
                        )
                        self.assertEqual(normalized.roots(prime), expected)
                        # Enumerate complete powers separately, including
                        # singular lifts where derivative inversion fails.
                        modulus = prime**2
                        lifted = {
                            x
                            for x in range(modulus)
                            if normalized.value(x) % modulus == 0
                        }
                        self.assertTrue(
                            all(x % prime in expected for x in lifted)
                        )
                    cases += 1
        self.assertGreater(cases, 100)

    def test_invalid_discriminants_and_zero_recovery(self):
        for args in ((7, 3, 1), (5, 2, 1), (5, 3, 2), (5, 3, 1)):
            with self.assertRaises(ValueError):
                DiscriminantPolynomial(*args)
        with self.assertRaises(ValueError):
            divide_over_base(0, (2, 3))
        with self.assertRaises(ValueError):
            DiscriminantPolynomial(5, 1, 1).roots(9)

    def test_actual_product_quality_does_not_prove_diversity(self):
        families = ((3, 5, 7), (3, 5, 11), (3, 7, 11))
        census = product_census(families, 165)
        self.assertEqual(census["unique_products"], 3)
        self.assertEqual(census["error_sum"], 126)
        self.assertEqual(census["pair_overlap"], {2: 3})

    def test_threshold_loss_and_residual_certainty_are_separate(self):
        # Residual 49 fits the report bound but is composite; it cannot be
        # labeled a prime endpoint just because it survived this census.
        result = threshold_losses(
            (12, -15, 49, 77),
            (8, 7, 6, 5),
            (6, 6, 6, 6),
            (2, 3, 5),
            49,
            (0, 2),
        )
        self.assertEqual(result[0], dict(reports=3, eligible=3, lost=0))
        self.assertEqual(result[2], dict(reports=1, eligible=3, lost=2))
        self.assertEqual(multiplier_divisor(15, 6), 3)
        self.assertIsNone(multiplier_divisor(15, 15))
        self.assertIsNone(multiplier_divisor(15, 2))


if __name__ == "__main__":
    unittest.main()
