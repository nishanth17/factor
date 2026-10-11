"""Check prepared confirmation construction using revealed fixtures only."""

import copy
import json
import unittest
from types import SimpleNamespace
from unittest.mock import Mock, patch

from v2.benchmarks.ecm.c3 import round2_confirmation_inputs as inputs
from v2.benchmarks.ecm.c3.round2_corpus import TRAINING


class ConfirmationInputTests(unittest.TestCase):
    def setUp(self):
        # These are already revealed training integers, never fresh evidence.
        corpus = json.loads(TRAINING.read_text())
        corpus["counts"] = dict(balanced_count=2, uneven_count=2)
        corpus["fixtures"] = [
            f for f in corpus["fixtures"] if f["id"].endswith(("_0", "_1"))
        ]
        self.corpus = corpus

    def test_existing_certificates_products_and_exact_strata_validate(self):
        inputs.validate(self.corpus, [])
        self.assertEqual(len(self.corpus["fixtures"]), 20)

    def test_wrong_stratum_and_product_are_rejected(self):
        for key, value in (("small_digits", 8), ("n", 123)):
            corpus = copy.deepcopy(self.corpus)
            corpus["fixtures"][0][key] = value
            with self.assertRaisesRegex(ValueError, "factorization"):
                inputs.validate(corpus, [])

    def test_revealed_modulus_or_reused_top_level_prime_is_rejected(self):
        fixture = self.corpus["fixtures"][0]
        for excluded in (
            [fixture],
            [dict(n=fixture["n"] + 1, factors=[fixture["factors"][0]])],
        ):
            with self.assertRaisesRegex(ValueError, "overlaps"):
                inputs.validate(self.corpus, excluded)

    def test_missing_stratum_and_duplicate_input_are_rejected(self):
        corpus = copy.deepcopy(self.corpus)
        corpus["fixtures"].pop()
        with self.assertRaisesRegex(ValueError, "incomplete"):
            inputs.validate(corpus, [])
        corpus["fixtures"].append(corpus["fixtures"][0])
        with self.assertRaisesRegex(ValueError, "overlaps"):
            inputs.validate(corpus, [])

    def test_corrupt_certificate_cannot_be_accepted_from_factor_labels(self):
        corpus = copy.deepcopy(self.corpus)
        prime = corpus["fixtures"][0]["factors"][0][0]
        corpus["certificates"][str(prime)]["kind"] = "probable"
        with self.assertRaisesRegex(ValueError, "unknown"):
            inputs.validate(corpus, [])

    def test_unfinished_proof_aborts_without_resampling(self):
        source = SimpleNamespace(
            certificates={}, prime=Mock(side_effect=RuntimeError("unfinished"))
        )
        with patch.object(inputs, "UniformPrimeSource", return_value=source):
            with self.assertRaisesRegex(RuntimeError, "unfinished"):
                inputs.generate(17, [], balanced_count=1, uneven_count=1)
        source.prime.assert_called_once_with(15)

    def test_joint_decimal_rejection_has_a_finite_attempt_limit(self):
        source = SimpleNamespace(certificates={}, prime=Mock(return_value=7))
        with patch.object(inputs, "UniformPrimeSource", return_value=source):
            with self.assertRaisesRegex(RuntimeError, "pair rejection"):
                inputs.generate(17, [], balanced_count=1, uneven_count=1)
        self.assertEqual(source.prime.call_count, 2000)

    def test_count_recipe_rejects_unbounded_or_boolean_requests(self):
        for count in (0, 9, True):
            with self.assertRaises(ValueError):
                tuple(inputs.strata(count, 2))
