"""Calibration evidence cannot turn invalid or censored results into wins."""

import unittest
from copy import deepcopy

from v2.benchmarks.b1_40d import configurations
from v2.benchmarks.b1_calibration import rank, training_configs, validate
from v2.benchmarks.b1_qs_width import configurations as qs_configurations


class CalibrationTests(unittest.TestCase):
    def setUp(self):
        self.fixture = {"n": 77, "factors": [7, 11]}
        self.row = dict(
            factors=[7, 11],
            remaining=[],
            complete=True,
            divisor=7,
            certainty=["proven_prime", "proven_prime"],
            work=1,
            stats={},
        )

    def test_unresolved_reconstruction_and_certainty(self):
        validate(self.row, self.fixture)
        unresolved = dict(
            self.row,
            factors=[],
            remaining=[77],
            complete=False,
            divisor=None,
            certainty=[],
        )
        validate(unresolved, self.fixture)

        for field, value in (
            ("remaining", [7]),
            ("certainty", []),
            ("divisor", 77),
            ("factors", [7, 7]),
        ):
            corrupt = deepcopy(self.row)
            corrupt[field] = value
            with self.assertRaises(AssertionError):
                validate(corrupt, self.fixture)

    def test_proof_corpus_does_not_upgrade_runtime_label(self):
        large = 2**127 - 1
        row = dict(
            self.row,
            factors=[large],
            divisor=None,
            certainty=["probable_prime"],
        )
        fixture = {"n": large, "factors": [large]}
        validate(row, fixture)

        row["certainty"] = ["proven_prime"]
        with self.assertRaises(AssertionError):
            validate(row, fixture)

    def test_fresh_a10_deterministic_labels_above_word_domain(self):
        # This prime lies in A10's wider domain, as do B1's 40-digit children.
        prime = 2**64 + 13
        fixture = {"n": prime, "factors": [prime]}
        row = dict(
            self.row, factors=[prime], divisor=None, certainty=["proven_prime"]
        )

        validate(row, fixture)
        row["certainty"] = ["probable_prime"]
        with self.assertRaises(AssertionError):
            validate(row, fixture)

    def test_fast_refusal_cannot_outrank_completion(self):
        complete = dict(complete=True, seconds=2, stats={})
        refusal = dict(complete=False, seconds=0.01, stats={"relations": 1})
        useful = dict(refusal, seconds=5, stats={"relations": 100})
        self.assertEqual(
            rank(
                {
                    "complete": [complete],
                    "refusal": [refusal],
                    "useful": [useful],
                }
            ),
            ["complete", "useful", "refusal"],
        )

    def test_sample_extension_does_not_inflate_completion(self):
        row = dict(complete=True, seconds=1, stats={})
        self.assertEqual(
            rank(
                {"nine": [row] * 9, "fifteen": [row] * 15},
                {"nine": 1, "fifteen": 2},
            ),
            ["nine", "fifteen"],
        )

    def test_bundles_preserve_exact_conservative_collector(self):
        configs = training_configs()
        self.assertEqual(
            {c.mode for c in configs.values()}, {"qs", "mpqs", "siqs"}
        )
        for config in configs.values():
            self.assertEqual(config.backend, "python-int")
            self.assertEqual(config.collector.threshold_extra, 0)
            self.assertGreater(config.polynomial_limit, 0)

    def test_upper_extension_has_separate_joint_controls(self):
        configs = configurations()
        self.assertEqual(len(configs), 8)
        self.assertEqual(configs["siqs_10k_c4"].factor_count, 4)
        self.assertEqual(configs["siqs_10k_gray1"].gray_limit, 1)

        for config in configs.values():
            self.assertLessEqual(config.memory_bytes, 256 * 2**20)
            self.assertEqual(config.collector.threshold_extra, 0)

    def test_wider_qs_retains_one_finite_fixed_polynomial(self):
        for config in qs_configurations().values():
            self.assertEqual(config.mode, "qs")
            self.assertEqual(config.polynomial_limit, 1)
            self.assertLess(2 * config.half_width, 1_000_000)


if __name__ == "__main__":
    unittest.main()
