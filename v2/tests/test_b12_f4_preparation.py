"""Queued runner/corpus/analyser QA; not the narrow scheduler-only harness."""

import unittest

from v2.benchmarks.ecm.b12_f4.inputs import random_integers, retained
from v2.benchmarks.ecm.b12_f4.runner import (
    clustered_interval,
    load_freeze,
    paired_survival_interval,
    quantile_with_censoring,
    validate_output,
)
from v2.benchmarks.support.paths import REPOSITORY_ROOT


class PreparationTests(unittest.TestCase):
    def test_random_panel_keeps_all_exact_bands_without_survival_filter(self):
        rows = random_integers(7919, 12)
        self.assertEqual(len(rows), 36)
        self.assertEqual({row["digits"] for row in rows}, {60, 80, 100})
        self.assertTrue(any(row["n"] % 2 == 0 for row in rows))
        self.assertTrue(any(row["n"] % 2 == 1 for row in rows))
        self.assertTrue(
            all(len(str(row["n"])) == row["digits"] for row in rows)
        )

    def test_censored_quantiles_do_not_report_successes_only(self):
        self.assertIsNone(
            quantile_with_censoring([(1, True), (10, False), (10, False)], 0.5)
        )
        self.assertEqual(
            quantile_with_censoring([(1, True), (2, True), (10, False)], 0.5),
            2,
        )
        self.assertIsNone(
            quantile_with_censoring([(1, True), (2, True), (10, False)], 0.9)
        )

    def test_tail_uncertainty_and_modulus_cluster_independence(self):
        rows = {
            i: dict(control=[(10, True)] * 9, nominee=[(5, True)] * 9)
            for i in range(12)
        }
        report = paired_survival_interval(rows, repetitions=40)
        self.assertEqual(report["intervals"]["median"], [5, 5])
        self.assertGreater(report["one_sided_loss_upper_95"], 0)
        rows[0]["nominee"] = [(10, False)] * 9
        report = paired_survival_interval(rows, repetitions=40)
        self.assertGreater(report["modulus_completion_losses"], 0)
        self.assertGreater(report["unidentified_replicates"]["p90"], 0)
        with self.assertRaises(ValueError):
            clustered_interval({0: [(10, 5)] * 100}, repetitions=40)

    def test_retained_exposure_selector_and_prospective_freeze_refusal(self):
        fixture, checksum = retained(
            REPOSITORY_ROOT, "b3_corpus.json", "small10_60d"
        )
        self.assertEqual(len(str(fixture["n"])), 60)
        self.assertEqual(len(checksum), 64)
        with self.assertRaises(ValueError):
            load_freeze(
                REPOSITORY_ROOT
                / "v2/benchmarks/inputs/controls/b12_f4_preparation.json"
            )


class IndependentOutputTests(unittest.TestCase):
    def test_raw_factors_unresolved_and_certainty_are_checked(self):
        good = dict(
            n=35,
            reconstruct=35,
            factors=[[5, 1, "proven_prime"]],
            remaining=[7],
            factorization_complete=False,
        )
        validate_output(good, 35, [[5, 1], [7, 1]])
        for changed in (
            dict(good, remaining=[]),
            dict(good, factors=[[6, 1, "proven_prime"]]),
            dict(good, factors=[[5, 1, "composite"]]),
        ):
            with self.assertRaises(ValueError):
                validate_output(changed, 35)
