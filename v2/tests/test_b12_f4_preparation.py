"""Queued runner/corpus/analyser QA; not the narrow scheduler-only harness."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from v2.benchmarks.ecm.b12_f4 import runner
from v2.benchmarks.ecm.b12_f4.inputs import random_integers, retained
from v2.benchmarks.ecm.b12_f4.runner import (
    clustered_interval,
    forced_arm,
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

    def test_option_catalog_is_prospective_and_arm_counts_are_explicit(self):
        path = (
            REPOSITORY_ROOT
            / "v2/benchmarks/inputs/controls/b12_f4_options_preparation.json"
        )
        data = json.loads(path.read_text())
        self.assertEqual(data["mechanism_comparisons"]["wall"], 720)
        self.assertEqual(data["mechanism_comparisons"]["cpu"], 2880)
        self.assertFalse(data["wider_counts"]["funded"])
        with self.assertRaises(ValueError):
            load_freeze(path)
        self.assertTrue(forced_arm("P4", 4))
        self.assertTrue(forced_arm("P8", 8))
        self.assertFalse(forced_arm("N", 4))
        for arm, workers in (("S0", 4), ("S1", 4), ("P4", 2), ("N", 0)):
            with self.assertRaises(ValueError):
                forced_arm(arm, workers)

    def test_cold_endpoint_includes_external_result_validation(self):
        clock = [0.0]
        output = dict(
            n=35,
            reconstruct=35,
            factors=[[5, 1, "proven_prime"], [7, 1, "proven_prime"]],
            remaining=[],
            factorization_complete=True,
            aggregate_cpu=0,
        )

        def launch(*args, **kwargs):
            kwargs["stdout"].write(json.dumps(output).encode())
            kwargs["stdout"].flush()
            return Mock(returncode=0)

        original = runner.validate_output

        def validate(*args):
            clock[0] += 2
            original(*args)

        cell = dict(
            id="endpoint",
            arm="S1",
            mode="cold",
            workers=1,
            wall=10,
            cpu=10,
            startup_warmup_wall=1,
            pool="fresh",
            input="one",
        )
        freeze = dict(
            identity="fake",
            wall_limit=20,
            cpu_limit=20,
            settings={},
            seed=7,
            inputs={"one": dict(n=35)},
            work_limit=100,
            pypy="unused",
            candidate_root="unused",
        )
        with (
            tempfile.TemporaryDirectory() as directory,
            (
                patch.object(runner.fcntl, "flock"),
                patch.object(runner.subprocess, "Popen", side_effect=launch),
                patch.object(
                    runner.time, "monotonic", side_effect=lambda: clock[0]
                ),
                patch.object(runner, "validate_output", side_effect=validate),
            ),
        ):
            record = runner.run_cell(freeze, cell, Path(directory))
        self.assertEqual(record["process_wall"], 0)
        self.assertEqual(record["validation_seconds"], 2)
        self.assertEqual(record["wall"], 2)

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
