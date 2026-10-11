"""Matched upstream-preset screening and bounded uncertainty extensions."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.ecm.c3 import round2_compare as matched
from v2.benchmarks.ecm.c3 import source_preset_study as study


def groups(candidate_cost=0.8, candidate_complete=True):
    result = {}
    for band, kind in ((30, "balanced"), (40, "uneven_16")):
        samples = {}
        for sample in range(9):
            rows = {}
            for arm in study.ARMS:
                candidate = arm in study.PRESETS
                rows[arm] = dict(
                    cpu=candidate_cost if candidate else 1.0,
                    wall=candidate_cost if candidate else 1.0,
                    complete=candidate_complete if candidate else True,
                    cap_seconds=5 if band == 30 else 30,
                )
            samples[sample] = rows
        result[(band, kind)] = {str(band): samples}
    return result


class SourcePresetStudyTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        control = Path(temporary.name) / "control.json"
        control.write_text("{}")
        replacement = patch.object(study, "CONTROL", control)
        replacement.start()
        self.addCleanup(replacement.stop)

    def test_whole_group_reuse_restores_prior_comparison_identity(self):
        identity = matched.CONTROL, matched.ARMS
        with self.assertRaisesRegex(RuntimeError, "stop"):
            with study.layout():
                self.assertEqual(matched.CONTROL, study.CONTROL)
                self.assertEqual(matched.ARMS, study.ARMS)
                raise RuntimeError("stop")
        self.assertEqual((matched.CONTROL, matched.ARMS), identity)

    def test_assignment_covers_declared_inputs_seeds_and_rotated_arms(self):
        fixtures = [{"id": name} for name in study.CASE_IDS]
        assignments = list(study.assignment(fixtures, 27))
        self.assertEqual(len(assignments), 162)
        self.assertEqual(len(set(study.SEEDS)), 27)
        self.assertEqual(
            len({(f["id"], sample) for f, sample, _, _ in assignments}), 162
        )
        for fixture, sample, seed, order in assignments:
            self.assertEqual(seed, study.SEEDS[sample])
            self.assertEqual(set(order), set(study.ARMS))
        self.assertNotEqual(assignments[0][3], assignments[1][3])

    def test_clear_cost_reduction_selects_screen_candidate_without_extension(
        self,
    ):
        result = study.assess(groups(), 9, repetitions=25)
        self.assertIn(result["selected"], study.PRESETS)
        self.assertIsNone(result["next_samples"])
        self.assertTrue(
            all(r["stable"] for r in result["candidates"].values())
        )

    def test_unresolved_early_return_pays_full_service_cap(self):
        result = study.assess(
            groups(candidate_cost=0.001, candidate_complete=False),
            9,
            repetitions=25,
        )
        self.assertIsNone(result["selected"])
        self.assertIsNone(result["next_samples"])
        self.assertTrue(
            all(
                not r["completion_gate"] for r in result["candidates"].values()
            )
        )

    def test_promising_uncertainty_extends_only_to_predeclared_limit(self):
        uncertain = dict(relative_improvement=0.1, interval95=[-0.1, 0.3])
        with patch.object(
            study.reporter, "improvement_interval", return_value=uncertain
        ):
            first = study.assess(groups(), 9)
            second = study.assess(groups(), 18)
            last = study.assess(groups(), 27)
        self.assertEqual(first["next_samples"], 18)
        self.assertEqual(second["next_samples"], 27)
        self.assertIsNone(last["next_samples"])
        self.assertIsNone(last["selected"])

    def test_wall_regression_blocks_selection_and_extra_samples(self):
        observations = groups()
        for subjects in observations.values():
            for samples in subjects.values():
                for rows in samples.values():
                    for arm in study.PRESETS:
                        rows[arm]["wall"] = 1.06
        result = study.assess(observations, 9, repetitions=25)
        self.assertIsNone(result["selected"])
        self.assertIsNone(result["next_samples"])

    def test_changed_or_unrecorded_assessment_cannot_authorize_extension(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            self.assertEqual(study.target_samples(output), 9)
            path = output / "assessment-09.json"
            path.write_text("{}")
            with self.assertRaisesRegex(ValueError, "unrecorded"):
                study.target_samples(output)

    def test_recorded_extension_and_terminal_stop_are_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            control = output / "control.json"
            control.write_text("{}")
            path = output / "assessment-09.json"
            identity = study.first.digest(control)
            path.write_text(
                json.dumps(
                    dict(protocol_sha256=identity, samples=9, next_samples=18)
                )
            )
            ledger = dict(
                protocol_sha256=identity,
                assessments={"9": study.first.digest(path)},
            )
            (output / "ledger.json").write_text(json.dumps(ledger))
            with patch.object(study, "CONTROL", control):
                self.assertEqual(study.target_samples(output), 18)
                path.write_text("{}")
                with self.assertRaisesRegex(ValueError, "assessment changed"):
                    study.target_samples(output)

    def test_undeclared_extension_and_failed_phase_cannot_continue(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            control = output / "control.json"
            control.write_text("{}")
            identity = study.first.digest(control)
            path = output / "assessment-09.json"
            path.write_text(
                json.dumps(
                    dict(protocol_sha256=identity, samples=9, next_samples=100)
                )
            )
            ledger = dict(
                protocol_sha256=identity,
                assessments={"9": study.first.digest(path)},
            )
            (output / "ledger.json").write_text(json.dumps(ledger))
            with patch.object(study, "CONTROL", control):
                with self.assertRaisesRegex(ValueError, "undeclared"):
                    study.target_samples(output)
                ledger["failed"] = "interrupted"
                (output / "ledger.json").write_text(json.dumps(ledger))
                with self.assertRaisesRegex(ValueError, "failed/interrupted"):
                    study.target_samples(output)
