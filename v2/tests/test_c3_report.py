"""Hand-constructed paired outcomes for the training assessment."""

import unittest

from v2.benchmarks.ecm.c3.round2_report import (
    assess,
    capped_cost,
    weighted_cost,
)


def groups(candidate=1, control=2, fixed=1.5, complete=True):
    result = {}
    for band in (30, 40):
        subjects = {}
        for subject in ("a", "b", "c"):
            samples = {}
            for seed in range(9):
                samples[seed] = {
                    arm: dict(
                        cpu=value,
                        wall=value,
                        cap_seconds=10,
                        complete=complete if arm == "fitted" else True,
                    )
                    for arm, value in (
                        ("control", control),
                        ("fixed32", fixed),
                        ("fitted", candidate),
                    )
                }
            subjects[subject] = samples
        result[(band, "balanced")] = subjects
    return result


class ReportTests(unittest.TestCase):
    def test_clear_paired_improvement_is_eligible(self):
        report = assess(groups(), repetitions=100)
        self.assertEqual(
            report["verdict"], "eligible for new fresh confirmation"
        )
        self.assertEqual(
            report["comparisons"]["control"]["interval95"], [0.5, 0.5]
        )

    def test_beating_only_historical_control_does_not_pay_for_model(self):
        report = assess(groups(candidate=1.6), repetitions=100)
        self.assertIn("retain fixed", report["verdict"])

    def test_partial_fast_returns_are_penalized_and_fail_completion_gate(self):
        report = assess(
            groups(candidate=0.01, complete=False), repetitions=100
        )
        self.assertFalse(report["completion_gate"])
        self.assertIn("retain fixed", report["verdict"])
        self.assertEqual(
            capped_cost(dict(cpu=0.01, cap_seconds=10, complete=False), "cpu"),
            10,
        )

    def test_bands_are_equal_despite_unequal_subject_counts(self):
        example = groups()
        example[(40, "balanced")].pop("c")
        for samples in example[(40, "balanced")].values():
            for row in samples.values():
                row["control"]["cpu"] = 4
        self.assertEqual(weighted_cost(example, "control", "cpu"), 3)
