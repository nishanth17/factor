"""Repeated leases and interrupted matched groups retain evidence and costs."""

import json
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.qs.c10a05 import study
from v2.qs import SIQSConfig


class StudyTests(unittest.TestCase):
    def prepare(self, output):
        config = SIQSConfig()
        registry = {
            a: asdict(c) for a, c in study.configurations(config).items()
        }
        fixtures = [
            dict(id=f"fixture-{i}", n=77, factors=[7, 11], digits=30)
            for i in range(3)
        ]
        for arm in registry:
            (output / f"train-cold-30-{arm}.json").write_text("{}")
        return dict(fixtures=fixtures, configurations={"30": registry})

    def test_policy_is_the_only_configuration_difference(self):
        registry = study.configurations(SIQSConfig())
        control, candidate = (asdict(registry[a]) for a in study.ARMS)
        self.assertEqual(candidate.pop("solve_policy"), "collection")
        self.assertEqual(control.pop("solve_policy"), "window")
        self.assertEqual(candidate, control)

    def test_repeated_yield_leases_do_not_overwrite_or_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            control = self.prepare(output)
            with (
                patch.object(study, "load_control", return_value=control),
                patch.object(study, "digest", return_value="fixed"),
                patch.object(study, "warm"),
            ):
                study.measure(output, "train", 30, 0)
                first = {p: p.read_bytes() for p in output.iterdir()}
                study.measure(output, "train", 30, 0)
            self.assertEqual(
                len(list(output.glob("train-yield-30-*.json"))), 2
            )
            for path, content in first.items():
                self.assertEqual(path.read_bytes(), content)

    def test_interrupted_group_keeps_first_arm_and_failed_charge(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            control = self.prepare(output)
            calls = []

            def run(fixture, seed, config, arm, cap):
                calls.append(arm)
                if len(calls) == 2:
                    raise KeyboardInterrupt("injected interruption")
                return dict(
                    id=fixture["id"],
                    arm=arm,
                    seed=seed,
                    complete=True,
                    seconds=1,
                    capped_seconds=1,
                    reason="factor_found",
                )

            with (
                patch.object(study, "load_control", return_value=control),
                patch.object(study, "digest", return_value="fixed"),
                patch.object(study, "warm"),
                patch.object(study.c9, "run_one", side_effect=run),
            ):
                with self.assertRaises(KeyboardInterrupt):
                    study.measure(output, "train", 30, 1000)
            paths = list(output.glob("train-30-00-*.json"))
            self.assertEqual(len(paths), 1)
            saved = paths[0].read_bytes()
            expenses = study.expenses(output)
            self.assertGreaterEqual(study.phase_used(expenses, "train"), 12)
            with (
                patch.object(study, "load_control", return_value=control),
                patch.object(study, "digest", return_value="fixed"),
                patch.object(study, "warm"),
                patch.object(study.c9, "run_one", side_effect=run),
            ):
                study.measure(output, "train", 30, 1000)
            self.assertEqual(paths[0].read_bytes(), saved)
            self.assertEqual(
                len(list(output.glob("train-30-[0-9][0-9]-*.json"))), 18
            )
            self.assertEqual(len(calls), 19)
            data = json.loads((output / "expenses.json").read_text())
            self.assertEqual(
                sum(r["status"] == "failed" for r in data["calls"]), 1
            )


if __name__ == "__main__":
    unittest.main()
