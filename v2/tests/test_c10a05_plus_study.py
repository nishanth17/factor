"""Larger-band arm identity and durable interrupted study leases."""

import json
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.qs.c10a05_plus import study


class LargerBandStudyTests(unittest.TestCase):
    def test_candidate_arms_make_only_the_declared_deltas(self):
        configs = {
            name: asdict(value)
            for name, value in study.configurations().items()
        }
        c9 = configs["c9_window"]
        c10 = configs["c10_window"]
        paced = configs["collection"]
        self.assertEqual(configs["current"]["base_bound"], 1000)
        self.assertEqual(c9["collector"]["power_plan_bytes"], 0)
        self.assertEqual(c10["collector"]["power_plan_bytes"], 2 * 2**20)
        self.assertEqual(paced["solve_policy"], "collection")
        self.assertEqual(c9["solve_policy"], "window")
        self.assertEqual(c10["solve_policy"], "window")
        c9["collector"]["power_plan_bytes"] = 2 * 2**20
        self.assertEqual(c9, c10)
        paced["solve_policy"] = "window"
        self.assertEqual(paced, c10)

    def test_interrupted_group_keeps_first_arm_and_failed_charge(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            fixtures = [
                dict(id=f"test-{i}", n=77, factors=[7, 11], digits=51)
                for i in range(3)
            ]
            control = dict(
                fixtures=fixtures,
                configurations={
                    name: asdict(config)
                    for name, config in study.configurations().items()
                },
            )

            def cold(_fixture, _arm, _config, _seconds, path):
                path.write_text("{}\n")

            calls = []

            def sample(fixture, seed, _config, arm, cap):
                calls.append((fixture["id"], seed, arm))
                if len(calls) == 2:
                    raise RuntimeError("interrupt after the first saved arm")
                seconds = 100 if arm == "current" else 90
                return dict(
                    id=fixture["id"],
                    seed=seed,
                    arm=arm,
                    complete=True,
                    seconds=seconds,
                    cpu_seconds=seconds,
                    capped_seconds=seconds,
                    cap_seconds=cap,
                    reason="factor_found",
                )

            with (
                patch.object(study, "load", return_value=control),
                patch.object(study, "CONTROL", output / "control.json"),
                patch.object(study.prior, "cold", side_effect=cold),
                patch.object(study, "warm"),
                patch.object(study.c9, "run_one", side_effect=sample),
            ):
                (output / "control.json").write_text("{}\n")
                with self.assertRaisesRegex(RuntimeError, "interrupt"):
                    study.measure(output, "train", 7200)
                first = next(output.glob("train-51-00-*.json"))
                saved = first.read_bytes()
                expense = study.Expenses(output)
                self.assertEqual(expense.data["calls"][-1]["status"], "failed")
                self.assertGreaterEqual(expense.data["calls"][-1]["wall"], 602)

                study.measure(output, "train", 7200)
                self.assertEqual(first.read_bytes(), saved)
                self.assertEqual(len(study.read_rows(output, "train")), 36)
                self.assertEqual(len(calls), 37)
                self.assertTrue((output / "train-finished.json").exists())
                self.assertGreater(
                    study.Expenses(output).used(), study.PREPAID + 602
                )
                self.assertEqual(
                    json.loads((output / "train-finished.json").read_text())[
                        "samples"
                    ],
                    9,
                )


if __name__ == "__main__":
    unittest.main()
