"""Prepared arms preserve calibrated controls and refuse stale source pins."""

import json
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.qs.a7 import a7_r5
from v2.qs.sss import SSSJob


class PreparedArmTests(unittest.TestCase):
    def test_calibrated_controls_are_exact_and_budgets_are_common(self):
        plan = a7_r5.load_plan()
        for band, name in (
            ("30d", "b1_selected.json"),
            ("40d", "b1_40d_selected.json"),
        ):
            selected = json.loads(a7_r5.PLAN.with_name(name).read_text())
            configs = a7_r5.serial_configurations(band, plan=plan)
            for mode, label in selected["selected"].items():
                self.assertEqual(
                    asdict(configs[mode]), selected["configurations"][label]
                )
            for config in configs.values():
                self.assertEqual(config.memory_bytes, a7_r5.b1.MEMORY)
                self.assertEqual(config.backend, "python-int")
        self.assertNotIn("sss", a7_r5.serial_configurations("40d", plan=plan))
        self.assertEqual(plan["seeds"], [7, 29])

    def test_both_sssf_loss_policies_are_separately_labelled(self):
        configs = a7_r5.serial_configurations("30d")
        self.assertEqual(configs["sss"].mode, "sss")
        self.assertGreater(configs["sssf_filtered"].filter_bound, 0)
        self.assertEqual(configs["sssf_six_unfiltered"].filter_bound, 0)
        self.assertEqual(configs["sssf_six_unfiltered"].selection_size, 6)

    def test_current_adapter_uses_b1_classification_and_allowances(self):
        fixture = {"n": 10**29 + 1, "digits": 30}
        with patch.object(a7_r5.b1, "run_one", return_value="checked") as run:
            self.assertEqual(
                a7_r5.run_serial(fixture, 7, "30d", "sss"), "checked"
            )
        self.assertEqual(run.call_args.kwargs["job_type"], SSSJob)
        self.assertEqual(run.call_args.args[-1], 5)
        with self.assertRaises(ValueError):
            a7_r5.run_serial(fixture, 999, "30d", "sss")

    def test_worker_configuration_is_the_accepted_bounded_schedule(self):
        plan = a7_r5.load_plan()
        for band in ("small", "medium"):
            config = a7_r5.worker_configuration(band, plan=plan)
            self.assertEqual(config.batch_width, 256)
            self.assertEqual(config.poll_interval, 64)
            self.assertEqual(config.memory_bytes, 512 * 2**20)
            self.assertLessEqual(config.half_width, 8192)

    def test_changed_runtime_pin_cannot_load_as_the_frozen_control(self):
        with tempfile.TemporaryDirectory(prefix="a7-pin-") as directory:
            root = Path(directory)
            (root / "runtime.py").write_text("changed")
            plan = root / "plan.json"
            plan.write_text(
                json.dumps(dict(schema=1, source_sha256={"runtime.py": "old"}))
            )
            with (
                patch.object(a7_r5, "ROOT", root),
                patch.object(a7_r5, "PLAN", plan),
            ):
                with self.assertRaisesRegex(ValueError, "pin changed"):
                    a7_r5.load_plan()


if __name__ == "__main__":
    unittest.main()
