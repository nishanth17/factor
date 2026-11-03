"""Frozen nominee and unrun pilot resource invariants."""

import subprocess
import unittest
from tempfile import TemporaryDirectory

from v2.benchmarks.qs.c10a01 import nominee_v4, pilots_v4, plan_v4


class PilotV4FreezeTests(unittest.TestCase):
    def test_worker_command_imports_entrypoint(self):
        command = pilots_v4.worker_command("A", 0)
        self.assertEqual(command[2], pilots_v4.WORKER_MODULE)
        result = subprocess.run(
            command + ["--help"], capture_output=True, text=True, timeout=10
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--cell-index", result.stdout)

    def test_nominee_configs_and_receipt_bindings(self):
        frozen = nominee_v4.verify()
        self.assertEqual(len(frozen["calls"]), 8)
        self.assertEqual(frozen["pilot_status"], "unrun")
        self.assertEqual(frozen["control_commit"], plan_v4.START)
        self.assertEqual(frozen["policies"]["root_update_policy"], "batch")
        self.assertEqual(
            frozen["policies"]["multiplier_policy"], "parity2_wide"
        )
        for call in frozen["calls"]:
            config = call["configuration"]
            if call["arm"] == "control":
                self.assertNotIn("polynomial_policy", config)
                self.assertNotIn("multiplier_policy", config)
            else:
                self.assertEqual(config["polynomial_policy"], "parity2")
                self.assertEqual(config["multiplier_policy"], "parity2_wide")
                self.assertEqual(config["assignment_policy"], "windowed_pair")

    def test_exact_pilot_reservations(self):
        frozen = nominee_v4.verify()
        self.assertEqual(
            [frozen["stages"][stage]["total"]["wall"] for stage in "ABC"],
            [9200, 7600, 14800],
        )
        self.assertEqual(
            sum(frozen["stages"][stage]["total"]["wall"] for stage in "ABC"),
            31600,
        )
        for stage in "ABC":
            limits = frozen["stages"][stage]
            self.assertEqual(limits["total"]["wall"], limits["total"]["cpu"])
            self.assertEqual(
                limits["total"]["work"],
                sum(p["work"] for p in limits["phases"].values()),
            )

    def test_execution_ledger_accepts_first_frozen_reservation(self):
        plan = plan_v4.build_plan()
        call = plan["pilots"]["A"]["calls"][0]
        limits = pilots_v4.ledger_limits(plan, "A")
        self.assertEqual(limits["disk_bytes"], 16 * 2**20)
        with TemporaryDirectory() as output:
            ledger = pilots_v4.BoundedLedger(output, "test", limits)
            self.assertTrue(ledger.room("warmup", 35, 290_000_000_000))
            ledger.call(
                "warmup",
                "dry",
                35,
                290_000_000_000,
                {},
                lambda: {"work": 0},
            )
            self.assertTrue(
                ledger.room("measured", call["seconds"] + 5, call["work"])
            )


if __name__ == "__main__":
    unittest.main()
