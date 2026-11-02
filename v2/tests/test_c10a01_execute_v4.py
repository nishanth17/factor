"""Preflight the v4 worker source boundary and finite reservation math."""

import json
import subprocess
import sys
import time
import unittest

from v2.benchmarks.qs.c10a01 import execute_v4, execute_v4_stage, plan_v4


class ExecutionPreflightTests(unittest.TestCase):
    def test_primary_reservations_fit_frozen_calibration(self):
        plan = json.loads(plan_v4.PLAN.read_text())
        limits = execute_v4_stage.limits(plan)
        self.assertEqual(limits["total"]["wall"], 3400)
        self.assertEqual(limits["primary"]["wall"], 980 + 16 * 25)
        self.assertLessEqual(
            sum(
                v["wall"]
                for k, v in limits.items()
                if k not in ("total", "disk_bytes")
            ),
            3400,
        )
        self.assertLessEqual(
            sum(
                v["work"]
                for k, v in limits.items()
                if k not in ("total", "disk_bytes")
            ),
            limits["total"]["work"],
        )

    def test_owned_control_and_challenger_validate_tiny_split(self):
        plan = json.loads(plan_v4.PLAN.read_text())
        fixture = dict(n=10403, factors=[101, 103])
        for row in plan["calibration"]["calls"][:2]:
            modules = execute_v4.runtime_modules(row)
            result = execute_v4.sieve(
                row, fixture, modules, seconds=1, work=1000000, blocks=1
            )
            self.assertEqual(sorted(result["factors"]), [101, 103])
            self.assertEqual(result["reason"], "factor_found")

    def test_worker_output_exhaustion_is_harness_failure(self):
        process = subprocess.Popen(
            [sys.executable, "-c", "print('x' * 65537)"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        with self.assertRaisesRegex(RuntimeError, "recorder_exhaustion"):
            execute_v4_stage.read_worker(process, time.monotonic() + 3)
        self.assertIsNotNone(process.poll())
