"""Preflight the v4 worker source boundary and finite reservation math."""

import json
import subprocess
import sys
import tempfile
import time
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.qs.c10a01 import execute_v4, execute_v4_stage, plan_v4
from v2.benchmarks.qs.c10a04.accounting import Ledger


class ExecutionPreflightTests(unittest.TestCase):
    def test_isolated_control_budget_exception_keeps_unresolved(self):
        class ControlBudgetExhaustedError(Exception):
            pass

        class ControlBudget:
            def consume(self, amount):
                raise ControlBudgetExhaustedError()

        fixture = dict(n=10403, factors=[101, 103])
        result = SimpleNamespace(divisor=101, cofactor=103)
        validated = execute_v4.validate_result(
            fixture,
            result,
            ControlBudget(),
            None,
            ControlBudgetExhaustedError,
        )
        self.assertEqual(validated["factors"], [101, 103])
        self.assertEqual(validated["certainty"], ["unclassified"] * 2)
        self.assertEqual(validated["unresolved"], [101, 103])

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
        try:
            for row in plan["calibration"]["calls"][:2]:
                modules = execute_v4.runtime_modules(row)
                result = execute_v4.sieve(
                    row,
                    fixture,
                    modules,
                    seconds=1,
                    work=1000000,
                    blocks=1,
                )
                self.assertEqual(sorted(result["factors"]), [101, 103])
                self.assertEqual(result["reason"], "factor_found")
        finally:
            for name in tuple(sys.modules):
                if name == execute_v4.control_v4.PACKAGE or name.startswith(
                    execute_v4.control_v4.PACKAGE + "."
                ):
                    del sys.modules[name]

    def test_worker_work_is_settled_by_durable_ledger(self):
        plan = json.loads(plan_v4.PLAN.read_text())
        row = plan["calibration"]["calls"][0]
        result = dict(
            before=dict(wall=4, cpu=4, work=100),
            result=dict(
                seconds=2,
                cpu_seconds=2,
                work=300,
                reason="paused",
                positions=row["positions"],
            ),
            after=dict(wall=1, cpu=1, work=200),
            total_work=600,
        )
        process = type("FinishedProcess", (), {"poll": lambda self: 0})()
        with patch.object(
            execute_v4_stage.subprocess, "Popen", return_value=process
        ):
            with patch.object(
                execute_v4_stage, "read_worker", return_value=result
            ):
                reported = execute_v4_stage.run_worker(row)
        self.assertEqual(reported["work"], 600)
        with tempfile.TemporaryDirectory() as directory:
            limits = dict(
                primary=dict(wall=60, cpu=60, work=1000),
                total=dict(wall=60, cpu=60, work=1000),
                disk_bytes=1024 * 1024,
            )
            ledger = Ledger(directory, "settlement-test", limits)
            ledger.call("primary", "cell", 60, 1000, {}, lambda: reported)
            self.assertEqual(ledger.used("primary", "work"), 600)

    def test_worker_output_exhaustion_is_harness_failure(self):
        process = subprocess.Popen(
            [sys.executable, "-c", "print('x' * 65537)"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        with self.assertRaisesRegex(RuntimeError, "recorder_exhaustion"):
            execute_v4_stage.read_worker(process, time.monotonic() + 3)
        self.assertIsNotNone(process.poll())
