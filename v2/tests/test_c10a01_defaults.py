"""All-size automatic policies, exact prefixes and legacy resume isolation."""

import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.qs import SIQSConfig, SIQSJob, default_siqs_config
from v2.qs.collector_kernels import components
from v2.qs.defaults import default_siqs_budget
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_siqs_cli import cli


def allowance():
    return Budget(work_limit=10**10, seconds=30, cpu_seconds=30)


class CoupledDefaultTests(unittest.TestCase):
    def test_every_size_and_backend_selects_the_coupled_policies(self):
        for digits in (
            1,
            10,
            30,
            39,
            40,
            41,
            49,
            50,
            51,
            52,
            55,
            56,
            60,
            61,
            62,
            65,
            70,
            80,
            100,
            150,
        ):
            for backend in ("python-int", "gmpy2-mpz"):
                for sign in (-1, 1):
                    n = sign * (10 ** (digits - 1) + 7)
                    config = default_siqs_config(n, backend=backend)
                    self.assertEqual(config.multiplier, 0)
                    self.assertEqual(config.polynomial_policy, "parity2")
                    self.assertEqual(config.multiplier_policy, "parity2_wide")
                    self.assertEqual(config.root_update_policy, "batch")
                    self.assertEqual(config.assignment_policy, "windowed_pair")
                    self.assertIn(
                        "parity_two", components(config.collector.marking)
                    )
                    self.assertEqual(config.backend, backend)
        self.assertEqual(default_siqs_budget(10**69 + 7).work_limit, 2000000)
        self.assertEqual(default_siqs_budget(10**69 + 7).seconds, 30)

    def test_existing_marking_memory_and_explicit_controls_are_preserved(self):
        config = default_siqs_config(10**49 + 123)
        self.assertEqual(
            config.collector.marking, "density+presieve+scan_fused+parity_two"
        )
        self.assertEqual(config.collector.sparse_queue_bytes, 32 * 2**20)
        limited = adapt_resource_configuration(
            default_siqs_config,
            10**49 + 123,
            recipe="calibrated"
            if 128 * 2**20 in (None, 2**30)
            else "reference",
        )
        self.assertEqual(limited.memory_bytes, 128 * 2**20)
        self.assertEqual(limited.collector.power_plan_bytes, 0)
        explicit = SIQSConfig()
        self.assertIs(SIQSJob(10403, config=explicit).config, explicit)
        self.assertEqual(explicit.polynomial_policy, "standard")
        self.assertEqual(explicit.multiplier, 1)

    def test_actual_default_roots_match_direct_integer_evaluation_and_resume(
        self,
    ):
        n = 1000003 * 1000033
        job = SIQSJob(n, seed=7, budget=allowance())
        result = job.run(max_blocks=2)
        self.assertEqual((result.divisor or 1) * result.cofactor, n)
        self.assertIsNotNone(job.family)
        step = job.family.current
        for roots in step.roots:
            if roots.prime > 97:
                break
            expected = {
                x
                for x in range(roots.prime)
                if step.polynomial.value(x) % roots.prime == 0
            }
            actual = (
                set(range(roots.prime))
                if roots.all_positions
                else set(roots.roots)
            )
            self.assertEqual(actual, expected)
        saved = json.loads(json.dumps(job.checkpoint()))
        restored = SIQSJob.from_checkpoint(saved, budget=allowance())
        self.assertEqual(restored.config, job.config)
        self.assertGreaterEqual(restored.budget.used, job.budget.used)
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                saved,
                budget=allowance(),
                config=replace(job.config, root_update_policy="standard"),
                allow_extension=True,
            )

    def test_legacy_checkpoint_does_not_acquire_automatic_policies(self):
        config = SIQSConfig(base_bound=97, half_width=64)
        job = SIQSJob(1009 * 1013, config=config, budget=allowance())
        job.run(max_blocks=1)
        saved = json.loads(json.dumps(job.checkpoint()))
        with patch(
            "v2.qs.defaults.default_siqs_config",
            side_effect=AssertionError("automatic selection on restore"),
        ):
            restored = SIQSJob.from_checkpoint(saved, budget=allowance())
        self.assertEqual(restored.config, config)
        self.assertEqual(restored.config.root_update_policy, "standard")

    def test_cli_uses_coupled_defaults_outside_calibrated_bands(self):
        for digits in (30, 45, 61, 70, 100):
            with tempfile.TemporaryDirectory() as folder:
                path = Path(folder) / "checkpoint.json"
                result = cli(
                    10 ** (digits - 1) + 123,
                    "--method",
                    "siqs",
                    "--work-limit",
                    0,
                    "--checkpoint",
                    path,
                )
                self.assertEqual(result.returncode, 1, result.stderr)
                config = json.loads(path.read_text())["payload"]["config"][
                    "siqs"
                ]
                self.assertEqual(config["polynomial_policy"], "parity2")
                self.assertEqual(config["multiplier_policy"], "parity2_wide")
                self.assertEqual(config["root_update_policy"], "batch")
                self.assertEqual(config["assignment_policy"], "windowed_pair")


if __name__ == "__main__":
    unittest.main()
