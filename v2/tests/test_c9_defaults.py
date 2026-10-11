"""C9 defaults match the confirmed bundle and preserve saved policy."""

import json
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import (
    GrowingDoubleLargeSieveConfig,
    SIQSConfig,
    SIQSJob,
    default_siqs_config,
)
from v2.qs.defaults import (
    C9_MEMORY_BYTES,
    C9_PORTFOLIO_MEMORY_BYTES,
    C9_SECONDS,
    C9_WORK_LIMIT,
    C10_POWER_PLAN_BYTES,
    uses_c9_defaults,
)
from v2.tests.test_siqs_cli import cli

ROOT = Path(__file__).resolve().parents[2]
NUMBER = 10**49 + 123


def no_work():
    return Budget(work_limit=0, seconds=None, cpu_seconds=None)


class C9DefaultTests(unittest.TestCase):
    def test_default_is_exact_confirmed_configuration(self):
        path = ROOT / "v2/benchmarks/inputs/controls/c9_long_selected.json"
        accepted = json.loads(path.read_text())["configurations"]["50"][
            "challenger"
        ]
        selected = asdict(default_siqs_config(NUMBER))
        self.assertEqual(selected.pop("solve_policy"), "collection")
        self.assertEqual(
            selected["collector"]["power_plan_bytes"], C10_POWER_PLAN_BYTES
        )
        # C10/C10a-05 add only their confirmed plan/pacing to this bundle.
        selected["collector"]["power_plan_bytes"] = 0
        self.assertEqual(selected, accepted)
        self.assertEqual(SIQSJob(NUMBER).config, default_siqs_config(NUMBER))
        explicit = SIQSConfig()
        self.assertIs(SIQSJob(NUMBER, config=explicit).config, explicit)

    def test_omitted_grants_use_confirmation_limits_and_explicit_grants_win(
        self,
    ):
        job = SIQSJob(NUMBER)
        self.assertEqual(job.budget.work_limit, C9_WORK_LIMIT)
        self.assertEqual(job.budget.seconds, C9_SECONDS)
        self.assertEqual(job.budget.cpu_seconds, C9_SECONDS)
        explicit = no_work()
        self.assertIs(SIQSJob(NUMBER, budget=explicit).budget, explicit)

        with patch(
            "v2.portfolio._advance",
            side_effect=BudgetExhaustedError("work_limit"),
        ):
            run = factorize_bounded(NUMBER)
        limits = run.checkpoint["payload"]["limits"]
        self.assertEqual(limits["work"], C9_WORK_LIMIT)
        self.assertEqual(limits["seconds"], C9_SECONDS)
        self.assertEqual(limits["cpu_seconds"], C9_SECONDS)

    def test_exact_size_backend_and_sign_scope(self):
        for number, enabled in (
            (10**49 - 1, False),
            (10**49, True),
            (10**50 - 1, True),
            (10**50, False),
            (10**29, False),
            (10**39, False),
        ):
            for signed in (number, -number):
                self.assertEqual(uses_c9_defaults(signed), enabled)
        self.assertFalse(uses_c9_defaults(NUMBER, backend="gmpy2-mpz"))

    def test_parent_default_and_explicit_memory_are_finite(self):
        for original in (
            PortfolioConfig(),
            PortfolioConfig(memory_bytes=128 * 2**20),
        ):
            selected = original.for_input(NUMBER)
            self.assertIsInstance(
                selected.siqs.collector, GrowingDoubleLargeSieveConfig
            )
            expected_memory = (
                128 * 2**20
                if not original._automatic_memory
                else C9_PORTFOLIO_MEMORY_BYTES
            )
            self.assertEqual(selected.memory_bytes, expected_memory)
            self.assertEqual(
                selected.siqs.memory_bytes,
                selected.memory_bytes - selected.workspace_reserve - 8192
                if not original._automatic_memory
                else C9_MEMORY_BYTES,
            )
            self.assertIs(selected.for_input(10**29), selected)
            self.assertIsNone(selected.siqs.collector.max_partials)
        self.assertIsNone(PortfolioConfig(siqs=None).for_input(NUMBER).siqs)
        explicit = SIQSConfig()
        self.assertIs(
            PortfolioConfig(siqs=explicit, memory_bytes=80 * 2**20)
            .for_input(NUMBER)
            .siqs,
            explicit,
        )

    def test_default_and_legacy_checkpoints_keep_concrete_policy(self):
        for config in (None, PortfolioConfig(siqs=None)):
            first = factorize_bounded(NUMBER, config=config, budget=no_work())
            for _ in range(2):
                restored = factorize_bounded(
                    NUMBER, checkpoint=first.checkpoint, budget=no_work()
                )
                self.assertEqual(restored.result.reconstruct(), NUMBER)
                self.assertEqual(
                    restored.checkpoint["payload"]["config"],
                    first.checkpoint["payload"]["config"],
                )

    def test_cli_default_explicit_override_and_implicit_resume(self):
        for options, growing, memory in (
            ((), True, 1040),
            (("--method", "siqs"), True, 1040),
            (("--method", "siqs", "--memory-mib", 128), True, 128),
            (("--method", "siqs", "--qs-base-bound", 1000), False, 80),
        ):
            with (
                self.subTest(options=options),
                tempfile.TemporaryDirectory() as directory,
            ):
                path = Path(directory) / "checkpoint.json"
                first = cli(
                    NUMBER, *options, "--work-limit", 0, "--checkpoint", path
                )
                self.assertEqual(first.returncode, 1, first.stderr)
                self.assertNotIn("Traceback", first.stderr)
                saved = json.loads(path.read_text())["payload"]["config"]
                self.assertEqual(saved["memory_bytes"], memory * 2**20)
                self.assertEqual(
                    saved["siqs"]["collector"].get("storage_policy") == "grow",
                    growing,
                )
                resumed = cli(
                    "--resume", path, "--work-limit", 0, "--checkpoint", path
                )
                self.assertEqual(resumed.returncode, 1, resumed.stderr)
                self.assertNotIn("Traceback", resumed.stderr)
                self.assertEqual(
                    json.loads(path.read_text())["payload"]["config"], saved
                )
