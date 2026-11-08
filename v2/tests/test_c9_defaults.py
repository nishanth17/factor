"""C9 defaults match the confirmed bundle and preserve saved policy."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.execution.budget import BudgetExhaustedError
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import (
    GrowingDoubleLargeSieveConfig,
    SIQSConfig,
    SIQSJob,
    default_siqs_config,
)
from v2.qs.defaults import (
    C10_POWER_PLAN_BYTES,
    uses_c9_defaults,
)
from v2.tests.qs_default_controls import previous_siqs_config
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.resource_stops import recipe_values as asdict
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
        selected = asdict(previous_siqs_config(NUMBER))
        self.assertEqual(selected.pop("sign_policy"), "swap_high")
        self.assertEqual(selected.pop("polynomial_policy"), "standard")
        self.assertEqual(selected.pop("multiplier_policy"), "standard")
        self.assertEqual(selected.pop("root_update_policy"), "standard")
        self.assertEqual(
            selected["collector"].pop("sparse_queue_bytes"), 32 * 2**20
        )
        self.assertEqual(selected.pop("solve_policy"), "collection")
        self.assertEqual(
            selected["collector"]["power_plan_bytes"], C10_POWER_PLAN_BYTES
        )
        self.assertEqual(
            selected["collector"]["marking"], "density+presieve+scan_fused"
        )
        # Remove the later A1a marking and C10 plan/pacing deltas.
        selected["collector"]["marking"] = "sparse"
        selected["collector"]["power_plan_bytes"] = 0
        self.assertEqual(selected, asdict(accepted))
        self.assertEqual(SIQSJob(NUMBER).config, default_siqs_config(NUMBER))
        explicit = SIQSConfig()
        self.assertIs(SIQSJob(NUMBER, config=explicit).config, explicit)

    def test_omitted_grants_use_confirmation_limits_and_explicit_grants_win(
        self,
    ):
        job = SIQSJob(NUMBER)
        self.assertFalse(hasattr(job.budget, "work_limit"))
        self.assertIsNone(job.budget.seconds)
        self.assertIsNone(job.budget.cpu_seconds)
        explicit = no_work()
        self.assertIs(SIQSJob(NUMBER, budget=explicit).budget, explicit)

        with patch(
            "v2.portfolio._advance",
            side_effect=BudgetExhaustedError("work_limit"),
        ):
            run = factorize_bounded(NUMBER)
        limits = run.checkpoint["payload"]["limits"]
        self.assertIsNone(limits["work"])
        self.assertIsNone(limits["seconds"])
        self.assertIsNone(limits["cpu_seconds"])

    def test_exact_size_backend_and_sign_scope(self):
        for number, enabled in (
            (10**49 - 1, False),
            (10**49, True),
            (10**50 - 1, True),
            (10**50, True),
            (10**51 - 1, True),
            (10**51, False),
            (10**29, False),
            (10**39, False),
        ):
            for signed in (number, -number):
                self.assertEqual(uses_c9_defaults(signed), enabled)
        self.assertFalse(uses_c9_defaults(NUMBER, backend="gmpy2-mpz"))

    def test_parent_default_and_explicit_memory_are_finite(self):
        for original in (
            PortfolioConfig(),
            PortfolioConfig(),
        ):
            selected = original.for_input(NUMBER)
            self.assertIsInstance(
                selected.siqs.collector, GrowingDoubleLargeSieveConfig
            )
            self.assertFalse(hasattr(selected, "memory_bytes"))
            self.assertFalse(hasattr(selected.siqs, "memory_bytes"))
            self.assertIs(selected.for_input(10**29), selected)
            self.assertIsNone(selected.siqs.collector.max_partials)
        self.assertIsNone(PortfolioConfig(siqs=None).for_input(NUMBER).siqs)
        explicit = SIQSConfig()
        self.assertIs(
            PortfolioConfig(siqs=explicit).for_input(NUMBER).siqs,
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
            (("--method", "siqs", "--rss-bytes", 134217728), True, 128),
            (("--method", "siqs", "--qs-base-bound", 1000), False, 80),
        ):
            with (
                self.subTest(options=options),
                tempfile.TemporaryDirectory() as directory,
            ):
                path = Path(directory) / "checkpoint.json"
                first = cli(
                    NUMBER, *options, "--seconds", "0", "--checkpoint", path
                )
                self.assertEqual(first.returncode, 1, first.stderr)
                self.assertNotIn("Traceback", first.stderr)
                saved = json.loads(path.read_text())["payload"]["config"]
                self.assertNotIn("memory_bytes", saved)
                self.assertEqual(
                    saved["siqs"]["collector"].get("storage_policy") == "grow",
                    growing,
                )
                resumed = cli(
                    "--resume", path, "--seconds", "0", "--checkpoint", path
                )
                self.assertEqual(resumed.returncode, 1, resumed.stderr)
                self.assertNotIn("Traceback", resumed.stderr)
                self.assertEqual(
                    json.loads(path.read_text())["payload"]["config"], saved
                )
