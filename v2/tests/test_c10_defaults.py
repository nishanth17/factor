"""C10 automatic plan routing matches frozen bundles and saved policy."""

import json
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path

from v2.execution.budget import Budget
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import SIQSConfig, SIQSJob, default_siqs_config
from v2.qs.defaults import (
    C10_POWER_PLAN_BYTES,
    default_siqs_budget,
    uses_siqs_defaults,
)
from v2.tests.test_siqs_cli import cli

ROOT = Path(__file__).resolve().parents[2]
NUMBER_40 = 10**39 + 123
NUMBER_50 = 10**49 + 123
NUMBER_51 = 10**50 + 123


def no_work():
    return Budget(work_limit=0, seconds=None, cpu_seconds=None)


class C10DefaultTests(unittest.TestCase):
    def test_selected_bundles_match_frozen_plan_and_control_delta(self):
        path = ROOT / "v2/benchmarks/inputs/controls/c10_training.json"
        configurations = json.loads(path.read_text())["configurations"]

        for digits, number in (("40", NUMBER_40), ("50", NUMBER_50)):
            with self.subTest(digits=digits):
                selected = asdict(default_siqs_config(number))
                self.assertEqual(
                    selected.pop("solve_policy"),
                    "collection" if digits == "50" else "window",
                )
                self.assertEqual(selected, configurations[digits]["plan"])
                selected["collector"]["power_plan_bytes"] = 0
                self.assertEqual(selected, configurations[digits]["control"])
                self.assertEqual(
                    SIQSJob(number).config, default_siqs_config(number)
                )

    def test_exact_native_bands_and_resource_scope(self):
        for number, enabled in (
            (10**39 - 1, False),
            (10**39, True),
            (10**40 - 1, True),
            (10**40, False),
            (10**49 - 1, False),
            (10**49, True),
            (10**50 - 1, True),
            (10**50, True),
            (10**51 - 1, True),
            (10**51, False),
        ):
            for signed in (number, -number):
                self.assertEqual(uses_siqs_defaults(signed), enabled)
        self.assertFalse(uses_siqs_defaults(NUMBER_40, backend="gmpy2-mpz"))
        self.assertEqual(
            default_siqs_config(10**40).collector.power_plan_bytes, 0
        )
        self.assertEqual(SIQSConfig().collector.power_plan_bytes, 0)

        for number in (NUMBER_40, NUMBER_50, NUMBER_51):
            selected = default_siqs_config(number, memory_bytes=128 * 2**20)
            self.assertEqual(selected.memory_bytes, 128 * 2**20)
            self.assertEqual(selected.collector.power_plan_bytes, 0)
            self.assertEqual(default_siqs_budget(number).work_limit, 10**13)
            explicit = SIQSConfig()
            self.assertIs(SIQSJob(number, config=explicit).config, explicit)
            granted = no_work()
            self.assertIs(SIQSJob(number, budget=granted).budget, granted)

    def test_portfolio_and_checkpoint_keep_concrete_selection(self):
        for number in (NUMBER_40, NUMBER_50, NUMBER_51):
            selected = PortfolioConfig().for_input(number)
            self.assertEqual(selected.memory_bytes, 1040 * 2**20)
            self.assertEqual(
                selected.siqs.collector.power_plan_bytes,
                C10_POWER_PLAN_BYTES,
            )
            self.assertEqual(
                PortfolioConfig(**asdict(selected)).siqs, selected.siqs
            )
            self.assertEqual(
                PortfolioConfig(memory_bytes=128 * 2**20)
                .for_input(number)
                .siqs.collector.power_plan_bytes,
                0,
            )
            self.assertIsNone(
                PortfolioConfig(siqs=None).for_input(number).siqs
            )

            first = factorize_bounded(number, budget=no_work())
            self.assertEqual(first.result.reconstruct(), number)
            saved = first.checkpoint["payload"]["config"]
            self.assertEqual(
                saved["siqs"]["collector"]["power_plan_bytes"],
                C10_POWER_PLAN_BYTES,
            )
            restored = factorize_bounded(
                number, checkpoint=first.checkpoint, budget=no_work()
            )
            self.assertEqual(restored.result.reconstruct(), number)
            self.assertEqual(restored.checkpoint["payload"]["config"], saved)

    def test_cli_default_and_resume_at_40_digits(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.json"
            first = cli(
                NUMBER_40,
                "--method",
                "siqs",
                "--work-limit",
                0,
                "--checkpoint",
                path,
            )
            self.assertEqual(first.returncode, 1, first.stderr)
            saved = json.loads(path.read_text())["payload"]["config"]
            self.assertEqual(saved["memory_bytes"], 1040 * 2**20)
            self.assertEqual(
                saved["siqs"]["collector"]["power_plan_bytes"],
                C10_POWER_PLAN_BYTES,
            )

            restored = cli(
                "--resume", path, "--work-limit", 0, "--checkpoint", path
            )
            self.assertEqual(restored.returncode, 1, restored.stderr)
            self.assertEqual(
                json.loads(path.read_text())["payload"]["config"], saved
            )


if __name__ == "__main__":
    unittest.main()
