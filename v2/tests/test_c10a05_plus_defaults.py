"""The confirmed 51-digit bundle routes without changing saved policies."""

import json
import unittest
from dataclasses import asdict
from pathlib import Path

from v2.execution.budget import Budget
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import SIQSConfig, SIQSJob, default_siqs_config
from v2.qs.defaults import C9_MEMORY_BYTES, default_siqs_budget

ROOT = Path(__file__).resolve().parents[2]
NUMBER_51 = 10**50 + 123


def no_work():
    return Budget(work_limit=0, seconds=None, cpu_seconds=None)


class LargerBandDefaultTests(unittest.TestCase):
    def test_exact_51_digit_selected_configuration_and_explicit_scope(self):
        path = ROOT / "v2/benchmarks/inputs/controls/c10a05_plus_training.json"
        control = json.loads(path.read_text())
        expected = control["configurations"]["collection"]
        for number in (10**50, NUMBER_51, 10**51 - 1):
            for signed in (number, -number):
                self.assertEqual(asdict(default_siqs_config(signed)), expected)
                if signed > 0 and signed % 2:
                    self.assertEqual(
                        SIQSJob(signed).config, default_siqs_config(signed)
                    )

        self.assertEqual(default_siqs_budget(NUMBER_51).work_limit, 10**13)
        self.assertEqual(default_siqs_budget(NUMBER_51).seconds, 1200)
        self.assertEqual(default_siqs_config(10**51), SIQSConfig())
        self.assertEqual(
            default_siqs_config(NUMBER_51, backend="gmpy2-mpz").solve_policy,
            "window",
        )
        self.assertEqual(
            default_siqs_config(
                NUMBER_51, memory_bytes=128 * 2**20
            ).solve_policy,
            "window",
        )
        explicit = SIQSConfig()
        self.assertIs(
            SIQSJob(NUMBER_51, config=explicit, budget=no_work()).config,
            explicit,
        )

    def test_portfolio_checkpoint_repeats_concrete_51_digit_choice(self):
        selected = PortfolioConfig().for_input(NUMBER_51)
        self.assertEqual(selected.memory_bytes, 1040 * 2**20)
        self.assertEqual(selected.siqs.memory_bytes, C9_MEMORY_BYTES)
        self.assertEqual(selected.siqs.solve_policy, "collection")

        first = factorize_bounded(NUMBER_51, budget=no_work())
        saved = first.checkpoint["payload"]["config"]
        self.assertEqual(saved["siqs"]["solve_policy"], "collection")
        self.assertEqual(first.result.reconstruct(), NUMBER_51)
        for _ in range(2):
            first = factorize_bounded(
                NUMBER_51, checkpoint=first.checkpoint, budget=no_work()
            )
            self.assertEqual(first.result.reconstruct(), NUMBER_51)
            self.assertEqual(first.checkpoint["payload"]["config"], saved)


if __name__ == "__main__":
    unittest.main()
