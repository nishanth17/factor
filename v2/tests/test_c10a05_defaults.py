"""Scoped production selection and concrete legacy checkpoint policy."""

import json
import unittest
from dataclasses import asdict, replace
from pathlib import Path

from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import SIQSConfig, SIQSJob, default_siqs_config
from v2.qs.defaults import C9_MEMORY_BYTES
from v2.tests.test_siqs import allowance

ROOT = Path(__file__).resolve().parents[2]


class PacingDefaultTests(unittest.TestCase):
    def test_production_matches_only_the_confirmed_fixed_arm(self):
        control = json.loads(
            (
                ROOT / "v2/benchmarks/inputs/controls/c10a05_training.json"
            ).read_text()
        )["configurations"]
        for digits, arm in ((40, "control"), (50, "collection")):
            number = 10 ** (digits - 1) + 123
            for signed in (number, -number):
                self.assertEqual(
                    asdict(default_siqs_config(signed)),
                    control[str(digits)][arm],
                )
            self.assertEqual(
                PortfolioConfig().for_input(number).siqs,
                default_siqs_config(number),
            )

    def test_other_resources_backends_and_explicit_configs_keep_window(self):
        number = 10**49 + 123
        for memory in (128 * 2**20, C9_MEMORY_BYTES - 1, C9_MEMORY_BYTES + 1):
            self.assertEqual(
                default_siqs_config(number, memory_bytes=memory).solve_policy,
                "window",
            )
        for other in (10**29, 10**40, 10**49 - 1, 10**50):
            self.assertEqual(default_siqs_config(other).solve_policy, "window")
        self.assertEqual(
            default_siqs_config(number, backend="gmpy2-mpz").solve_policy,
            "window",
        )
        explicit = SIQSConfig()
        self.assertIs(SIQSJob(number, config=explicit).config, explicit)
        self.assertEqual(explicit.solve_policy, "window")

    def test_old_50_digit_snapshots_do_not_acquire_new_default(self):
        number = 10**49 + 123
        legacy = replace(default_siqs_config(number), solve_policy="window")
        job = SIQSJob(number, config=legacy, budget=allowance())
        checkpoint = json.loads(json.dumps(job.checkpoint()))
        self.assertEqual(checkpoint["version"], 5)
        self.assertNotIn(
            "solve_policy", json.loads(checkpoint["blob"])["config"]
        )
        restored = SIQSJob.from_checkpoint(checkpoint, budget=allowance())
        self.assertEqual(restored.config, legacy)

        config = PortfolioConfig(siqs=legacy, memory_bytes=1040 * 2**20)
        first = factorize_bounded(number, config=config, budget=allowance(0))
        saved = first.checkpoint["payload"]["config"]
        self.assertNotIn("solve_policy", saved["siqs"])
        resumed = factorize_bounded(
            number, checkpoint=first.checkpoint, budget=allowance(0)
        )
        self.assertEqual(resumed.checkpoint["payload"]["config"], saved)
        self.assertEqual(resumed.result.reconstruct(), number)


if __name__ == "__main__":
    unittest.main()
