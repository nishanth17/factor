"""Scoped production selection and concrete legacy checkpoint policy."""

import json
import unittest
from dataclasses import asdict, replace
from pathlib import Path

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import SIQSConfig, SIQSJob, default_siqs_config
from v2.qs.defaults import C9_MEMORY_BYTES
from v2.tests.qs_default_controls import previous_siqs_config
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
                selected = asdict(previous_siqs_config(signed))
                self.assertEqual(selected.pop("sign_policy"), "swap_high")
                self.assertEqual(selected.pop("polynomial_policy"), "standard")
                self.assertEqual(selected.pop("multiplier_policy"), "standard")
                self.assertEqual(
                    selected.pop("root_update_policy"), "standard"
                )
                self.assertEqual(
                    selected["collector"].pop("sparse_queue_bytes"), 32 * 2**20
                )
                self.assertEqual(
                    selected["collector"]["marking"],
                    "density+presieve+scan_fused"
                    if digits == 50
                    else "sparse",
                )
                selected["collector"]["marking"] = "sparse"
                self.assertEqual(
                    selected,
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
                adapt_resource_configuration(
                    default_siqs_config,
                    number,
                    recipe="calibrated"
                    if memory in (None, 2**30)
                    else "reference",
                ).solve_policy,
                "window",
            )
        for other in (10**29, 10**40, 10**49 - 1, 10**60):
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
        current = previous_siqs_config(number)
        legacy = replace(
            current,
            solve_policy="window",
            sign_policy="canonical",
            collector=replace(current.collector, sparse_queue_bytes=0),
        )
        job = SIQSJob(number, config=legacy, budget=allowance())
        checkpoint = json.loads(json.dumps(job.checkpoint()))
        self.assertEqual(checkpoint["version"], 5)
        self.assertNotIn(
            "solve_policy", json.loads(checkpoint["blob"])["config"]
        )
        restored = SIQSJob.from_checkpoint(checkpoint, budget=allowance())
        self.assertEqual(restored.config, legacy)

        config = PortfolioConfig(siqs=legacy)
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
