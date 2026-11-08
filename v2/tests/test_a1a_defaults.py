"""User-directed A1a routing preserves explicit and saved configurations."""

import json
import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import SIQSConfig, SIQSJob, default_siqs_config
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_siqs_cli import cli

SELECTED = "density+presieve+scan_fused"
DEFAULT_SELECTED = SELECTED + "+parity_two"
EXTENSION = "density+presieve_powers+unroll4+scan_fused"


class A1aDefaultTests(unittest.TestCase):
    def test_both_edges_signs_and_owned_memory_scope(self):
        for digits in (40, 49, *range(50, 62)):
            for n in (10 ** (digits - 1), 10**digits - 1):
                for signed in (n, -n):
                    for memory in (
                        None,
                        128 * 2**20,
                        2**30 - 1,
                        2**30,
                        2**30 + 1,
                    ):
                        config = adapt_resource_configuration(
                            default_siqs_config,
                            signed,
                            recipe="calibrated"
                            if memory in (None, 2**30)
                            else "reference",
                        )
                        selected = 50 <= digits <= 55 and memory in (
                            None,
                            2**30,
                        )
                        self.assertEqual(
                            config.collector.marking,
                            DEFAULT_SELECTED if selected else "parity_two",
                        )
        for digits in range(50, 56):
            config = default_siqs_config(
                10 ** (digits - 1), backend="gmpy2-mpz"
            )
            self.assertEqual(config.collector.marking, "parity_two")

    def test_explicit_and_parent_configs_take_precedence(self):
        n = 10**54 + 123
        explicit = SIQSConfig()
        self.assertIs(SIQSJob(n, config=explicit).config, explicit)
        self.assertEqual(explicit.collector.marking, "sparse")
        self.assertIsNone(PortfolioConfig(siqs=None).for_input(n).siqs)
        self.assertEqual(
            PortfolioConfig().for_input(n).siqs.collector.marking,
            DEFAULT_SELECTED,
        )

    def test_default_cli_and_saved_selector_survive_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.json"
            first = cli(
                10**54 + 123,
                "--method",
                "siqs",
                "--work-limit",
                0,
                "--checkpoint",
                path,
            )
            self.assertEqual(first.returncode, 1, first.stderr)
            saved = json.loads(path.read_text())["payload"]["config"]
            self.assertEqual(
                saved["siqs"]["collector"]["marking"], DEFAULT_SELECTED
            )
            resumed = cli(
                "--resume", path, "--work-limit", 0, "--checkpoint", path
            )
            self.assertEqual(resumed.returncode, 1, resumed.stderr)
            self.assertEqual(
                json.loads(path.read_text())["payload"]["config"], saved
            )

    def test_extension_stays_opt_in_and_retains_saved_policy(self):
        n = 10**54 + 123
        default = default_siqs_config(n)
        extended = replace(
            default, collector=replace(default.collector, marking=EXTENSION)
        )
        self.assertEqual(default.collector.marking, DEFAULT_SELECTED)
        config = PortfolioConfig(siqs=extended)
        self.assertEqual(
            adapt_resource_configuration(
                PortfolioConfig, **asdict(config)
            ).siqs,
            extended,
        )

        first = factorize_bounded(
            n, config=config, budget=Budget(work_limit=0)
        )
        saved = first.checkpoint["payload"]["config"]
        self.assertEqual(saved["siqs"]["collector"]["marking"], EXTENSION)
        resumed = factorize_bounded(
            n, checkpoint=first.checkpoint, budget=Budget(work_limit=0)
        )
        self.assertEqual(resumed.checkpoint["payload"]["config"], saved)
        self.assertEqual(resumed.result.reconstruct(), n)
