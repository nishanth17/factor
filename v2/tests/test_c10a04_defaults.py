"""User-directed larger defaults retain exact setup and saved policies."""

import json
import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.benchmarks.qs.c10a04 import configuration, short_decisions
from v2.execution.budget import BudgetExhaustedError
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import CountedSIQSConfig, SIQSConfig, SIQSJob, default_siqs_config
from v2.qs.count_policy import CountPolicyRefusedError, endpoint_for_count
from v2.qs.defaults import default_siqs_budget, uses_siqs_defaults
from v2.tests.qs_default_controls import previous_siqs_config
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_siqs_cli import cli


def allowance(work=10**10):
    return Budget(work_limit=work, seconds=30, cpu_seconds=30)


class LargerDefaultTests(unittest.TestCase):
    def test_exact_bands_counts_and_explicit_limits(self):
        for digits in range(50, 61):
            for n in (10 ** (digits - 1), 10**digits - 1):
                for signed in (n, -n):
                    self.assertTrue(uses_siqs_defaults(signed))
                    config = default_siqs_config(signed)
                    if digits < 52:
                        self.assertIs(type(config), SIQSConfig)
                        self.assertEqual(config.factor_count, 6)
                    else:
                        self.assertIs(type(config), CountedSIQSConfig)
                        self.assertEqual(
                            config.base_count, 2500 if digits <= 55 else 5000
                        )
                        self.assertEqual(
                            config.factor_count, 7 if digits <= 55 else 8
                        )
                        self.assertEqual(config.multiplier, 0)
                        self.assertEqual(config.base_bound, 250000)
                        self.assertEqual(config.gray_limit, 64)
                    self.assertEqual(config.half_width, 65536)
                    self.assertEqual(config.memory_bytes, 2**30)
                    self.assertEqual(
                        default_siqs_budget(signed).seconds,
                        1200 if digits <= 55 else 1800,
                    )
        self.assertFalse(uses_siqs_defaults(10**49 - 1))
        self.assertFalse(uses_siqs_defaults(10**60))
        self.assertEqual(
            previous_siqs_config(10**60), SIQSConfig(sign_policy="swap_high")
        )
        self.assertFalse(uses_siqs_defaults(10**59, backend="gmpy2-mpz"))
        n = 10**59 + 123
        explicit = SIQSConfig()
        self.assertIs(SIQSJob(n, config=explicit).config, explicit)
        grant = allowance(0)
        self.assertIs(SIQSJob(n, budget=grant).budget, grant)
        limited = adapt_resource_configuration(
            default_siqs_config,
            n,
            recipe="calibrated"
            if 128 * 2**20 in (None, 2**30)
            else "reference",
        )
        self.assertEqual(limited.memory_bytes, 128 * 2**20)
        self.assertEqual(limited.collector.power_plan_bytes, 0)
        self.assertEqual(limited.solve_policy, "window")

    def test_explicit_previous_60_digit_bundle_matches_frozen_control(self):
        corpus = json.loads(short_decisions.CORPUS.read_text())
        for fixture in corpus["fixtures"]:
            expected, divisor = configuration.resolve(
                fixture["n"], 60, "compact", allowance()
            )
            self.assertIsNone(divisor)
            job = SIQSJob(
                fixture["n"],
                config=previous_siqs_config(fixture["n"]),
                budget=allowance(),
            )
            job.run(max_blocks=1)
            self.assertEqual(len(job.base.entries), 5000)
            self.assertEqual(job.base.bound, expected.base_bound)
            self.assertEqual(job.multiplier, expected.multiplier)
            values = asdict(job.config)
            values.pop("base_count")
            values.update(base_bound=job.base.bound, multiplier=job.multiplier)
            self.assertEqual(values.pop("sign_policy"), "swap_high")
            self.assertEqual(
                values["collector"].pop("sparse_queue_bytes"), 32 * 2**20
            )
            legacy = asdict(expected)
            legacy.pop("sign_policy")
            legacy["collector"].pop("sparse_queue_bytes")
            self.assertEqual(values, legacy)
            self.assertLess(job.budget.used, 10**10)

    def test_cli_checkpoint_keeps_counted_policy(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.json"
            first = cli(
                10**59 + 123,
                "--method",
                "siqs",
                "--work-limit",
                0,
                "--checkpoint",
                path,
            )
            self.assertEqual(first.returncode, 1, first.stderr)
            saved = json.loads(path.read_text())["payload"]["config"]
            self.assertEqual(saved["siqs"]["base_count"], 5000)
            resumed = cli(
                "--resume", path, "--work-limit", 0, "--checkpoint", path
            )
            self.assertEqual(resumed.returncode, 1, resumed.stderr)
            self.assertEqual(
                json.loads(path.read_text())["payload"]["config"], saved
            )

    def test_counted_resume_rebuilds_same_base_and_charges_replay(self):
        config = CountedSIQSConfig(
            base_bound=500, base_count=10, factor_count=3, pool_size=16
        )
        job = SIQSJob(1009 * 1013, config=config, budget=allowance())
        job.run(max_blocks=1)
        checkpoint = json.loads(json.dumps(job.checkpoint()))
        used = job.budget.used
        restored = SIQSJob.from_checkpoint(checkpoint, budget=allowance())
        self.assertIs(type(restored.config), CountedSIQSConfig)
        self.assertEqual(restored.config, config)
        self.assertEqual(restored.base.primes, job.base.primes)
        self.assertGreater(restored.budget.used, used)
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                checkpoint,
                budget=allowance(),
                config=replace(config, base_count=11),
                allow_extension=True,
            )

    def test_portfolio_serialization_and_zero_work_resume_keep_policy(self):
        n = 10**59 + 123
        selected = PortfolioConfig().for_input(n)
        self.assertEqual(selected.memory_bytes, 1040 * 2**20)
        self.assertIs(type(selected.siqs), CountedSIQSConfig)
        self.assertEqual(
            adapt_resource_configuration(
                PortfolioConfig, **asdict(selected)
            ).siqs,
            selected.siqs,
        )
        self.assertIsNone(PortfolioConfig(siqs=None).for_input(n).siqs)
        first = factorize_bounded(n, budget=allowance(0))
        saved = first.checkpoint["payload"]["config"]
        self.assertEqual(saved["siqs"]["base_count"], 5000)
        resumed = factorize_bounded(
            n, checkpoint=first.checkpoint, budget=allowance(0)
        )
        self.assertEqual(resumed.checkpoint["payload"]["config"], saved)
        self.assertEqual(resumed.result.reconstruct(), n)

    def test_count_resolution_memory_budget_and_refusal(self):
        n = 1009 * 1013
        primes = [p for p in range(2, 500) if all(p % d for d in range(2, p))]
        eligible = [
            p for p in primes if any(x * x % p == n % p for x in range(p))
        ]
        bound, divisor = endpoint_for_count(n, 1, 10, allowance(), ceiling=500)
        self.assertIsNone(divisor)
        self.assertEqual(bound, eligible[9] + 1)
        self.assertEqual(
            endpoint_for_count(77, 7, 10, allowance(), ceiling=500),
            (None, 7),
        )
        with self.assertRaises(TypeError):
            endpoint_for_count(
                n, 1, 10, allowance(), ceiling=500, memory_bytes=0
            )
        with self.assertRaises(BudgetExhaustedError):
            endpoint_for_count(n, 1, 10, allowance(0), ceiling=500)
        config = CountedSIQSConfig(base_bound=10, base_count=9)
        result = SIQSJob(n, config=config, budget=allowance()).run()
        self.assertEqual(result.reason, "configuration_refusal")
        self.assertEqual((result.divisor or 1) * result.cofactor, n)
        with patch(
            "v2.qs.count_policy.prime_sieve.prime_sieve", return_value=[]
        ):
            with self.assertRaises(CountPolicyRefusedError):
                endpoint_for_count(n, 1, 10, allowance(), ceiling=500)
