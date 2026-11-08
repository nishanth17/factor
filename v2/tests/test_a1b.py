"""Current fixed-score fallback, exact coverage and saved selector checks."""

import json
import unittest
from dataclasses import asdict, replace
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.benchmarks.infrastructure.resource_control import (
    decode_historical_siqs as _decode_siqs_config,
)
from v2.benchmarks.qs.a1b import CONTROL, ready
from v2.qs import Polynomial, SieveCollector, SIQSJob, build_factor_base
from v2.qs.sieve_collector import SieveConfig
from v2.tests.test_qs import reference_positions, unlimited_budget
from v2.tests.test_qs_sieve import signature
from v2.tests.test_siqs import allowance, configuration


class A1bTests(unittest.TestCase):
    def test_literal_arms_change_only_fixed_policy(self):
        protocol = json.loads(CONTROL.read_text())
        self.assertEqual(len({f["n"] for f in protocol["fixtures"]}), 9)
        self.assertEqual(len(protocol["schedule"]), 30)
        for arms in protocol["configurations"].values():
            control = _decode_siqs_config(arms["control"])
            fixed = _decode_siqs_config(arms["fixed"])
            self.assertEqual(
                fixed,
                replace(
                    control,
                    collector=replace(control.collector, score_policy="fixed"),
                ),
            )
            self.assertEqual(control.collector.power_plan_bytes, 2**21)
            self.assertEqual(fixed.collector.sparse_queue_bytes, 32 * 2**20)
            self.assertEqual(_decode_siqs_config(asdict(fixed)), fixed)

    def test_fixed_fallback_presieve_scan_against_independent_division(self):
        # R2 covers clipping, high powers, zeros, p=2 and refusal. This adds
        # today's queue selector plus current kernels, including multiplier/A.
        for multiplier in (1, 3):
            base = build_factor_base(
                10403, bound=100, multiplier=multiplier
            ).factor_base
            coefficients = [(1, 102)]
            for entry in base.entries:
                if entry.prime != 2 and len(coefficients) < 3:
                    coefficients.append(
                        (
                            entry.prime,
                            next(
                                x
                                for x in range(entry.prime)
                                if x * x % entry.prime
                                == base.n_prime % entry.prime
                            ),
                        )
                    )
            for a, b in coefficients:
                poly = Polynomial(10403, multiplier, a, b)
                expected = reference_positions(poly, base, -71, 80, 97)
                for marking in ("sparse", "density+presieve+scan_fused"):
                    cfg = SieveConfig(
                        score_policy="fixed",
                        score_backend="list",
                        division="bucket",
                        power_plan_bytes=2**21,
                        sparse_queue_bytes=32 * 2**20,
                        block_width=17,
                        residual_bound=97,
                        marking=marking,
                        max_partials=4096,
                        preparation_cache_bytes=min(2 * 2**20, 2048 * 4096),
                        target_rows=4096,
                    )
                    worker = SieveCollector(
                        poly, base, config=cfg, budget=unlimited_budget()
                    )
                    with patch(
                        "v2.qs.sparse_queue.build",
                        side_effect=AssertionError("fixed queued"),
                    ):
                        result = worker.collect(-71, 80)
                    self.assertEqual(result.reason, "complete")
                    self.assertEqual(signature(result), expected)
                    self.assertNotIn("queue_builds", result.stats)
                    self.assertGreaterEqual(result.workspace_bytes, 0)

    def test_fixed_saved_policy_restores_with_charges_and_rejects_powers(self):
        cfg = adapt_resource_configuration(
            configuration,
            sign_policy="swap_high",
            checkpoint_bytes=2**20,
            memory_bytes=128 * 2**20,
            collector=SieveConfig(
                score_policy="fixed",
                division="bucket",
                power_plan_bytes=2**20,
                sparse_queue_bytes=2**20,
                residual_bound=500,
                marking="density+presieve+scan_fused",
                max_partials=4096,
                preparation_cache_bytes=min(2 * 2**20, 2048 * 4096),
                target_rows=4096,
            ),
        )
        job = SIQSJob(4001 * 5003, config=cfg, seed=73, budget=allowance())
        job.run(max_blocks=1)
        saved = json.loads(json.dumps(job.checkpoint()))
        restored = SIQSJob.from_checkpoint(saved, budget=allowance())
        self.assertEqual(restored.config, cfg)
        self.assertGreater(
            restored.budget.used, saved["resources"]["work_used"]
        )
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                saved,
                budget=allowance(),
                config=replace(
                    cfg,
                    collector=replace(cfg.collector, score_policy="powers"),
                ),
            )
        left, right = job.run(), restored.run()
        self.assertEqual(
            (left.divisor, left.cofactor), (right.divisor, right.cofactor)
        )
        self.assertEqual(left.divisor * left.cofactor, job.n)

    def test_readiness_accumulates_both_clocks_and_requires_complete_prefix(
        self,
    ):
        rules = json.loads(CONTROL.read_text())["warm"]
        rows = [dict(seconds=1.1, cpu_seconds=1.05, reason="paused")] * 3
        self.assertTrue(ready(rows, rules))
        self.assertFalse(ready(rows[:2], rules))
        self.assertFalse(
            ready(
                rows
                + [
                    dict(seconds=1.1, cpu_seconds=1.05, reason="stalled_yield")
                ],
                rules,
            )
        )


if __name__ == "__main__":
    unittest.main()
