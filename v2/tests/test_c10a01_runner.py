"""Recorder passivity, direct coordinate replay and frozen accounting."""

import json
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.benchmarks.qs.c10a01 import runner
from v2.benchmarks.qs.c10a01.configuration import configuration
from v2.benchmarks.qs.c10a01.inputs import CORPUS, validate


class CalibrationContractTests(unittest.TestCase):
    def test_exact_budget_sums_and_exposed_input_proofs(self):
        plan = json.loads(runner.BOOTSTRAP.read_text())
        phases = [
            v
            for k, v in plan["limits"].items()
            if k not in ("total", "disk_bytes")
        ]
        self.assertEqual(
            {k: sum(p[k] for p in phases) for k in ("wall", "cpu", "work")},
            plan["limits"]["total"],
        )
        self.assertEqual(3600 + 9200 + 7600 + 14800, 35200)
        validate(json.loads(CORPUS.read_text()))

    def test_capacity_fixture_kernels_and_reservation(self):
        from v2.benchmarks.qs.c10a01.diagnostics import matrix_check
        from v2.execution.budget import Budget

        for kind in ("cycle", "sparse35"):
            budget = Budget(work_limit=10**9, seconds=None, cpu_seconds=None)
            result = matrix_check(
                dict(rows=128, columns=37, kind=kind), budget
            )
            self.assertEqual(result["reason"], "complete")
            self.assertTrue(result["refused_below_reserve"])
            self.assertGreaterEqual(result["elimination"]["dependencies"], 64)

    def test_pilot_accounts_and_censoring_rules_without_execution(self):
        from v2.benchmarks.qs.c10a01 import pilots

        plan = json.loads(pilots.PLAN.read_text())
        self.assertTrue(pilots.verify_plan(plan))
        self.assertEqual(
            {
                c["digits"]
                for stage in plan["stages"].values()
                for c in stage["calls"]
            },
            {60, 62, 65, 70},
        )
        rows = [
            dict(
                call=dict(digits=70, arm=arm),
                result=dict(complete=done, seconds=100, cpu_seconds=99),
            )
            for arm, done in (("control", False), ("parity2", True))
        ]
        result = pilots.analyze(rows)
        self.assertNotIn("descriptive_wall_ratio", result["pairs"][0])
        self.assertFalse(result["adoption"])

    def test_lower_band_explicit_control_and_promoted_default_boundaries(self):
        from v2.execution.budget import Budget
        from v2.qs import SIQSJob, default_siqs_config
        from v2.tests.qs_default_controls import previous_siqs_config

        files = ("c10_fresh_corpus.json", "c10a05_plus_feasibility.json")
        for digits, name in zip((50, 55), files):
            path = CORPUS.with_name(name)
            fixture = next(
                r
                for r in json.loads(path.read_text())["fixtures"]
                if r["digits"] == digits
            )
            n = fixture["n"]
            default = previous_siqs_config(n)
            self.assertEqual(default.polynomial_policy, "standard")
            config = replace(
                default,
                polynomial_policy="parity2",
                multiplier_policy="parity2",
            )
            job = SIQSJob(
                n,
                seed=7,
                config=config,
                budget=Budget(work_limit=10**10, seconds=10, cpu_seconds=10),
            )
            result = job.run(max_blocks=1)
            self.assertEqual(result.reason, "paused")
            polynomial = job.engine.collector.polynomial
            for x in (-65536, -1, 0, 65536):
                self.assertEqual(
                    (polynomial.a * x + polynomial.b) ** 2
                    - polynomial.n_prime,
                    polynomial.a * polynomial.value(x),
                )
        for digits in (49, 50, 51, 52, 55, 56, 60, 61):
            for n in (10 ** (digits - 1), 10**digits - 1):
                cfg = default_siqs_config(n)
                self.assertEqual(cfg.polynomial_policy, "parity2")
                self.assertEqual(cfg.multiplier_policy, "parity2_wide")

    def test_observer_and_direct_coordinate_replay(self):
        for arm in ("control", "parity2"):
            cfg = configuration(60, arm, multiplier=3)
            cfg = replace(
                cfg,
                base_count=10,
                base_bound=97,
                factor_count=3,
                pool_size=4,
                half_width=32,
                max_half_width=32,
                batch_width=16,
                polynomials_per_family=4,
                collector=replace(cfg.collector, block_width=16),
            )
            with patch.object(runner, "configuration", return_value=cfg):
                off = runner.run_one(4001 * 5003, 60, arm, 20, 10**9, 10)
                first = runner.describe(runner.LAST_JOB)
                on = runner.run_one(
                    4001 * 5003, 60, arm, 20, 10**9, 10, observer=True
                )
                second = runner.describe(runner.LAST_JOB)
            self.assertEqual(off["state_digest"], on["state_digest"])
            self.assertEqual(first["records"], second["records"])
            self.assertTrue(on["observation_records"])
            self.assertEqual(
                sum(r["hi"] - r["lo"] for r in first["records"]),
                off["positions"],
            )


if __name__ == "__main__":
    unittest.main()
