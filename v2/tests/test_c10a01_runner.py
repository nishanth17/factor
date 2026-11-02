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
