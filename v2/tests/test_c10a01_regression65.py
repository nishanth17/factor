"""Diagnostic continuation, semantic-work and reservation checks."""

import copy
import json
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.qs.c10a01 import regression65


class RegressionDiagnosticTests(unittest.TestCase):
    def test_observation_and_segmentation_preserve_work_and_atoms(self):
        def make_job(*args, **kwargs):
            budget = SimpleNamespace(used=0)
            engine = SimpleNamespace(
                stats={"collected_positions": 0},
                collector=SimpleNamespace(_atoms={}),
            )
            job = SimpleNamespace(
                engine=engine,
                multiplier=3,
                family_index=0,
                stats={"polynomials": 0, "stage_seconds": {"roots": 0}},
            )

            def run(*, max_blocks):
                self.assertEqual(max_blocks % 33, 0)
                count = max_blocks // 33
                for _ in range(count):
                    index = job.stats["polynomials"]
                    engine.collector._atoms[(index, "atom")] = index
                    job.stats["polynomials"] += 1
                job.family_index = job.stats["polynomials"] // 64
                engine.stats["collected_positions"] += count * 131073
                budget.used += count * 77
                return SimpleNamespace(reason="paused", stats=job.stats)

            job.run = run
            return job, budget, None

        base = dict(polynomials=1024, seconds=90, work=10**20)
        modules = (None, SimpleNamespace(BudgetExhaustedError=Exception), None)
        with (
            patch.object(regression65.execute_v4, "instantiate", make_job),
            patch.object(
                regression65.execute_v4, "validate_result", return_value={}
            ),
            patch.object(regression65, "snapshot", return_value={"gc": 0}),
        ):
            direct = regression65.prefix(
                dict(base, observed=False, chunk_polynomials=1024), {}, modules
            )
            observed = regression65.prefix(
                dict(base, observed=True, chunk_polynomials=128), {}, modules
            )

        self.assertEqual(direct["digest"], observed["digest"])
        self.assertEqual(direct["work"], 1024 * 77)
        self.assertEqual(observed["work"], direct["work"])
        self.assertEqual(observed["positions"], 1024 * 131073)
        self.assertEqual(len(observed["chunks"]), 8)
        observed["stats"]["stage_seconds"]["roots"] = 999
        self.assertEqual(observed["chunks"][0]["stages"]["roots"], 0)

    def test_frozen_phase_totals_and_literal_policy_match(self):
        plan = json.loads(regression65.PLAN.read_text())
        for field in ("wall", "cpu", "work"):
            expected = sum(
                plan["limits"][phase][field]
                for phase in (
                    "screen",
                    "long",
                    "warmup",
                    "analysis",
                    "orchestration",
                )
            )
            self.assertEqual(plan["limits"]["total"][field], expected)
        self.assertEqual(plan["limits"]["total"]["wall"], 2190)
        for arm in ("control", "nominee"):
            rows = [r for r in plan["rows"] if r["arm"] == arm]
            expected = copy.deepcopy(rows[0]["configuration"])
            self.assertTrue(all(r["configuration"] == expected for r in rows))

    def test_retired_runtime_guard_rejects_promoted_defaults(self):
        with self.assertRaisesRegex(RuntimeError, "source drift"):
            regression65.check()


if __name__ == "__main__":
    unittest.main()
