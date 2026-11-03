"""Check diagnostic reservations, shared roots and observer passivity."""

import json
import unittest
from pathlib import Path

from v2.benchmarks.qs.c10a01 import (
    execute_v4,
    execute_v4_diagnostics,
    execute_v4_group,
    plan_v4,
)


class DiagnosticAdapterTests(unittest.TestCase):
    def test_all_declared_cells_fit_remaining_phase_reservations(self):
        plan = json.loads(plan_v4.PLAN.read_text())
        calls = execute_v4_diagnostics.schedule(plan)
        self.assertEqual(sum(len(c["rows"]) for c in calls), 22)
        self.assertEqual(len(calls), 16)
        limits = execute_v4_diagnostics.limits(plan)
        self.assertEqual(limits["total"]["wall"], 2020)
        self.assertEqual(
            sum(
                v["wall"]
                for k, v in limits.items()
                if k not in ("total", "disk_bytes")
            ),
            1530,
        )
        self.assertLessEqual(
            sum(
                v["work"]
                for k, v in limits.items()
                if k not in ("total", "disk_bytes")
            ),
            limits["total"]["work"],
        )
        self.assertEqual(
            [
                r["configuration"]["assignment_policy"]
                for c in calls
                for r in c["rows"]
                if r["phase"] == "population" and r["digits"] == 70
            ],
            ["paired_flyer", "windowed_pair"],
        )

    def test_root_policies_share_exact_family_endpoints(self):
        corpus = json.loads(
            (
                Path(__file__).resolve().parents[1]
                / "benchmarks/inputs/corpora/a1a_integration_fresh.json"
            ).read_text()
        )
        sample = corpus["fixtures"][0]
        fixture = dict(n=sample["n"], factors=sample["factors"])
        rows = [
            r
            for r in json.loads(plan_v4.PLAN.read_text())["calibration"][
                "calls"
            ]
            if r["phase"] == "roots" and r["digits"] == 62
        ]
        shared, _, _ = execute_v4.instantiate(
            rows[0],
            fixture,
            execute_v4.runtime_modules(rows[0]),
            seconds=45,
            work=rows[0]["work"],
        )
        shared._setup()
        self.assertIsNone(shared.divisor)
        assignments = tuple(shared.assignments[i] for i in range(8))
        receipts = [
            execute_v4_group.sample_root(
                row, fixture, shared.base, shared.multiplier, assignments
            )
            for row in rows
        ]
        self.assertEqual({r["endpoint_count"] for r in receipts}, {8})
        self.assertEqual(len({r["fingerprint"] for r in receipts}), 1)

    def test_short_observation_preserves_semantic_work(self):
        corpus = json.loads(
            (
                Path(__file__).resolve().parents[1]
                / "benchmarks/inputs/corpora/a1a_integration_fresh.json"
            ).read_text()
        )
        sample = corpus["fixtures"][0]
        fixture = dict(n=sample["n"], factors=sample["factors"])
        row = next(
            r
            for r in json.loads(plan_v4.PLAN.read_text())["calibration"][
                "calls"
            ]
            if r["phase"] == "observer" and r["digits"] == 60
        )
        modules = execute_v4.runtime_modules(row)
        off = execute_v4.sieve(
            row, fixture, modules, seconds=15, work=row["work"], blocks=4
        )
        on = execute_v4.sieve(
            row,
            fixture,
            modules,
            seconds=15,
            work=row["work"],
            blocks=4,
            observer=True,
        )
        self.assertEqual(off["reason"], on["reason"])
        self.assertEqual(off["positions"], on["positions"])
        self.assertEqual(off["work"], on["work"])
        self.assertEqual(off["semantic_digest"], on["semantic_digest"])
        self.assertEqual(on["observer_records"], 4)
