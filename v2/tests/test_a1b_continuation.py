"""Preflight guards for the separately identified A1b continuation."""

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from v2.benchmarks.qs import a1b
from v2.benchmarks.qs import a1b_continuation as study
from v2.benchmarks.qs.c10a04.accounting import Ledger


class ContinuationTests(unittest.TestCase):
    def test_endpoint_reservation_and_generator_regression(self):
        self.assertTrue(study.can_warm(7, 2, 10))
        self.assertFalse(study.can_warm(8, 2, 10))
        self.assertTrue(study.can_warm(8, 1, 10))
        self.assertAlmostEqual(study.geometric_mean(v for v in (1, 4, 16)), 4)

    def test_helper_sort_and_graph_refuse_before_unbounded_allocation(self):
        self.assertEqual(study.bounded_sorted([2, 1]), [1, 2])
        with self.assertRaises(MemoryError):
            study.bounded_sorted(range(8193))
        with self.assertRaises(MemoryError):
            study.helper_guard([], {"oversized": "x" * 2**20})
        self.assertGreater(study.deep_bytes({"k": [1, 2]}), 1024)

    def test_prospective_scope_and_accounts_fit_without_refunds(self):
        control = json.loads(study.CONTROL.read_text())
        parent = json.loads(a1b.CONTROL.read_text())
        self.assertEqual(
            control["schedule"],
            [c for c in parent["schedule"] if c["digits"] == 55],
        )
        self.assertEqual(
            [c["cell"] for c in control["schedule"]], list(range(24, 30))
        )
        limits = control["limits"]
        self.assertEqual(
            sum(
                limits[k]["wall"]
                for k in limits
                if k not in ("total", "disk_bytes")
            ),
            230,
        )
        self.assertEqual(1240 + limits["total"]["wall"], 1470)
        self.assertEqual(limits["measured"]["wall"], 6 * 25)
        self.assertEqual(limits["warm"]["wall"], 10 * 5)
        self.assertEqual(limits["legacy_analysis"]["wall"], 10)

    def test_actual_literal_warm_prefix_outcomes_before_measurement(self):
        parent = json.loads(a1b.CONTROL.read_text())
        protocol = json.loads(study.CONTROL.read_text())
        fixture = next(f for f in parent["fixtures"] if f["digits"] == 55)
        for arm in ("control", "fixed"):
            cfg = a1b._decode_siqs_config(parent["configurations"]["55"][arm])
            row = study.call_prefix(
                protocol,
                SimpleNamespace(data=[]),
                fixture,
                cfg,
                1056,
                4,
                10**10,
                retained=[parent],
            )
            self.assertEqual(row["reason"], "paused")
            self.assertEqual(row["physical"]["blocks"], 1056)
            self.assertEqual(row["cofactor"], fixture["n"])
            self.assertLessEqual(row["helper_owned_bytes"], 8 * 2**20)

    def test_recovery_censors_lost_cell_and_never_invokes_it_again(self):
        limits = json.loads(study.CONTROL.read_text())["limits"]
        with tempfile.TemporaryDirectory() as directory:
            ledger = Ledger(Path(directory), "test", limits)
            identity = dict(
                sample=24,
                fixture="test55",
                n=91,
                digits=55,
                seed=7,
                arm="control",
                cap=20,
            )
            with self.assertRaises(KeyboardInterrupt):
                ledger.call(
                    "measured",
                    "cell-24",
                    25,
                    2 * 10**10,
                    identity,
                    lambda: (_ for _ in ()).throw(KeyboardInterrupt()),
                )
            recovered = Ledger(Path(directory), "test", limits)
            result = recovered.call(
                "measured",
                "cell-24",
                25,
                2 * 10**10,
                identity,
                lambda: self.fail("lost cell was selectively repeated"),
            )
            self.assertEqual(result["reason"], "interrupted")
            self.assertEqual(recovered.used("measured"), 25)
            self.assertEqual(recovered.used("measured", "work"), 2 * 10**10)


if __name__ == "__main__":
    unittest.main()
