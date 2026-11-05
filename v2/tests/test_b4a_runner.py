"""Finite reservations, censoring and raw-output validation."""

import json
import tempfile
import unittest
from pathlib import Path

from v2.benchmarks.ecm.b4a.analysis import interval
from v2.benchmarks.ecm.b4a.plan import (
    downstream_reservation,
    stage_a,
    uniform_panel,
)
from v2.benchmarks.ecm.b4a.runner import check_cell, reserve, validate_child


class B4aRunnerTests(unittest.TestCase):
    def test_early_incomplete_search_cannot_create_a_fast_tail(self):
        clusters = {
            n: dict(control=[(0.01, False)] * 9, nominee=[(2, True)] * 9)
            for n in range(12)
        }
        result = interval(clusters, grant=5, draws=100)
        self.assertIsNone(result["estimates"]["control"]["median"]["seconds"])
        self.assertIsNone(result["intervals"]["p90"])
        self.assertEqual(result["independent_moduli"], 12)
        clusters[0]["control"] = [(0.01, False)]
        with self.assertRaises(ValueError):
            interval(clusters, grant=5, draws=100)

    def test_literal_stage_a_fits_and_has_no_thread_impersonation(self):
        cells = stage_a()
        self.assertEqual(len(cells), 93)
        self.assertLessEqual(sum(c["wall"] for c in cells), 1800)
        self.assertLessEqual(sum(c["cpu"] for c in cells), 7200)
        schedules = {
            c["schedule"] for c in cells if c["id"].startswith("schedule-")
        }
        self.assertEqual(schedules, {"static", "dynamic", "guided"})
        self.assertTrue(any(c["arm"] == "Native-T4" for c in cells))
        with self.assertRaises(ValueError):
            check_cell(dict(cells[0], schedule="weighted"))

    def test_reservations_survive_failure_and_reject_changed_freeze(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.json"
            freeze = dict(identity="test", wall=50, cpu=200)
            cell = stage_a()[0]
            reserve(path, freeze, cell)
            self.assertEqual(
                json.loads(path.read_text())["wall"], cell["wall"]
            )
            with self.assertRaises(ValueError):
                reserve(path, freeze, cell)
            with self.assertRaises(ValueError):
                reserve(
                    path,
                    dict(freeze, identity="changed"),
                    dict(cell, id="another"),
                )

    def test_calibration_cannot_enlarge_an_uninformative_campaign(self):
        with self.assertRaises(ValueError):
            downstream_reservation("B", [], cells=1, curves=8)
        clean = [dict(status="complete", wall=1, work=10**7) for _ in range(3)]
        funded = downstream_reservation("B", clean, cells=9, curves=8)
        self.assertGreaterEqual(funded["work"], 32_000_000 * 8)
        with self.assertRaises(ValueError):
            downstream_reservation("B", clean, cells=1000, curves=8)

    def test_uniform_panels_are_unconditioned_exact_size_and_disjoint(self):
        training, confirmation = (
            uniform_panel(2026110417, 36),
            uniform_panel(2026110419, 36),
        )
        self.assertFalse(
            {r["n"] for r in training} & {r["n"] for r in confirmation}
        )
        for rows in (training, confirmation):
            self.assertTrue(any(r["n"] % 2 == 0 for r in rows))
            self.assertTrue(all(len(str(r["n"])) == r["digits"] for r in rows))

    def test_incomplete_outputs_reconstruct_and_are_not_fast_completions(self):
        cell = dict(stage_a()[0], endpoint="complete")
        data = dict(
            n=15,
            arm=cell["arm"],
            mode="cold",
            endpoint="complete",
            complete=False,
            total_work=5,
            factors=[[3, 1, "proven_prime"]],
            remaining=[5],
            reconstruct=15,
        )
        validate_child(data, cell, dict(n=15))
        with self.assertRaises(ValueError):
            validate_child(dict(data, complete=True), cell, dict(n=15))
        with self.assertRaises(ValueError):
            validate_child(dict(data, remaining=[7]), cell, dict(n=15))


if __name__ == "__main__":
    unittest.main()
