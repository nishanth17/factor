"""Independent incidence and domain checks for the C9 study harness."""

import unittest
from dataclasses import asdict
from types import SimpleNamespace

from v2.benchmarks.qs.b1.b1_calibration import common_config
from v2.benchmarks.qs.c9.c9_study import configuration, graph_diagnostics


class C9StudyTests(unittest.TestCase):
    def test_alignment_preserves_dlp_candidate_and_product_domains(self):
        base = common_config(30000, 65536, 6, rows=8192, partials=8192)
        control = asdict(configuration(base, "dlp"))
        aligned = asdict(configuration(base, "aligned"))
        self.assertEqual(control["collector"]["residual_bound"], 900000000)
        self.assertEqual(aligned["collector"]["residual_bound"], 3000000)
        aligned["collector"]["residual_bound"] = 900000000
        self.assertEqual(control, aligned)

    def test_multigraph_core_and_disconnected_cycles(self):
        # A triangle, a parallel pair, a loop and a two-edge hanging path.
        pairs = [
            (3, 5),
            (5, 7),
            (7, 3),
            (11, 13),
            (11, 13),
            (17, 17),
            (7, 19),
            (19, 23),
        ]
        collector = SimpleNamespace(
            _atoms={
                i: SimpleNamespace(large_primes=pair, residual=1)
                for i, pair in enumerate(pairs)
            },
            config=SimpleNamespace(large_prime_bound=100),
        )
        result = graph_diagnostics(collector)
        # Enumerate every even-incidence subset independently of peeling.
        kernel = []
        for mask in range(1 << len(pairs)):
            odd = set()
            for i, pair in enumerate(pairs):
                if mask & (1 << i):
                    for vertex in pair:
                        odd.symmetric_difference_update((vertex,))
            if not odd:
                kernel.append(mask)
        self.assertEqual(len(kernel), 1 << result["retained_cycle_rank"])
        self.assertEqual(result["retained_cycle_rank"], 3)
        self.assertEqual(result["core_edges"], 6)
        self.assertEqual(result["components"], 3)
        self.assertEqual(result["core_vertices"], 6)


class C9ConfirmationTests(unittest.TestCase):
    def test_extension_preserves_negative_and_censored_pairs(self):
        from v2.benchmarks.qs.c9.c9_confirmation import unstable

        rows = []
        for sample, candidate in enumerate((8, 12, 10)):
            for arm, seconds in (
                ("slp", 10),
                ("dlp", 10),
                ("challenger", candidate),
            ):
                rows.append(
                    dict(
                        id="case",
                        sample=sample,
                        seed=sample,
                        arm=arm,
                        seconds=seconds,
                        capped_seconds=seconds,
                    )
                )
        self.assertEqual(len(unstable(rows)), 2)
        for row in rows:
            row["capped_seconds"] = 240
        self.assertEqual(unstable(rows), [])

    def test_training_selection_charges_fast_refusals_at_full_cap(self):
        from v2.benchmarks.qs.c9.c9_selection import ranked

        rows = [
            dict(arm=arm, complete=complete, capped_seconds=cost)
            for arm, complete, cost in (
                ("refused", False, 240),
                ("refused", False, 240),
                ("completed", True, 120),
                ("completed", True, 100),
                ("mixed", True, 1),
                ("mixed", False, 240),
            )
        ]
        result = ranked(rows, ("refused", "completed", "mixed"))
        self.assertEqual(
            [r["arm"] for r in result], ["completed", "mixed", "refused"]
        )
        with self.assertRaises(ValueError):
            ranked(rows[:-1], ("mixed",))

    def test_input_heterogeneity_alone_is_not_temporal_instability(self):
        from v2.benchmarks.qs.c9.c9_confirmation import unstable

        rows = [
            dict(
                id=identity,
                sample=0,
                seed=7,
                arm=arm,
                seconds=cost,
                capped_seconds=cost,
            )
            for identity, candidate in (("easy", 8), ("hard", 12))
            for arm, cost in (
                ("slp", 10),
                ("dlp", 10),
                ("challenger", candidate),
            )
        ]
        self.assertEqual(unstable(rows), [])

    def test_expense_refusal_and_failed_call_remain_charged(self):
        from pathlib import Path
        from tempfile import TemporaryDirectory
        from unittest.mock import patch

        from v2.benchmarks.qs.c9.c9_confirmation import Expenses

        with TemporaryDirectory() as directory:
            selection = {
                "spent_before_confirmation": {"wall": 14398, "cpu": 14398}
            }
            expense = Expenses(
                Path(directory), selection, {"generation_seconds": 0}
            )
            with self.assertRaises(RuntimeError):
                expense.call("over cap", 3, lambda: self.fail("ran"))
            self.assertEqual(expense.data["calls"], [])
            with (
                patch("time.monotonic", side_effect=[10, 11]),
                patch("time.process_time", side_effect=[20, 21]),
            ):
                with self.assertRaises(ZeroDivisionError):
                    expense.call("failed", 2, lambda: 1 / 0)
            self.assertEqual(expense.data["calls"][0]["status"], "failed")
            self.assertGreaterEqual(expense.used(), 14399)


if __name__ == "__main__":
    unittest.main()
