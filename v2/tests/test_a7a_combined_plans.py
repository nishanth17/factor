"""Amended stage membership, ordering and finite confirmation allowances."""

import unittest
from collections import Counter
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.qs.a7a.combined import (
    _capture,
    corpus_identity,
    fixtures,
)
from v2.benchmarks.qs.a7a.combined_plans import (
    LIMITS,
    ablation_cells,
    stage1a_cells,
    stage1b_cells,
)


class CombinedPlanTests(unittest.TestCase):
    def test_broad_panel_preserves_old_identity_and_binds_both_sources(self):
        import hashlib

        from v2.benchmarks.qs.a7a.protocol import INPUTS

        self.assertEqual(
            corpus_identity(),
            hashlib.sha256(
                (INPUTS / "combined_corpus.json").read_bytes()
            ).hexdigest(),
        )
        self.assertNotEqual(corpus_identity("all"), corpus_identity())
        self.assertEqual(len(fixtures()), 5)
        broad = fixtures("all")
        for digits in (30, 40, 50, 55, 60, 62, 65, 70, 75, 80, 100):
            self.assertTrue(any(r["digits"] == digits for r in broad.values()))
        for name, row in fixtures().items():
            self.assertEqual(broad[name], row)
        with self.assertRaises(ValueError):
            corpus_identity("other")

    def test_cpu_only_capture_never_requests_live_heaps(self):
        events = {
            pid: [
                {"kind": "ready", "role": role},
                {"kind": "cpu_started"},
                {"kind": "cpu_stopped"},
            ]
            for pid, role in ((101, "parent"), (102, "worker"))
        }
        state = dict(
            targets=[101, 102],
            needs_worker=True,
            started={101, 102},
            stopped={101, 102},
            heaps={},
        )
        with (
            patch(
                "v2.benchmarks.qs.a7a.combined._events", return_value=events
            ),
            patch("v2.benchmarks.qs.a7a.combined.os.kill") as send,
        ):
            _capture(
                SimpleNamespace(pid=101),
                Path("unused"),
                5,
                [1, 3, 4],
                state,
                live_heaps=False,
            )
        send.assert_not_called()
        self.assertEqual(state["heaps"], {})

    def test_all60_cells_precede_diagnostics_and_retain_c6(self):
        rows = stage1a_cells()
        clean, profiles = rows[:36], rows[36:]
        counts = Counter((r["fixture"], r["arm"]) for r in clean)
        self.assertEqual(len(counts), 12)
        self.assertEqual(set(counts.values()), {3})
        self.assertEqual(
            Counter((r["fixture"], r["arm"]) for r in profiles),
            Counter({key: 1 for key in counts}),
        )
        self.assertTrue(all(r.get("mode") != "profile" for r in clean))
        self.assertTrue(all(r["mode"] == "profile" for r in profiles))
        self.assertTrue(all(r["tier_a"] for r in rows))
        self.assertEqual(LIMITS["60"], (7200, 72000, 8 * 2**30))

    def test70_has_no_c6_or_unconditional_profiles_or_rollover(self):
        for candidate, count in (("C4", 6), (None, 4)):
            rows = stage1b_cells(candidate)
            self.assertEqual(len(rows), count)
            self.assertTrue(all(r["arm"] != "C6" for r in rows))
            self.assertTrue(all(r.get("mode") != "profile" for r in rows))
            self.assertEqual(
                set(Counter((r["fixture"], r["arm"]) for r in rows).values()),
                {1},
            )
        with self.assertRaises(ValueError):
            stage1b_cells("C6")
        self.assertEqual(LIMITS["70"], (9000, 90000, 8 * 2**30))

    def test_ablations_keep_fixed_order_and_common_grants(self):
        rows = ablation_cells()
        self.assertEqual(
            [r["variant"] for r in rows],
            [
                "all",
                "without_L",
                "without_R",
                "without_S",
                "without_S",
                "without_R",
                "without_L",
                "all",
            ],
        )
        for row in rows:
            self.assertFalse(row["features"].get("arrival_order", False))
            self.assertEqual(row["polynomials_per_task"], 64)
            self.assertEqual(row["memory"], 20 * 2**30)
        self.assertLess(sum(r["wall"] for r in rows), LIMITS["ablation"][0])


if __name__ == "__main__":
    unittest.main()
