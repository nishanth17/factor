"""Coverage and conditional-validity guards for joint experiment recipes."""

import unittest
from collections import Counter
from itertools import combinations

from v2.benchmarks.qs.a7a.cell import settings
from v2.benchmarks.qs.a7a.joint_search import (
    initial_design,
    local_design,
    pair_coverage,
)
from v2.qs.collector_kernels import components


class JointSearchTests(unittest.TestCase):
    def test_every_raw_factor_pair_is_covered_in_the_full_design(self):
        rows = initial_design()
        self.assertEqual(len(rows), 32)
        coverage = pair_coverage(rows)
        self.assertEqual(len(coverage), 253)
        self.assertTrue(
            all(
                pair["observed"] == pair["possible"]
                for pair in coverage.values()
            )
        )
        self.assertEqual(
            Counter(row["arm"] for row in rows),
            {"C4": 8, "C6": 8, "C8": 8, "C10": 8},
        )

    def test_conditional_settings_are_valid_without_hiding_the_raw_design(
        self,
    ):
        for row in initial_design():
            worker = settings(
                row["arm"],
                22 * 2**30,
                False,
                row["polynomials_per_task"],
                2048,
                row["features"],
            )
            components(row["collector_overrides"]["marking"])
            self.assertLessEqual(worker.lookahead, 4 * worker.workers)
            if not worker.lease_budget:
                self.assertFalse(worker.private_stop)
                self.assertFalse(worker.cpu_snapshot)
                self.assertEqual(worker.stop_backend, "event")
            self.assertEqual(
                worker.whole_family, row["polynomials_per_task"] == 64
            )

    def test_calls_do_not_share_mutable_configuration(self):
        first = initial_design()
        first[0]["features"]["reuse_state"] = True
        self.assertFalse(initial_design()[0]["features"]["reuse_state"])

    def test_small_family_geometry_does_not_collapse_the_packet_design(self):
        for gray in (16, 32, 64):
            rows = initial_design(gray_limit=gray)
            self.assertEqual(
                {row["polynomials_per_task"] for row in rows},
                {gray // 16, gray // 8, gray // 4, gray},
            )
            self.assertTrue(
                all(
                    pair["observed"] == pair["possible"]
                    for pair in pair_coverage(rows, gray_limit=gray).values()
                )
            )
            for row in rows:
                self.assertEqual(
                    row["features"]["whole_family"],
                    row["polynomials_per_task"] == gray,
                )
                self.assertEqual(gray % row["polynomials_per_task"], 0)
        with self.assertRaises(ValueError):
            pair_coverage(initial_design(gray_limit=16))
        with self.assertRaises(ValueError):
            initial_design(gray_limit=8)

    def test_local_refinement_covers_all_pairs_without_mutating_anchor(self):
        anchor = next(
            row
            for row in initial_design()
            if row["features"]["reuse_state"]
            and row["features"]["stop_backend"] == "semaphore"
        )
        before = repr(anchor)
        rows = local_design(anchor)
        self.assertEqual(len(rows), 8)
        self.assertEqual(repr(anchor), before)
        names = tuple(rows[0]["design_levels"])
        for left, right in combinations(names, 2):
            frequencies = Counter(
                (row["design_levels"][left], row["design_levels"][right])
                for row in rows
            )
            self.assertEqual(sorted(frequencies.values()), [2, 2, 2, 2])
        for row in rows:
            settings(
                row["arm"],
                22 * 2**30,
                False,
                row["polynomials_per_task"],
                2048,
                row["features"],
            )
            pieces = components(row["collector_overrides"]["marking"])
            self.assertNotEqual(
                "presieve" in pieces, "presieve_powers" in pieces
            )
        with self.assertRaises(ValueError):
            local_design(initial_design()[0])


if __name__ == "__main__":
    unittest.main()
