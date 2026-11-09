"""Independent packed predicates and first-uncommitted scan positions."""

import unittest
from array import array
from unittest.mock import patch

from v2.qs import Polynomial, build_factor_base
from v2.qs.collector_kernels import word_candidates
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_a1a_density import worker
from v2.tests.test_qs import unlimited_budget


class WordScanTests(unittest.TestCase):
    def test_every_byte_threshold_and_wide_lane_boundaries(self):
        scores = bytearray(range(256))
        for threshold in range(258):
            self.assertEqual(
                list(word_candidates(scores, threshold, len(scores))),
                [i for i, score in enumerate(scores) if score >= threshold],
            )

        for threshold in (1, 127, 128, 255, 256, 32768, 65536, 2**32 + 1):
            values = [0, threshold - 1, threshold, threshold + 1, 2**40]
            # Alternating low/high lanes catch cross-lane borrow and tail bugs.
            scores = (values + values[::-1]) * 14
            for width in (1, 63, 64, 65, 129, len(scores)):
                expected = [
                    i
                    for i, value in enumerate(scores[:width])
                    if value >= threshold
                ]
                self.assertEqual(
                    list(word_candidates(scores, threshold, width)), expected
                )
                if max(scores) < 2**64:
                    self.assertEqual(
                        list(
                            word_candidates(
                                array("Q", scores), threshold, width
                            )
                        ),
                        expected,
                    )

    def scan_worker(self, selected):
        base = build_factor_base(10403, bound=100).factor_base
        collector = worker(
            Polynomial(10403, 1, 49, 8), base, "scan_fused", block_width=130
        )
        visited = []

        def sieve(lo, hi, stats):
            collector._scores[: hi - lo] = [
                int(x in selected) for x in range(lo, hi)
            ]
            stats["blocks"] += 1
            return 1

        def divide(position, offset, stats):
            collector.budget.consume(1)
            visited.append(position)
            return None, None

        return collector, visited, sieve, divide

    def test_word_boundary_refusal_does_not_replay_committed_candidate(self):
        collector, visited, sieve, divide = self.scan_worker({63, 128})
        collector.budget = Budget(
            work_limit=10**9, cancelled=lambda: bool(visited)
        )
        with (
            patch.object(collector, "_sieve", sieve),
            patch.object(collector, "_divide", divide),
        ):
            stopped = collector.collect(0, 130)
            self.assertEqual(stopped.reason, "cancelled")
            self.assertEqual(stopped.next_position, 64)
            collector.budget = unlimited_budget()
            resumed = collector.collect(stopped.next_position, 130)
        self.assertEqual(resumed.next_position, 130)
        self.assertEqual(visited, [63, 128])

    def test_reject_span_reserves_each_position_and_polls_within_64(self):
        for work, cursor in ((131 + 17, 17), (131 + 64, 64)):
            collector, _, sieve, divide = self.scan_worker(set())
            collector.budget = Budget(work_limit=work)
            with (
                patch.object(collector, "_sieve", sieve),
                patch.object(collector, "_divide", divide),
            ):
                result = collector.collect(0, 130)
            self.assertEqual(result.reason, "work_limit")
            self.assertEqual(result.next_position, cursor)
            self.assertEqual(result.stats["scanned"], cursor)
            self.assertEqual(collector.budget.used, work)
