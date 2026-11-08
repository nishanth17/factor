"""Disjoint-work accounting checks without spawning benchmark workers."""

import unittest

from v2.benchmarks.qs.a7a.independent_packets import partition, summarize


class IndependentPacketTests(unittest.TestCase):
    def test_every_concurrency_preserves_the_same_disjoint_work(self):
        tasks = list(range(16))
        for workers in (1, 4, 8, 10):
            groups = partition(tasks, workers)
            self.assertEqual(sorted(t for g in groups for t in g), tasks)
            self.assertEqual(max(map(len, groups)), (16 + workers - 1) // workers)
            self.assertLessEqual(max(map(len, groups)) - min(map(len, groups)), 1)

    def test_summary_rejects_missing_and_duplicate_work(self):
        def stream(index):
            return dict(
                cpu=1,
                wall=2,
                reports=[
                    dict(
                        index=index,
                        sha256=str(index),
                        scanned=100,
                        work=10,
                        family_reuses=0,
                    )
                ],
            )

        expected = summarize([stream(0), stream(1)], 2)
        self.assertEqual(expected["scanned"], 200)
        self.assertEqual(expected["worker_cpu"], 2)
        self.assertEqual(
            summarize([stream(1), stream(0)], 2)["signature"],
            expected["signature"],
        )
        for rows in ([stream(0)], [stream(0), stream(0)]):
            with self.assertRaises(ValueError):
                summarize(rows, 2)
