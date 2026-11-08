"""Grouped marking retains exact buffers, poll prefixes and refusal charges."""

import unittest
from collections import Counter
from contextlib import nullcontext
from dataclasses import replace
from unittest.mock import patch

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.qs.current_parallel import CurrentParallelSIQSJob, QSWorkerConfig
from v2.qs.hit_words import WordHits
from v2.qs.packet_lease import PacketLeaseBudget
from v2.qs.parallel import CollectionPool
from v2.tests.test_current_qs_workers import (
    allowance,
    configuration,
    settings,
)
from v2.tests.test_qs_hit_words import make_worker


class GroupedLeaseTests(unittest.TestCase):
    def test_pair_never_crosses_a_poll_or_work_boundary(self):
        for remaining in (0, 1, 2, 63):
            for available in (0, 6, 7, 12):
                view = PacketLeaseBudget(
                    Budget(
                        work_limit=available, seconds=None, cpu_seconds=None
                    )
                )
                view._actions_left = remaining
                before = (view.used, view._actions_left, view.marking_pairs)
                grouped = view._precharge_marking_pair(3, 4)
                self.assertEqual(grouped, remaining >= 2 and available >= 7)
                expected = (7, remaining - 2, 1) if grouped else before
                self.assertEqual(
                    (view.used, view._actions_left, view.marking_pairs),
                    expected,
                )

    def sieve(self, grouped, kernels, work_limit=10**8, cancel=False):
        collector = make_worker(64, kernels)
        collector._set_plan_interval(-79, 82)
        trace = []
        stopped = [False]
        raw = Budget(
            work_limit=work_limit,
            seconds=None,
            cpu_seconds=None,
            cancelled=lambda: stopped[0],
        )
        view = PacketLeaseBudget(raw, poll=lambda b: trace.append(b.used))
        collector.budget = view
        collector._group_lease = view if grouped else None
        original_add = WordHits.add

        def mark_and_stop(*args, **kwargs):
            result = original_add(*args, **kwargs)
            stopped[0] = True
            return result

        stats = Counter()
        reason = None
        try:
            # Raise cancellation only after marking has started, so both
            # paths must finish the same already-polled action prefix.
            context = (
                patch.object(WordHits, "add", mark_and_stop)
                if cancel
                else nullcontext()
            )
            with context:
                collector._sieve(-79, -15, stats)
        except BudgetExhaustedError:
            reason = view.reason
        view.reconcile()
        return dict(
            scores=list(collector._scores),
            hits=collector._hits[:],
            work=view.used,
            remaining=view._remaining,
            polls=trace,
            stats=dict(stats),
            reason=reason,
        ), view.marking_pairs

    def test_exact_buffers_and_poll_prefix_with_presieve_and_fallbacks(self):
        for kernels in (
            ("unroll4",),
            ("presieve", "unroll4", "parity_two"),
            ("presieve_powers", "unroll4"),
            (),
        ):
            expected, _ = self.sieve(False, kernels)
            actual, pairs = self.sieve(True, kernels)
            self.assertEqual(expected, actual)
            if "unroll4" in kernels and "presieve_powers" not in kernels:
                self.assertGreater(pairs, 0)
            else:
                self.assertEqual(pairs, 0)

    def test_near_limit_and_cancellation_keep_the_original_paid_prefix(self):
        kernels = ("presieve", "unroll4", "parity_two")
        complete, _ = self.sieve(False, kernels)
        limits = {1, 64, 512, complete["work"] - 1, complete["work"]}
        # Poll-prefix neighborhoods exercise scalar fallback as well as
        # previously successful pairs before a later refusal.
        for value in complete["polls"][:12]:
            limits.update((max(0, value - 1), value, value + 1))
        for limit in sorted(limits):
            self.assertEqual(
                self.sieve(False, kernels, limit)[0],
                self.sieve(True, kernels, limit)[0],
            )
        self.assertEqual(
            self.sieve(False, kernels, cancel=True)[0],
            self.sieve(True, kernels, cancel=True)[0],
        )

    def test_worker_flag_validates_and_reports_exercised_groups(self):
        self.assertNotIn("grouped_lease_charging", QSWorkerConfig().record())
        with self.assertRaises(ValueError):
            QSWorkerConfig(grouped_lease_charging=True)
        cfg = configuration()
        cfg = replace(
            cfg,
            base_bound=2000,
            collector=replace(
                cfg.collector, marking="presieve+unroll4+parity_two"
            ),
        )
        jobs = []
        for grouped in (False, True):
            job = CurrentParallelSIQSJob(
                1000003 * 1000033,
                seed=29,
                config=cfg,
                workers=replace(
                    settings(),
                    lease_budget=True,
                    reuse_state=True,
                    grouped_lease_charging=grouped,
                ),
                budget=allowance(),
            )
            with CollectionPool() as pool:
                job.run(pool=pool, fixed_work=True)
                self.assertIsNone(
                    pool.state["_qs_context"].collector._group_lease
                )
            jobs.append(job)
        self.assertEqual(jobs[0].stats["scanned"], jobs[1].stats["scanned"])
        self.assertEqual(
            dict(jobs[0].engine.collector._atoms),
            dict(jobs[1].engine.collector._atoms),
        )
        self.assertGreater(jobs[1].stats["grouped_marking_pairs"], 0)

        restored = CurrentParallelSIQSJob.from_checkpoint(
            jobs[1].checkpoint(), budget=allowance()
        )
        self.assertTrue(restored.workers.grouped_lease_charging)
        self.assertEqual(
            restored.stats["grouped_marking_pairs"],
            jobs[1].stats["grouped_marking_pairs"],
        )
