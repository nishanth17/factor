"""Pinned family placement preserves exact rows, leases and pool retirement."""

import multiprocessing
import unittest
from dataclasses import replace
from unittest.mock import Mock, patch

from v2.qs.affinity_executor import FamilyAffinityExecutor
from v2.qs.current_parallel import CurrentParallelSIQSJob, QSWorkerConfig
from v2.qs.parallel import CollectionPool
from v2.tests.test_current_qs_workers import (
    allowance,
    crash_packet,
    prepared,
    settings,
)


def affinity_settings():
    return replace(
        settings("process", 2),
        reuse_state=True,
        lease_budget=True,
        packed_transport=True,
        verification="worker",
        family_affinity=True,
        lookahead=2,
    )


class AffinityPlacementTests(unittest.TestCase):
    def test_broadcast_and_explicit_packet_placement_visit_correct_slots(self):
        with patch("v2.qs.affinity_executor.ProcessPoolExecutor") as factory:
            factory.side_effect = [Mock() for _ in range(3)]
            executor = FamilyAffinityExecutor(
                3, mp_context=None, initializer=None, initargs=()
            )
            for _ in range(2):
                for slot in range(3):
                    executor.submit(str, slot)
            self.assertEqual(factory.call_count, 3)
            for worker in executor.executors:
                self.assertEqual(worker.submit.call_count, 2)
            executor.submit_to(2, str, "family")
            executor.executors[2].submit.assert_called_with(str, "family")
            self.assertTrue(
                all(
                    c.kwargs["max_workers"] == 1
                    for c in factory.call_args_list
                )
            )

    def test_off_identity_and_conditional_admission(self):
        self.assertNotIn("family_affinity", QSWorkerConfig().record())
        with self.assertRaises(ValueError):
            QSWorkerConfig(family_affinity=True)
        with self.assertRaises(ValueError):
            QSWorkerConfig(
                mode="serial", reuse_state=True, family_affinity=True
            )
        self.assertTrue(affinity_settings().record()["family_affinity"])

    def test_exact_packets_reuse_families_and_resume_in_a_fresh_pool(self):
        reference = prepared()
        reference.run(fixed_work=True)
        job = prepared(worker_settings=affinity_settings())

        with CollectionPool("process", 2, family_affinity=True) as pool:
            paused = job.run(pool=pool, max_assignments=2, fixed_work=True)
            self.assertEqual(paused.reason, "paused")
            self.assertEqual(len(pool._owned_processes), 2)
            checkpoint = job.checkpoint()
        restored = CurrentParallelSIQSJob.from_checkpoint(
            checkpoint, budget=allowance()
        )
        result = restored.run(fixed_work=True)

        self.assertTrue(result.stats["schedule_complete"])
        self.assertEqual(
            dict(reference.engine.collector._atoms),
            dict(restored.engine.collector._atoms),
        )
        self.assertEqual(reference.stats["scanned"], restored.stats["scanned"])
        self.assertGreater(restored.stats["state_family_reuses"], 0)
        self.assertTrue(
            all(r["status"] != "issued" for r in restored.receipts)
        )

    def test_crash_retains_paid_leases_and_retires_every_owned_slot(self):
        before = {p.pid for p in multiprocessing.active_children()}
        job = prepared(worker_settings=affinity_settings())

        with CollectionPool("process", 2, family_affinity=True) as pool:
            with patch("v2.qs.current_parallel.execute_packet", crash_packet):
                result = job.run(pool=pool, fixed_work=True)

        self.assertEqual(result.reason, "worker_failure")
        self.assertTrue(pool.closed)
        self.assertGreater(job.stats["lost_work"], 0)
        self.assertGreaterEqual(job.budget.used, job.stats["lost_work"])
        self.assertTrue(all(not p.is_alive() for p in pool._owned_processes))
        self.assertEqual(
            before, {p.pid for p in multiprocessing.active_children()}
        )
