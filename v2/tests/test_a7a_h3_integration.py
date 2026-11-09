"""Observed resource ownership for the current-collector QS integration."""

import unittest
from unittest.mock import patch

from v2.execution.budget import Budget
from v2.execution.resources import (
    AllocationFailureError,
    UnsupportedResourceVersionError,
)
from v2.qs.current_parallel import CurrentParallelSIQSJob
from v2.qs.parallel import CollectionPool
from v2.qs.storage import RelationPublicationError
from v2.tests.test_current_qs_workers import (
    N,
    configuration,
    prepared,
    settings,
)


def allocation_failure(task):
    raise AllocationFailureError("injected private allocation", "qs_worker")


class H3IntegrationTests(unittest.TestCase):
    def test_omitted_wall_cpu_and_rss_limits_do_not_cancel_packets(self):
        signatures = []
        for mode, count in (("serial", 1), ("process", 2)):
            job = CurrentParallelSIQSJob(
                N,
                seed=29,
                config=configuration(),
                workers=settings(mode, count),
                budget=Budget(seconds=None, cpu_seconds=None, rss_bytes=None),
            )
            result = job.run(fixed_work=True)
            self.assertTrue(result.stats["schedule_complete"])
            signatures.append(
                (dict(job.engine.collector._atoms), result.stats["scanned"])
            )
        self.assertEqual(signatures[0], signatures[1])

    def test_snapshots_reject_the_old_resource_contract(self):
        job = prepared()
        snapshot = job.checkpoint()
        snapshot.pop("resource_version")
        with self.assertRaises(UnsupportedResourceVersionError):
            CurrentParallelSIQSJob.from_checkpoint(snapshot, budget=Budget())

    def test_interrupted_parent_publication_cannot_resume_or_checkpoint(self):
        job = prepared()
        original = job.engine.collector._admit

        def interrupted(atom, stats):
            original(atom, stats)
            raise MemoryError("injected after insertion")

        with patch.object(job.engine.collector, "_admit", interrupted):
            with self.assertRaises(AllocationFailureError) as caught:
                job.run(fixed_work=True)
        self.assertIsInstance(caught.exception.__cause__, MemoryError)
        with self.assertRaises(RelationPublicationError):
            job.checkpoint()
        with self.assertRaises(RelationPublicationError):
            job.run(fixed_work=True)

    def test_worker_allocation_failure_retains_its_cause_and_drains_pool(self):
        job = prepared(worker_settings=settings("process", 2))
        with CollectionPool("process", 2) as pool:
            with patch(
                "v2.qs.current_parallel.execute_packet", allocation_failure
            ):
                with self.assertRaises(AllocationFailureError) as caught:
                    job.run(pool=pool, fixed_work=True)
        self.assertEqual(
            caught.exception.payload,
            dict(
                reason="allocation_failure",
                operation="injected private allocation",
                owner="qs_worker",
            ),
        )
        self.assertTrue(pool.closed)
        self.assertTrue(all(not p.is_alive() for p in pool._owned_processes))
