"""Partition, lease, verifier and restart contracts for coarse collection."""

import copy
import io
import json
import unittest
from contextlib import redirect_stdout
from dataclasses import replace
from unittest.mock import patch

from v2.budget import Budget
from v2.qs.families import _checksum
from v2.qs.parallel import (
    CollectionPool,
    ParallelConfig,
    ParallelSIQSJob,
    _collect,
)
from v2.qs.relations import verify_atomic


def allowance(work=200_000_000, **options):
    return Budget(
        work_limit=work,
        seconds=options.get("seconds"),
        cpu_seconds=options.get("cpu_seconds"),
    )


def configuration(**options):
    return replace(ParallelConfig(), family_count=4, pool_size=8, **options)


def mutate(checkpoint, change):
    result = copy.deepcopy(checkpoint)
    payload = json.loads(result["blob"])
    change(payload)
    result["blob"] = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    result["sha256"] = _checksum(payload)
    return result


class ParallelTests(unittest.TestCase):
    def test_fixed_work_keeps_scanning_after_a_direct_residual_split(self):
        config = replace(configuration(), pool_size=16, max_batch_atoms=512)
        expected = None

        for mode, workers in (("serial", 1), ("process", 2)):
            job = ParallelSIQSJob(
                4513 * 7717, config=config, budget=allowance()
            )
            with CollectionPool(mode, workers) as pool:
                result = job.run(pool=pool, fixed_work=True)

            self.assertEqual(result.reason, "factor_found")
            self.assertEqual(job.next_assignment, len(job.assignments))
            self.assertEqual(job.stats["scanned"], len(job.assignments) * 513)
            self.assertTrue(result.stats["schedule_complete"])
            self.assertIsNotNone(job.deferred_divisor)
            signature = (
                job.assignments,
                job.engine.collector._atoms,
                job.budget.used,
            )
            if expected is None:
                expected = signature
            else:
                self.assertEqual(signature, expected)

        job = ParallelSIQSJob(4513 * 7717, config=config, budget=allowance())
        job.run(max_assignments=3, fixed_work=True)

        restored = ParallelSIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )

        self.assertEqual(restored.deferred_divisor, job.deferred_divisor)

        result = restored.run(fixed_work=True)

        self.assertTrue(result.stats["schedule_complete"])
        used = restored.budget.used

        self.assertEqual(restored.run(fixed_work=True).divisor, result.divisor)
        self.assertEqual(restored.budget.used, used)

    def test_failed_worker_accounting_releases_job_and_pool_guards(self):
        with CollectionPool("thread", 2) as pool:
            job = ParallelSIQSJob(
                4001 * 5003, config=configuration(), budget=allowance()
            )
            with patch.object(
                pool, "finish", side_effect=RuntimeError("snapshot")
            ):
                with self.assertRaisesRegex(RuntimeError, "snapshot"):
                    job.run(pool=pool, fixed_work=True)

            self.assertFalse(pool.lock.locked())
            self.assertFalse(job._running)
            self.assertIs(job.engine.budget, job.budget)
            self.assertEqual(
                job.run(pool=pool, fixed_work=True).reason, "factor_found"
            )

    def test_unpublished_family_refusal_is_not_overwritten(self):
        job = ParallelSIQSJob(
            4001 * 5003, config=configuration(), budget=allowance()
        )

        def selective_refusal(task, state):
            if task[2] == 0:
                task = (
                    *task[:3],
                    replace(task[3], max_batch_atoms=1),
                    *task[4:],
                )
            return _collect(task, state)

        with CollectionPool("thread", 4) as pool:
            with patch(
                "v2.qs.parallel._collect", side_effect=selective_refusal
            ):
                result = job.run(pool=pool)

        self.assertEqual(result.reason, "batch_limit")
        self.assertEqual(result.next_position, 0)
        self.assertFalse(job.engine.collector._atoms)
        self.assertTrue(job.pending)
        self.assertEqual(result.cofactor, job.n)
        used = job.budget.used

        restored = ParallelSIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )

        self.assertEqual(restored.run().reason, "factor_found")
        self.assertGreater(restored.budget.used, used)

    def test_fixed_assignments_match_serial_threads_and_spawned_processes(
        self,
    ):
        expected = None

        for mode, workers in (
            ("serial", 1),
            ("thread", 2),
            ("process", 2),
            ("process", 4),
        ):
            with CollectionPool(mode, workers) as pool:
                job = ParallelSIQSJob(
                    4001 * 5003,
                    seed=29,
                    config=configuration(),
                    budget=allowance(),
                )

                result = job.run(pool=pool, fixed_work=True)

            self.assertEqual(result.reason, "factor_found")
            self.assertEqual(result.divisor * result.cofactor, job.n)
            self.assertEqual(job.next_assignment, len(job.assignments))
            self.assertFalse(job.pending)
            store = job.engine.collector._atoms

            for atom in store.values():
                self.assertTrue(
                    verify_atomic(
                        atom,
                        job.base,
                        residual_bound=10000,
                        budget=allowance(),
                    )
                )

            signature = (job.assignments, store, job.budget.used)
            if expected is None:
                expected = signature
            else:
                self.assertEqual(signature, expected)

            self.assertLessEqual(
                result.stats["peak_owned_bytes"], job.config.memory_bytes
            )

    def test_pause_checkpoint_and_resume_with_different_worker_count(self):
        job = ParallelSIQSJob(
            4001 * 5003, config=configuration(), budget=allowance()
        )
        with CollectionPool("process", 2) as pool:
            result = job.run(pool=pool, max_assignments=1, fixed_work=True)

        self.assertEqual(result.reason, "paused")
        checkpoint = json.loads(json.dumps(job.checkpoint()))

        restored = ParallelSIQSJob.from_checkpoint(
            checkpoint, budget=allowance()
        )

        self.assertEqual(restored.next_assignment, 1)
        self.assertEqual(
            restored.engine.collector._atoms, job.engine.collector._atoms
        )
        self.assertGreater(
            restored.budget.used, checkpoint["resources"]["work_used"]
        )
        self.assertGreaterEqual(
            restored.budget.prior_cpu, checkpoint["resources"]["cpu_used"]
        )
        with CollectionPool("thread", 4) as pool:
            result = restored.run(pool=pool, fixed_work=True)

        self.assertEqual(result.divisor * result.cofactor, job.n)
        original = ParallelSIQSJob(
            job.n, config=configuration(), budget=allowance()
        )
        original.run(fixed_work=True)

        self.assertEqual(restored.assignments, original.assignments)
        self.assertEqual(
            restored.engine.collector._atoms, original.engine.collector._atoms
        )

    def test_work_reservation_refusal_and_cumulative_resume(self):
        job = ParallelSIQSJob(
            4001 * 5003, config=configuration(), budget=allowance(100_000)
        )

        result = job.run()

        self.assertEqual(result.reason, "work_limit")
        # A partial lease may start work; refusal still retains exact costs
        # and the unpublished assignment for cumulative replay.
        self.assertGreater(job.stats["attempts"], 0)
        self.assertLessEqual(job.budget.used, 100_000)
        used = job.budget.used
        job.budget.work_limit = 200_000_000

        self.assertEqual(job.run().reason, "factor_found")
        self.assertGreater(job.budget.used, used)
        job.budget = allowance()
        with self.assertRaises(ValueError):
            job.run()

    def test_private_prefix_refusal_replays_without_publishing_rows(self):
        job = ParallelSIQSJob(
            4001 * 5003,
            config=configuration(assignment_work=10000),
            budget=allowance(),
        )

        result = job.run()

        self.assertEqual(result.reason, "work_limit")
        self.assertEqual(job.next_assignment, 0)
        self.assertFalse(job.engine.collector._atoms)
        self.assertGreater(job.stats["cancelled_work"], 0)
        used, assignments = job.budget.used, job.assignments
        job.config = replace(job.config, assignment_work=10_000_000)

        self.assertEqual(job.run().reason, "factor_found")
        self.assertEqual(job.assignments, assignments)
        self.assertGreater(job.budget.used, used)

    def test_pending_admission_and_solver_interruption_keep_provenance(self):
        job = ParallelSIQSJob(
            4001 * 5003, config=configuration(), budget=allowance()
        )
        original = job._merge

        def refuse(result):
            job.budget.work_limit = job.budget.used
            return original(result)

        with patch.object(job, "_merge", side_effect=refuse):
            self.assertEqual(job.run().reason, "work_limit")

        self.assertIn(0, job.pending)
        job.budget.work_limit = 200_000_000

        restored = ParallelSIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )

        self.assertEqual(restored.run().reason, "factor_found")
        self.assertEqual(restored.stats["attempts"], 1)

    def test_corrupt_and_rehashed_checkpoint_math_and_caps_are_rejected(self):
        job = ParallelSIQSJob(
            4001 * 5003, config=configuration(), budget=allowance()
        )
        job.run(max_assignments=1, fixed_work=True)
        checkpoint = job.checkpoint()

        for change in (
            lambda p: p["store"]["atoms"][0].__setitem__(2, 0),
            lambda p: p.__setitem__("next_assignment", 100),
            lambda p: p.__setitem__("cursor", 10000),
            lambda p: p.__setitem__("divisor", 1),
        ):
            with self.assertRaises((ValueError, TypeError)):
                ParallelSIQSJob.from_checkpoint(
                    mutate(checkpoint, change), budget=allowance()
                )

        with self.assertRaises(ValueError):
            ParallelSIQSJob.from_checkpoint(checkpoint, budget=allowance(0))
        job.config = replace(job.config, checkpoint_bytes=4096)
        with self.assertRaises(MemoryError):
            job.checkpoint()

    def test_aggregate_cpu_wall_and_owned_memory_refusals(self):
        for options, reason in (
            ({"seconds": 0}, "wall_limit"),
            ({"cpu_seconds": 0}, "cpu_limit"),
        ):
            with CollectionPool("process", 2) as pool:
                job = ParallelSIQSJob(4001 * 5003, budget=allowance(**options))

                self.assertEqual(job.run(pool=pool).reason, reason)
                self.assertFalse(pool.lock.locked())

        job = ParallelSIQSJob(
            4001 * 5003,
            config=configuration(memory_bytes=4 * 2**20),
            budget=allowance(),
        )

        self.assertEqual(job.run().reason, "memory_limit")
        self.assertEqual(job.stats["attempts"], 0)
        with CollectionPool("process", 4) as pool:
            job = ParallelSIQSJob(
                4001 * 5003, budget=allowance(cpu_seconds=0.03)
            )

            result = job.run(pool=pool)

            self.assertEqual(result.reason, "cpu_limit")
            self.assertGreater(job.budget.cpu_used, 0.03)
            self.assertLessEqual(job.budget.used, job.budget.work_limit)

    def test_early_split_cancellation_pool_reuse_and_quiet_calls(self):
        with CollectionPool("process", 4) as pool:
            with redirect_stdout(io.StringIO()) as output:
                for seed in (7, 29):
                    job = ParallelSIQSJob(
                        4001 * 5003,
                        seed=seed,
                        config=configuration(),
                        budget=allowance(),
                    )

                    result = job.run(pool=pool)

                    self.assertEqual(result.divisor * result.cofactor, job.n)
                    self.assertLessEqual(len(job.pending), 4)
                    self.assertGreaterEqual(job.stats["attempts"], 1)
                    self.assertFalse(pool.lock.locked())

            self.assertEqual(output.getvalue(), "")


if __name__ == "__main__":
    unittest.main()
