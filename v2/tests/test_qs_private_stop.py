"""Per-worker cancellation retains exact collection and pool ownership."""

import os
import unittest
from concurrent.futures.process import BrokenProcessPool

from v2.benchmarks.infrastructure.resource_control import Budget
from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration as resource_configuration,
)
from v2.benchmarks.infrastructure.resource_control import (
    historical_replace as replace,
)
from v2.qs.current_parallel import CurrentParallelSIQSJob, QSWorkerConfig
from v2.qs.packet_lease import ParentCPUSnapshot
from v2.qs.parallel import CollectionPool
from v2.qs.worker_stop import configure_stop
from v2.tests.test_current_qs_workers import (
    N,
    allowance,
    configuration,
    prepared,
    settings,
)


def observe_worker_stop():
    from v2.qs.parallel import _PROCESS_STATE

    state, slot = _PROCESS_STATE
    state["barrier"].wait(timeout=10)
    return slot, state["stop"].worker_event(slot).is_set()


def crash_holding_private_stop():
    from v2.qs.parallel import _PROCESS_STATE

    state, slot = _PROCESS_STATE
    state["stop"].worker_event(slot)._cond.acquire()
    os._exit(23)


class PrivateStopTests(unittest.TestCase):
    def test_settings_require_process_leases_and_keep_disabled_identity(self):
        self.assertNotIn("private_stop", QSWorkerConfig().record())
        for kwargs in (
            dict(private_stop=1),
            dict(private_stop=True),
            dict(private_stop=True, lease_budget=True, mode="serial"),
        ):
            with self.assertRaises(ValueError):
                resource_configuration(QSWorkerConfig, **kwargs)
        config = QSWorkerConfig(private_stop=True, lease_budget=True)
        self.assertTrue(config.record()["private_stop"])
        self.assertNotIn("cpu_snapshot", QSWorkerConfig().record())
        for kwargs in (dict(cpu_snapshot=1), dict(cpu_snapshot=True)):
            with self.assertRaises(ValueError):
                resource_configuration(QSWorkerConfig, **kwargs)

    def test_parent_snapshot_charges_updates_once_and_rebinds_baseline(self):
        with CollectionPool("process", 2) as pool:
            first = allowance()
            pool.begin(first)
            snapshot = ParentCPUSnapshot(pool)
            pool.state["cpu"][:] = [1.0, 2.0]
            snapshot(first)
            self.assertEqual(first.prior_cpu, 0.0)
            snapshot(first)
            self.assertEqual(first.prior_cpu, 0.0)
            pool.state["cpu"][1] = 3.0
            snapshot(first)
            pool.poll(first)
            self.assertEqual(first.prior_cpu, 0.0)
            pool.lock.release()

            second = allowance()
            pool.begin(second)
            snapshot = ParentCPUSnapshot(pool)
            pool.state["cpu"][:] = [2.0, 4.0]
            snapshot(second)
            self.assertEqual(second.prior_cpu, 0.0)
            pool.lock.release()

    def test_broadcast_and_quiescent_clear_reach_every_worker(self):
        with CollectionPool("process", 4) as pool:
            configure_stop(pool, True)
            pool.state["stop"].events[0].set()
            rows = [
                pool.executor.submit(observe_worker_stop) for _ in range(4)
            ]
            self.assertEqual(
                sorted(f.result() for f in rows),
                [(slot, slot == 0) for slot in range(4)],
            )
            pool.state["stop"].set()
            rows = [
                pool.executor.submit(observe_worker_stop) for _ in range(4)
            ]
            self.assertEqual(
                sorted(f.result() for f in rows),
                [(slot, True) for slot in range(4)],
            )
            pool.begin(Budget(work_limit=1000, seconds=10, cpu_seconds=60))
            try:
                rows = [
                    pool.executor.submit(observe_worker_stop) for _ in range(4)
                ]
                self.assertEqual(
                    sorted(f.result() for f in rows),
                    [(slot, False) for slot in range(4)],
                )
            finally:
                pool.lock.release()
            with self.assertRaisesRegex(ValueError, "fresh process pool"):
                configure_stop(pool, False)

    def test_poisoned_private_flag_does_not_block_owned_abort(self):
        with CollectionPool("process", 1) as pool:
            configure_stop(pool, True)
            future = pool.executor.submit(crash_holding_private_stop)
            with self.assertRaises(BrokenProcessPool):
                future.result(timeout=10)
            pool.abort()
            self.assertTrue(pool.closed)
            self.assertFalse(any(p.is_alive() for p in pool._owned_processes))

    def test_shared_and_private_flags_preserve_exact_atoms_and_resume(self):
        reference = prepared()
        reference.run(fixed_work=True)
        expected = dict(reference.engine.collector._atoms)

        worker_settings = settings(
            "process",
            2,
            lease_budget=True,
            private_stop=True,
            reuse_state=True,
            packed_transport=True,
            verification="worker",
            cpu_snapshot=True,
        )
        job = CurrentParallelSIQSJob(
            N,
            seed=29,
            config=configuration(),
            workers=worker_settings,
            budget=allowance(),
        )
        with CollectionPool("process", 2) as pool:
            job.run(pool=pool, fixed_work=True, max_assignments=1)
            saved = job.checkpoint()
            restored = CurrentParallelSIQSJob.from_checkpoint(
                saved,
                budget=allowance(),
                workers=worker_settings,
            )
            self.assertTrue(restored.workers.private_stop)
            restored.run(pool=pool, fixed_work=True)
            self.assertEqual(dict(restored.engine.collector._atoms), expected)
            self.assertEqual(restored.stats["lost_work"], 0)
            incompatible = replace(worker_settings, private_stop=False)
            with self.assertRaises(ValueError):
                CurrentParallelSIQSJob.from_checkpoint(
                    saved,
                    budget=allowance(),
                    workers=incompatible,
                )


if __name__ == "__main__":
    unittest.main()
