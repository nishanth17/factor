"""Single-reader stop flags preserve cancellation, resume and owned cleanup."""

import multiprocessing
import os
import unittest
from concurrent.futures.process import BrokenProcessPool

from v2.benchmarks.infrastructure.resource_control import Budget
from v2.benchmarks.infrastructure.resource_control import (
    historical_replace as replace,
)
from v2.execution.budget import BudgetExhaustedError
from v2.qs.current_parallel import CurrentParallelSIQSJob, QSWorkerConfig
from v2.qs.packet_lease import PacketLeaseBudget
from v2.qs.parallel import CollectionPool
from v2.qs.worker_stop import configure_stop
from v2.tests.test_current_qs_workers import (
    N,
    allowance,
    configuration,
    prepared,
    settings,
)


def observe_setup_stop():
    from v2.qs.parallel import _PROCESS_STATE

    state, slot = _PROCESS_STATE
    state["barrier"].wait(timeout=10)
    return (
        slot,
        state["stop"].worker_event(slot).is_set(),
        state["stop"].is_set(),
    )


def exercise_paid_stop():
    from v2.qs.parallel import _PROCESS_STATE

    state, slot = _PROCESS_STATE
    raw = Budget(
        work_limit=1000,
        seconds=None,
        cpu_seconds=None,
        cancelled=state["stop"].worker_event(slot).is_set,
    )
    view = PacketLeaseBudget(raw)
    view.consume()
    state["test_ready"].set()
    if not state["test_continue"].wait(timeout=10):
        raise AssertionError("parent did not publish the test stop")
    for _ in range(63):
        view.consume()
    try:
        view.consume()
    except BudgetExhaustedError:
        return view.used, raw.used, view.reason
    raise AssertionError("cancelled packet consumed another polling interval")


def crash_with_stop_token():
    from v2.qs.parallel import _PROCESS_STATE

    state, slot = _PROCESS_STATE
    if not state["stop"].worker_event(slot).flag.acquire(False):
        raise AssertionError("missing published stop token")
    os._exit(23)


class SemaphoreStopTests(unittest.TestCase):
    def test_validation_and_default_private_event_identity(self):
        original = QSWorkerConfig(private_stop=True, lease_budget=True)
        self.assertNotIn("stop_backend", original.record())
        for value in ("unknown", None, "semaphore"):
            with self.assertRaises(ValueError):
                QSWorkerConfig(stop_backend=value)
        enabled = replace(original, stop_backend="semaphore")
        self.assertEqual(enabled.record()["stop_backend"], "semaphore")

    def test_setup_and_packet_checks_have_the_same_exclusive_reader(self):
        with CollectionPool("process", 2) as pool:
            configure_stop(pool, True, "semaphore")
            pool.state["stop"].events[0].set()
            self.assertTrue(pool.state["stop"].is_set())
            futures = [
                pool.executor.submit(observe_setup_stop) for _ in range(2)
            ]
            self.assertEqual(
                sorted(f.result() for f in futures),
                [(0, True, True), (1, False, False)],
            )
            for _ in range(3):
                pool.state["stop"].set()
            futures = [
                pool.executor.submit(observe_setup_stop) for _ in range(2)
            ]
            self.assertEqual(
                sorted(f.result() for f in futures),
                [(0, True, True), (1, True, True)],
            )
            pool.begin(allowance())
            try:
                futures = [
                    pool.executor.submit(observe_setup_stop) for _ in range(2)
                ]
                self.assertEqual(
                    sorted(f.result() for f in futures),
                    [(0, False, False), (1, False, False)],
                )
            finally:
                pool.lock.release()
            with self.assertRaises(ValueError):
                configure_stop(pool, True, "event")

    def test_parent_stop_retains_the_63_additional_action_bound(self):
        context = multiprocessing.get_context("spawn")
        with CollectionPool("process", 1) as pool:
            pool.state["test_ready"] = context.Event()
            pool.state["test_continue"] = context.Event()
            configure_stop(pool, True, "semaphore")
            future = pool.executor.submit(exercise_paid_stop)
            self.assertTrue(pool.state["test_ready"].wait(timeout=10))
            pool.state["stop"].set()
            pool.state["test_continue"].set()
            self.assertEqual(future.result(timeout=10), (64, 64, "cancelled"))

    def test_exact_atoms_survive_checkpoint_and_reused_pool(self):
        reference = prepared()
        reference.run(fixed_work=True)
        worker_settings = settings(
            "process",
            2,
            lease_budget=True,
            private_stop=True,
            stop_backend="semaphore",
            reuse_state=True,
            packed_transport=True,
            verification="worker",
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
            checkpoint = job.checkpoint()
            restored = CurrentParallelSIQSJob.from_checkpoint(
                checkpoint,
                budget=allowance(),
                workers=worker_settings,
            )
            restored.run(pool=pool, fixed_work=True)
            self.assertEqual(
                dict(restored.engine.collector._atoms),
                dict(reference.engine.collector._atoms),
            )
            self.assertEqual(restored.stats["lost_work"], 0)
        with self.assertRaises(ValueError):
            CurrentParallelSIQSJob.from_checkpoint(
                checkpoint,
                budget=allowance(),
                workers=replace(worker_settings, stop_backend="event"),
            )

    def test_crash_with_consumed_stop_token_retires_only_owned_pool(self):
        with CollectionPool("process", 1) as pool:
            configure_stop(pool, True, "semaphore")
            pool.state["stop"].set()
            future = pool.executor.submit(crash_with_stop_token)
            with self.assertRaises(BrokenProcessPool):
                future.result(timeout=10)
            pool.abort()
            self.assertTrue(pool.closed)
            self.assertFalse(any(p.is_alive() for p in pool._owned_processes))


if __name__ == "__main__":
    unittest.main()
