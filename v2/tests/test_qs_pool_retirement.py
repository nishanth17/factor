"""A crash can wake the owner while the executor manager is still reaping."""

import multiprocessing
import os
import signal
import unittest
from unittest.mock import patch

from v2.qs.parallel import CollectionPool


def _crash_holding_stop(task):
    from v2.qs.parallel import _PROCESS_STATE

    state, _ = _PROCESS_STATE
    state["stop"]._cond.acquire()
    os._exit(24)


def _poisoned_stop_scenario(connection):
    """Own a separate group so a regressed cleanup cannot hang the suite."""
    os.setsid()
    from v2.tests.test_current_qs_workers import prepared, settings

    try:
        job = prepared(worker_settings=settings("process", 2))
        with CollectionPool("process", 2) as pool:
            with patch(
                "v2.qs.current_parallel.execute_packet", _crash_holding_stop
            ):
                result = job.run(pool=pool, fixed_work=True)
            connection.send(
                dict(
                    reason=result.reason,
                    lost_work=job.stats["lost_work"],
                    closed=pool.closed,
                    retired=all(
                        not p.is_alive() for p in pool._owned_processes
                    ),
                )
            )
    finally:
        connection.close()


class RetirementTests(unittest.TestCase):
    def test_hard_abort_never_signals_a_possibly_poisoned_stop_event(self):
        class PoisonedEvent:
            def set(self):
                raise AssertionError("acquired a dead worker's stop lock")

        with CollectionPool("process", 1) as pool:
            pool.state["stop"] = PoisonedEvent()
            pool.abort()
            self.assertTrue(pool.closed)

    @unittest.skipUnless(os.name == "posix", "requires an owned process group")
    def test_detected_crash_while_polling_stop_has_finite_retirement(self):
        context = multiprocessing.get_context("spawn")
        receive, send = context.Pipe(duplex=False)
        process = context.Process(target=_poisoned_stop_scenario, args=(send,))
        process.start()
        send.close()
        try:
            process.join(timeout=10)
            if process.is_alive():
                # The helper creates its group before starting descendants.
                if os.getpgid(process.pid) == process.pid:
                    os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
                process.join(timeout=2)
                self.fail(
                    "detected worker crash did not retire within 10 seconds"
                )
            self.assertEqual(process.exitcode, 0)
            self.assertTrue(receive.poll(1))
            result = receive.recv()
            self.assertEqual(result["reason"], "worker_failure")
            self.assertGreater(result["lost_work"], 0)
            self.assertTrue(result["closed"])
            self.assertTrue(result["retired"])
        finally:
            receive.close()

    def test_failed_future_owner_does_not_join_concurrently_with_manager(self):
        class Manager:
            alive = True

            def join(self, timeout):
                self.alive = False

            def is_alive(self):
                return self.alive

        manager = Manager()

        class Executor:
            _executor_manager_thread = manager

            def shutdown(self, **kwargs):
                self._executor_manager_thread = None

        class Process:
            pid = 1  # A fake handle: no OS signal is sent by this oracle.
            alive = True

            def is_alive(self):
                return self.alive

            def terminate(self):
                self.alive = False

            def kill(self):
                self.alive = False

            def join(self, timeout):
                if manager.alive:
                    raise AssertionError("concurrent child reaping")

        with CollectionPool("process", 1) as pool:
            pool.executor.shutdown(wait=True, cancel_futures=True)
            pool.executor = Executor()
            pool._owned_processes = [Process()]
            pool.abort()
        self.assertTrue(pool.closed)
        self.assertFalse(manager.alive)


if __name__ == "__main__":
    unittest.main()
