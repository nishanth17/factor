"""The diagnostic changes guard placement, never private grants or charging."""

import multiprocessing
import unittest
from unittest.mock import patch

from v2.benchmarks.qs.a7a.parent_guard_probe import execute_packet
from v2.qs import parallel
from v2.tests.test_current_qs_workers import prepared


def state():
    context = multiprocessing.get_context("spawn")
    return dict(
        stop=context.Event(),
        cpu=context.Array("d", (0.0,)),
        rss=context.Array("q", 1),
    )


class ParentGuardProbeTests(unittest.TestCase):
    def test_probe_receives_identical_task_and_retains_stop_and_publication(
        self,
    ):
        shared, task = state(), {"private_grants": "unchanged"}

        def collect(received, *, cancelled, publish_cpu):
            self.assertIs(received, task)
            self.assertFalse(cancelled())
            publish_cpu(None)
            self.assertGreater(shared["cpu"][0], 0)
            shared["stop"].set()
            self.assertTrue(cancelled())
            return "checked"

        with (
            patch.object(parallel, "_PROCESS_STATE", (shared, 0)),
            patch("v2.qs.worker_packet.collect_packet", side_effect=collect),
        ):
            self.assertEqual(execute_packet(task), "checked")
        self.assertGreater(shared["rss"][0], 0)

    def test_real_packet_keeps_private_work_cpu_and_cancel_refusal(self):
        job = prepared()
        for reason in ("work_limit", "cpu_limit", "cancelled"):
            with self.subTest(reason=reason):
                shared = state()
                task = job._task(0, 10**8, 10, 10)
                if reason == "work_limit":
                    task["work"] = 0
                elif reason == "cpu_limit":
                    task["cpu"] = 0
                else:
                    shared["stop"].set()

                with patch.object(parallel, "_PROCESS_STATE", (shared, 0)):
                    result = execute_packet(task)

                self.assertEqual(
                    result["reason"],
                    "cancelled" if reason == "cancelled" else "complete",
                )
                self.assertGreaterEqual(result["work"], 0)
                self.assertGreater(shared["cpu"][0], 0)


if __name__ == "__main__":
    unittest.main()
