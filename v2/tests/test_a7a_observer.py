"""Owned RSS ABI, recorder failures and distinct observation allowances."""

import ctypes
import errno
import os
import select
import subprocess
import sys
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.qs.a7a.observer import (
    DarwinRSS,
    OwnedRSS,
    _TaskInfo,
    validate_observation,
)


class ObserverTests(unittest.TestCase):
    def test_task_info_reads_bytes_and_distinguishes_a_gone_pid(self):
        reader = DarwinRSS.__new__(DarwinRSS)

        def query(pid, flavor, arg, buffer, size):
            self.assertEqual((flavor, arg, size), (4, 0, 96))
            if pid == 13:
                ctypes.set_errno(errno.ESRCH)
                return 0
            buffer.contents.resident_size = 123456
            return size

        reader.query = query
        self.assertEqual(reader([12, 13]), ({12: 123456}, [13]))
        self.assertEqual(ctypes.sizeof(_TaskInfo), 96)

        reader.query = lambda *args: 48
        with self.assertRaisesRegex(OSError, "48 task-info bytes"):
            reader([12])

    @unittest.skipUnless(sys.platform == "darwin", "Darwin RSS API")
    def test_current_process_has_positive_native_resident_bytes(self):
        rss, missing = DarwinRSS()([os.getpid()])
        self.assertFalse(missing)
        self.assertGreater(rss[os.getpid()], 0)

    @unittest.skipUnless(sys.platform == "darwin", "Darwin RSS API")
    def test_owned_child_rss_matches_ps_for_a_quiescent_process(self):
        command = (
            "import sys; print('ready', flush=True); sys.stdin.readline()"
        )
        with subprocess.Popen(
            [sys.executable, "-B", "-c", command],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
        ) as child:
            try:
                self.assertTrue(select.select([child.stdout], [], [], 5)[0])
                self.assertEqual(child.stdout.readline(), "ready\n")
                pool = SimpleNamespace(
                    mode="process",
                    executor=SimpleNamespace(_processes={child.pid: child}),
                )
                observer = OwnedRSS(pool, limit=2**40, samples=1)
                observer._run()
                rss = observer.rows[0]["processes"]
                self.assertIn(os.getpid(), rss)
                self.assertIn(child.pid, rss)
                reported = subprocess.check_output(
                    ["ps", "-o", "rss=", "-p", str(child.pid)],
                    text=True,
                    timeout=3,
                )
                self.assertLessEqual(
                    abs(rss[child.pid] - 1024 * int(reported)), 2**20
                )
            finally:
                try:
                    child.communicate("done\n", timeout=5)
                except subprocess.TimeoutExpired:
                    child.kill()
                    child.communicate(timeout=5)

    def test_recorder_cap_is_separate_from_memory_exhaustion(self):
        for limit, exceeded in ((100, False), (98, True)):
            with (
                patch("sys.platform", "darwin"),
                patch(
                    "v2.benchmarks.qs.a7a.observer.DarwinRSS",
                    return_value=lambda pids: ({os.getpid(): 99}, []),
                ),
            ):
                observer = OwnedRSS(None, limit=limit, samples=1)
            observer._run()
            self.assertEqual(observer.exceeded, exceeded)
            self.assertEqual(observer.rows[0]["backend"], "darwin-libproc")
            if exceeded:
                self.assertIsNone(observer.error)
            else:
                self.assertIn("sample allowance", observer.error)

    def test_live_guard_counts_the_external_controller(self):
        parent = os.getpid()
        controller = parent + 100000

        def read(pids):
            self.assertEqual(set(pids), {parent, controller})
            return {parent: 99, controller: 2}, []

        with (
            patch("sys.platform", "darwin"),
            patch(
                "v2.benchmarks.qs.a7a.observer.DarwinRSS", return_value=read
            ),
            patch(
                "multiprocessing.resource_tracker._resource_tracker._pid", None
            ),
        ):
            observer = OwnedRSS(None, limit=100, owner_pid=controller)
        observer._run()

        self.assertTrue(observer.exceeded)
        self.assertEqual(observer.rows[0]["bytes"], 101)
        self.assertIsNone(observer.error)

    def test_failed_observation_never_passes_as_factorer_cancellation(self):
        result = dict(
            cell=dict(observe=True),
            observer_failure=None,
            observer_exhausted=False,
            rss_samples=[dict(bytes=1)],
        )
        validate_observation(result)
        for change, message in (
            (dict(observer_failure="timeout"), "RSS recorder failure"),
            (dict(observer_exhausted=True), "RSS exceeded"),
            (dict(rss_samples=[]), "no samples"),
        ):
            with self.assertRaisesRegex(RuntimeError, message):
                validate_observation(dict(result, **change))


if __name__ == "__main__":
    unittest.main()
