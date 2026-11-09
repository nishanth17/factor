"""Future Stage Q owned-group lifecycle checks using tiny synthetic work."""

import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.ecm.b4a.ownership import GroupObserver, guard
from v2.benchmarks.ecm.b4a.runner import ROOT, run_owned


@unittest.skipUnless(
    os.environ.get("FACTOR_B4A_LIFECYCLE_QA") == "1",
    "future Stage Q lifecycle lane",
)
class OwnershipTests(unittest.TestCase):
    def limits(self, wall=8):
        return dict(
            wall=wall, cpu=12, shutdown=1, rss=2**30, output=2**20, processes=4
        )

    def test_early_exit_and_reaped_grandchild_cpu(self):
        workload = (
            "import time; t=time.process_time();\n"
            "while time.process_time()-t<0.25: pass"
        )
        script = (
            "import subprocess,sys,resource,json;"
            f"subprocess.run([sys.executable,'-c',{workload!r}],check=True);"
            "a=resource.getrusage(resource.RUSAGE_SELF);"
            "b=resource.getrusage(resource.RUSAGE_CHILDREN);"
            "print(json.dumps({'cpu':a.ru_utime+a.ru_stime"
            "+b.ru_utime+b.ru_stime}))"
        )
        with tempfile.TemporaryDirectory() as directory:
            record = run_owned(
                [sys.executable, "-B", "-c", script],
                Path(directory) / "cell",
                self.limits(),
            )
            reported = json.loads(
                (Path(directory) / "cell/stdout").read_text()
            )["cpu"]
            self.assertEqual(record["exit"], 0)
            self.assertTrue(record["group_empty"])
            self.assertGreaterEqual(record["observed_cpu"], reported - 0.01)
            self.assertGreater(reported, 0.24)

    def test_deadline_kills_and_drains_a_blocked_child(self):
        with tempfile.TemporaryDirectory() as directory:
            record = run_owned(
                [sys.executable, "-B", "-c", "import time; time.sleep(30)"],
                Path(directory) / "cell",
                self.limits(wall=3),
            )
            self.assertNotEqual(record["exit"], 0)
            self.assertTrue(record["group_empty"])
            self.assertLess(record["wall"], 5)

    def test_controller_eof_stops_owned_group(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            read_pipe, write_pipe = os.pipe()
            child = subprocess.Popen(
                [sys.executable, "-B", "-c", "import time; time.sleep(30)"],
                start_new_session=True,
            )
            guardian = None
            try:
                request = dict(
                    pgid=child.pid,
                    deadline=time.monotonic() + 10,
                    cpu=10,
                    rss=2**30,
                    processes=2,
                    outputs=[],
                    output=2**20,
                    shutdown=1,
                )
                path = directory / "request.json"
                receipt = directory / "receipt.json"
                path.write_text(json.dumps(request))
                guardian = subprocess.Popen(
                    [
                        sys.executable,
                        "-B",
                        "-m",
                        "v2.benchmarks.ecm.b4a.runner",
                        "guard",
                        str(path),
                        str(read_pipe),
                        str(receipt),
                    ],
                    cwd=ROOT,
                    pass_fds=(read_pipe,),
                )
                os.close(read_pipe)
                read_pipe = None
                os.close(write_pipe)
                write_pipe = None
                child.wait(timeout=4)
                guardian.wait(timeout=4)
                report = json.loads(receipt.read_text())
                self.assertEqual(report["reason"], "controller_closed")
                self.assertFalse(GroupObserver().sample(child.pid)["pids"])
            finally:
                for descriptor in (read_pipe, write_pipe):
                    if descriptor is not None:
                        os.close(descriptor)
                for process in (child, guardian):
                    if process is not None and process.poll() is None:
                        process.kill()
                        process.wait(timeout=2)

    def test_denied_recorder_is_a_harness_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            receipt = Path(directory) / "receipt.json"
            read_pipe, write_pipe = os.pipe()
            try:
                with (
                    patch.object(
                        GroupObserver,
                        "sample",
                        side_effect=RuntimeError("denied"),
                    ),
                    patch("v2.benchmarks.ecm.b4a.ownership.terminate"),
                    patch("v2.benchmarks.ecm.b4a.ownership.os.killpg"),
                ):
                    guard(dict(pgid=999999, shutdown=0), read_pipe, receipt)
                self.assertIn(
                    "recorder_failure",
                    json.loads(receipt.read_text())["reason"],
                )
            finally:
                os.close(read_pipe)
                os.close(write_pipe)


if __name__ == "__main__":
    unittest.main()


class CpuAccountingTests(unittest.TestCase):
    def test_absolute_ticks_convert_process_and_reaped_cpu_once(self):
        from v2.benchmarks.ecm.b4a.ownership import Usage, cpu_seconds

        parent = Usage(
            user=24_000_000,
            system=12_000_000,
            child_user=6_000_000,
            child_system=6_000_000,
        )
        self.assertEqual(cpu_seconds([parent], 125, 3), 2)
        # A live child's process total is separate from already-reaped work.
        child = Usage(user=24_000_000)
        self.assertEqual(cpu_seconds([parent, child], 125, 3), 3)
        # Intel's 1:1 clock is supported without an arm64 frequency constant.
        self.assertEqual(cpu_seconds([Usage(user=1_000_000_000)], 1, 1), 1)
        for numer, denom in ((0, 1), (1, 0)):
            with self.assertRaisesRegex(RuntimeError, "timebase"):
                cpu_seconds([parent], numer, denom)

    def test_observer_does_not_add_native_thread_cpu_twice(self):
        import ctypes
        from types import SimpleNamespace

        from v2.benchmarks.ecm.b4a.ownership import TaskInfo, Usage

        def pids(kind, group, buffer, capacity):
            buffer[0] = 1234
            return ctypes.sizeof(ctypes.c_int)

        def usage(pid, version, pointer):
            result = ctypes.cast(pointer, ctypes.POINTER(Usage)).contents
            result.user = 24_000_000
            return 0

        def task(pid, kind, argument, pointer, size):
            result = ctypes.cast(pointer, ctypes.POINTER(TaskInfo)).contents
            result.user = result.threads_user = 24_000_000
            result.thread_count = 4
            return ctypes.sizeof(TaskInfo)

        observer = object.__new__(GroupObserver)
        observer.numer, observer.denom = 125, 3
        observer.lib = SimpleNamespace(
            proc_listpids=pids, proc_pid_rusage=usage, proc_pidinfo=task
        )
        for _ in range(2):
            sample = observer.sample(1234)
            self.assertEqual(sample["cpu"], 1)
            self.assertEqual(sample["threads"], 4)
