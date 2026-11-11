"""Same-counter RSS observation, portable fallback and native failures."""

import ctypes
import errno
import os
import resource
import sys
import unittest
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.support.resources import PeakRSSReader
from v2.execution.resources import (
    MacResources,
    ProcessObservationDeniedError,
    ProcessSample,
)


class ResourceReaderTests(unittest.TestCase):
    def test_portable_reads_are_not_cached(self):
        reader = PeakRSSReader(native=False)
        with patch(
            "v2.benchmarks.support.resources.resource.getrusage",
            side_effect=[
                SimpleNamespace(ru_maxrss=11),
                SimpleNamespace(ru_maxrss=13),
            ],
        ) as query:
            scale = 1 if sys.platform == "darwin" else 1024
            self.assertEqual(reader(), 11 * scale)
            self.assertEqual(reader(), 13 * scale)
            self.assertEqual(query.call_count, 2)

    def test_unsupported_native_provider_falls_back_quietly(self):
        with (
            patch("v2.benchmarks.support.resources.sys.platform", "darwin"),
            patch(
                "v2.benchmarks.support.resources.platform.machine",
                return_value="arm64",
            ),
            patch(
                "v2.benchmarks.support.resources._darwin_reader",
                side_effect=ImportError,
            ),
        ):
            reader = PeakRSSReader()
        self.assertEqual(reader.identity, "resource-getrusage")
        self.assertGreater(reader(), 0)

    def test_native_reader_brackets_actual_rss_and_uses_thread_local_buffers(
        self,
    ):
        reader = PeakRSSReader()
        scale = 1 if sys.platform == "darwin" else 1024
        for _ in range(8):
            before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale
            observed = reader()
            after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale
            self.assertLessEqual(before, observed)
            self.assertLessEqual(observed, after)
        with ThreadPoolExecutor(max_workers=2) as pool:
            results = list(pool.map(lambda _: reader(), range(16)))
        self.assertTrue(all(value > 0 for value in results))

    def test_query_failure_is_never_hidden_after_native_admission(self):
        reader = PeakRSSReader(native=False)

        def fail():
            raise OSError("resource query failed")

        reader._native = fail
        with self.assertRaisesRegex(OSError, "resource query failed"):
            reader()


class MacCpuReaderTests(unittest.TestCase):
    def setUp(self):
        self.pid = os.getpid()
        self.reader = MacResources.__new__(MacResources)
        self.reader._cpu_numer = 125
        self.reader._cpu_denom = 3_000_000_000
        self.reader.system = SimpleNamespace(waitid=self.waitid)
        self.reader.lib = SimpleNamespace(
            proc_pid_rusage=self.rusage,
            proc_listchildpids=lambda *args: 0,
        )

    def waitid(self, kind, pid, info, flags):
        self.assertEqual(
            (kind, pid, flags), (1, self.pid, 4 | 32 | os.WNOHANG)
        )
        info._obj.pid = pid
        return 0

    def rusage(self, pid, version, pointer):
        self.assertEqual((pid, version), (self.pid, 2))
        usage = pointer._obj
        usage.start = 42
        usage.resident = 2048
        usage.user, usage.system = 3_000_000, 3_000_000
        usage.child_user, usage.child_system = 12_000_000, 0
        return 0

    def test_live_and_final_use_same_seconds_conversion(self):
        self.assertEqual(
            self.reader.sample({self.pid: 42}),
            [ProcessSample(self.pid, 42, 2048, 0.25, 0.5)],
        )
        self.assertEqual(
            self.reader.peek_exit(self.pid, os.WNOHANG),
            ProcessSample(self.pid, 42, 0, 0.25, 0.5),
        )
        self.assertEqual(self.reader.sample({self.pid: 43}), [])

    def test_live_preflight_and_final_report_access_denial(self):
        def denied(*args):
            ctypes.set_errno(errno.EPERM)
            return -1

        self.reader.lib.proc_pid_rusage = denied
        for observe in (
            lambda: self.reader.sample([self.pid]),
            lambda: self.reader.peek_exit(self.pid, os.WNOHANG),
        ):
            with self.assertRaisesRegex(
                ProcessObservationDeniedError, f"proc_pid_rusage.*{self.pid}"
            ):
                observe()

    def test_empty_child_list_cannot_hide_access_denial(self):
        def denied(*args):
            ctypes.set_errno(errno.EPERM)
            return 0

        self.reader.lib.proc_listchildpids = denied
        with self.assertRaisesRegex(
            ProcessObservationDeniedError, f"proc_listchildpids.*{self.pid}"
        ):
            self.reader.sample([self.pid])
