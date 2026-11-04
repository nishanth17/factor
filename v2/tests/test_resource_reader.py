"""Same-counter RSS observation, portable fallback and native failures."""

import resource
import sys
import unittest
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.support.resources import PeakRSSReader


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
