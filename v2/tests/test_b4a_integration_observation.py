"""Canonical QS evidence survives packet partitioning and refuses overflow."""

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from v2.benchmarks.ecm.b4a.integration_observation import IntegrationObserver


class CanonicalQSEvidenceTests(unittest.TestCase):
    def test_scalar_handoff_releases_context_without_thread_lifecycle(self):
        from v2.benchmarks.ecm.b4a.integration_comparison import (
            native_handoff,
        )

        state = dict(
            native=True,
            topology="scalar",
            closed=True,
            contexts=0,
            live_contexts=0,
            shared_native_blocks=0,
            lifecycle=None,
        )
        evidence = dict(
            snapshots=[
                dict(label=label, native=state, native_threads=None)
                for label in ("qs-entry", "after-ecm-close")
            ]
        )
        result = native_handoff(evidence)
        self.assertTrue(result["reached_native_QS"])
        self.assertEqual(result["topology"], "scalar")
        state["live_contexts"] = 1
        with self.assertRaisesRegex(RuntimeError, "remain live"):
            native_handoff(evidence)

    def test_partitioning_preserves_every_common_prefix(self):
        rows = [
            [index, index - 17, 1, [[2, 1]], 1, [], None, False]
            for index in range(37)
        ]
        job = SimpleNamespace(identity="same-mathematical-job")
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory) / str(i) for i in range(2)]
            observers = [IntegrationObserver(None, path) for path in paths]
            for observer, chunks in zip(
                observers, ((rows,), (rows[:3], rows[3:19], rows[19:]))
            ):
                with observer.prefix_path.open("xb") as stream:
                    observer.stream = stream
                    for chunk in chunks:
                        observer._prefix(job, 0, chunk)
            self.assertEqual(paths[0].read_bytes(), paths[1].read_bytes())
            self.assertEqual(len(paths[0].read_bytes()), 37 * 32)
            self.assertEqual(
                paths[0].read_bytes()[: 19 * 32],
                paths[1].read_bytes()[: 19 * 32],
            )

    def test_capture_refusal_preserves_the_accepted_prefix(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "prefixes"
            observer = IntegrationObserver(None, path, prefix_bytes=64)
            job = SimpleNamespace(identity="job")
            with path.open("xb") as stream:
                observer.stream = stream
                observer._prefix(job, 0, [[1], [2]])
                before = stream.tell()
                with self.assertRaisesRegex(RuntimeError, "recorder"):
                    observer._prefix(job, 0, [[3]])
                self.assertEqual(stream.tell(), before)
                self.assertEqual(observer.assignments["job/0"]["rows"], 2)

    def test_extra_speculation_compares_common_work_and_detects_tampering(
        self,
    ):
        from v2.benchmarks.ecm.b4a.integration_comparison import (
            compare_prefixes,
        )

        job = SimpleNamespace(identity="same-job")
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory) / str(i) for i in range(2)]
            evidence = []
            for path, count in zip(paths, (5, 7)):
                observer = IntegrationObserver(None, path)
                with path.open("xb") as stream:
                    observer.stream = stream
                    observer._prefix(job, 0, [[i] for i in range(count)])
                evidence.append(
                    dict(
                        jobs={"same-job": dict(seed=73)},
                        assignments=list(observer.assignments.values()),
                    )
                )
            result = compare_prefixes(*evidence, *paths)
            self.assertEqual(result["common_published_rows"], 5)
            self.assertEqual(result["right_extra_published_rows"], 2)
            value = bytearray(paths[1].read_bytes())
            value[31] ^= 1
            paths[1].write_bytes(value)
            with self.assertRaises(ArithmeticError):
                compare_prefixes(*evidence, *paths)

    def test_publication_clock_follows_the_real_checked_split(self):
        from v2 import portfolio
        from v2.benchmarks.ecm.b4a.confirmation_runtime import measured_call

        state = dict(current=dict(n=15, mult=1), pending=[])

        class Session:
            def call(self, request):
                portfolio._split(state, request["divisor"])
                return dict(complete=True)

        original = portfolio._split
        result = measured_call(Session(), dict(divisor=3))
        self.assertEqual(state["pending"], [[3, 1], [5, 1]])
        self.assertEqual(result["publications"][0]["divisor"], 3)
        self.assertGreaterEqual(
            result["publications"][0]["at"], result["dispatch_started"]
        )
        self.assertIs(portfolio._split, original)
