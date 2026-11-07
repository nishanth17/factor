"""Bounded diagnostics preserve returns, exceptions and nested spans."""

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.qs.a7a import performance_lock
from v2.benchmarks.qs.a7a.packet_observer import (
    PacketCounters,
    Recorder,
    initialize_worker,
    validate_capture,
)


class PacketObserverTests(unittest.TestCase):
    def test_spawned_observer_preserves_exact_admitted_atoms(self):
        from v2.qs import parallel
        from v2.tests.test_current_qs_workers import prepared, settings

        control = prepared(extended=True)
        control.run(fixed_work=True)
        observed = prepared(
            extended=True, worker_settings=settings("process", 2)
        )
        with tempfile.TemporaryDirectory() as temporary:
            with patch.dict(os.environ, FACTOR_A7A_PACKET_DIR=temporary):
                recorder = Recorder("parent")
                with patch.object(
                    parallel, "_initialize_process", initialize_worker
                ):
                    with parallel.CollectionPool("process", 2) as pool:
                        result = observed.run(pool=pool, fixed_work=True)
                recorder.stream.close()
            validate_capture(temporary, 2)
            self.assertNotEqual(result.reason, "worker_failure")
            self.assertEqual(
                set(control.engine.collector._atoms),
                set(observed.engine.collector._atoms),
            )
            for key, atom in control.engine.collector._atoms.items():
                other = observed.engine.collector._atoms[key]
                self.assertEqual(atom.exponents, other.exponents)
                self.assertEqual(atom.residual_primes, other.residual_primes)

    def test_nested_spans_are_counted_once_and_exceptions_escape(self):
        counters = PacketCounters()
        counters.begin()
        inner = counters.wrap("split", lambda: 7)
        outer = counters.wrap("split", inner)
        self.assertEqual(outer(), 7)

        def fail():
            raise ValueError("original")

        with self.assertRaisesRegex(ValueError, "original"):
            counters.wrap("recovery", fail)()
        result = counters.finish()
        self.assertEqual(result["split"][0], 1)
        self.assertEqual(result["recovery"][0], 1)
        self.assertFalse(counters.active)

    def test_sampling_has_fixed_space_and_never_changes_stop(self):
        counters = PacketCounters(interval=4)
        counters.begin()
        for _ in range(100):
            self.assertTrue(counters.stop_poll(lambda: True))
        result = counters.finish()
        self.assertEqual(result["poll_calls"], 100)
        self.assertEqual(result["poll_samples"], 25)
        self.assertEqual(sum(result["poll_histogram"]), 25)
        self.assertEqual(len(result["poll_histogram"]), 64)
        with self.assertRaises(ValueError):
            PacketCounters(interval=3)

    def test_recorder_ownership_and_exhaustion_are_explicit(self):
        with tempfile.TemporaryDirectory() as temporary:
            with patch.dict(os.environ, FACTOR_A7A_PACKET_DIR=temporary):
                recorder = Recorder("parent")
            recorder.rows = 16384
            with self.assertRaisesRegex(RuntimeError, "recorder allowance"):
                recorder.write("extra")
            recorder.stream.close()
            with self.assertRaisesRegex(RuntimeError, "exhausted"):
                validate_capture(temporary, 4)
            (Path(temporary) / (str(os.getpid()) + ".failure.json")).unlink()
            with self.assertRaisesRegex(RuntimeError, "all owned workers"):
                validate_capture(temporary, 4)
            first = json.loads(recorder.path.read_text())
            self.assertEqual(first["role"], "parent")

    def test_inherited_lock_duplicates_only_the_real_lock_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "lock"
            path.touch()
            with path.open("r+") as original:
                with patch.object(performance_lock, "LOCK_PATH", path):
                    with patch.dict(
                        os.environ,
                        FACTOR_A7A_PERFORMANCE_LOCK_FD=str(original.fileno()),
                    ):
                        with performance_lock.open_lock() as copied:
                            self.assertNotEqual(
                                copied.fileno(), original.fileno()
                            )
                    self.assertFalse(original.closed)
                    other = Path(temporary) / "other"
                    other.touch()
                    with patch.object(performance_lock, "LOCK_PATH", other):
                        with patch.dict(
                            os.environ,
                            FACTOR_A7A_PERFORMANCE_LOCK_FD=str(
                                original.fileno()
                            ),
                        ):
                            with self.assertRaises(ValueError):
                                performance_lock.open_lock()
