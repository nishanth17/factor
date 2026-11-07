"""Recorder framing, finite reservations and process ownership admission."""

import inspect
import json
import os
import tempfile
import unittest
from array import array
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.qs.a7a.diagnostic_analysis import heap_census
from v2.benchmarks.qs.a7a.diagnostic_hook import _logical_owners
from v2.benchmarks.qs.a7a.diagnostics import _events, plan


class DiagnosticTests(unittest.TestCase):
    def test_optional_worker_counters_compose_after_cpu_heap_installation(
        self,
    ):
        from v2.benchmarks.qs.a7a.diagnostic_hook import initialize_worker

        for enabled in (False, True):
            calls = []
            with (
                patch.dict(
                    os.environ,
                    FACTOR_A7A_PACKET_DIR="capture" if enabled else "",
                ),
                patch(
                    "v2.qs.parallel._initialize_process",
                    side_effect=lambda state: calls.append("base"),
                ),
                patch(
                    "v2.benchmarks.qs.a7a.diagnostic_hook.install",
                    side_effect=lambda role: calls.append("cpu_heap"),
                ),
                patch(
                    "v2.benchmarks.qs.a7a.packet_observer.install_worker",
                    side_effect=lambda state, **kwargs: calls.append(
                        ("counters", kwargs)
                    ),
                ),
            ):
                initialize_worker({})
            self.assertEqual(
                calls,
                ["base", "cpu_heap"]
                + ([("counters", {"stop_samples": False})] if enabled else []),
            )

    def test_collector_snapshot_exposes_actual_numerical_and_memory_domain(
        self,
    ):
        from v2.tests.test_qs_hit_words import make_worker

        worker = make_worker(64)
        worker._opaque_staging = SimpleNamespace(
            atoms={1: None, 2: None}, used=4096
        )
        frame = inspect.currentframe()
        try:
            owners = _logical_owners(frame)
        finally:
            del frame
        snapshot = next(
            owner["derived_buffers"]
            for owner in owners
            if owner["owner"] == "worker"
        )
        self.assertEqual(
            snapshot["collector_memory_limit"], worker.config.memory_limit
        )
        self.assertEqual(
            snapshot["factor_base_bound"], worker.factor_base.bound
        )
        self.assertEqual(
            snapshot["factor_base_count"], len(worker.factor_base.entries)
        )
        self.assertEqual(snapshot["candidate_bound"], worker._candidate_bound)
        self.assertEqual(snapshot["split_calls"], worker._split_calls)
        self.assertEqual(snapshot["opaque_staging_entries"], 2)
        self.assertEqual(snapshot["opaque_staging_accounted_bytes"], 4096)

    def test_heap_census_rejects_truncation_and_duplicate_root_marker(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "live.heap"
            for words in ([123, 2, 48, 456], [0, 0, 0, -1] * 2):
                path.write_bytes(array("q", words).tobytes())
                with self.assertRaises(ValueError):
                    heap_census(path)

    def test_heap_census_counts_roots_sizes_and_references(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "live.heap"
            words = [123, 2, 48, 456, -1, 0, 0, 0, -1, 456, 2, 24, -1]
            path.write_bytes(array("q", words).tobytes())
            result = heap_census(path)
            self.assertEqual(result["objects"], 2)
            self.assertEqual(result["root_objects"], 1)
            self.assertEqual(result["object_bytes"], 72)
            self.assertEqual(result["references"], 1)

    def test_recorder_only_signals_direct_owned_processes(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            child = folder / "12"
            child.mkdir()
            path = child / "events.json"
            ready = dict(kind="ready", pid=12, ppid=11, role="worker")
            path.write_text(json.dumps([ready]))
            self.assertEqual(list(_events(folder, 11)), [12])
            with self.assertRaises(RuntimeError):
                _events(folder, 10)

    def test_plan_prices_all_attempts_and_matches_observer_work(self):
        frozen = plan()
        cells = frozen["cells"]
        for key in ("wall", "cpu"):
            self.assertEqual(
                sum(row["cell"][key] for row in cells)
                + frozen[f"reserve_{key}"],
                frozen[key],
            )
        self.assertEqual(
            [r["mode"] for r in cells[:4]], ["off", "cpu", "cpu", "off"]
        )
        self.assertEqual({r["cell"]["polynomials"] for r in cells[:5]}, {768})
        self.assertEqual(cells[-1]["cell"]["endpoint"], "complete")


if __name__ == "__main__":
    unittest.main()
