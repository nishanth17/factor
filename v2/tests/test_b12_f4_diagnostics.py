"""Independent diagnostic-format and factor-publication checks."""

import array
import contextlib
import io
import itertools
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.ecm.b12_f4 import profile_runner
from v2.benchmarks.ecm.b12_f4.dump_analysis import (
    heap_summary,
    sample_summary,
)
from v2.benchmarks.ecm.b12_f4.profile_runner import component_arm, signature


class DiagnosticTests(unittest.TestCase):
    def test_heap_live_sizes_exclude_the_root_marker(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sample.heap"
            data = array.array(
                "l", [8, 2, 24, 16, -1, 0, 0, 0, -1, 16, 3, 40, -1]
            )
            path.write_bytes(data.tobytes())
            summary = heap_summary(path)
            self.assertEqual(summary["objects"], 2)
            self.assertEqual(summary["references"], 1)
            self.assertEqual(summary["gc_object_bytes"], 64)

    def test_truncated_oversized_and_missing_marker_heaps_refuse(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "bad.heap"
            for data in (
                [8, 2, 24],
                [8, 2, 24, -1],
                [0, 0, 0, -1, 8, 2, 24, 16],
            ):
                path.write_bytes(array.array("l", data).tobytes())
                with self.assertRaises(ValueError):
                    heap_summary(path)
            path.write_bytes(array.array("l", [0, 0, 0, -1]).tobytes())
            with self.assertRaises(ValueError):
                heap_summary(path, maximum=1)

    def test_result_oracle_rejects_nondivisors_and_missing_curves(self):
        result = dict(
            id="curve",
            sigma=17,
            factor=5,
            reason="complete",
            phase="stage_two",
            chain_chunks=0,
        )
        checkpoint = {"search": {"payload": {"completed": {"0": result}}}}
        self.assertEqual(signature(checkpoint, 35, 1)[0][2], 5)
        result["factor"] = 6
        with self.assertRaises(ValueError):
            signature(checkpoint, 35, 1)
        with self.assertRaises(ValueError):
            signature(checkpoint, 35, 2)

    def test_live_integer_keys_and_json_string_keys_have_the_same_identity(
        self,
    ):
        result = dict(
            id="curve",
            sigma=17,
            factor=5,
            reason="complete",
            phase="stage_two",
            chain_chunks=0,
        )
        live = {"search": {"payload": {"completed": {0: result}}}}
        encoded = {"search": {"payload": {"completed": {"0": result}}}}
        self.assertEqual(signature(live, 35, 1), signature(encoded, 35, 1))

    def test_alias_positions_cannot_hide_a_duplicate_result(self):
        result = dict(
            id="curve",
            sigma=17,
            factor=5,
            reason="complete",
            phase="stage_two",
            chain_chunks=0,
        )
        for rows in ({0: result, "0": result}, {"00": result}):
            checkpoint = {"search": {"payload": {"completed": rows}}}
            with self.assertRaises(ValueError):
                signature(checkpoint, 35, 1)

    def test_jit_machine_pc_does_not_hide_its_logical_python_frame(self):
        class Assembler(int):
            pass

        profiles = [([1, 2, Assembler(99)], 4, 0, 0), ([1, 1], 1, 0, 0)]
        summary = sample_summary(
            profiles, {1: "caller", 2: "kernel"}, Assembler
        )
        self.assertEqual(
            summary["leaf"][0], dict(function="kernel", percent=80)
        )
        self.assertEqual(summary["jit_compiled_percent"], 80)
        self.assertEqual(summary["unmapped_logical_leaf_percent"], 0)
        self.assertEqual(
            summary["inclusive"][0], dict(function="caller", percent=100)
        )

    def test_unmapped_logical_frames_stay_visible(self):
        class Assembler(int):
            pass

        summary = sample_summary(
            [([12, Assembler(99)], 1, 0, 0)], {}, Assembler
        )
        self.assertEqual(summary["unmapped_logical_leaf_percent"], 100)

    def test_component_count_labels_cannot_silently_select_serial(self):
        self.assertFalse(component_arm("M0", 1))
        for workers in (4, 6, 8, 10, 12, 14, 16):
            self.assertTrue(component_arm("P" + str(workers), workers))
        for arm, workers in (("S0", 1), ("P8", 4), ("M0", 8)):
            with self.assertRaises(ValueError):
                component_arm(arm, workers)

    def test_refused_first_call_records_elapsed_time_and_closes(self):
        closed = []

        class RefusedExecutor:
            lifetime_stats = {}

            def __init__(self, **kwargs):
                pass

            def run(self, *args, **kwargs):
                return SimpleNamespace(
                    stats={"termination_reason": "cpu_limit"},
                    reason="cpu_limit",
                    checkpoint={},
                )

            def close(self):
                closed.append(True)

        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            request = dict(
                observer="none",
                capture=str(root / "capture"),
                config={"ecm_tiers": [[2, 3, 1]]},
                settings={"workers": 1},
                seed=7,
                arm="S1",
                n=35,
                wall=100,
                cpu=100,
                entry=1_000_000,
                mode="cold",
                repetitions=1,
            )
            path = root / "request.json"
            path.write_text(json.dumps(request))
            ticks = itertools.count(1_000_000, 0.001)
            original_stat = profile_runner.os.stat
            lock_identity = SimpleNamespace(st_dev=1, st_ino=2)

            def stat(name, *args, **kwargs):
                if str(name) == "/private/tmp/factor-performance.lock":
                    return lock_identity
                return original_stat(name, *args, **kwargs)

            # This is a fake resource refusal: no curve or process is run.
            with contextlib.ExitStack() as stack:
                stack.enter_context(
                    patch(
                        "v2.execution.ecm_workers.ECMExecutor", RefusedExecutor
                    )
                )
                stack.enter_context(
                    patch.dict("os.environ", FACTOR_DIAGNOSTIC_LOCK_FD="1")
                )
                stack.enter_context(
                    patch.object(profile_runner.os, "stat", side_effect=stat)
                )
                stack.enter_context(
                    patch.object(
                        profile_runner.os, "fstat", return_value=lock_identity
                    )
                )
                stack.enter_context(
                    patch.object(profile_runner.fcntl, "flock")
                )
                stack.enter_context(
                    patch.object(profile_runner.resource, "setrlimit")
                )
                stack.enter_context(
                    patch.object(
                        profile_runner.time,
                        "monotonic",
                        side_effect=lambda: next(ticks),
                    )
                )
                stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
                profile_runner.child(path)

            data = json.loads((root / "capture/result.json").read_text())
            self.assertFalse(data["complete"])
            self.assertEqual(data["stop_reason"], "cpu_limit")
            self.assertEqual(data["calls"], [])
            self.assertGreater(data["measured_wall"], 0)
            self.assertLess(data["measured_wall"], 1)
            self.assertEqual(closed, [True])

    def test_full_heap_histogram_reconciles_the_displayed_tail(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "many-types.heap"
            words = [0, 0, 0, -1]
            for kind in range(30):
                words.extend([kind + 1, kind, 32 + kind, -1])
            path.write_bytes(array.array("l", words).tobytes())

            short = heap_summary(path)
            complete = heap_summary(path, limit=None)

            self.assertEqual(short["type_count"], 30)
            self.assertEqual(short["omitted_gc_bytes"], 170)
            self.assertEqual(len(complete["types"]), 30)
            self.assertEqual(complete["omitted_gc_bytes"], 0)
            self.assertEqual(complete["gc_object_bytes"], 1395)
            self.assertEqual(
                sum(row["bytes"] for row in complete["types"]), 1395
            )

    def test_full_cpu_histogram_reconciles_samples_outside_top_entries(self):
        class Assembler(int):
            pass

        profiles = [([kind], kind, 0, 0) for kind in range(1, 31)]
        symbols = {kind: str(kind) for kind in range(1, 31)}

        short = sample_summary(profiles, symbols, Assembler)
        complete = sample_summary(profiles, symbols, Assembler, limit=None)

        self.assertEqual(short["leaf_count"], 30)
        self.assertEqual(short["omitted_leaf_samples"], 15)
        self.assertEqual(len(complete["leaf"]), 30)
        self.assertEqual(complete["weighted_samples"], 465)
        self.assertEqual(complete["omitted_leaf_samples"], 0)
        self.assertAlmostEqual(
            sum(row["percent"] for row in complete["leaf"]), 100
        )
