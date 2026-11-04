"""Independent diagnostic-format and factor-publication checks."""

import array
import tempfile
import unittest
from pathlib import Path

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
