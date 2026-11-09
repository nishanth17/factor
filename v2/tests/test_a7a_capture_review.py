"""Bounded sampling does not imply an all-worker capture contract."""

import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.qs.a7a.capture_review import review


class CaptureReviewTests(unittest.TestCase):
    def test_explicit_scope_preserves_inventory_and_raw_capture(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            capture = root / "capture"
            for pid, role in ((11, "parent"), (12, "worker"), (13, "worker")):
                folder = capture / str(pid)
                folder.mkdir(parents=True)
                events = [dict(kind="ready", pid=pid, role=role)]
                if pid != 13:
                    events.extend(
                        dict(kind=kind)
                        for kind in ("cpu_started", "cpu_stopped")
                    )
                    (folder / "cpu.vmprof").write_bytes(b"profile")
                (folder / "events.json").write_text(json.dumps(events))
            stats = SimpleNamespace(
                profiles=[([1, 2], 3)],
                adr_dict={1: "caller", 2: "leaf"},
                n_samples=3,
            )
            reader = SimpleNamespace(read_profile=lambda path: stats)
            with patch.dict(sys.modules, vmprof=reader):
                result = review(
                    capture,
                    {},
                    source="source",
                    input_id="input",
                    expected_pids=[11, 12],
                    output=root / "scoped-review.json",
                )
                self.assertEqual(result["unreviewed_pids"], [13])
                self.assertEqual(result["expected_pids"], [11, 12])
                self.assertEqual(len(result["processes"]), 2)
                self.assertFalse((capture / "analysis.json").exists())

                # The default/full contract still rejects a missing stream.
                with self.assertRaisesRegex(ValueError, "valid flush"):
                    review(capture, {}, source="source", input_id="input")
                with self.assertRaisesRegex(ValueError, "capture processes"):
                    review(
                        capture,
                        {},
                        source="s",
                        input_id="i",
                        expected_pids=[11, 99],
                    )

                (capture / "12/cpu.vmprof").unlink()
                with self.assertRaisesRegex(ValueError, "CPU sample"):
                    review(
                        capture,
                        {},
                        source="s",
                        input_id="i",
                        expected_pids=[11, 12],
                    )


if __name__ == "__main__":
    unittest.main()
