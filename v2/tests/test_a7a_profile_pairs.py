"""Paired admission and immutable, pause-separated CPU/heap captures."""

import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.qs.a7a import diagnostic_hook
from v2.benchmarks.qs.a7a.combined import cell, corpus_identity, run
from v2.benchmarks.qs.a7a.profile_pairs import paired_cells


class ProfilePairTests(unittest.TestCase):
    def test_every_clean_configuration_has_an_unchanged_diagnostic_peer(self):
        original = cell("C4", variant="case", samples=9)
        clean, diagnostic = paired_cells([original])
        self.assertTrue(clean["profile_required"])
        self.assertEqual(clean["pair"], diagnostic["diagnostic_for"])
        for key in (
            "fixture",
            "arm",
            "features",
            "seed",
            "work",
            "memory",
            "lifetime",
        ):
            self.assertEqual(clean[key], diagnostic[key])
        self.assertEqual(diagnostic["mode"], "profile")
        self.assertTrue(diagnostic["full_profile"])
        self.assertTrue(diagnostic["observe"])
        self.assertEqual((clean["samples"], diagnostic["samples"]), (9, 9))
        self.assertNotIn("pair", original)

    def test_clean_run_is_not_issued_without_the_whole_pair_reservation(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            cells = paired_cells([cell("C4", wall=3, cpu=30)])
            plan = dict(
                identity="reservation",
                source="test",
                wall=9,
                cpu=1000,
                storage=4 * 2**30,
                prior_charges=[],
                cells=cells,
                corpus_sha256=corpus_identity(),
                profile_folder=str(folder / "profiles"),
            )
            path = folder / "plan.json"
            path.write_text(json.dumps(plan))
            lock = (folder / "lock").open("w+")

            def archive(command, **kwargs):
                if command[:2] == ["git", "archive"]:
                    Path(command[4]).write_bytes(b"source")

            with (
                patch(
                    "v2.benchmarks.qs.a7a.combined.source", return_value="test"
                ),
                patch(
                    "v2.benchmarks.qs.a7a.combined.open_lock",
                    return_value=lock,
                ),
                patch(
                    "v2.benchmarks.qs.a7a.combined.subprocess.run",
                    side_effect=archive,
                ),
                patch(
                    "v2.benchmarks.qs.a7a.combined.subprocess.Popen"
                ) as launch,
            ):
                run(path, folder / "results")
            launch.assert_not_called()
            ledger = json.loads((folder / "results/ledger.json").read_text())
            self.assertEqual(ledger["attempts"], [])
            self.assertIn("reservation", ledger["stopped"])

    def test_heap_pause_splits_cpu_streams_without_overwriting_prior_data(
        self,
    ):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            profiler = SimpleNamespace(
                enable=lambda *a: None,
                disable=lambda: None,
                is_enabled=lambda: False,
            )

            def dump(path):
                Path(path).write_bytes(b"heap")

            with (
                patch.multiple(
                    diagnostic_hook,
                    _folder=folder,
                    _role="parent",
                    _profile=None,
                    _cpu_file_limit=None,
                    _events=[],
                    _cpu_segments=0,
                    _closed=False,
                ),
                patch.dict(sys.modules, {"_vmprof": profiler}),
                patch.dict(
                    "os.environ",
                    {
                        "FACTOR_A7A_PROFILE_ALL": "1",
                        "FACTOR_A7A_LIVE_HEAPS": "0",
                    },
                ),
                patch.object(
                    diagnostic_hook.resource,
                    "getrlimit",
                    return_value=(2**31, 2**31),
                ),
                patch.object(diagnostic_hook.resource, "setrlimit"),
                patch.object(
                    diagnostic_hook.gc,
                    "dump_rpy_heap",
                    side_effect=dump,
                    create=True,
                ),
                patch.object(
                    diagnostic_hook.gc,
                    "get_stats",
                    return_value="GC",
                    create=True,
                ),
            ):
                diagnostic_hook._cpu_signal(None, None)
                diagnostic_hook._heap_signal(None, None)
                diagnostic_hook._cpu_signal(None, None)
                kinds = [row["kind"] for row in diagnostic_hook._events]
                self.assertEqual(
                    kinds,
                    [
                        "cpu_started",
                        "cpu_stopped",
                        "heap_finished",
                        "cpu_started",
                        "cpu_stopped",
                    ],
                )
            self.assertTrue((folder / "cpu.vmprof").exists())
            self.assertTrue((folder / "cpu-002.vmprof").exists())
            self.assertEqual((folder / "live.heap").read_bytes(), b"heap")


if __name__ == "__main__":
    unittest.main()
