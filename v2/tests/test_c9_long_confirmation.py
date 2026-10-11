"""Duration amendments must reach cold workers without leaking study state."""

import unittest
from pathlib import Path
from subprocess import CompletedProcess
from unittest.mock import patch

from v2.benchmarks.qs.c9 import c9_growth_confirmation as original
from v2.benchmarks.qs.c9 import c9_long_confirmation as longer


class LongerConfirmationTests(unittest.TestCase):
    def test_scoped_amendment_restores_original_identity(self):
        before = (original.SELECTED, original.CORPUS, original.CAPS)
        with self.assertRaisesRegex(RuntimeError, "interrupt"):
            with longer.amended_protocol():
                self.assertEqual(original.CAPS[50], 1200)
                self.assertEqual(original.SELECTED, longer.SELECTED)
                self.assertEqual(original.CORPUS, longer.CORPUS)
                raise RuntimeError("interrupt")
        self.assertEqual(
            (original.SELECTED, original.CORPUS, original.CAPS), before
        )
        self.assertEqual(original.CAPS[50], 360)

    def test_cold_child_uses_amended_entry_point_and_allowance(self):
        completed = CompletedProcess([], 0, stdout="", stderr="")
        with longer.amended_protocol():
            with patch.object(
                original.subprocess, "run", return_value=completed
            ) as run:
                with patch.object(original, "save"):
                    original.cold(
                        {"id": "fixture", "digits": 50}, "dlp", Path("unused")
                    )
        args, kwargs = run.call_args
        self.assertIn(longer.__spec__.name, args[0])
        self.assertEqual(kwargs["timeout"], 1230)
