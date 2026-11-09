"""Observation requirements must fail before a window or actor is started."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.qs.a7a import cell, combined
from v2.benchmarks.qs.a7a.protocol import validate_worker_timeline


class TimelinePolicyTests(unittest.TestCase):
    def test_required_stages_and_serial_exemptions(self):
        for stage in ("G2H", "G3"):
            for arm in ("C1", "C4", "P4"):
                for missing in (None, "", " ", False):
                    with self.subTest(stage=stage, arm=arm, path=missing):
                        with self.assertRaises(ValueError):
                            validate_worker_timeline(
                                stage, dict(arm=arm, timeline_path=missing)
                            )
                validate_worker_timeline(
                    stage, dict(arm=arm, timeline_path="timelines/cell")
                )
            for arm in ("S0", "M0", "MAINLINE"):
                validate_worker_timeline(stage, dict(arm=arm))

        for arm in ("C4", "S0"):
            validate_worker_timeline("G4", dict(arm=arm))
            with self.assertRaises(ValueError):
                validate_worker_timeline(
                    "G4", dict(arm=arm, timeline_path="timelines/cell")
                )
        validate_worker_timeline("G2", dict(arm="C4"))

    def test_direct_cell_refuses_before_configuration_or_pool_setup(self):
        with patch.object(cell, "configuration") as configure:
            for stage, timeline in (("G3", None), ("G4", "timelines")):
                with self.assertRaises(ValueError):
                    cell.execute(
                        dict(
                            n=15,
                            cell=dict(
                                stage=stage, arm="C4", timeline_path=timeline
                            ),
                        )
                    )
            configure.assert_not_called()

    def test_combined_refuses_before_creating_output_or_opening_lock(self):
        with tempfile.TemporaryDirectory() as folder:
            plan = Path(folder) / "plan.json"
            output = Path(folder) / "receipts"
            plan.write_text(
                json.dumps(
                    dict(
                        stage="G3",
                        source="pin",
                        corpus_kind="regimes",
                        corpus_sha256="corpus",
                        cells=[dict(arm="C4")],
                    )
                )
            )
            with (
                patch.object(combined, "source", return_value="pin"),
                patch.object(
                    combined, "corpus_identity", return_value="corpus"
                ),
                patch.object(combined, "open_lock") as open_lock,
            ):
                with self.assertRaises(ValueError):
                    combined.run(plan, output)
                open_lock.assert_not_called()
            self.assertFalse(output.exists())
