"""Host-audit charges, native refusals and fixed-work forwarding."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.qs.c10a01 import regression65_recheck as recheck


class HostAuditedRunnerTests(unittest.TestCase):
    def test_canary_is_charged_separately_from_siqs_work(self):
        with (
            patch.object(
                recheck, "ORIGINAL_SENTINEL", return_value={"work": 70}
            ),
            patch.object(
                recheck.host_state,
                "power_state",
                return_value={"low_power": 1},
            ),
            patch.object(
                recheck.host_state, "runtime_state", return_value={"qos": 33}
            ),
            patch.object(
                recheck.host_state,
                "native_canary",
                return_value={"work": 1800000},
            ),
        ):
            value = recheck.audited_sentinel({}, {}, ())
        self.assertEqual(value["work"], 1800070)
        self.assertEqual(value["native_canary"]["work"], 1800000)
        self.assertEqual(value["host"]["runtime"]["qos"], 33)

    def test_power_mismatch_refuses_before_job_and_blocks_are_forwarded(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "freeze.json"
            path.write_text(json.dumps(dict(expected_power={"low_power": 1})))
            with (
                patch.object(recheck, "FREEZE", path),
                patch.object(
                    recheck.host_state,
                    "power_state",
                    return_value={"low_power": 0},
                ),
                patch.object(recheck, "ORIGINAL_MEASURED") as run,
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "host power changed"
                ):
                    recheck.measured({"blocks": 4224}, {}, (), "identity")
                run.assert_not_called()
            with (
                patch.object(recheck, "FREEZE", path),
                patch.object(
                    recheck.host_state,
                    "power_state",
                    return_value={"low_power": 1},
                ),
                patch.object(
                    recheck, "ORIGINAL_MEASURED", return_value={}
                ) as run,
            ):
                row = {"blocks": 4224}
                recheck.measured(row, {}, (), "identity")
                run.assert_called_once_with(
                    row, {}, (), "identity", blocks=4224
                )


if __name__ == "__main__":
    unittest.main()
