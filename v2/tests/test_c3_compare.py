"""Substantive admission, interruption and matched-capture invariants."""

import contextlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.ecm.c3 import round2_compare as compare


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name)
        self.control = self.output / "control.json"
        self.control.write_text("{}")
        control_patch = patch.object(compare, "CONTROL", self.control)
        control_patch.start()
        self.addCleanup(control_patch.stop)
        quiet_patch = patch.object(compare, "check_quiet")
        quiet_patch.start()
        self.addCleanup(quiet_patch.stop)
        self.fixture = dict(id="independent", n=10**39 + 1)
        self.calls = []

    def run_arm(self, fixture, seed, arm):
        self.calls.append(arm)
        return dict(arm=arm, instrumented=False, reason="complete", cpu=0.1)

    def capture(self, allowance):
        return compare.capture_group(
            self.output,
            self.fixture,
            0,
            17,
            compare.ARMS,
            allowance,
            self.run_arm,
        )

    def test_fifty_seconds_admits_no_arm_of_ninety_three_second_group(self):
        allowance = compare.Allowance(self.output / "ledger.json", 50)
        self.assertFalse(self.capture(allowance))
        self.assertEqual(self.calls, [])
        self.assertFalse((self.output / "00-independent").exists())

    def test_all_arms_share_lease_and_completed_group_never_reruns(self):
        allowance = compare.Allowance(self.output / "ledger.json", 120)
        self.assertTrue(self.capture(allowance))
        self.assertTrue(self.capture(allowance))
        self.assertEqual(self.calls, list(compare.ARMS))
        rows = compare.read_group(
            self.output / "00-independent", self.fixture, 0, 17
        )
        self.assertEqual({r["lease"] for r in rows}, {allowance.lease})
        self.assertNotIn("pending", allowance.ledger)

    def test_interrupted_group_keeps_grant_and_cannot_resume(self):
        allowance = compare.Allowance(self.output / "ledger.json", 120)
        self.assertTrue(allowance.reserve("group", 93, "00-independent"))
        value = json.loads((self.output / "ledger.json").read_text())
        self.assertGreaterEqual(value["active_seconds"], 93)
        with self.assertRaisesRegex(ValueError, "interrupted"):
            compare.Allowance(self.output / "ledger.json", 120)

    def test_interrupted_warmup_is_charged_before_reacquisition(self):
        allowance = compare.Allowance(self.output / "ledger.json", 120)
        allowance.reserve("warmup", 60, "40:fitted")
        resumed = compare.Allowance(self.output / "ledger.json", 120)
        self.assertGreaterEqual(resumed.spent, 60)
        self.assertTrue(resumed.reserve("warmup", 60, "40:fitted"))
        self.assertGreaterEqual(resumed.ledger["active_seconds"], 120)

    def test_partial_capture_is_inconclusive_even_if_ledger_was_lost(self):
        (self.output / "00-independent").mkdir()
        with self.assertRaisesRegex(ValueError, "interrupted"):
            self.capture(compare.Allowance(self.output / "ledger.json", 120))
        self.assertEqual(self.calls, [])

    def test_row_and_assignment_tampering_are_rejected(self):
        self.capture(compare.Allowance(self.output / "ledger.json", 120))
        path = self.output / "00-independent"
        with self.assertRaisesRegex(ValueError, "assignment"):
            compare.read_group(path, self.fixture, 0, 43)
        (path / "fitted.json").write_text("{}")
        with self.assertRaisesRegex(ValueError, "row changed"):
            compare.read_group(path, self.fixture, 0, 17)

    def test_failure_cannot_be_discarded_into_acceptance(self):
        allowance = compare.Allowance(self.output / "ledger.json", 120)
        allowance.reserve("group", 93, "00-independent")
        allowance.fail(ValueError("factor reconstruction failed"))
        with self.assertRaisesRegex(ValueError, "failed comparison"):
            compare.Allowance(self.output / "ledger.json", 120)

    def test_failed_grant_is_retained_without_double_charging_elapsed(self):
        clock = [0.0]
        with patch.object(compare.time, "monotonic", lambda: clock[0]):
            allowance = compare.Allowance(self.output / "ledger.json", 120)
            allowance.reserve("group", 93, "00-independent")
            clock[0] = 20
            allowance.fail(ValueError("interrupted"))
            allowance.finish()
        self.assertEqual(allowance.ledger["active_seconds"], 93)
        self.assertEqual(allowance.ledger["leases"][0]["wall"], 20)

    def test_measurement_preserves_full_report_grant(self):
        ledger = dict(
            protocol_sha256=compare.first.digest(self.control),
            active_seconds=3350,
            leases=[],
            warmups=[],
        )
        (self.output / "ledger.json").write_text(json.dumps(ledger))
        allowance = compare.Allowance(
            self.output / "ledger.json",
            1200,
            phase_reserve=compare.REPORT_GRANT + 3,
        )
        self.assertFalse(allowance.fits(93))
        self.assertTrue(allowance.fits(60))

    def test_both_bands_and_forced_relation_paths_are_warmed(self):
        fixtures = [
            dict(id=f"r2_{band}_balanced_0", n=10 ** (band - 1) + 1)
            for band in (30, 40)
        ]
        allowance = compare.Allowance(self.output / "ledger.json", 1200)
        clock = [0.0]

        def tick():
            return clock[0]

        def run(*_args, **_kwargs):
            clock[0] += 4
            return dict(complete=True)

        with (
            patch.object(compare.time, "monotonic", tick),
            patch.object(compare.time, "process_time", tick),
            patch.object(compare.training, "run_one", run),
            patch.object(
                compare, "deadline", lambda _s: contextlib.nullcontext()
            ),
        ):
            # Start the synthetic clock at the lease's own origin.
            allowance.started = 0
            self.assertTrue(compare.warmup(fixtures, allowance, run, None))
        self.assertEqual(
            {(r["band"], r["arm"]) for r in allowance.ledger["warmups"]},
            {
                (band, arm)
                for band in (30, 40)
                for arm in (*compare.ARMS, "old_no_ecm", "no_ecm")
            },
        )
        self.assertTrue(
            all(r["validated"] >= 1 for r in allowance.ledger["warmups"])
        )
