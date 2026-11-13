"""Untimed protocol, matching, allocation and analysis checks for NEON v3."""

import copy
import unittest

from v2.benchmarks.ecm.b4a.neon_plan import amend
from v2.benchmarks.ecm.b4a.neon_screen import (
    ARMS,
    IDENTITY,
    WIDTHS,
    check_screen_cell,
    panel,
    screen_cell,
    summarize,
    validate_report,
)
from v2.benchmarks.ecm.b4a.plan import stage_a
from v2.benchmarks.ecm.b4a.runner import check_cell


def observations(neon_complete_ratio=1.2):
    result = []
    for width in WIDTHS:
        for operation in ("multiply", "square"):
            for cycle in range(9):
                for index, arm in enumerate(ARMS):
                    wall = 0.02
                    if index == 2:
                        wall = 0.005  # Arithmetic-only win is insufficient.
                    if index == 5:
                        wall *= neon_complete_ratio
                    result.append(
                        dict(
                            width=width,
                            operation=operation,
                            cycle=cycle,
                            arm=arm,
                            wall=wall,
                            cpu=wall,
                            independent_operations=0 if index == 6 else 8192,
                            pairs=4096,
                            repeats=64,
                            offset=cycle,
                            status=0,
                            validated=True,
                            qos_set=0,
                            qos_get=0,
                            requested_qos=0x21,
                            observed_qos=0x21,
                            relative_priority=0,
                            thread_id=19,
                            thread_id_status=0,
                        )
                    )
    return result


class NEONProtocolTests(unittest.TestCase):
    def test_replacement_preserves_old_freeze_and_stage_ceiling(self):
        base = dict(
            identity="old",
            commit="source",
            cells=stage_a(),
            wall=1800,
            cpu=7200,
        )
        saved = copy.deepcopy(base)
        admission = dict(
            stage="NEON-correctness-admission",
            passed=True,
            commit="source",
            timing_executed=False,
        )
        revised = amend(base, admission)
        self.assertEqual(base, saved)
        self.assertEqual(len(revised["removed_cells"]), 12)
        self.assertEqual(len(revised["cells"]), 82)
        self.assertEqual(sum(c["wall"] for c in revised["cells"]), 1695)
        self.assertEqual(sum(c["cpu"] for c in revised["cells"]), 6420)
        check_cell(screen_cell())
        for changed in (dict(wall=181), dict(cpu=361), dict(warmup=2)):
            with self.assertRaises(ValueError):
                check_screen_cell(dict(screen_cell(), **changed))

    def test_timing_requires_exact_untimed_admission(self):
        base = dict(identity="old", commit="source", cells=stage_a())
        for admission in (
            {},
            dict(passed=True, commit="other"),
            dict(passed=True, commit="source", timing_executed=True),
        ):
            with self.assertRaises(ValueError):
                amend(base, admission)

    def test_arithmetic_only_win_does_not_advance_a_slow_complete_path(self):
        report = summarize(observations(), 17)
        self.assertTrue(all(r["arithmetic_neon_over_mpn"] < 1 for r in report))
        self.assertTrue(
            all(r["complete_neon_over_production"] > 1 for r in report)
        )
        self.assertFalse(
            any(r["supports_further_integration"] for r in report)
        )
        favorable = summarize(observations(0.8), 17)
        self.assertTrue(
            all(r["supports_further_integration"] for r in favorable)
        )

    def test_mismatched_counts_or_qos_cannot_supply_clean_simd_evidence(self):
        rows = observations(0.8)
        rows[0]["independent_operations"] += 2
        with self.assertRaises(ValueError):
            summarize(rows, 17)
        rows = observations(0.8)
        rows[0]["qos_set"] = 1
        result = summarize(rows, 17)
        self.assertFalse(result[0]["supports_further_integration"])
        rows = observations(0.8)[:-1]
        with self.assertRaises(ValueError):
            summarize(rows, 17)

    def test_complete_report_reconstructs_and_requires_all_warmup_cases(self):
        rows = observations(0.8)
        warmup = [
            dict(width=w, operation=o, arm=a, validated_driver_seconds=3)
            for w in WIDTHS
            for o in ("multiply", "square")
            for a in ARMS
        ]
        report = dict(
            identity=IDENTITY,
            complete=True,
            measurements=rows,
            warmup=warmup,
            summary=summarize(rows, screen_cell()["seed"]),
        )
        validate_report(report, screen_cell())
        report["warmup"][-1] = dict(report["warmup"][0])
        with self.assertRaises(ValueError):
            validate_report(report, screen_cell())

    def test_panel_is_reproducible_varied_and_full_modulus(self):
        for width in WIDTHS:
            n, values = panel(width, 17)
            self.assertEqual((n, values), panel(width, 17))
            self.assertEqual(n.bit_length(), width)
            self.assertEqual(len(values), 64)
            self.assertEqual(len(set(values)), 64)
            self.assertTrue(
                all(0 <= value < n for row in values for value in row)
            )


if __name__ == "__main__":
    unittest.main()
