"""Relocation must preserve strict source pins and isolated benchmark arms."""

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.qs.phase_three import phase_three_audit
from v2.benchmarks.support import paths
from v2.execution.budget import Budget


class LayoutTests(unittest.TestCase):
    def test_migration_pin_requires_both_recorded_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "moved.py"
            source.write_text("current bytes")
            actual = hashlib.sha256(source.read_bytes()).hexdigest()
            manifest = root / "inputs/controls/layout_migration.json"
            manifest.parent.mkdir(parents=True)
            manifest.write_text(
                json.dumps(
                    {
                        "sources": {
                            "moved.py": {
                                "before_sha256": "original hash",
                                "after_sha256": actual,
                            }
                        }
                    }
                )
            )

            with (
                patch.object(paths, "REPOSITORY_ROOT", root),
                patch.object(paths, "BENCHMARK_ROOT", root),
            ):
                self.assertTrue(
                    paths.matches_source_pin(source, "original hash")
                )
                self.assertFalse(paths.matches_source_pin(source, "unknown"))

                source.write_text("tampered bytes")

                self.assertFalse(
                    paths.matches_source_pin(source, "original hash")
                )

    def test_existing_historical_file_takes_precedence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            original = root / "v2/utils.py"
            original.parent.mkdir(parents=True)
            original.write_text("historical source")

            with patch.object(paths, "REPOSITORY_ROOT", root):
                self.assertEqual(paths.source_path(original), original)

                original.unlink()

                self.assertEqual(
                    paths.source_path(original), root / "v2/common/utils.py"
                )

    def test_current_private_qs_arm_keeps_its_own_relation_types(self):
        arm = phase_three_audit.load_arm("_layout_test_current", current=True)
        budget = Budget(work_limit=10**8)

        result = arm.qs.factor_base.build_factor_base(
            101 * 137, bound=50, budget=budget
        )

        self.assertIsInstance(
            result.factor_base, arm.qs.factor_base.FactorBase
        )
        self.assertEqual(result.factor_base.n, 101 * 137)
        self.assertEqual(
            arm.qs.relations.Polynomial, arm.qs.polynomial.Polynomial
        )

    def test_standalone_drivers_do_not_import_a_runtime_early(self):
        for name in (
            "infrastructure/performance/performance_audit.py",
            "infrastructure/performance/performance_followup.py",
            "qs/p38/p38_r1_regression.py",
        ):
            with self.subTest(driver=name):
                source = paths.BENCHMARK_ROOT / name
                command = (
                    "import runpy, sys; "
                    "runpy.run_path(sys.argv[1], run_name='layout_probe'); "
                    "assert 'v2' not in sys.modules"
                )

                subprocess.run(
                    [sys.executable, "-B", "-c", command, str(source)],
                    check=True,
                    capture_output=True,
                    text=True,
                )


if __name__ == "__main__":
    unittest.main()
