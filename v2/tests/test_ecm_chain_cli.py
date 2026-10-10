"""Default/optional bounded ECM flags and preserved saved policy."""

import json
import tempfile
import unittest
from pathlib import Path

from v2.tests.test_siqs_cli import cli


class ChainCLITests(unittest.TestCase):
    def test_default_and_optional_flags_pin_policy(self):
        for choice, version in (
            (None, 10),
            ("prac", 10),
            ("lucas", 11),
            ("cf", 11),
            ("off", 9),
        ):
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "progress.json"
                arguments = [] if choice is None else ["--ecm-chain", choice]
                run = cli(
                    1009 * 1013, "--bounded", "--checkpoint", path, *arguments
                )
                self.assertEqual(run.returncode, 0, run.stderr)
                payload = json.loads(path.read_text())["payload"]
                self.assertEqual(payload["version"], version)
                self.assertIn("1009^1 * 1013^1", run.stdout)
                if choice in ("lucas", "cf"):
                    self.assertEqual(
                        payload["config"]["ecm_chain_family"], choice
                    )

    def test_resume_keeps_new_and_old_memory_and_execution(self):
        for choice, memory, version in (
            (None, 16, 10),
            ("off", 8, 9),
            ("cf", 16, 11),
        ):
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "progress.json"
                arguments = [] if choice is None else ["--ecm-chain", choice]
                if choice == "off":
                    arguments += ["--memory-mib", 8]
                first = cli(
                    1009 * 1013,
                    "--bounded",
                    "--work-limit",
                    0,
                    "--checkpoint",
                    path,
                    *arguments,
                )
                self.assertEqual(first.returncode, 1, first.stderr)
                resumed = cli("--resume", path, "--checkpoint", path)
                self.assertEqual(resumed.returncode, 0, resumed.stderr)
                payload = json.loads(path.read_text())["payload"]
                self.assertEqual(payload["version"], version)
                self.assertEqual(
                    payload["config"]["memory_bytes"], memory * 1024**2
                )
                self.assertIn("1009^1 * 1013^1", resumed.stdout)

    def test_explicit_small_cap_and_unsupported_schedule_keep_ladder(self):
        for arguments in (("--memory-mib", 8), ("--ecm-curves", 2)):
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "progress.json"
                run = cli(101, "--bounded", "--checkpoint", path, *arguments)
                self.assertEqual(run.returncode, 0, run.stderr)
                payload = json.loads(path.read_text())["payload"]
                self.assertEqual(payload["version"], 9)
                self.assertEqual(
                    payload["config"]["memory_bytes"], 8 * 1024**2
                )
                self.assertNotIn("chains", payload)

    def test_unrelated_engine_and_native_prac_backend_are_rejected(self):
        for arguments in (
            ("--method", "siqs", "--ecm-chain", "cf"),
            ("--backend", "gmpy2-mpz", "--ecm-chain", "prac"),
        ):
            run = cli(101, *arguments)
            self.assertEqual(run.returncode, 2)
            self.assertIn("--ecm-chain", run.stderr)


if __name__ == "__main__":
    unittest.main()
