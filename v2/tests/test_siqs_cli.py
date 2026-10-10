"""QS/MPQS/SIQS CLI selection, recursion and shared checkpoint limits."""

import ast
import io
import json
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

from v2.factor import main

ROOT = Path(__file__).resolve().parents[2]
NUMBER = 100003 * 100019


def cli(*arguments):
    return subprocess.run(
        [sys.executable, "-B", "-m", "v2.factor", *map(str, arguments)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=30,
    )


class SIQSCLITests(unittest.TestCase):
    def test_interactive_and_direct_script_entries_use_siqs(self):
        commands = (
            [sys.executable, "-B", "-m", "v2.factor", "--method", "siqs"],
            [
                sys.executable,
                "-B",
                str(ROOT / "v2/factor.py"),
                str(NUMBER),
                "--method",
                "siqs",
            ],
        )

        for command in commands:
            with self.subTest(command=command):
                run = subprocess.run(
                    command,
                    input=f"{NUMBER}\n",
                    cwd=ROOT,
                    capture_output=True,
                    text=True,
                    timeout=30,
                )

                self.assertEqual(run.returncode, 0, run.stderr)
                self.assertIn(f"{NUMBER} = 100003^1 * 100019^1", run.stdout)

    def test_selected_qs_modes_retain_sign_and_multiplicity(self):
        number = -8 * NUMBER**2

        for mode in ("qs", "mpqs", "siqs"):
            with self.subTest(mode=mode):
                run = cli(
                    number,
                    "--method",
                    mode,
                    "--qs-base-bound",
                    400,
                    "--seed",
                    7,
                )

                self.assertEqual(run.returncode, 0, run.stderr)
                self.assertEqual(
                    run.stdout.strip(),
                    f"{number} = -1 * 2^3 * 100003^2 * 100019^2",
                )

    def test_selected_qs_modes_factor_a_composite_child(self):
        number = NUMBER * 100043
        for mode in ("qs", "mpqs", "siqs"):
            with (
                self.subTest(mode=mode),
                tempfile.TemporaryDirectory() as directory,
            ):
                path = Path(directory) / "progress.json"
                # Fixed-polynomial QS needs more positions for this fixture;
                # family engines can collect across their default schedule.
                options = ["--qs-half-width", 8192] if mode == "qs" else []

                run = cli(
                    number,
                    "--method",
                    mode,
                    *options,
                    "--seed",
                    7,
                    "--checkpoint",
                    path,
                )
                events = json.loads(path.read_text())["payload"]["state"][
                    "events"
                ]

                self.assertEqual(run.returncode, 0, run.stderr)
                self.assertEqual(
                    run.stdout.strip(),
                    f"{number} = 100003^1 * 100019^1 * 100043^1",
                )
                qs_events = [
                    event for event in events if event["stage"] == "siqs"
                ]
                self.assertGreaterEqual(len(qs_events), 2)
                self.assertTrue(
                    any(event["n"] != number for event in qs_events)
                )

    def test_siqs_fallback_follows_exhausted_pretests(self):
        # Force pretest exhaustion without replacing SIQS or child factoring.
        def exhaust(job, budget, context, config):
            budget.consume()
            job.update(done=True, factor=None)

        arguments = [
            "factor",
            str(NUMBER),
            "--siqs",
            "--ecm-curves",
            "1",
            "--siqs-base-bound",
            "400",
            "--seed",
            "7",
            "--verbose",
        ]
        output = io.StringIO()

        with patch.object(sys, "argv", arguments), redirect_stdout(output):
            with patch("v2.portfolio.advance_job", side_effect=exhaust):
                status = main()

        lines = output.getvalue().splitlines()
        events = [
            ast.literal_eval(line) for line in lines if line.startswith("{")
        ]
        stages = [event["stage"] for event in events]
        self.assertEqual(status, 0)
        self.assertEqual(stages, ["rho"] * 4 + ["pm1", "ecm", "siqs"])
        self.assertEqual(events[-1]["outcome"], "factor")
        self.assertEqual(lines[-1], f"{NUMBER} = 100003^1 * 100019^1")

    def test_checkpoint_restores_each_qs_mode_and_rejects_changed_search(self):
        for mode in ("qs", "mpqs", "siqs"):
            with (
                self.subTest(mode=mode),
                tempfile.TemporaryDirectory() as directory,
            ):
                path = Path(directory) / "progress.json"
                selection = [
                    "--method",
                    mode,
                    "--qs-base-bound",
                    400,
                    "--seed",
                    7,
                ]

                first = cli(
                    NUMBER,
                    *selection,
                    "--work-limit",
                    130000,
                    "--checkpoint",
                    path,
                )
                checkpoint = json.loads(path.read_text())
                mismatched = cli(
                    "--method",
                    mode,
                    "--qs-base-bound",
                    401,
                    "--resume",
                    path,
                )
                resumed = cli(
                    *selection, "--resume", path, "--checkpoint", path
                )
                restored = json.loads(path.read_text())

                self.assertEqual(first.returncode, 1, first.stderr)
                self.assertIn("unresolved", first.stdout)
                self.assertEqual(
                    checkpoint["payload"]["state"]["current"]["stage"], "siqs"
                )
                self.assertEqual(
                    checkpoint["payload"]["config"]["siqs"]["mode"], mode
                )
                self.assertEqual(mismatched.returncode, 2)
                self.assertEqual(resumed.returncode, 0, resumed.stderr)
                self.assertIn("100003^1 * 100019^1", resumed.stdout)
                self.assertGreater(
                    restored["payload"]["work_used"],
                    checkpoint["payload"]["work_used"],
                )
                self.assertEqual(
                    restored["payload"]["config"],
                    checkpoint["payload"]["config"],
                )

    def test_explicit_qs_modes_report_the_selected_engine(self):
        for mode in ("qs", "mpqs", "siqs"):
            with self.subTest(mode=mode):
                run = cli(
                    NUMBER,
                    "--method",
                    mode,
                    "--qs-base-bound",
                    400,
                    "--seed",
                    7,
                    "--verbose",
                )

                self.assertEqual(run.returncode, 0, run.stderr)
                self.assertIn(f"'stage': '{mode}'", run.stdout)
                self.assertIn(f"{NUMBER} = 100003^1 * 100019^1", run.stdout)

    def test_help_explains_explicit_methods_and_resume_controls(self):
        run = cli("--help")

        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertIn("{auto,qs,mpqs,siqs,sss,sssf}", run.stdout)
        self.assertIn("--qs-family-count", run.stdout)
        self.assertIn("--qs-dlp", run.stdout)
        self.assertIn("--checkpoint", run.stdout)
        self.assertIn("--resume", run.stdout)
        self.assertIn("QS uses one polynomial", " ".join(run.stdout.split()))

    def test_invalid_or_unused_siqs_settings_are_explicit_errors(self):
        for options in (
            ("--siqs", "--method", "sss"),
            ("--siqs-base-bound", 400),
            ("--method", "siqs", "--siqs-half-width", 0),
            ("--method", "siqs", "--memory-mib", 8),
            ("--method", "qs", "--qs-assignment-policy", "nearest"),
            ("--method", "mpqs", "--qs-polynomials-per-family", 2),
            ("--method", "siqs", "--qs-dlp"),
            ("--method", "siqs", "--qs-large-prime-bound", 10000),
            (
                "--method",
                "qs",
                "--qs-dlp",
                "--qs-large-prime-bound",
                10000,
                "--qs-large-product-bound",
                50000000,
            ),
        ):
            with self.subTest(options=options):
                run = cli(NUMBER, *options)

                self.assertEqual(run.returncode, 2)
                self.assertIn("error:", run.stderr)
                self.assertNotIn("Traceback", run.stderr)

    def test_dlp_flag_serializes_explicit_bounds_and_replays(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "dlp.json"
            options = (
                "--method",
                "siqs",
                "--qs-dlp",
                "--qs-residual-bound",
                1000,
                "--qs-large-prime-bound",
                10000,
                "--qs-large-product-bound",
                50000000,
                "--qs-dlp-candidate-bound",
                50000000,
                "--qs-dlp-split-call-limit",
                64,
                "--work-limit",
                0,
            )
            first = cli(NUMBER, *options, "--checkpoint", path)

            self.assertEqual(first.returncode, 1, first.stderr)
            checkpoint = json.loads(path.read_text())
            collector = checkpoint["payload"]["config"]["siqs"]["collector"]
            self.assertEqual(collector["residual_bound"], 1000)
            self.assertEqual(collector["large_prime_bound"], 10000)
            self.assertEqual(collector["large_product_bound"], 50000000)
            self.assertEqual(collector["candidate_bound"], 50000000)
            self.assertEqual(collector["split_call_limit"], 64)

            resumed = cli(*options, "--resume", path)
            self.assertEqual(resumed.returncode, 1, resumed.stderr)
            self.assertNotIn("Traceback", resumed.stderr)

    def test_bounded_prac_default_keeps_bounds_and_optional_qs_disabled(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "progress.json"

            run = cli(
                NUMBER,
                "--bounded",
                "--work-limit",
                0,
                "--checkpoint",
                path,
            )
            config = json.loads(path.read_text())["payload"]["config"]

            self.assertEqual(run.returncode, 1, run.stderr)
            self.assertIsNone(config["siqs"])
            self.assertIsNone(config["sss"])
            self.assertEqual(config["memory_bytes"], 16 * 1024 * 1024)
            self.assertEqual(config["ecm_chain_mode"], "reuse")
            self.assertEqual(config["ecm_tiers"], [[2000, 147396, 32]])


if __name__ == "__main__":
    unittest.main()
