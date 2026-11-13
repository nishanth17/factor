"""Actual default entry and replay with admitted scalar/pthread builds."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.ecm.auto import AutomaticECMOwner, _runtime_supported
from v2.ecm.execution_identity import require_native
from v2.ecm.native_admission import build_present
from v2.ecm.native_auto import SOURCE
from v2.execution.budget import Budget
from v2.portfolio import PortfolioConfig, factorize_bounded

ROOT = Path(__file__).resolve().parents[2]
N = 595063068176486113232605383931148599031685550719226328338409
AVAILABLE = _runtime_supported() and all(
    build_present(SOURCE, variant) for variant in ("scalar", "pthread")
)


@unittest.skipUnless(
    AVAILABLE, "compatible local scalar/pthread builds required"
)
class ActualAutomaticNativeTests(unittest.TestCase):
    def test_explicit_scalar_and_pthread_receipts_count_actual_entry(self):
        from v2.ecm.workers import ECMExecutor
        from v2.execution.ecm_native import NativeECMExecutor
        from v2.execution.worker_contract import WorkerSettings

        for workers in (1, 4, 8):
            settings = WorkerSettings(workers=workers, allow_wide_workers=True)
            owner = (
                ECMExecutor(seed=7919, settings=settings)
                if workers == 1
                else NativeECMExecutor(
                    seed=7919, settings=settings, capacity=32
                )
            )
            with owner:
                ledger = Budget(seconds=5, cpu_seconds=30)
                try:
                    run = factorize_bounded(
                        N,
                        seed=7919,
                        ecm_executor=owner,
                        budget=ledger,
                        config=PortfolioConfig(ecm_kernel="c-gmp-full"),
                        stop_after_split=True,
                    )
                    require_native(run, workers=workers)
                    self.assertTrue(run.ecm_execution["jobs"][0]["build"])
                finally:
                    ledger.close()

    def test_normal_default_library_actually_executes_whole_native(self):
        ledger = Budget(seconds=5, cpu_seconds=30)
        try:
            run = factorize_bounded(
                N,
                seed=7919,
                budget=ledger,
                stop_after_split=True,
            )
            self.assertEqual(run.result.reconstruct(), N)
            self.assertEqual(require_native(run)["requested"], "auto-full")
            self.assertEqual(run.ecm_execution["jobs"][0]["workers"], 8)
        finally:
            ledger.close()

    def test_scalar_only_build_preserves_manifest_and_parent_rng(self):
        real_presence = build_present
        with patch(
            "v2.ecm.native_admission.build_present",
            side_effect=lambda source, variant: (
                variant == "scalar" and real_presence(source, variant)
            ),
        ):
            ledger = Budget(seconds=5, cpu_seconds=30)
            try:
                run = factorize_bounded(
                    N,
                    seed=7919,
                    budget=ledger,
                    stop_after_split=True,
                )
                require_native(run, workers=1)
                self.assertEqual(
                    run.ecm_execution["jobs"][0]["mathematical_stream"],
                    "b12-manifest",
                )
                self.assertEqual(run.result.reconstruct(), N)
            finally:
                ledger.close()

    def test_interrupted_native_resume_keeps_build_and_search(self):
        with AutomaticECMOwner(seed=7919) as owner:

            def after_native_entry():
                selected = owner.native_context
                return bool(
                    selected is not None
                    and selected.team is not None
                    and selected.team.slices
                )

            ledger = Budget(
                seconds=5, cpu_seconds=30, cancelled=after_native_entry
            )
            try:
                first = factorize_bounded(N, budget=ledger, ecm_executor=owner)
            finally:
                ledger.close()
            self.assertFalse(owner.closed)
            self.assertEqual(first.reason, "cancelled")
            require_native(first)
        saved = first.checkpoint
        with patch(
            "v2.ecm.native_admission.admit_packet",
            side_effect=ValueError("saved build disappeared"),
        ):
            with self.assertRaisesRegex(ValueError, "disappeared"):
                factorize_bounded(N, checkpoint=saved)
        ledger = Budget(seconds=5, cpu_seconds=30)
        try:
            resumed = factorize_bounded(
                N,
                seed=999,
                budget=ledger,
                checkpoint=saved,
                stop_after_split=True,
            )
            require_native(resumed)
            self.assertEqual(resumed.result.reconstruct(), N)
            self.assertEqual(
                resumed.checkpoint["payload"]["state"]["ecm_auto"]["seed"],
                7919,
            )
        finally:
            ledger.close()

    def test_normal_cli_without_native_flag_records_actual_native_entry(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.json"
            call = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    "-m",
                    "v2.factor",
                    str(N),
                    "--seed",
                    "7919",
                    "--seconds",
                    "5",
                    "--checkpoint",
                    str(path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=10,
            )
            self.assertIn(call.returncode, (0, 1), call.stderr)
            state = json.loads(path.read_text())["payload"]["state"]
            records = state["ecm_execution"]["jobs"]
            self.assertTrue(records)
            self.assertTrue(all(r["backend"] == "c-gmp-full" for r in records))
            self.assertGreater(sum(r["native_entries"] for r in records), 0)


if __name__ == "__main__":
    unittest.main()
