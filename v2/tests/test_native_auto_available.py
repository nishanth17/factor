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
    def test_borrowed_80_digit_owner_retains_unadopted_python_route(self):
        corpus = ROOT / (
            "v2/benchmarks/inputs/corpora/b4a_production_3p_v1.json"
        )
        n = next(
            f["n"]
            for f in json.loads(corpus.read_text())["fixtures"]
            if f["digits"] == 80 and f["population"] == "failed"
        )
        with AutomaticECMOwner(seed=73) as owner:
            ledger = Budget(seconds=10, cpu_seconds=70)
            try:
                run = factorize_bounded(
                    n,
                    config=PortfolioConfig(siqs=None, sss=None),
                    budget=ledger,
                    ecm_executor=owner,
                    stop_after_split=True,
                )
            finally:
                ledger.close()
            self.assertEqual(run.result.reconstruct(), n)
            self.assertEqual(run.ecm_execution["backend"], "python")
            self.assertEqual(run.ecm_execution["native_entries"], 0)
            self.assertEqual(
                run.ecm_execution["jobs"][0]["reason"],
                "unconfirmed_retained_native_size",
            )
            self.assertIsNone(owner.native_context)
            with self.assertRaisesRegex(ValueError, "require-native"):
                require_native(run)

    def test_expected_unavailability_falls_back_before_native_arithmetic(self):
        from v2.ecm.native_admission import NativeUnavailableError

        with patch(
            "v2.ecm.native_admission.admit_team",
            side_effect=NativeUnavailableError(
                "native_unavailable", "missing build"
            ),
        ):
            ledger = Budget(seconds=5, cpu_seconds=40)
            try:
                run = factorize_bounded(
                    N, seed=7919, budget=ledger, stop_after_split=True
                )
                self.assertEqual(run.result.reconstruct(), N)
                self.assertEqual(run.ecm_execution["backend"], "python")
                self.assertEqual(run.ecm_execution["native_entries"], 0)
                with self.assertRaisesRegex(ValueError, "require-native"):
                    require_native(run)
            finally:
                ledger.close()

    def test_allocation_fault_does_not_rerun_the_python_backend(self):
        from v2.execution.resources import AllocationFailureError

        ledger = Budget(seconds=5)
        try:
            with (
                patch(
                    "v2.ecm.native_admission.admit_team",
                    side_effect=MemoryError("partial native allocation"),
                ),
                patch.object(
                    AutomaticECMOwner,
                    "_run_incumbent",
                    side_effect=AssertionError("silent backend rerun"),
                ),
            ):
                with self.assertRaises(AllocationFailureError):
                    factorize_bounded(N, seed=7919, budget=ledger)
        finally:
            ledger.close()

    def test_automatic_production_runner_retains_and_closes_the_actual_owner(
        self,
    ):
        from v2.benchmarks.ecm.b4a.production_batch import ProductionSession

        request = dict(
            role="automatic",
            n=N,
            seed=7919,
            seconds=5,
            cpu_seconds=30,
            ecm_tiers=[
                [200, 7700, 10],
                [2000, 81000, 30],
                [10000, 554000, 100],
            ],
            production_defaults=True,
            stop_after_split=True,
            require_native=True,
            expected_workers=8,
        )
        session = ProductionSession()
        try:
            first = session.call(request)
            second = session.call(request)
            for row in (first, second):
                self.assertNotIn("native_assertion_failure", row)
                self.assertEqual(row["ecm_execution"]["backend"], "c-gmp-full")
            self.assertGreater(second["before_state"]["live_contexts"], 0)
        finally:
            session.close()
        self.assertTrue(session.owner.closed)
        self.assertEqual(session.state()["live_contexts"], 0)

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
