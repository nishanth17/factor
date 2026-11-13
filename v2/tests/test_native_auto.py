"""Automatic whole-native intent, execution receipts and legacy semantics."""

import copy
import random
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from v2.common.arithmetic import BackendUnavailableError
from v2.ecm.auto import AutomaticECMOwner
from v2.ecm.execution_identity import require_native
from v2.ecm.native_admission import NativeAdmission, _fingerprint, admit_packet
from v2.ecm.native_auto import eligibility
from v2.execution.budget import Budget
from v2.portfolio import PortfolioConfig, factorize_bounded


class AutomaticNativeIntentTests(unittest.TestCase):
    def test_default_intent_and_explicit_modes_are_distinct(self):
        self.assertEqual(PortfolioConfig().ecm_kernel, "auto-full")
        for mode in ("python", "auto", "c-gmp", "c-gmp-full", "auto-full"):
            self.assertEqual(PortfolioConfig(ecm_kernel=mode).ecm_kernel, mode)

    def test_small_input_and_unsupported_formula_skip_runtime_probe(self):
        config = PortfolioConfig()
        with patch(
            "v2.ecm.auto._runtime_supported", side_effect=AssertionError
        ):
            self.assertEqual(
                eligibility(1013, config), "unmeasured_native_size"
            )
            unsupported = replace(config, ecm_chain_family="lucas")
            self.assertEqual(
                eligibility(10**59 + 7, unsupported),
                "unsupported_native_formula",
            )

    def test_new_selector_keeps_the_incumbent_seed_domain(self):
        generator = random.Random(911)
        before = generator.getstate()
        native = AutomaticECMOwner(seed=None)
        python = AutomaticECMOwner(seed=None)
        native.configure(PortfolioConfig(), generator, None)
        python.configure(PortfolioConfig(ecm_kernel="python"), generator, None)
        self.assertEqual(native.seed, python.seed)
        self.assertNotEqual(native.identity, python.identity)
        self.assertEqual(generator.getstate(), before)

    def test_custom_bounds_outside_native_manifest_range_fall_back_cheaply(
        self,
    ):
        config = PortfolioConfig(
            ecm_tiers=((2000, 2**32, 8),),
            ecm_chain_mode="off",
            ecm_program_bytes=0,
        )
        with patch(
            "v2.ecm.auto._runtime_supported", side_effect=AssertionError
        ):
            self.assertEqual(
                eligibility(10**59 + 7, config),
                "unsupported_native_work_range",
            )

    def test_checkpoint_before_ecm_does_not_admit_or_hash_native(self):
        with (
            patch(
                "v2.ecm.native_admission.admit_packet",
                side_effect=AssertionError,
            ),
            patch(
                "v2.ecm.native_packet.source_identity",
                side_effect=AssertionError,
            ),
        ):
            run = factorize_bounded(1009 * 1013, budget=Budget(seconds=0))
        self.assertEqual(run.checkpoint["payload"]["version"], 14)
        self.assertNotIn("native_packet", run.checkpoint["payload"])
        self.assertEqual(run.ecm_execution["backend"], "not-engaged")

    def test_historical_omitted_kernel_keeps_python(self):
        first = factorize_bounded(
            1009 * 1013,
            config=PortfolioConfig(ecm_kernel="python"),
            budget=Budget(seconds=0),
        )
        self.assertNotIn("ecm_kernel", first.checkpoint["payload"]["config"])
        with patch(
            "v2.ecm.native_admission.admit_packet", side_effect=AssertionError
        ):
            resumed = factorize_bounded(
                1009 * 1013,
                checkpoint=first.checkpoint,
                budget=Budget(seconds=0),
            )
        self.assertNotIn("ecm_kernel", resumed.checkpoint["payload"]["config"])
        self.assertEqual(resumed.checkpoint["payload"]["version"], 13)

    def test_legacy_auto_checkpoint_keeps_stage_one_semantics(self):
        first = factorize_bounded(
            1009 * 1013,
            config=PortfolioConfig(ecm_kernel="auto"),
            budget=Budget(seconds=0),
        )
        with patch(
            "v2.ecm.native_admission.admit_packet", side_effect=AssertionError
        ):
            resumed = factorize_bounded(
                1009 * 1013,
                checkpoint=first.checkpoint,
                budget=Budget(seconds=0),
            )
        self.assertEqual(
            resumed.checkpoint["payload"]["config"]["ecm_kernel"], "auto"
        )
        self.assertEqual(resumed.checkpoint["payload"]["version"], 13)

    def test_cheap_completed_inputs_never_admit_native(self):
        with patch(
            "v2.ecm.native_admission.admit_packet", side_effect=AssertionError
        ):
            for n in (1, 101, 2**80, 1009**2):
                run = factorize_bounded(n)
                self.assertTrue(run.result.complete)
                self.assertEqual(run.result.reconstruct(), n)
                self.assertEqual(run.ecm_execution["backend"], "not-engaged")


class RequireNativeTests(unittest.TestCase):
    def native(self):
        return dict(
            backend="c-gmp-full",
            native_entries=3,
            dropped=0,
            jobs=[
                dict(
                    backend="c-gmp-full",
                    implementation="whole-native",
                    workers=4,
                    native_entries=3,
                )
            ],
        )

    def test_assertion_accepts_only_actual_expected_execution(self):
        run = SimpleNamespace(ecm_execution=self.native())
        self.assertIs(require_native(run, workers=4), run.ecm_execution)
        with self.assertRaisesRegex(ValueError, "recipe"):
            require_native(run, workers=8)

    def test_fallback_mixed_stage_one_and_no_entry_are_rejected(self):
        for changes in (
            dict(backend="python"),
            dict(backend="mixed"),
            dict(backend="c-gmp"),
            dict(native_entries=0),
            dict(dropped=1),
        ):
            with self.subTest(changes=changes):
                record = dict(self.native(), **changes)
                with self.assertRaisesRegex(ValueError, "require-native"):
                    require_native(SimpleNamespace(ecm_execution=record))
        record = copy.deepcopy(self.native())
        record["jobs"][0]["implementation"] = "stage-one-native"
        with self.assertRaisesRegex(ValueError, "implementation"):
            require_native(SimpleNamespace(ecm_execution=record))


class NativeAdmissionTests(unittest.TestCase):
    def test_retained_identity_refuses_changed_or_removed_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "binding"
            path.write_bytes(b"original")
            admission = NativeAdmission({}, ((str(path), _fingerprint(path)),))
            admission.check()
            path.write_bytes(b"changed bytes")
            with self.assertRaisesRegex(ValueError, "changed"):
                admission.check()
            path.unlink()
            with self.assertRaisesRegex(ValueError, "disappeared"):
                admission.check()

    def test_resource_provider_error_cannot_become_optional_fallback(self):
        original = OSError("observer denied")

        class FailingBudget:
            def consume(self, amount=1):
                raise original

        class Packet:
            full_ecm, calls, variant = True, 0, "scalar"

            def available(self, n, ledger):
                try:
                    ledger.consume(1)
                except OSError as error:
                    raise BackendUnavailableError(
                        "wrapped by loader"
                    ) from error

        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)
            (source / "build").mkdir()
            for name in ("manifest.json", "libfactor_ecm.dylib"):
                (source / "build" / name).write_bytes(b"test")
            with patch("v2.ecm.native_packet.SOURCE", source):
                with self.assertRaises(OSError) as caught:
                    admit_packet(Packet(), 1013, FailingBudget())
        self.assertIs(caught.exception, original)


if __name__ == "__main__":
    unittest.main()
