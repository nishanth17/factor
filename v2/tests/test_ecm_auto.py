"""Automatic adoption, resource refusal, exact restart and owned lifecycle."""

import copy
import hashlib
import json
import multiprocessing
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
    recording_budget,
)
from v2.ecm.selection import COMPATIBILITY
from v2.execution.ecm_auto import (
    DEFAULT_MEMORY_BYTES,
    POLICY_PATH,
    POLICY_SHA256,
    SELECTOR_OWNERSHIP,
    AutomaticECMOwner,
)
from v2.execution.ecm_defaults import ESCALATING140
from v2.execution.ecm_workers import WorkerRun
from v2.portfolio import PortfolioConfig, _canonical, factorize_bounded
from v2.tests.resource_stops import InterruptBudget as Budget

ROOT = Path(__file__).resolve().parents[2]


def config(**changes):
    values = dict(
        trial_bound=100,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=((50, 500, 8),),
        siqs=None,
        sss=None,
    )
    values.update(changes)
    return adapt_resource_configuration(PortfolioConfig, **values)


def budget(**changes):
    values = dict(work_limit=10**9, seconds=20, cpu_seconds=4096)
    values.update(changes)
    return Budget(**values)


class FakeSelected:
    """Checked divisor supplied independently of the production ECM path."""

    instances = []
    failure = None
    load_failure = None
    divisor = 1009

    def __init__(self, *args, **kwargs):
        self.arguments = kwargs
        self.closed = False
        self.close_count = 0
        self.last_stats = dict(engaged=True, worker_dispatches=1)
        self.ledger = None
        type(self).instances.append(self)

    def _load(self, ledger):
        if self.load_failure is not None:
            raise self.load_failure

    def run_portfolio(self, n, *, budget, **kwargs):
        self.ledger = getattr(budget, "budget", budget)
        if self.failure is not None:
            raise self.failure
        budget.consume(13)
        return WorkerRun(self.divisor, "factor_found", {}, self.last_stats)

    def close(self):
        self.close_count += 1
        self.closed = True
        if self.ledger is not None:
            # A final snapshot must include shutdown charges, not just curves.
            self.ledger.consume(17)
            self.ledger.prior_cpu += 0.125


class AutomaticTests(unittest.TestCase):
    def setUp(self):
        FakeSelected.instances = []
        FakeSelected.failure = FakeSelected.load_failure = None
        FakeSelected.divisor = 1009
        self.runtime = patch(
            "v2.execution.ecm_auto._runtime_supported", return_value=True
        )
        self.runtime.start()
        self.addCleanup(self.runtime.stop)
        # The fake catalogue covers our small independent factor oracle.
        self.sizes = patch(
            "v2.execution.ecm_auto._size_supported", return_value=True
        )
        self.sizes.start()
        self.addCleanup(self.sizes.stop)

    def test_default_dispatch_validates_factor_and_closes_before_snapshot(
        self,
    ):
        with patch(
            "v2.execution.ecm_selection.SelectedECMContext", FakeSelected
        ):
            run = factorize_bounded(
                1009 * 1013, seed=7, config=config(), budget=budget()
            )
        self.assertTrue(run.result.complete)
        self.assertEqual(run.result.reconstruct(), 1009 * 1013)
        self.assertEqual([f.value for f in run.result.factors], [1009, 1013])
        owner = FakeSelected.instances[0]
        self.assertEqual(owner.close_count, 1)
        self.assertEqual(owner.arguments["seed"], 7)
        self.assertGreaterEqual(run.cpu_seconds, 0.125)
        self.assertGreaterEqual(run.checkpoint["payload"]["cpu_used"], 0.125)
        self.assertEqual(run.work_used, owner.ledger.used)
        self.assertEqual(owner.arguments["memory_bytes"], DEFAULT_MEMORY_BYTES)

    def test_small_or_unbounded_grants_stay_serial_without_profile_load(self):
        for options in (
            dict(work_limit=2_000_000),
            dict(cpu_seconds=30),
            dict(cpu_seconds=None),
            dict(seconds=None),
        ):
            with (
                self.subTest(options=options),
                patch(
                    "v2.execution.ecm_selection.SelectedECMContext",
                    side_effect=AssertionError("unfunded profile load"),
                ),
            ):
                run = factorize_bounded(
                    1009 * 1013, config=config(), budget=budget(**options)
                )
                self.assertEqual(run.result.reconstruct(), 1009 * 1013)
                self.assertTrue(
                    any(e["stage"] == "ecm_selection" for e in run.events)
                )

    def test_explicit_parent_memory_or_zero_aggregate_cap_stays_serial(self):
        for settings, cap in (
            (config(memory_bytes=16 * 2**20), None),
            (config(), 0),
        ):
            with patch(
                "v2.execution.ecm_selection.SelectedECMContext",
                side_effect=AssertionError("unfunded pool"),
            ):
                run = factorize_bounded(
                    1009 * 1013, config=settings, budget=budget()
                )
            self.assertEqual(run.result.reconstruct(), 1009 * 1013)

    def test_parent_limit_includes_selector_coexistence_before_loading(self):
        # A 4 GiB fallback parent remains valid, but cannot also own a selector
        # inside the worker contract's separate 4 GiB parent ceiling.
        maximum = 4 * 2**30 - SELECTOR_OWNERSHIP
        for parent in (maximum + 1, 4 * 2**30):
            with (
                self.subTest(parent=parent),
                patch(
                    "v2.execution.ecm_selection.SelectedECMContext",
                    side_effect=AssertionError("unadmitted selector"),
                ),
            ):
                settings = config(memory_bytes=parent)
                generator = random.Random(7)
                before = generator.getstate()
                with AutomaticECMOwner(seed=7) as owner:
                    owner.configure(
                        settings, generator, None, automatic_memory=False
                    )
                    result = owner.run_portfolio(
                        1009 * 1013, config=settings, budget=budget()
                    )
                    self.assertIsNone(result)
                    self.assertTrue(owner.last_stats["legacy_serial"])
                    self.assertEqual(
                        owner.last_stats["engagement"]["reason"],
                        "worker_ownership_refused",
                    )
                self.assertEqual(generator.getstate(), before)

        with patch(
            "v2.execution.ecm_selection.SelectedECMContext", FakeSelected
        ):
            settings = config(memory_bytes=maximum)
            with AutomaticECMOwner(seed=7) as owner:
                owner.configure(
                    settings, random.Random(7), None, automatic_memory=False
                )
                result = owner.run_portfolio(
                    1009 * 1013, config=settings, budget=budget()
                )
                self.assertEqual(result.divisor, 1009)

    def test_unmeasured_sizes_do_not_import_or_load_the_selector(self):
        with (
            patch("v2.execution.ecm_auto._size_supported", return_value=False),
            patch(
                "v2.execution.ecm_selection.SelectedECMContext",
                side_effect=AssertionError("unmeasured catalogue import"),
            ),
        ):
            run = factorize_bounded(
                1009 * 1013, config=config(), budget=budget()
            )
        self.assertEqual(run.result.reconstruct(), 1009 * 1013)
        selection = next(
            e for e in run.events if e["stage"] == "ecm_selection"
        )
        self.assertEqual(
            selection["stats"]["engagement"]["reason"],
            "unmeasured_selection_cost",
        )

    def test_easy_exit_does_not_load_or_construct_a_worker_context(self):
        with patch(
            "v2.execution.ecm_selection.SelectedECMContext",
            side_effect=AssertionError("easy-exit worker setup"),
        ):
            for n in (1, 101, 2**80, 1009**2):
                run = factorize_bounded(n)
                self.assertTrue(run.result.complete)
                self.assertEqual(run.result.reconstruct(), n)

    def test_explicit_serial_preserves_rng_and_algorithm_result(self):
        settings = config()
        serial = factorize_bounded(
            1009 * 1013, seed=91, config=settings, ecm_executor=None
        )
        declined = factorize_bounded(
            1009 * 1013,
            seed=91,
            config=settings,
            budget=budget(cpu_seconds=30),
        )
        self.assertEqual(serial.result, declined.result)

        def curve_receipts(run):
            return [
                (event["seed"], event.get("sigma"), event["outcome"])
                for event in run.events
                if event["stage"] == "ecm"
            ]

        self.assertTrue(curve_receipts(serial))
        self.assertEqual(curve_receipts(serial), curve_receipts(declined))
        self.assertEqual(
            serial.checkpoint["payload"]["rng"],
            declined.checkpoint["payload"]["rng"],
        )
        self.assertNotIn("ecm_auto", serial.checkpoint["payload"]["state"])

    def test_failure_closes_owned_context(self):
        FakeSelected.failure = MemoryError("injected result failure")
        with patch(
            "v2.execution.ecm_selection.SelectedECMContext", FakeSelected
        ):
            with self.assertRaisesRegex(MemoryError, "injected"):
                factorize_bounded(
                    1009 * 1013, config=config(), budget=budget()
                )
        self.assertTrue(FakeSelected.instances[0].closed)
        self.assertEqual(FakeSelected.instances[0].close_count, 1)

    def test_corrupt_factor_publication_closes_owned_context(self):
        FakeSelected.divisor = 17
        with (
            patch(
                "v2.execution.ecm_selection.SelectedECMContext", FakeSelected
            ),
            self.assertRaises(ValueError),
        ):
            factorize_bounded(1009 * 1013, config=config(), budget=budget())
        self.assertTrue(FakeSelected.instances[0].closed)

    def test_result_construction_failure_still_closes_once(self):
        with (
            patch(
                "v2.execution.ecm_selection.SelectedECMContext", FakeSelected
            ),
            patch("v2.portfolio._pack", side_effect=MemoryError("checkpoint")),
        ):
            with self.assertRaises(MemoryError):
                factorize_bounded(
                    1009 * 1013, config=config(), budget=budget()
                )
        self.assertEqual(FakeSelected.instances[0].close_count, 1)

    def test_refusal_before_dispatch_has_no_worker_side_effect(self):
        for options in (
            dict(work_limit=0),
            dict(seconds=0),
            dict(cancelled=lambda: True),
        ):
            with patch(
                "v2.execution.ecm_selection.SelectedECMContext",
                side_effect=AssertionError("refused preparation"),
            ):
                run = factorize_bounded(
                    1009 * 1013, config=config(), budget=budget(**options)
                )
            self.assertEqual(run.result.reconstruct(), 1009 * 1013)
            self.assertFalse(FakeSelected.instances)

    def test_selected_pool_must_fit_explicit_aggregate_memory(self):
        from v2.execution.ecm_selection import SelectedECMContext

        n = 10**79 + 7
        settings = PortfolioConfig().for_input(n)
        with SelectedECMContext(POLICY_PATH, POLICY_SHA256) as owner:
            owner._load(budget())
            with self.assertRaises(MemoryError):
                owner._select(n, settings, None, None, portfolio=True)

    def test_retained_parent_coexistence_refuses_with_accounted_serial_path(
        self,
    ):
        from v2.execution.ecm_selection import SelectedECMContext

        n = 10**59 + 7
        settings = PortfolioConfig().for_input(n)
        ledger = budget()
        with (
            SelectedECMContext(POLICY_PATH, POLICY_SHA256) as owner,
            patch(
                "v2.execution.ecm_workers.ECMExecutor",
                side_effect=AssertionError("unadmitted worker construction"),
            ),
        ):
            result = owner.run_portfolio(n, config=settings, budget=ledger)
            self.assertIsNone(result)
            self.assertIsNone(owner.inner)
            self.assertEqual(
                owner.last_stats["engagement"]["reason"],
                "worker_ownership_refused",
            )
            self.assertTrue(owner.last_stats["legacy_serial"])
            # Loading and validating metadata was real work, despite refusal.
            self.assertGreater(ledger.used, 0)
        # Component and portfolio modes intentionally require distinct owners.
        with SelectedECMContext(POLICY_PATH, POLICY_SHA256) as component:
            with self.assertRaisesRegex(MemoryError, "parent ownership"):
                component.run(n, config=settings, budget=budget())

    def test_borrowed_context_is_not_closed(self):
        class Borrowed:
            identity = "borrowed"
            settings = type("Settings", (), dict(parent_memory=2**30))()
            last_stats = dict(engaged=True)
            closed = False

            def run(self, n, **kwargs):
                return WorkerRun(1009, "factor_found", {}, self.last_stats)

            def close(self):
                self.closed = True

        owner = Borrowed()
        run = factorize_bounded(
            1009 * 1013, config=config(), budget=budget(), ecm_executor=owner
        )
        self.assertTrue(run.result.complete)
        self.assertFalse(owner.closed)

    def test_auto_resume_retains_seed_policy_and_memory(self):
        settings = config()
        first = factorize_bounded(
            1009 * 1013, seed=991, config=settings, budget=budget(work_limit=0)
        )
        header = first.checkpoint["payload"]["state"]["ecm_auto"]
        with patch(
            "v2.execution.ecm_selection.SelectedECMContext", FakeSelected
        ):
            resumed = factorize_bounded(
                1009 * 1013,
                seed=123,
                config=settings,
                budget=budget(),
                checkpoint=first.checkpoint,
            )
        self.assertTrue(resumed.result.complete)
        self.assertEqual(
            resumed.checkpoint["payload"]["state"]["ecm_auto"], header
        )
        self.assertEqual(
            FakeSelected.instances[0].arguments["seed"], header["seed"]
        )
        for options in (dict(ecm_memory_bytes=0), dict(ecm_executor=None)):
            with self.assertRaises(ValueError):
                factorize_bounded(
                    1009 * 1013,
                    config=settings,
                    budget=budget(),
                    checkpoint=first.checkpoint,
                    **options,
                )

    def test_automatic_envelope_rejects_downgrade_or_recursive_schema(self):
        first = factorize_bounded(
            1009 * 1013, config=config(), budget=budget(work_limit=0)
        )
        self.assertEqual(first.checkpoint["payload"]["version"], 13)
        for downgrade in (True, False):
            checkpoint = copy.deepcopy(first.checkpoint)
            if downgrade:
                checkpoint["payload"]["version"] = checkpoint["payload"].pop(
                    "portfolio_version"
                )
            else:
                checkpoint["payload"]["portfolio_version"] = 13
            checkpoint["sha256"] = hashlib.sha256(
                _canonical(checkpoint["payload"]).encode()
            ).hexdigest()
            with self.assertRaises(ValueError):
                factorize_bounded(
                    1009 * 1013,
                    config=config(),
                    budget=budget(),
                    checkpoint=checkpoint,
                )

    def test_legacy_reader_rejects_new_automatic_checkpoint(self):
        from dataclasses import asdict

        from v2.benchmarks.ecm.b3.b3_production import baseline

        settings = config(
            pm1_gap_mode="cached",
            pm1_chunk_size=None,
            ecm_chain_mode="off",
            ecm_chain_bytes=0,
        )
        first = factorize_bounded(
            1009 * 1013, config=settings, budget=budget(work_limit=0)
        )
        old = baseline()
        options = asdict(settings)
        for key in (
            "allocation",
            "ecm_chain_mode",
            "ecm_chain_bytes",
            "ecm_chain_family",
        ):
            options.pop(key)
        with self.assertRaises(ValueError):
            old.factorize_bounded(
                1009 * 1013,
                config=adapt_resource_configuration(
                    old.PortfolioConfig, **options
                ),
                checkpoint=first.checkpoint,
                budget=recording_budget(
                    old.Budget, work_limit=1000, seconds=10, cpu_seconds=10
                ),
            )

    def test_legacy_checkpoint_does_not_adopt_workers(self):
        first = factorize_bounded(
            1009 * 1013,
            config=config(),
            budget=budget(work_limit=0),
            ecm_executor=None,
        )
        with patch(
            "v2.execution.ecm_selection.SelectedECMContext",
            side_effect=AssertionError("legacy adoption"),
        ):
            resumed = factorize_bounded(
                1009 * 1013,
                config=config(),
                budget=budget(),
                checkpoint=first.checkpoint,
            )
        self.assertNotIn("ecm_auto", resumed.checkpoint["payload"]["state"])
        self.assertEqual(resumed.result.reconstruct(), 1009 * 1013)

    def test_tampered_automatic_header_is_rejected(self):
        first = factorize_bounded(
            1009 * 1013, config=config(), budget=budget(work_limit=0)
        )
        for field, value in (
            ("version", "unknown"),
            ("seed", -1),
            ("profile", "0" * 64),
        ):
            checkpoint = copy.deepcopy(first.checkpoint)
            checkpoint["payload"]["state"]["ecm_auto"][field] = value
            checkpoint["sha256"] = hashlib.sha256(
                _canonical(checkpoint["payload"]).encode()
            ).hexdigest()
            with self.assertRaises(ValueError):
                factorize_bounded(
                    1009 * 1013,
                    config=config(),
                    budget=budget(),
                    checkpoint=checkpoint,
                )

    def test_nested_process_is_serial_without_loading_policy(self):
        with (
            patch("multiprocessing.parent_process", return_value=object()),
            patch(
                "v2.execution.ecm_selection.SelectedECMContext",
                side_effect=AssertionError("nested pool"),
            ),
        ):
            run = factorize_bounded(
                1009 * 1013, config=config(), budget=budget()
            )
        self.assertEqual(run.result.reconstruct(), 1009 * 1013)
        selection = next(
            e for e in run.events if e["stage"] == "ecm_selection"
        )
        self.assertEqual(
            selection["stats"]["engagement"]["reason"], "nested_process"
        )

    def test_missing_or_stale_policy_falls_back_before_search(self):
        FakeSelected.load_failure = ValueError("stale source")
        with patch(
            "v2.execution.ecm_selection.SelectedECMContext", FakeSelected
        ):
            run = factorize_bounded(
                1009 * 1013, config=config(), budget=budget()
            )
        self.assertTrue(FakeSelected.instances[0].closed)
        self.assertEqual(run.result.reconstruct(), 1009 * 1013)
        selection = next(
            e for e in run.events if e["stage"] == "ecm_selection"
        )
        self.assertEqual(
            selection["stats"]["engagement"]["reason"],
            "unavailable_automatic_policy",
        )

    def test_method_seed_derivation_does_not_consume_random_draws(self):
        generator = random.Random(71)
        before = generator.getstate()
        owner = AutomaticECMOwner()
        owner.configure(config(), generator, None, automatic_memory=True)
        self.assertEqual(generator.getstate(), before)
        other = AutomaticECMOwner()
        other.configure(config(), generator, None, automatic_memory=True)
        self.assertEqual(owner.identity, other.identity)

    def test_packaged_policy_preserves_measured_rows_and_compatibility(
        self,
    ):
        from v2.ecm.worker_policy import policy_data

        policy = policy_data()
        original = json.loads(
            (
                ROOT
                / "v2/benchmarks/inputs/controls"
                / "b12_f4_targeted_selected_v1.json"
            ).read_text()
        )
        self.assertEqual(policy["rows"], original["rows"])
        self.assertEqual(
            policy["serial_settings"], original["serial_settings"]
        )
        self.assertEqual(policy["compatibility"], COMPATIBILITY)

    def test_policy_packaging_preserves_costs_and_rejects_wrong_reference(
        self,
    ):
        from v2.benchmarks.ecm.b12_f4.adoption_profile import build_profile
        from v2.execution.ecm_auto import MEASURED_PROFILE_SHA256

        reference = (
            ROOT
            / "v2/benchmarks/inputs/controls"
            / "b12_f4_targeted_selected_v1.json"
        )
        data = build_profile(reference, MEASURED_PROFILE_SHA256)
        from v2.ecm.worker_policy import policy_data

        self.assertEqual(data, policy_data())
        with self.assertRaises(ValueError):
            build_profile(reference, "0" * 64)

    def test_policy_repinning_changes_only_one_bounded_hash(self):
        from v2.benchmarks.ecm.b12_f4.adoption_profile import repin_adapter

        text = 'POLICY_SHA256 = (\n    "' + "a" * 64 + '"\n)'
        expected = text.replace("a" * 64, "b" * 64)
        self.assertEqual(repin_adapter(text, "b" * 64), expected)
        for malformed in ("", text + text):
            with self.assertRaises(ValueError):
                repin_adapter(malformed, "b" * 64)

    def test_new_numerical_defaults_and_explicit_plans(self):
        for digits in (60, 80, 100):
            self.assertEqual(
                PortfolioConfig().for_input(10 ** (digits - 1)).ecm_tiers,
                ESCALATING140,
            )
            explicit = PortfolioConfig(ecm_tiers=((50, 500, 8),))
            self.assertIs(explicit.for_input(10 ** (digits - 1)), explicit)
        self.assertEqual(
            PortfolioConfig().for_input(10**69).ecm_tiers,
            ((2000, 147396, 32),),
        )

    def test_default_auto_cli_and_serial_override_checkpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "state.json"
            for flags, automatic in (
                ([], True),
                (["--ecm-workers", "serial"], False),
            ):
                call = subprocess.run(
                    [
                        sys.executable,
                        "-B",
                        "-m",
                        "v2.factor",
                        str(10**79 + 7),
                        "--work-limit",
                        "0",
                        "--checkpoint",
                        str(path),
                        *flags,
                    ],
                    cwd=ROOT,
                    capture_output=True,
                    text=True,
                    timeout=15,
                )
                self.assertEqual(call.returncode, 1, call.stderr)
                payload = json.loads(path.read_text())["payload"]
                self.assertEqual(
                    payload["config"]["ecm_tiers"],
                    [list(t) for t in ESCALATING140],
                )
                self.assertEqual("ecm_auto" in payload["state"], automatic)
            call = subprocess.run(
                [sys.executable, "-B", "-m", "v2.factor", "101", "--verbose"],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=15,
            )
            self.assertEqual(call.returncode, 0, call.stderr)
            self.assertIn("Portfolio: complete", call.stdout)

    def test_real_default_spawn_dispatch_is_drained_without_orphans(self):
        before = {child.pid for child in multiprocessing.active_children()}
        corpus = json.loads(
            (
                ROOT
                / "v2/benchmarks/inputs/corpora"
                / "b12_f4_targeted_confirmation_v1.json"
            ).read_text()
        )
        n = next(
            row["n"] for row in corpus["fixtures"] if row["id"] == "ecm-80-0"
        )
        run = factorize_bounded(
            n,
            seed=7919,
            budget=budget(seconds=6),
            stop_after_split=True,
        )
        self.assertEqual(run.result.reconstruct(), n)
        events = [e for e in run.events if e["stage"] == "ecm_workers"]
        self.assertTrue(events, run.reason)
        from v2.execution.ecm_workers import arithmetic_options
        from v2.execution.worker_contract import Manifest, WorkerSettings

        manifest = Manifest(
            n,
            7919,
            f"{n}/{n}/1",
            arithmetic_options(
                PortfolioConfig().for_input(n),
                WorkerSettings(),
            ),
        )
        self.assertEqual(events[0]["stats"]["manifest"], manifest.identity)
        self.assertTrue(
            any(e["stats"]["worker_dispatches"] > 0 for e in events)
        )
        self.assertEqual(
            {child.pid for child in multiprocessing.active_children()}, before
        )


class ColdSerialAdmissionTests(unittest.TestCase):
    def settings(self, digits, **changes):
        from dataclasses import replace

        n = 10 ** (digits - 1) + 7
        config = PortfolioConfig().for_input(n)
        return n, replace(config, **changes) if changes else config

    def test_default_cold_plans_decline_before_context_or_runtime_loading(
        self,
    ):
        for digits in (20, 30, 50, 70):
            n, settings = self.settings(digits)
            generator = random.Random(7)
            before = generator.getstate()
            fallback = settings.siqs
            with (
                self.subTest(digits=digits),
                patch(
                    "v2.execution.ecm_auto.selected_context",
                    side_effect=AssertionError("unexpected selector"),
                ),
                patch(
                    "v2.execution.ecm_auto._runtime_supported",
                    side_effect=AssertionError("unnecessary platform load"),
                ),
                AutomaticECMOwner(seed=7) as owner,
            ):
                owner.configure(
                    settings, generator, None, automatic_memory=True
                )
                ledger = budget()
                self.assertIsNone(
                    owner.run_portfolio(n, config=settings, budget=ledger)
                )
                self.assertIsNone(owner.context)
                self.assertEqual(ledger.used, 0)
                self.assertTrue(owner.last_stats["legacy_serial"])
            self.assertEqual(generator.getstate(), before)
            self.assertIs(settings.siqs, fallback)
            if digits == 50:
                self.assertIsNotNone(fallback)

    def test_borrowed_packaged_cold_context_skips_profile_io(self):
        from v2.execution.ecm_auto import selected_context

        for digits in (30, 50, 70):
            n, settings = self.settings(digits)
            with (
                self.subTest(digits=digits),
                selected_context() as owner,
                patch.object(
                    owner, "_load", side_effect=AssertionError("profile I/O")
                ),
            ):
                ledger = budget()
                self.assertIsNone(
                    owner.run_portfolio(n, config=settings, budget=ledger)
                )
                self.assertEqual(ledger.used, 0)
                self.assertIsNone(owner.inner)
                self.assertEqual(
                    owner.last_stats["engagement"]["reason"],
                    "cold_serial_default",
                )

    def test_warm_custom_or_larger_plan_keeps_regular_selection(self):
        from dataclasses import replace

        from v2.execution.ecm_auto import selected_context
        from v2.execution.ecm_selection import SelectedECMContext

        n, settings = self.settings(30)
        larger = replace(settings, ecm_tiers=((2000, 50000, 65),))
        owners = [
            (selected_context(prefer_warm=True), settings),
            (selected_context(), larger),
            (
                SelectedECMContext(ROOT / "custom.json", "0" * 64),
                settings,
            ),
        ]
        for owner, config in owners:
            with (
                owner,
                patch.object(
                    owner,
                    "_load",
                    side_effect=LookupError("regular selection"),
                ),
            ):
                with self.assertRaisesRegex(LookupError, "regular selection"):
                    owner.run_portfolio(n, config=config, budget=budget())

    def test_expired_cancelled_or_corrupt_resume_does_not_use_bypass(self):
        from v2.execution.budget import BudgetExhaustedError
        from v2.execution.ecm_auto import selected_context

        n, settings = self.settings(50)
        for limits in (dict(seconds=0), dict(cancelled=lambda: True)):
            with selected_context() as owner:
                with self.assertRaises(BudgetExhaustedError):
                    owner.run_portfolio(
                        n, config=settings, budget=budget(**limits)
                    )
        with selected_context() as owner:
            with self.assertRaises(ValueError):
                owner.run_portfolio(
                    n, config=settings, budget=budget(), checkpoint={}
                )


if __name__ == "__main__":
    unittest.main()
