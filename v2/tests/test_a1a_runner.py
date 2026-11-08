"""Frozen replay, counted arithmetic bridge and synthetic acceptance gates."""

import copy
import json
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    historical_asdict as asdict,
)
from v2.benchmarks.infrastructure.resource_control import (
    historical_replace as replace,
)
from v2.benchmarks.qs.a1a_integration import analysis, inputs, runtime, study
from v2.benchmarks.qs.c10a04 import configuration
from v2.benchmarks.qs.c10a04.accounting import Ledger
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import SIQSConfig, SIQSJob, default_siqs_config
from v2.tests.qs_default_controls import previous_siqs_config
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_qs import unlimited_budget
from v2.tests.test_siqs_cli import cli


def synthetic_rows():
    return [
        dict(
            id=str(index),
            arm=arm,
            repeat=repeat,
            seconds=100 if arm == "control" else 95,
            cpu_seconds=100 if arm == "control" else 95,
            complete=True,
            stats=dict(blocks=7, polynomials=2, relations=3),
        )
        for index in range(9)
        for repeat in range(2)
        for arm in study.ARMS
    ]


class FrozenRunnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.TemporaryDirectory()
        cls.control = runtime.load_control(Path(cls.directory.name))

    @classmethod
    def tearDownClass(cls):
        for name in tuple(sys.modules):
            if name == runtime.PACKAGE or name.startswith(
                runtime.PACKAGE + "."
            ):
                del sys.modules[name]
        cls.directory.cleanup()

    def test_exact_configuration_bridge_and_private_control_types(self):
        for digits in (50, 55):
            old = study.configuration(
                10 ** (digits - 1), "control", self.control
            )
            new = study.configuration(
                10 ** (digits - 1), "selected", self.control
            )
            values = asdict(new)
            self.assertEqual(values.pop("polynomial_policy"), "standard")
            self.assertEqual(values.pop("multiplier_policy"), "standard")
            self.assertEqual(values.pop("root_update_policy"), "standard")
            self.assertEqual(values.pop("sign_policy"), "swap_high")
            self.assertEqual(
                values["collector"].pop("sparse_queue_bytes"), 32 * 2**20
            )
            values["collector"]["marking"] = old.collector.marking
            self.assertEqual(values, asdict(old))
            self.assertIsNot(type(old), type(new))
        for digits in (40, 49, 51, 52, 54, 56, 60, 65):
            with self.assertRaises(ValueError):
                study.configuration(
                    10 ** (digits - 1), "selected", self.control
                )
        with self.assertRaises(ValueError):
            study.configuration(10**49, "inverse", self.control)

    def test_revealed_55_count_setup_matches_external_resolution(self):
        fixtures = json.loads(study.WARM.read_text())["fixtures"]
        n = next(f["n"] for f in fixtures if f["digits"] == 55)
        budget = unlimited_budget()
        previous, divisor = configuration.resolve(n, 55, "compact", budget)
        self.assertIsNone(divisor)
        prior_work = budget.used
        job = SIQSJob(
            n, config=previous_siqs_config(n), budget=unlimited_budget()
        )
        job._setup()
        self.assertEqual(len(job.base.entries), 2500)
        self.assertEqual(job.base.bound, previous.base_bound)
        self.assertEqual(job.multiplier, previous.multiplier)
        values = asdict(job.config)
        self.assertEqual(values.pop("sign_policy"), "swap_high")
        self.assertEqual(
            values["collector"].pop("sparse_queue_bytes"), 32 * 2**20
        )
        values.pop("base_count")
        values.update(base_bound=job.base.bound, multiplier=job.multiplier)
        self.assertEqual(values["collector"]["marking"], study.SELECTED)
        values["collector"]["marking"] = "sparse"
        previous_values = asdict(previous)
        self.assertEqual(previous_values.pop("sign_policy"), "swap_high")
        self.assertEqual(
            previous_values["collector"].pop("sparse_queue_bytes"), 32 * 2**20
        )
        from v2.tests.resource_stops import recipe_values

        self.assertEqual(recipe_values(values), recipe_values(previous_values))
        # Count-stage reservations are identical; the new boundary additionally
        # preflights temporary sieve memory and then constructs the same base.
        count_budget = unlimited_budget()
        from v2.qs.count_policy import endpoint_for_count
        from v2.qs.multiplier import select_multiplier

        choice = select_multiplier(n, budget=count_budget)
        endpoint_for_count(
            n, choice.multiplier, 2500, count_budget, ceiling=250000
        )
        self.assertEqual(count_budget.used, prior_work)

    def test_small_correctness_call_and_teardown_without_accepted_measurement(
        self,
    ):
        fixture = dict(id="qa-only", n=10403, factors=[101, 103], digits=5)
        for arm in study.ARMS:
            cls = (
                self.control[0].SIQSConfig if arm == "control" else SIQSConfig
            )
            cfg = cls(base_bound=50, half_width=64, max_half_width=64)
            if arm == "selected":
                cfg = replace(
                    cfg,
                    collector=replace(cfg.collector, marking=study.SELECTED),
                )
            with patch.object(study, "configuration", return_value=cfg):
                row = study.run_one(fixture, arm, 3, self.control)
            study.validate(row, fixture)
            self.assertTrue(row["complete"])
            self.assertEqual(row["certainty"], ["proven_prime"] * 2)
            self.assertEqual(sorted(row["factors"]), [101, 103])

    def test_legacy_checkpoint_uses_saved_policy_without_default_resolution(
        self,
    ):
        cfg = default_siqs_config(10**54 + 123)
        cfg = replace(cfg, collector=replace(cfg.collector, marking="sparse"))
        job = SIQSJob(10**54 + 123, config=cfg, budget=Budget(work_limit=0))
        job.run()
        saved = job.checkpoint()
        with patch(
            "v2.qs.defaults.default_siqs_config", side_effect=AssertionError
        ):
            restored = SIQSJob.from_checkpoint(
                saved, budget=Budget(work_limit=0)
            )
        self.assertEqual(restored.config, cfg)
        self.assertEqual(restored.config.collector.marking, "sparse")

    def test_selected_counted_portfolio_checkpoint_survives_cli_resume(self):
        n = 10**54 + 123
        selected = study.configuration(n, "selected", self.control)
        cfg = PortfolioConfig(siqs=selected)
        first = factorize_bounded(n, config=cfg, budget=Budget(work_limit=0))
        saved = first.checkpoint["payload"]["config"]
        self.assertEqual(saved["siqs"]["collector"]["marking"], study.SELECTED)
        self.assertEqual(saved["siqs"]["base_count"], 2500)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.json"
            path.write_text(json.dumps(first.checkpoint))
            resumed = cli(
                "--resume", path, "--seconds", "0", "--checkpoint", path
            )
            self.assertEqual(resumed.returncode, 1, resumed.stderr)
            self.assertEqual(
                json.loads(path.read_text())["payload"]["config"], saved
            )

    def test_owned_baseline_corruption_and_frozen_live_hashes_fail_closed(
        self,
    ):
        data = runtime.checked_baseline()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "baseline.json"
            data["source"]["v2/qs/sieve_collector.py"] += "\n# corrupt\n"
            path.write_text(json.dumps(data))
            with patch.object(runtime, "BASELINE", path):
                with self.assertRaises(ValueError):
                    runtime.checked_baseline()
        if study.CONTROL.exists():
            with patch.object(study, "sources", return_value={}):
                with self.assertRaises(ValueError):
                    study.load()

    def test_per_band_cost_gate_requires_independent_complete_stable_inputs(
        self,
    ):
        rows = synthetic_rows()
        report = analysis.summarize(rows, stable_warmup=True)
        self.assertTrue(report["cost_gate_passed"])
        self.assertEqual(report["independent_inputs"], 9)
        self.assertEqual(report["exact_sign_flip_p"], 1 / 512)
        self.assertFalse(
            analysis.summarize(rows, stable_warmup=False)["cost_gate_passed"]
        )
        incomplete = copy.deepcopy(rows)
        incomplete[0]["complete"] = False
        self.assertFalse(
            analysis.summarize(incomplete, stable_warmup=True)[
                "cost_gate_passed"
            ]
        )
        drifting = copy.deepcopy(rows)
        drifting[0]["seconds"] = 150
        self.assertFalse(
            analysis.summarize(drifting, stable_warmup=True)[
                "cost_gate_passed"
            ]
        )
        for row in rows:
            row["repeat"] = 0
        with self.assertRaises(ValueError):
            analysis.summarize(rows, stable_warmup=True)

    def test_missing_measurement_is_censored_with_its_full_reservation(self):
        limits = dict(
            total=dict(wall=10, cpu=10, work=100),
            base50=dict(wall=10, cpu=10, work=100),
            disk_bytes=2**20 * 8,
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            ledger = Ledger(output, "qa-freeze", limits)

            def interrupt():
                raise KeyboardInterrupt

            identity = dict(
                sample=0,
                repeat=0,
                fixture="qa-only",
                n=10403,
                digits=50,
                seed=7,
                arm="selected",
                cap=3,
            )
            with self.assertRaises(KeyboardInterrupt):
                ledger.call("base50", "lost", 4, 100, identity, interrupt)
            recovered = Ledger(output, "qa-freeze", limits)
            wrapper = json.loads((output / "lost.json").read_text())
            self.assertEqual(wrapper["result"]["reason"], "interrupted")
            self.assertFalse(wrapper["result"]["complete"])
            self.assertEqual(recovered.used(field="work"), 100)
            self.assertGreaterEqual(recovered.used(field="wall"), 4)

    def test_generated_corpus_proofs_and_binding_when_present(self):
        if inputs.CORPUS.exists():
            inputs.check(
                json.loads(inputs.CORPUS.read_text()),
                study.digest(study.CONTROL),
            )

    def test_warmup_stability_is_bounded_and_stays_in_the_stage_process(self):
        for unstable in (False, True):
            with tempfile.TemporaryDirectory() as directory:
                ledger = Ledger(Path(directory), "qa-only", study.LIMITS)
                calls = 0

                def prefix(fixture, arm, seconds, control, *, blocks):
                    nonlocal calls
                    duration = 1 + int(unstable and (calls // 2) % 2 == 1)
                    calls += 1
                    return dict(
                        divisor=None,
                        factors=[],
                        remaining=[fixture["n"]],
                        certainty=[],
                        base=dict(count=2500),
                        reason="paused",
                        seconds=duration,
                        cpu_seconds=duration,
                        work=0,
                        stats=dict(blocks=blocks, polynomials=32, relations=3),
                    )

                with patch.object(study, "run_one", prefix):
                    stable = study.warmup(ledger, 50, "base", self.control)
                self.assertEqual(stable, not unstable)
                self.assertEqual(calls, 40 if unstable else 6)

    def test_failed_warmup_records_stage_and_refuses_cold_restart(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            ledger = Ledger(output, "qa-only", study.LIMITS)
            with patch.object(study, "warmup", return_value=False):
                study.execute(
                    ledger,
                    50,
                    "base",
                    self.control,
                    time.monotonic(),
                    time.process_time(),
                )
            marker = json.loads((output / "base-finished-50.json").read_text())
            self.assertFalse(marker["stable_warmup"])
            self.assertGreaterEqual(marker["stage_wall"], 0)
            self.assertEqual(ledger.used(field="wall"), 60)
            with self.assertRaises(ValueError):
                study.execute(
                    ledger,
                    50,
                    "base",
                    self.control,
                    time.monotonic(),
                    time.process_time(),
                )

    def test_orchestration_overrun_is_charged_before_stopping(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            ledger = Ledger(output, "qa-only", study.LIMITS)
            ledger.call(
                "orchestration",
                "base-overhead-50",
                60,
                0,
                dict(stage="base", digits=50),
                lambda: {"work": 0},
            )
            with (
                patch.object(study.time, "monotonic", return_value=100),
                patch.object(study.time, "process_time", return_value=100),
            ):
                with self.assertRaises(RuntimeError):
                    study.finish_stage(ledger, 50, "base", False, 35, 36)
            self.assertEqual(ledger.used(field="wall"), 65)
            self.assertEqual(ledger.used(field="cpu"), 64)
            self.assertEqual(ledger.data["calls"][-1]["status"], "failed")
            marker = json.loads((output / "base-finished-50.json").read_text())
            self.assertFalse(marker["overhead_within_prepayment"])
