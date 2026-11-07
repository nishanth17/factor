"""Admission, durable restart and corpus/progress contracts for G1."""

import copy
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.qs.a7a import confirmation, stage_driver
from v2.benchmarks.qs.a7a.combined import corpus_identity, fixtures
from v2.benchmarks.qs.a7a.progress import Progress
from v2.benchmarks.qs.a7a.tier_a import TierA


class G1Tests(unittest.TestCase):
    def plan(self):
        batches = []
        for digits, lifetime in sorted(stage_driver.ROWS):
            cell = dict(
                digits=digits,
                lifetime=lifetime,
                samples=1,
                endpoint="complete",
                arm="S0",
                memory=1024,
                wall=1,
                cpu=10,
                work=100,
                seed=7,
                fixture=f"regime_calibration_{digits}_0",
            )
            batches.append(
                dict(
                    identity=f"r{digits}-{lifetime}",
                    source="pin",
                    corpus_kind="regimes",
                    corpus_sha256=corpus_identity("regimes"),
                    wall=6,
                    cpu=10,
                    cells=[cell],
                )
            )
        return dict(
            stage="G2",
            source="pin",
            corpus_kind="regimes",
            corpus_sha256=corpus_identity("regimes"),
            batches=batches,
            wall=100,
            cpu=100,
            recipe_caps={f"{d}:{life}": 3 for d, life in stage_driver.ROWS},
        )

    def approve(self, plan):
        return dict(
            stage=plan["stage"],
            plan_sha256=stage_driver.digest(plan),
            approved=True,
        )

    def test_coverage_grants_and_exact_approval(self):
        plan = self.plan()
        with patch.object(stage_driver, "source", return_value="pin"):
            self.assertEqual(
                len(stage_driver.validate(plan, self.approve(plan))), 8
            )
            for changed in (
                "coverage",
                "fixture",
                "grant",
                "source",
                "cap",
                "nan",
                "inf",
                "unknown",
            ):
                bad = copy.deepcopy(plan)
                if changed == "coverage":
                    bad["batches"].pop()
                elif changed == "fixture":
                    bad["batches"][0]["cells"][0]["digits"] = 70
                elif changed == "grant":
                    bad["wall"] = 1
                elif changed == "source":
                    bad["source"] = "other"
                elif changed == "unknown":
                    bad["batches"][0]["cells"][0]["unsupported_flag"] = True
                elif changed in ("nan", "inf"):
                    bad["wall"] = float(changed)
                else:
                    bad["recipe_caps"]["40:cold"] = 0
                with (
                    self.subTest(changed=changed),
                    self.assertRaises(ValueError),
                ):
                    stage_driver.validate(bad, self.approve(bad))
            with self.assertRaises(ValueError):
                stage_driver.validate(plan, {})

    def test_diagnostics_never_supply_complete_coverage(self):
        plan = self.plan()
        for batch in plan["batches"]:
            batch["cells"][0].update(endpoint="prefix", purpose="diagnostic")
        with patch.object(stage_driver, "source", return_value="pin"):
            with self.assertRaises(ValueError):
                stage_driver.validate(plan, self.approve(plan))

    def test_restart_skips_completed_batches(self):
        calls = []

        def run(plan, output):
            calls.append(plan)
            output.mkdir()
            (output / "ledger.json").write_text(
                json.dumps(
                    dict(wall=0.1, cpu=0.1, attempts=[dict(status="reported")])
                )
            )

        plan = self.plan()
        with (
            tempfile.TemporaryDirectory() as folder,
            patch.object(stage_driver, "source", return_value="pin"),
            patch.object(stage_driver, "run", side_effect=run),
        ):
            stage_driver.execute(plan, self.approve(plan), folder)
            stage_driver.execute(plan, self.approve(plan), folder)
            self.assertEqual(len(calls), 8)
            changed = copy.deepcopy(plan)
            changed["wall"] += 1
            with self.assertRaises(ValueError):
                stage_driver.execute(changed, self.approve(changed), folder)

    def test_unfinished_attempt_is_not_reissued(self):
        plan = self.plan()
        with (
            tempfile.TemporaryDirectory() as folder,
            patch.object(stage_driver, "source", return_value="pin"),
            patch.object(
                stage_driver, "run", side_effect=RuntimeError("crash")
            ) as run,
        ):
            with self.assertRaises(RuntimeError):
                stage_driver.execute(plan, self.approve(plan), folder)
            with self.assertRaises(ValueError):
                stage_driver.execute(plan, self.approve(plan), folder)
            self.assertEqual(run.call_count, 1)

    def test_reference_honors_declared_checkpoint_allowance(self):
        from dataclasses import asdict

        from v2.benchmarks.qs.a7a.mainline_reference import (
            reference_configuration,
        )
        from v2.qs.defaults import default_siqs_config

        n = fixtures("regimes")["regime_calibration_60_0"]["n"]
        original = default_siqs_config(n, memory_bytes=4 * 2**30)
        self.assertEqual(original.checkpoint_bytes, 2**20)
        declared = reference_configuration(
            n, dict(qs_memory_bytes=4 * 2**30, checkpoint_bytes=64 * 2**20)
        )
        self.assertEqual(declared.checkpoint_bytes, 64 * 2**20)
        before, after = asdict(original), asdict(declared)
        before.pop("checkpoint_bytes")
        after.pop("checkpoint_bytes")
        self.assertEqual(before, after)

    def test_confirmation_proofs_and_disjointness(self):
        data = json.loads(confirmation.PATH.read_text())
        old = [r["n"] for r in fixtures("regimes").values() if r["n"] > 1]
        confirmation.validate(data, old)
        with self.assertRaises(ValueError):
            confirmation.validate(data, old + [data["fixtures"][0]["n"]])
        with self.assertRaises(FileExistsError):
            confirmation.generate()

    def test_fixed_relation_target_tracks_initial_and_topup_gates(self):
        collector = SimpleNamespace(
            _rows=[1, 2],
            _full=[1, 2],
            _combined=[],
            factor_base=SimpleNamespace(entries=[1] * 10),
        )
        engine = SimpleNamespace(
            collector=collector,
            _useful_target=None,
            row_excess=2,
            last_count=-1,
            filter_row_growth=2,
        )
        counters = TierA()
        counters.owner = SimpleNamespace(engine=engine, n=221, stats={})
        self.assertEqual(counters.snapshot()["target"], 13)
        self.assertEqual(counters.snapshot()["target_fraction"], 2 / 13)
        engine.last_count = 100
        self.assertEqual(counters.snapshot()["target"], 105)

    def test_progress_failure_is_observation_only(self):
        counters = TierA()
        with (
            tempfile.TemporaryDirectory() as folder,
            patch(
                "v2.benchmarks.qs.a7a.progress.thermal_status", return_value={}
            ),
        ):
            progress = Progress(
                Path(folder) / "progress.jsonl", counters, None
            )
            progress.sample()
            self.assertEqual(progress.rows, 1)
            progress.rows = 1000
            with (
                patch.object(progress.thread, "join"),
                patch.object(progress.thread, "is_alive", return_value=False),
            ):
                progress.__exit__()
            self.assertIsNotNone(progress.error)
            self.assertEqual(counters.snapshot()["stage"], "setup")
