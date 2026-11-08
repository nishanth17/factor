"""Observable selection, lazy ownership and adversarial drained resume."""

import copy
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.execution.budget import BudgetExhaustedError
from v2.execution.ecm_selection import (
    COMPATIBILITY,
    HEADER_BYTES,
    PROFILE_VERSION,
    SelectedECMContext,
    SelectionCosts,
    SelectionRule,
    _saved_budget,
)
from v2.execution.ecm_workers import WorkerRun, arithmetic_options
from v2.execution.portfolio_config import PortfolioConfig
from v2.execution.worker_contract import (
    SEARCH_VERSION,
    VERSION,
    Manifest,
    WorkerSettings,
    digest,
)
from v2.tests.resource_stops import InterruptBudget as Budget


class SelectionRuleTests(unittest.TestCase):
    def setUp(self):
        self.config = PortfolioConfig(
            ecm_tiers=((50, 500, 4),), siqs=None, sss=None
        )
        self.manifest = Manifest(
            1000003,
            7,
            "root",
            arithmetic_options(self.config, WorkerSettings()),
        )
        self.model = SelectionCosts(
            7, digest(self.manifest.options), (1.0,), 5.0, 2.0, None
        )
        self.rule = SelectionRule("0" * 64, "1" * 64, (self.model,))

    def decision(self, **kwargs):
        return self.rule.decision(self.manifest, execution="1" * 64, **kwargs)

    def test_cold_warm_and_unmeasured_metadata_switch(self):
        self.assertFalse(self.decision(ready=False)["engage"])
        self.assertTrue(self.decision(ready=True)["engage"])
        self.assertFalse(
            self.decision(ready=False, pool_started=True)["engage"]
        )
        self.assertFalse(self.decision(ready=True, remaining=(1,))["engage"])

    def test_invalid_estimates_are_serial_and_no_outcomes_are_inputs(self):
        for invalid in (0, -1, math.inf, math.nan):
            model = SelectionCosts(
                7, digest(self.manifest.options), (invalid,), 1, 1, 1
            )
            rule = SelectionRule("0" * 64, "1" * 64, (model,))
            self.assertFalse(
                rule.decision(self.manifest, ready=True, execution="1" * 64)[
                    "engage"
                ]
            )
        # Different moduli of the same declared regime yield the same choice;
        # no factor labels, selected sigma, or future factor result is read.
        other = Manifest(1000033, 99, "other", self.manifest.options)
        self.assertEqual(
            self.decision(ready=True),
            self.rule.decision(other, ready=True, execution="1" * 64),
        )

    def test_unmeasured_size_and_changed_execution_remain_serial(self):
        other = Manifest(10000019, 7, "root", self.manifest.options)
        self.assertFalse(self.rule.decision(other, ready=True)["engage"])
        self.assertFalse(
            self.rule.decision(self.manifest, ready=True, execution="2" * 64)[
                "engage"
            ]
        )


class FakeExecutor:
    """Deterministic boundary oracle; no arithmetic or worker startup."""

    def __init__(self, **kwargs):
        self.settings = kwargs["settings"]
        self.rule = kwargs["rule"]
        self.seed = kwargs["seed"]
        self.lifetime_stats = {}
        self.closed = False

    def run(self, n, *, budget, checkpoint, **kwargs):
        if checkpoint is not None:
            saved = checkpoint["budget"]
            assert budget.used >= saved["work"]
            assert budget.cpu_used >= saved["cpu"]
            assert budget.wall_used >= saved["wall"]
        budget.consume(7)
        base = getattr(budget, "budget", budget)
        inner = dict(
            budget=dict(
                work=base.used, cpu=base.cpu_used, wall=base.wall_used
            ),
            sha256="a" * 64,
        )
        self.lifetime_stats["charged_work"] = (
            self.lifetime_stats.get("charged_work", 0) + 7
        )
        return WorkerRun(None, "curves_exhausted", inner, dict(parent_cpu=0))

    def close(self):
        self.closed = True


class LazySelectionTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "model.json"
        self.config = PortfolioConfig(
            ecm_tiers=((50, 500, 4),), siqs=None, sss=None
        )
        self.settings = WorkerSettings(workers=2)
        options = digest(arithmetic_options(self.config, self.settings))
        costs = asdict(SelectionCosts(7, options, (1.0,), 2, 1, None))
        candidate = dict(settings=asdict(self.settings), costs=costs)
        self.model = dict(
            version=PROFILE_VERSION,
            status="training",
            executor_version=VERSION,
            search_version=SEARCH_VERSION,
            compatibility=COMPATIBILITY,
            serial_settings=asdict(WorkerSettings()),
            rows=[
                dict(
                    digits=7,
                    options=options,
                    tiers=[[50, 500, 4]],
                    cold=candidate,
                    warm=candidate,
                )
            ],
        )
        self.save()

    def save(self):
        self.path.write_text(json.dumps(self.model))
        self.sha = hashlib.sha256(self.path.read_bytes()).hexdigest()

    def owner(self, **kwargs):
        return adapt_resource_configuration(
            SelectedECMContext, self.path, self.sha, seed=7, **kwargs
        )

    def run_fake(self, owner, *, checkpoint=None, budget=None):
        with patch("v2.execution.ecm_workers.ECMExecutor", FakeExecutor):
            return owner.run(
                1000003,
                config=self.config,
                budget=budget or self.budget(),
                checkpoint=checkpoint,
            )

    def budget(self):
        return Budget(work_limit=10**9, seconds=20, cpu_seconds=20)

    def test_constructor_does_not_read_model_or_import_multiprocessing(self):
        code = (
            "import sys;"
            "from v2.execution.ecm_selection import SelectedECMContext;"
            "c=SelectedECMContext('/missing/model','0'*64);"
            "assert c.inner is None and c.catalogue is None;"
            "assert 'multiprocessing' not in sys.modules;c.close()"
        )
        result = subprocess.run(
            [sys.executable, "-B", "-c", code], capture_output=True, timeout=10
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_once_per_context_preparation_and_bounded_checkpoint(self):
        owner = self.owner()
        first = self.run_fake(owner)
        inner = owner.inner
        self.assertEqual(inner.settings.workers, 2)
        self.assertLess(inner.settings.checkpoint_bytes, 8 * 2**20)
        self.assertGreater(first.stats["selector_work"], HEADER_BYTES)
        self.path.unlink()
        self.run_fake(owner)
        self.assertIs(owner.inner, inner)
        owner.close()
        self.assertTrue(inner.closed)
        self.assertIsNone(owner.catalogue)
        with self.assertRaises(RuntimeError):
            self.run_fake(owner)

    def test_zero_cancelled_and_expired_allowances_do_no_preparation(self):
        for budget in (
            Budget(work_limit=0),
            Budget(cancelled=lambda: True),
            Budget(seconds=0),
        ):
            owner = self.owner()
            with self.assertRaises(BudgetExhaustedError):
                self.run_fake(owner, budget=budget)
            self.assertIsNone(owner.inner)
            self.assertIsNone(owner.catalogue)

    def test_drained_restart_preserves_charges_and_original_count(self):
        owner = self.owner()
        first_budget = self.budget()
        first = self.run_fake(owner, budget=first_budget)
        owner.close()
        restored = self.owner()
        budget = self.budget()
        second = self.run_fake(
            restored, checkpoint=first.checkpoint, budget=budget
        )
        self.assertGreater(budget.used, first_budget.used + HEADER_BYTES)
        self.assertEqual(restored.inner.settings.workers, 2)
        self.assertEqual(second.reason, first.reason)
        smaller = Budget(work_limit=1, seconds=20, cpu_seconds=20)
        with self.assertRaises(ValueError):
            self.run_fake(
                self.owner(), checkpoint=first.checkpoint, budget=smaller
            )

    def test_changed_header_identity_inner_or_settings_are_rejected(self):
        first = self.run_fake(self.owner())
        for field in ("identity", "inner"):
            damaged = copy.deepcopy(first.checkpoint)
            damaged["header"][field] = "b" * 64
            damaged["sha256"] = digest(damaged["header"])
            with self.assertRaises(ValueError):
                self.run_fake(self.owner(), checkpoint=damaged)
        changed = self.owner(prefer_warm=True)
        with self.assertRaises(ValueError):
            self.run_fake(changed, checkpoint=first.checkpoint)
        damaged = copy.deepcopy(first.checkpoint)
        damaged["header"]["settings"]["workers"] = 4
        damaged["sha256"] = digest(damaged["header"])
        with self.assertRaises(ValueError):
            self.run_fake(self.owner(), checkpoint=damaged)
        current = self.owner()
        self.run_fake(current)
        with self.assertRaises(ValueError):
            self.run_fake(current, checkpoint=damaged)

    def test_oversized_stale_duplicate_and_hidden_rows_fail_closed(self):
        for mutation in (
            lambda m: m.update(rows=m["rows"] * 257),
            lambda m: m["rows"].append(copy.deepcopy(m["rows"][0])),
            lambda m: m["rows"][0].update(known_factor=1009),
            lambda m: m.update(compatibility="unsupported"),
        ):
            original = copy.deepcopy(self.model)
            mutation(self.model)
            self.save()
            with self.assertRaises(ValueError):
                self.run_fake(self.owner())
            self.model = original
        self.save()
        owner = self.owner()
        self.path.write_text("{}")
        with self.assertRaises(ValueError):
            self.run_fake(owner)

    def test_unbounded_header_fields_reject_before_canonical_encoding(self):
        first = self.run_fake(self.owner())
        for mutate in (
            lambda h: h.update(origin=[0] * 100000),
            lambda h: h["settings"].update({"x" * 100000: 0}),
        ):
            damaged = copy.deepcopy(first.checkpoint)
            mutate(damaged["header"])
            with patch(
                "v2.execution.ecm_selection.canonical",
                side_effect=AssertionError("must reject before encoding"),
            ):
                with self.assertRaises(ValueError):
                    self.run_fake(self.owner(), checkpoint=damaged)

    def test_live_settings_cannot_resize_on_resume(self):
        owner = self.owner()
        first = self.run_fake(owner)
        owner.seed += 1
        with self.assertRaises(ValueError):
            self.run_fake(owner, checkpoint=first.checkpoint)

    def test_restore_never_reduces_already_owned_charges(self):
        budget = self.budget()
        budget.used = 100
        _saved_budget(budget, dict(work=50, cpu=0, wall=0))
        self.assertEqual(budget.used, 150)
        _saved_budget(budget, dict(work=101, cpu=0, wall=0))
        self.assertEqual(budget.used, 201)
        _saved_budget(budget, dict(work=101, cpu=0, wall=0))
        self.assertEqual(budget.used, 201)


class ArithmeticSelectionTests(unittest.TestCase):
    def test_selected_spawn_matches_serial_curves_and_checked_resume(self):
        from v2.benchmarks.ecm.b12_f4.profile_runner import signature
        from v2.benchmarks.ecm.b12_f4.selection_profile import build_profile
        from v2.execution.ecm_workers import ECMExecutor

        n = 1009 * 1013
        config = PortfolioConfig(
            ecm_tiers=((50, 500, 8),), siqs=None, sss=None
        )
        settings = WorkerSettings(workers=2)
        costs = dict(
            seconds_per_curve=(1,),
            cold_seconds=0.001,
            warm_seconds=0.001,
            pool_seconds=None,
        )
        reviewed = dict(
            serial_settings=asdict(WorkerSettings()),
            rows=[
                dict(
                    digits=len(str(n)),
                    config=asdict(config),
                    cold=dict(settings=asdict(settings), costs=costs),
                    warm=dict(settings=asdict(settings), costs=costs),
                )
            ],
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test-only-costs.json"
            path.write_text(json.dumps(build_profile(reviewed)))
            sha = hashlib.sha256(path.read_bytes()).hexdigest()
            with SelectedECMContext(path, sha, seed=7) as selected:
                result = selected.run(
                    n,
                    config=config,
                    budget=Budget(
                        work_limit=10**9, seconds=20, cpu_seconds=300
                    ),
                    fixed_work=True,
                )
                self.assertTrue(result.stats["engaged"])
                pids = [slot["process"].pid for slot in selected.inner.slots]
            self.assertEqual(selected.inner.slots, [])
            for pid in pids:
                with self.assertRaises(ChildProcessError):
                    os.waitpid(pid, os.WNOHANG)
            with ECMExecutor(seed=7, settings=WorkerSettings()) as serial:
                control = serial.run(
                    n,
                    config=config,
                    budget=Budget(
                        work_limit=10**9, seconds=20, cpu_seconds=300
                    ),
                    fixed_work=True,
                )
            rows = signature(result.checkpoint["inner"], n, 8)
            self.assertEqual(rows, signature(control.checkpoint, n, 8))
            # The independent trial-division oracle checks every discovered
            # divisor, including setup and stage-one/stage-two outcomes.
            for row in rows:
                if row[2] is not None:
                    self.assertIn(row[2], (1009, 1013))
            with SelectedECMContext(path, sha, seed=7) as restored:
                resumed = restored.run(
                    n,
                    config=config,
                    budget=Budget(
                        work_limit=10**9, seconds=30, cpu_seconds=600
                    ),
                    checkpoint=result.checkpoint,
                    fixed_work=True,
                )
                self.assertEqual(
                    rows, signature(resumed.checkpoint["inner"], n, 8)
                )
                self.assertGreater(
                    resumed.checkpoint["header"]["budget"]["work"],
                    result.checkpoint["header"]["budget"]["work"],
                )


if __name__ == "__main__":
    unittest.main()
