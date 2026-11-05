"""Cofactor cost coverage and incumbent serial/RNG integration oracles."""

import copy
import hashlib
import json
import os
import random
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.ecm.b12_f4.selection_profile import build_profile
from v2.common import utils
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.ecm_selection import (
    SelectedECMContext,
    SelectionCosts,
    SelectionRule,
)
from v2.execution.ecm_workers import WorkerRun, arithmetic_options
from v2.execution.portfolio_config import PortfolioConfig
from v2.execution.worker_contract import Manifest, WorkerSettings, digest
from v2.portfolio import _advance, _canonical, factorize_bounded


def allowance(**changes):
    values = dict(work_limit=10**9, seconds=30, cpu_seconds=300)
    values.update(changes)
    return Budget(**values)


def configuration(**changes):
    values = dict(
        trial_bound=100,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=((50, 500, 8),),
        siqs=None,
        sss=None,
    )
    values.update(changes)
    return PortfolioConfig(**values)


class InterpolatedCostTests(unittest.TestCase):
    def setUp(self):
        self.settings = WorkerSettings(workers=2)
        self.config = configuration()
        self.options = arithmetic_options(self.config, self.settings)
        self.key = digest(self.options)
        self.lower = SelectionCosts(6, self.key, (2,), 5, 4, None)
        self.upper = SelectionCosts(8, self.key, (4,), 9, 8, 7)
        self.rule = SelectionRule(
            "0" * 64,
            "1" * 64,
            (self.lower, self.upper),
            (self.lower, self.upper),
            workers=2,
        )

    def decision(self, n=1000003, **changes):
        values = dict(ready=False, execution="1" * 64)
        values.update(changes)
        return self.rule.decision(
            Manifest(n, 7, "root", self.options), **values
        )

    def test_known_linear_rate_and_bundle_crossover(self):
        result = self.decision()
        self.assertEqual(result["estimate"], 24)
        self.assertEqual(
            result["threshold"], 6
        )  # Two curves per dispatch wave.
        self.assertEqual(result["anchor_digits"], [6, 8])
        self.assertEqual(result["threshold_anchor"], 6)  # Lower wins a tie.
        self.assertEqual(result["size_model"], "bracket-linear-v1")
        self.assertTrue(result["engage"])
        self.assertFalse(self.decision(remaining=(1,))["engage"])
        self.assertFalse(self.decision(pool_started=True)["engage"])

    def test_exact_costs_keep_literal_thresholds(self):
        result = self.decision(n=100003)
        self.assertEqual(result["estimate"], 16)
        self.assertEqual(result["threshold"], 5)
        self.assertEqual(result["size_model"], "measured")

    def test_crossovers_cannot_be_borrowed_from_another_bundle(self):
        costly = SelectionCosts(6, self.key, (2,), 100, 100, 100)
        self.rule = SelectionRule(
            "0" * 64,
            "1" * 64,
            (costly,),
            (self.lower, self.upper),
            workers=2,
        )
        result = self.decision()
        self.assertEqual(result["estimate"], 24)
        self.assertEqual(result["threshold"], 100)
        self.assertFalse(result["engage"])

    def test_no_extrapolation_different_options_or_wide_gaps(self):
        for n in (10007, 100000007):
            self.assertFalse(self.decision(n)["engage"])
        distant = SelectionCosts(30, self.key, (4,), 9, 8, 7)
        self.rule = SelectionRule(
            "0" * 64,
            "1" * 64,
            (self.lower,),
            (self.lower, distant),
            workers=2,
        )
        self.assertFalse(self.decision()["engage"])
        changed = Manifest(1000003, 7, "root", dict(self.options, gcd_batch=7))
        self.assertFalse(self.rule.decision(changed, ready=True)["engage"])

    def test_invalid_rates_and_execution_fail_closed(self):
        for bad in ((), (0,), (float("nan"),), (1, 2)):
            upper = SelectionCosts(8, self.key, bad, 9, 8, 7)
            self.rule = SelectionRule(
                "0" * 64,
                "1" * 64,
                (self.lower,),
                (self.lower, upper),
                workers=2,
            )
            self.assertFalse(self.decision()["engage"])
        self.assertFalse(self.decision(execution="2" * 64)["engage"])

    def test_chain_eligibility_boundary_cannot_be_interpolated(self):
        for chain_mode, expected in (("auto", False), ("off", True)):
            config = configuration(
                ecm_tiers=((2000, 5000, 8),),
                ecm_chain_mode=chain_mode,
            )
            options = arithmetic_options(config, self.settings)
            key = digest(options)
            anchors = tuple(
                SelectionCosts(digits, key, (10,), 1, 1, 1)
                for digits in (80, 100)
            )
            rule = SelectionRule(
                "0" * 64,
                "1" * 64,
                anchors,
                anchors,
                workers=2,
            )
            result = rule.decision(
                Manifest(10**84 + 7, 7, "root", options),
                ready=True,
                execution="1" * 64,
            )
            self.assertEqual(result["engage"], expected)

    def test_cofactor_value_seed_and_lineage_do_not_predict_outcomes(self):
        first = self.decision()
        other = Manifest(1000033, 991, "other", self.options)
        self.assertEqual(
            first,
            self.rule.decision(other, ready=False, execution="1" * 64),
        )


class FakePool:
    """Independent lifecycle boundary, with no child or ECM arithmetic."""

    def __init__(self, *, seed, settings, rule, force_workers=False):
        self.seed, self.settings, self.rule = seed, settings, rule
        self.force_workers = force_workers
        self.helper_funded = False
        self.lifetime_stats = {}
        self.calls = []

    def pool_started(self):
        return False

    def ready(self, manifest):
        return False

    def run(self, n, *, budget, checkpoint, **kwargs):
        self.calls.append(n)
        if checkpoint is not None:
            assert budget.used >= checkpoint["budget"]["work"]
        budget.consume(7)
        base = getattr(budget, "budget", budget)
        snapshot = dict(
            sha256="a" * 64,
            budget=dict(
                work=base.used, cpu=base.cpu_used, wall=base.wall_used
            ),
        )
        return WorkerRun(
            None, "curves_exhausted", snapshot, dict(parent_cpu=0)
        )

    def close(self):
        pass


class PortfolioSelectionTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "profile.json"
        self.config = configuration()
        self.write_profile()

    def write_profile(self, *, threshold=1, rates=(1, 3)):
        rows = []
        for digits, workers, rate in zip((6, 8), (2, 4), rates):
            settings = WorkerSettings(workers=workers, memory_bytes=8 * 2**30)
            candidate = dict(
                settings=asdict(settings),
                costs=dict(
                    seconds_per_curve=(rate,),
                    cold_seconds=threshold,
                    warm_seconds=threshold,
                    pool_seconds=None,
                ),
            )
            rows.append(
                dict(
                    digits=digits,
                    config=asdict(self.config),
                    cold=candidate,
                    warm=candidate,
                )
            )
        profile = build_profile(
            dict(
                serial_settings=asdict(WorkerSettings(memory_bytes=8 * 2**30)),
                rows=rows,
            )
        )
        self.path.write_text(json.dumps(profile))
        self.sha = hashlib.sha256(self.path.read_bytes()).hexdigest()

    def owner(self):
        return SelectedECMContext(self.path, self.sha, seed=7)

    def test_catalogue_retains_cost_anchors_and_frozen_v1_rejects(self):
        root = (
            Path(__file__).resolve().parents[1] / "benchmarks/inputs/controls"
        )
        for version in (1, 2):
            path = root / f"b12_f4_selection_v{version}.json"
            sha = hashlib.sha256(path.read_bytes()).hexdigest()
            with SelectedECMContext(path, sha) as owner:
                if version == 1:
                    with self.assertRaises(ValueError):
                        owner._load(allowance())
                    self.assertEqual(
                        sha,
                        "1dc664ac2d8c2b22eee39d53efe9212a9efd2d8e98e8f8668ad38501d9615ecc",
                    )
                else:
                    data = json.loads(path.read_text())
                    original = json.loads(
                        (root / "b12_f4_selection_v1.json").read_text()
                    )
                    self.assertEqual(data["version"], "b12-f4-selection-v2")
                    self.assertEqual(data["status"], "training")
                    self.assertEqual(data["rows"], original["rows"])
                    self.assertEqual(
                        data["serial_settings"], original["serial_settings"]
                    )
                    self.assertIsNone(owner.inner)

    def test_unknown_first_call_does_not_pin_a_serial_pool(self):
        with self.owner() as owner:
            with patch("v2.execution.ecm_workers.ECMExecutor", FakePool):
                self.assertIsNone(
                    owner.run_portfolio(
                        10007,
                        config=self.config,
                        budget=allowance(),
                    )
                )
                self.assertIsNone(owner.inner)
                self.assertIsNone(owner.origin)
                result = owner.run_portfolio(
                    1000003,
                    config=self.config,
                    budget=allowance(),
                )
                self.assertEqual(result.reason, "curves_exhausted")
                self.assertEqual(owner.inner.settings.workers, 2)
                self.assertTrue(owner.inner.force_workers)
                self.assertEqual(owner.inner.calls, [1000003])

    def test_retained_pool_keeps_its_bundle_and_declines_unknown_options(self):
        with self.owner() as owner:
            with patch("v2.execution.ecm_workers.ECMExecutor", FakePool):
                owner.run_portfolio(
                    1000003, config=self.config, budget=allowance()
                )
                pool = owner.inner
                owner.run_portfolio(
                    10000019, config=self.config, budget=allowance()
                )
                self.assertIs(owner.inner, pool)
                self.assertEqual(pool.settings.workers, 2)
                changed = configuration(gcd_batch=7)
                self.assertIsNone(
                    owner.run_portfolio(
                        1000003,
                        config=changed,
                        budget=allowance(),
                    )
                )
                self.assertEqual(pool.calls, [1000003, 10000019])

    def test_unfunded_start_stays_incumbent_without_constructing_executor(
        self,
    ):
        with self.owner() as owner:
            with patch("v2.execution.ecm_workers.ECMExecutor") as factory:
                self.assertIsNone(
                    owner.run_portfolio(
                        1000003,
                        config=self.config,
                        budget=allowance(cpu_seconds=3),
                    )
                )
                factory.assert_not_called()
                self.assertFalse(owner.last_stats["engagement"]["funded"])
                self.assertTrue(owner.last_stats["legacy_serial"])

    def test_unknown_large_quota_declines_before_worker_manifest_limits(self):
        config = configuration(ecm_tiers=((2, 2, 4097),))
        with self.owner() as owner:
            self.assertIsNone(
                owner.run_portfolio(
                    1000003,
                    config=config,
                    budget=allowance(),
                )
            )
            self.assertIsNone(owner.inner)

    def test_worker_ownership_refusal_returns_to_incumbent(self):
        profile = json.loads(self.path.read_text())
        for row in profile["rows"]:
            for mode in ("cold", "warm"):
                row[mode]["settings"].update(
                    parent_memory=16 * 2**20,
                    memory_bytes=512 * 2**20,
                )
        self.path.write_text(json.dumps(profile))
        self.sha = hashlib.sha256(self.path.read_bytes()).hexdigest()
        with self.owner() as owner:
            self.assertIsNone(
                owner.run_portfolio(
                    1000003,
                    config=self.config,
                    budget=allowance(),
                )
            )
            self.assertIsNone(owner.inner)
            self.assertEqual(
                owner.last_stats["engagement"]["reason"],
                "worker_ownership_refused",
            )

    def test_zero_expired_and_cancelled_allowances_do_not_load_profile(self):
        for budget in (
            allowance(work_limit=0),
            allowance(seconds=0),
            allowance(cancelled=lambda: True),
        ):
            with self.owner() as owner:
                with self.assertRaises(BudgetExhaustedError):
                    owner.run_portfolio(
                        1000003, config=self.config, budget=budget
                    )
                self.assertIsNone(owner.catalogue)
                self.assertIsNone(owner.inner)

    def test_component_and_portfolio_modes_cannot_be_mixed(self):
        with self.owner() as owner:
            owner.run_portfolio(10007, config=self.config, budget=allowance())
            with self.assertRaises(ValueError):
                owner.run(1000003, config=self.config, budget=allowance())

    def test_engaged_resume_binds_mode_settings_and_prior_charges(self):
        with patch("v2.execution.ecm_workers.ECMExecutor", FakePool):
            with self.owner() as owner:
                first = owner.run_portfolio(
                    1000003,
                    config=self.config,
                    budget=allowance(),
                )
            with self.owner() as owner:
                resumed = owner.run_portfolio(
                    1000003,
                    config=self.config,
                    budget=allowance(),
                    checkpoint=first.checkpoint,
                )
                self.assertTrue(owner.inner.force_workers)
                self.assertGreater(
                    resumed.checkpoint["header"]["budget"]["work"],
                    first.checkpoint["header"]["budget"]["work"],
                )
            for mutate in (
                lambda h: h.update(mode="component"),
                lambda h: h.update(version="b12-f4-selection-v1"),
                lambda h: h["settings"].update(workers=4),
            ):
                damaged = copy.deepcopy(first.checkpoint)
                mutate(damaged["header"])
                damaged["sha256"] = digest(damaged["header"])
                with self.owner() as owner:
                    with self.assertRaises(ValueError):
                        owner.run_portfolio(
                            1000003,
                            config=self.config,
                            budget=allowance(),
                            checkpoint=damaged,
                        )

    def test_serial_curves_results_and_rng_match_incumbent(self):
        self.write_profile(threshold=100000)
        n = 1009 * 1013
        control = factorize_bounded(
            n, seed=29, config=self.config, budget=allowance()
        )
        with self.owner() as owner:
            candidate = factorize_bounded(
                n,
                seed=29,
                config=self.config,
                budget=allowance(),
                ecm_executor=owner,
            )
            self.assertIsNone(owner.inner)
        self.assertEqual(candidate.result, control.result)
        self.assertEqual(candidate.result.reconstruct(), n)
        self.assertEqual(
            candidate.checkpoint["payload"]["rng"],
            control.checkpoint["payload"]["rng"],
        )

        def curves(run):
            return [
                (e["n"], e["seed"], e["b1"], e["b2"], e["outcome"])
                for e in run.events
                if e["stage"] == "ecm"
            ]

        self.assertTrue(curves(control))
        self.assertEqual(curves(candidate), curves(control))

    def test_serial_resume_is_sticky_and_rejects_corrupt_branch_marker(self):
        self.write_profile(threshold=100000)
        n = 1009 * 1013
        with self.owner() as owner:
            paused = factorize_bounded(
                n,
                seed=29,
                config=self.config,
                budget=allowance(
                    cancelled=lambda: owner.last_stats.get(
                        "legacy_serial", False
                    )
                ),
                ecm_executor=owner,
            )
        self.assertEqual(paused.reason, "cancelled")
        self.assertTrue(
            paused.checkpoint["payload"]["state"]["current"]["ecm_serial"]
        )
        with self.owner() as owner:
            with patch.object(
                owner,
                "run_portfolio",
                side_effect=AssertionError("must stay serial"),
            ):
                resumed = factorize_bounded(
                    n,
                    config=self.config,
                    budget=allowance(),
                    checkpoint=paused.checkpoint,
                    ecm_executor=owner,
                )
        self.assertEqual(resumed.result.reconstruct(), n)
        for bad in (False, "yes", 1):
            damaged = copy.deepcopy(paused.checkpoint)
            damaged["payload"]["state"]["current"]["ecm_serial"] = bad
            damaged["sha256"] = hashlib.sha256(
                _canonical(damaged["payload"]).encode()
            ).hexdigest()
            with self.owner() as owner:
                with self.assertRaises(ValueError):
                    factorize_bounded(
                        n,
                        config=self.config,
                        budget=allowance(),
                        checkpoint=damaged,
                        ecm_executor=owner,
                    )

    def test_cancellation_after_selection_never_switches_to_serial_stream(
        self,
    ):
        with self.owner() as owner:
            cancelled = [False]
            original = owner._portfolio_decision

            def choose_then_cancel(*args):
                decision = original(*args)
                self.assertTrue(decision["engage"])
                cancelled[0] = True
                return decision

            with patch.object(
                owner, "_portfolio_decision", side_effect=choose_then_cancel
            ):
                result = owner.run_portfolio(
                    1000003,
                    config=self.config,
                    budget=allowance(cancelled=lambda: cancelled[0]),
                )
            self.assertEqual(result.reason, "cancelled")
            self.assertTrue(owner.inner.force_workers)
            self.assertIsNone(owner.inner.serial_engine)
            self.assertEqual(result.checkpoint["header"]["mode"], "portfolio")
            self.assertEqual(
                result.checkpoint["inner"]["search"]["payload"]["completed"],
                {},
            )

    def test_real_spawn_uses_actual_cofactor_and_cleans_up(self):
        n = 1009 * 1013
        with self.owner() as owner:
            run = factorize_bounded(
                n * 2**10,
                seed=29,
                config=self.config,
                budget=allowance(),
                ecm_executor=owner,
            )
            self.assertIsNotNone(owner.inner)
            self.assertEqual(owner.inner.settings.workers, 2)
            self.assertTrue(owner.inner.force_workers)
            workers = [slot["process"].pid for slot in owner.inner.slots]
            self.assertTrue(workers)
            events = [
                event
                for event in run.events
                if event["stage"] == "ecm_workers"
            ]
            self.assertTrue(events)
            self.assertTrue(all(event["n"] == n for event in events))
            self.assertEqual(
                events[0]["stats"]["engagement"]["anchor_digits"], [6, 8]
            )
        self.assertEqual(run.result.reconstruct(), n * 2**10)
        for factor in run.result.factors:
            self.assertIn(factor.value, (2, 1009, 1013))
        for pid in workers:
            with self.assertRaises(ChildProcessError):
                os.waitpid(pid, os.WNOHANG)


class FallbackStreamTests(unittest.TestCase):
    def test_serial_decline_preserves_curve_and_later_siqs_seed(self):
        from v2.qs import SIQSConfig

        config = configuration(
            ecm_tiers=((2, 2, 3), (3, 3, 2)),
            siqs=SIQSConfig(),
            memory_bytes=2**31,
        )

        def execute(selected):
            state = dict(
                original=1009 * 1013,
                current=dict(
                    n=1009 * 1013,
                    mult=1,
                    stage="ecm",
                    attempt=0,
                    tier=0,
                    job=None,
                ),
                remaining=[],
                events=[],
                dropped_events=0,
                stage_seconds={},
            )
            rng = random.Random(12345)
            budget = allowance()
            seeds = []

            def new_job(kind, n, seed, b1, b2):
                seeds.append(seed)
                return dict(
                    kind=kind,
                    n=n,
                    seed=seed,
                    b1=b1,
                    b2=b2,
                    phase="done",
                    done=True,
                    factor=None,
                )

            class FinishedFallback:
                finished_reason = "prefix_exhausted"

                def __init__(self, n, *, seed, **kwargs):
                    self.n, self.seed = n, seed
                    seeds.append(seed)

                def run(self, **kwargs):
                    return SimpleNamespace(
                        divisor=None, reason=self.finished_reason, stats={}
                    )

            with (
                patch("v2.portfolio.new_job", side_effect=new_job),
                patch("v2.portfolio.advance_job"),
                patch("v2.portfolio.SIQSJob", FinishedFallback),
            ):
                for _ in range(30):
                    if state["current"] is None:
                        break
                    _advance(
                        state,
                        config,
                        budget,
                        None,
                        rng,
                        {},
                        None,
                        utils.PRIMALITY_POLICY,
                        selected,
                    )
            self.assertIsNone(state["current"])
            self.assertEqual(state["remaining"], [1009 * 1013])
            return seeds, rng.getstate()

        selector = SimpleNamespace(
            settings=SimpleNamespace(parent_memory=4 * 2**30),
            run=lambda *a, **kw: self.fail("component stream must not run"),
            run_portfolio=lambda *a, **kw: None,
            last_stats={},
        )
        control = execute(None)
        self.assertEqual(len(control[0]), 6)  # Five curves, then SIQS seed.
        self.assertEqual(execute(selector), control)


if __name__ == "__main__":
    unittest.main()
