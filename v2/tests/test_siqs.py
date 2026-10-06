"""Full SIQS shared stores, finite controls, compact resume and dispatch."""

import copy
import hashlib
import json
import unittest
from unittest.mock import patch

from v2.budget import Budget, BudgetExhaustedError
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import SieveConfig, SIQSConfig, SIQSJob
from v2.qs.families import _checksum
from v2.qs.multiplier import select_multiplier
from v2.qs.pipeline import DependencyExtractor, DependencySolver


def allowance(work=200_000_000):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


def configuration(**options):
    values = dict(
        base_bound=200,
        half_width=256,
        collector=SieveConfig(
            residual_bound=500,
            max_atoms=4096,
            max_relations=4096,
            max_partials=4096,
        ),
    )
    values.update(options)
    return SIQSConfig(**values)


def mutate(checkpoint, change):
    result = copy.deepcopy(checkpoint)
    payload = json.loads(result["blob"])
    change(payload)
    result["blob"] = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    result["sha256"] = hashlib.sha256(result["blob"].encode()).hexdigest()
    return result


class SIQSTests(unittest.TestCase):
    def test_shared_families_complete_previously_exhausted_inputs(self):
        for n in (13643477, 16187767, 18513647, 4001 * 5003):
            job = SIQSJob(n, config=configuration(), budget=allowance())

            result = job.run()

            self.assertEqual(result.reason, "factor_found")
            self.assertGreater(result.divisor, 1)
            self.assertEqual(result.divisor * result.cofactor, n)
            self.assertGreater(result.stats["polynomials"], 1)
            polynomials = {
                a.polynomial.identity
                for a in job.engine.collector._atoms.values()
            }

            self.assertGreater(len(polynomials), 1)
            self.assertLessEqual(
                result.stats["workspace_bytes"], job.config.memory_bytes
            )
            used = job.budget.used

            self.assertEqual(job.run().divisor, result.divisor)
            self.assertEqual(job.budget.used, used)

    def test_block_pause_roundtrip_keeps_positions_rows_seed_and_resources(
        self,
    ):
        config = configuration()
        original = SIQSJob(
            4001 * 5003, seed=73, config=config, budget=allowance()
        )

        paused = original.run(max_blocks=1)

        self.assertEqual(paused.reason, "paused")
        checkpoint = json.loads(json.dumps(original.checkpoint()))

        self.assertLess(
            len(checkpoint["blob"]) + 1024, config.checkpoint_bytes
        )

        restored = SIQSJob.from_checkpoint(
            checkpoint, budget=allowance(), config=config
        )

        self.assertEqual(restored.seed, original.seed)
        self.assertEqual(
            restored.engine.next_position, original.engine.next_position
        )
        self.assertEqual(
            restored.engine.collector._atoms, original.engine.collector._atoms
        )
        self.assertGreater(
            restored.budget.used, checkpoint["resources"]["work_used"]
        )
        self.assertGreaterEqual(
            restored.budget.prior_wall, checkpoint["resources"]["wall_used"]
        )
        self.assertEqual(restored.run().divisor, original.run().divisor)

    def test_every_block_can_checkpoint_until_complete(self):
        job = SIQSJob(
            18513647, seed=81, config=configuration(), budget=allowance()
        )

        for _ in range(64):
            result = job.run(max_blocks=1)
            if result.divisor:
                self.assertEqual(result.divisor * result.cofactor, job.n)
                break
            prior = job.budget.used

            job = SIQSJob.from_checkpoint(job.checkpoint(), budget=allowance())

            self.assertGreaterEqual(job.budget.used, prior)
        else:
            self.fail("bounded small fixture did not complete through resumes")

    def test_work_refusal_is_reported_and_resumable(self):
        job = SIQSJob(4001 * 5003, config=configuration(), budget=allowance(0))

        result = job.run(max_blocks=1)

        self.assertEqual(result.reason, "work_limit")
        self.assertEqual(result.cofactor, job.n)

        restored = SIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )

        self.assertEqual(restored.run().reason, "factor_found")

    def test_cancelled_and_deadline_jobs_keep_checked_cofactor(self):
        for budget in (
            Budget(cancelled=lambda: True),
            Budget(seconds=0),
            Budget(cpu_seconds=0),
        ):
            job = SIQSJob(4001 * 5003, config=configuration(), budget=budget)

            result = job.run(max_blocks=1)

            self.assertIn(
                result.reason, ("cancelled", "wall_limit", "cpu_limit")
            )
            self.assertEqual(result.cofactor, job.n)

            restored = SIQSJob.from_checkpoint(
                job.checkpoint(), budget=allowance()
            )

            self.assertEqual(restored.run().reason, "factor_found")

    def test_compact_checkpoint_replays_pending_solver_and_extraction(self):
        for cls in (DependencySolver, DependencyExtractor):
            job = SIQSJob(
                4001 * 5003,
                config=configuration(mode="qs"),
                budget=allowance(),
            )

            def refuse(state):
                if isinstance(state, DependencySolver):
                    state.step()
                job.budget.reason = "work_limit"
                raise BudgetExhaustedError("work_limit")

            with patch.object(cls, "run", refuse):
                result = job.run()

            self.assertEqual(result.reason, "work_limit")
            self.assertIsNotNone(job.engine.solver)
            checkpoint = job.checkpoint()

            self.assertNotIn('"pivots"', checkpoint["blob"])

            restored = SIQSJob.from_checkpoint(checkpoint, budget=allowance())

            self.assertIsNotNone(restored.engine.solver)
            self.assertEqual(restored.run().reason, "factor_found")
            self.assertEqual(restored.stats["matrix_reconstructions"], 1)

    def test_corrupt_envelope_resource_and_store_identity_rejected(self):
        job = SIQSJob(4001 * 5003, config=configuration(), budget=allowance())
        job.run(max_blocks=1)
        checkpoint = job.checkpoint()
        for field in ("blob", "sha256", "resources_sha256"):
            corrupt = copy.deepcopy(checkpoint)
            corrupt[field] += "x"
            with self.assertRaises(ValueError):
                SIQSJob.from_checkpoint(corrupt, budget=allowance())

        cases = (
            lambda p: p.update(seed=-1),
            lambda p: p.update(epoch=999),
            lambda p: p.update(gray_index=999),
            lambda p: p.update(half_width=999),
            lambda p: p.update(base_identity="x" * 64),
            lambda p: p["engine"].update(next_position=9999),
            lambda p: p["store"]["atoms"][0].__setitem__(4, 499),
        )
        for change in cases:
            with self.assertRaises((ValueError, TypeError)):
                SIQSJob.from_checkpoint(
                    mutate(checkpoint, change), budget=allowance()
                )

    def test_recomputed_checksums_still_require_exact_relation_math(self):
        job = SIQSJob(4001 * 5003, config=configuration(), budget=allowance())
        job.run(max_blocks=1)

        def corrupt(payload):
            payload["store"]["atoms"][0][2] *= -1
            payload["store_identity"] = _checksum(payload["store"])

        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                mutate(job.checkpoint(), corrupt), budget=allowance()
            )

    def test_config_and_consumed_budget_cannot_be_reset(self):
        job = SIQSJob(4001 * 5003, config=configuration(), budget=allowance())
        job.run(max_blocks=1)
        checkpoint = job.checkpoint()
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(checkpoint, budget=allowance(0))
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                checkpoint,
                budget=allowance(),
                config=configuration(base_bound=300),
            )

        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                checkpoint, budget=Budget(work_limit=200_000_000, used=1)
            )

    def test_finite_stall_and_width_recovery_preserve_prime_cofactor(self):
        config = configuration(
            base_bound=30,
            half_width=8,
            max_half_width=16,
            factor_count=1,
            family_count=1,
            max_stalled=1,
            growth_steps=1,
        )
        job = SIQSJob(1000003, config=config, budget=allowance())

        result = job.run()

        self.assertIsNone(result.divisor)
        self.assertEqual(result.cofactor, job.n)
        self.assertIn(
            result.reason,
            (
                "families_exhausted",
                "stalled_yield",
                "trivial_dependency_limit",
            ),
        )
        self.assertEqual(result.stats["recoveries"], 1)
        used = job.budget.used

        self.assertEqual(job.run().reason, result.reason)
        self.assertEqual(used, job.budget.used)

        restored = SIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )

        self.assertEqual(restored.run().cofactor, job.n)

    def test_storage_caps_attempt_final_extraction(self):
        config = configuration(
            mode="qs",
            base_bound=100,
            collector=SieveConfig(residual_bound=1, max_relations=8),
        )
        job = SIQSJob(4001 * 5003, config=config, budget=allowance())

        result = job.run()

        self.assertEqual(result.reason, "factor_found")
        self.assertEqual(result.divisor * result.cofactor, job.n)

    def test_multiplier_control_gcd_and_exact_score_order(self):
        choice = select_multiplier(15, candidates=(1, 3), budget=allowance())

        self.assertEqual(choice.divisor, 3)
        for n in (4001 * 5003, 2003 * 8353, 10403):
            choice = select_multiplier(n, budget=allowance())

            self.assertIsNone(choice.divisor)
            self.assertIn(
                (1, next(s for h, s in choice.scores if h == 1)), choice.scores
            )
            self.assertEqual(
                choice.multiplier,
                max(choice.scores, key=lambda x: (x[1], -x[0]))[0],
            )
            job = SIQSJob(
                n, config=configuration(multiplier=0), budget=allowance()
            )

            result = job.run()

            self.assertEqual(result.reason, "factor_found")
            self.assertEqual(result.divisor * result.cofactor, n)

        for candidates in ((2,), (9,), (1, 1), (), (257,)):
            with self.assertRaises(ValueError):
                select_multiplier(10403, candidates=candidates)

    def test_qs_and_multi_polynomial_mpqs_controls(self):
        for mode in ("qs", "mpqs"):
            result = SIQSJob(
                13643477, config=configuration(mode=mode), budget=allowance()
            ).run()
            self.assertEqual(result.reason, "factor_found")
            self.assertEqual(result.divisor * result.cofactor, 13643477)

    def test_legacy_portfolio_checkpoint_upgrades_without_siqs(self):
        config = PortfolioConfig(
            trial_bound=2, rho_attempts=0, pm1_attempts=0, ecm_tiers=()
        )

        first = factorize_bounded(
            4001 * 5003, config=config, budget=allowance(0)
        )
        checkpoint = copy.deepcopy(first.checkpoint)
        checkpoint["payload"]["version"] = 2
        checkpoint["payload"]["config"].pop("siqs")
        checkpoint["sha256"] = _checksum(checkpoint["payload"])

        result = factorize_bounded(
            4001 * 5003,
            config=config,
            checkpoint=checkpoint,
            budget=allowance(),
        )

        self.assertEqual(result.result.reconstruct(), 4001 * 5003)
        self.assertEqual(result.checkpoint["payload"]["version"], 5)

    def test_polynomial_switch_refusal_preserves_verified_store(self):
        job = SIQSJob(4001 * 5003, config=configuration(), budget=allowance())
        job.run(max_blocks=1)
        collector = job.engine.collector
        before = (
            collector.polynomial,
            dict(collector._atoms),
            dict(collector._pending),
        )
        step = job.family.next()
        collector.budget = Budget(
            work_limit=job.budget.used,
            used=job.budget.used,
            seconds=None,
            cpu_seconds=None,
        )
        with self.assertRaises(BudgetExhaustedError):
            collector.set_polynomial(step.polynomial, step.roots)

        self.assertEqual(
            (collector.polynomial, collector._atoms, collector._pending),
            before,
        )
        collector.budget = allowance()
        collector.set_polynomial(step.polynomial, step.roots)

        self.assertEqual(collector._atoms, before[1])
        self.assertEqual(collector._pending, before[2])

    def test_invalid_controls_and_checkpoint_caps(self):
        for options in (
            dict(mode="unknown"),
            dict(factor_count=9),
            dict(growth_steps=5),
            dict(pool_size=2),
            dict(checkpoint_bytes=4095),
        ):
            with self.assertRaises(ValueError):
                configuration(**options)

        with self.assertRaises(MemoryError):
            configuration(memory_bytes=4 * 1024 * 1024)
        job = SIQSJob(
            4001 * 5003,
            config=configuration(checkpoint_bytes=4096),
            budget=allowance(),
        )
        job.run(max_blocks=1)
        with self.assertRaises(MemoryError):
            job.checkpoint()

        self.assertEqual(job.run().reason, "factor_found")
        fresh = SIQSJob(
            4001 * 5003, config=configuration(), budget=allowance()
        )
        checkpoint = fresh.checkpoint()
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                dict(checkpoint, blob="x" * (4 * 1024 * 1024 + 1)),
                budget=allowance(),
            )

    def test_dispatch_after_ecm_and_recursive_reconstruction(self):
        config = PortfolioConfig(
            trial_bound=2,
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=((5, 7, 1),),
            memory_bytes=64 * 1024 * 1024,
            siqs=configuration(),
        )

        for n in (4001 * 5003, -4001 * 5003, (4001 * 5003) ** 2):
            result = factorize_bounded(
                n, config=config, seed=31, budget=allowance()
            )

            self.assertTrue(result.result.complete)
            self.assertEqual(result.result.reconstruct(), n)
            stages = [event["stage"] for event in result.events]
            if "siqs" in stages:
                self.assertLess(stages.index("ecm"), stages.index("siqs"))

        self.assertIsNone(PortfolioConfig().siqs)

    def test_dispatch_checkpoint_restores_siqs_under_same_budget(self):
        config = PortfolioConfig(
            trial_bound=2,
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=(),
            memory_bytes=64 * 1024 * 1024,
            siqs=configuration(),
        )
        n = 4001 * 5003

        first = factorize_bounded(n, config=config, budget=allowance(150000))

        self.assertFalse(first.result.complete)
        self.assertEqual(first.result.reconstruct(), n)
        self.assertEqual(
            first.checkpoint["payload"]["state"]["current"]["stage"], "siqs"
        )

        resumed = factorize_bounded(
            n, config=config, checkpoint=first.checkpoint, budget=allowance()
        )

        self.assertTrue(resumed.result.complete)
        self.assertEqual(resumed.result.reconstruct(), n)
        self.assertGreater(resumed.work_used, first.work_used)
