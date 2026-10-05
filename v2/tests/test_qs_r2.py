"""Independent scoring bounds, capped plan coverage and refusal contracts."""

import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.budget import BudgetExhaustedError
from v2.qs import Polynomial, build_factor_base
from v2.qs.polynomial import polynomial_roots
from v2.qs.score import SCORE_SCALE, log_bounds
from v2.qs.sieve_collector import SieveCollector, SieveConfig
from v2.tests.test_qs import reference_positions, unlimited_budget
from v2.tests.test_qs_sieve import signature


class R2CollectorTests(unittest.TestCase):
    def test_legacy_positional_configuration_remains_unchanged(self):
        config = SieveConfig(
            17,
            64,
            "list",
            "sparse",
            "bucket",
            0,
            0,
            "powers",
            97,
            4096,
            4096,
            4096,
            64 * 2**20,
        )
        self.assertEqual(config.residual_bound, 97)
        self.assertEqual(config.max_atoms, 4096)
        self.assertEqual(config.memory_bytes, 64 * 2**20)
        self.assertEqual(config.power_plan_bytes, 0)

    def test_plan_bound_reserves_large_coefficients_before_evaluation(self):
        base = build_factor_base(10403, bound=100).factor_base
        polynomial = Polynomial(10403, 1, 1, 2**200)
        budget = unlimited_budget()
        worker = SieveCollector(
            polynomial,
            base,
            config=SieveConfig(
                score_policy="fixed",
                power_plan_bytes=2**20,
                memory_bytes=64 * 2**20,
            ),
            budget=budget,
        )
        budget.work_limit = budget.used + 100
        with patch.object(Polynomial, "value", side_effect=AssertionError):
            with self.assertRaises(BudgetExhaustedError):
                worker._set_plan_interval(0, 1)
        self.assertIsNone(worker._plan_window)

    def test_experimental_tiny_batch_and_chunks_preserve_exact_coverage(self):
        from v2.benchmarks.p38_r2_experiments import experiment_collector
        from v2.qs import sieve_collector

        base = build_factor_base(10403, bound=100).factor_base
        polynomial = Polynomial(10403, 1, 49, 8)
        expected = reference_positions(polynomial, base, -71, 80, 97)
        work = {}
        for variant in ("current", "tiny", "batch", "chunks"):
            budget = unlimited_budget()
            cls = experiment_collector(sieve_collector, variant)
            worker = cls(
                polynomial,
                base,
                config=SieveConfig(
                    score_policy="powers",
                    division="bucket",
                    residual_bound=97,
                    block_width=17,
                    memory_bytes=64 * 2**20,
                    max_atoms=4096,
                    max_partials=4096,
                    max_relations=4096,
                ),
                budget=budget,
            )
            result = worker.collect(-71, 80)
            self.assertEqual(result.reason, "complete")
            self.assertEqual(signature(result), expected)
            if variant == "tiny":
                self.assertGreater(result.stats["tiny_evaluations"], 0)
            if variant == "batch":
                self.assertGreater(result.stats["tree_calls"], 0)
            work[variant] = budget.used
        self.assertEqual(work["current"], work["chunks"])

    def test_fixed_score_bounds_against_independent_integer_powers(self):
        values = set(range(1, 4097))
        for bits in (20, 53, 333, 4096):
            values.update((2**bits - 1, 2**bits, 2**bits + 1))
            values.update((1537 << bits, (1537 << bits) + 123))
        for value in values:
            lower, upper = log_bounds(value)
            powered = value**SCORE_SCALE
            self.assertLessEqual(1 << lower, powered)
            self.assertLessEqual(powered, 1 << upper)
            self.assertLessEqual(upper - lower, 2)
        with self.assertRaises(ValueError):
            log_bounds(0)

    def test_cached_singular_and_ramified_marks_cover_high_valuations(self):
        from v2.qs.factor_base import FactorBase, FactorBaseEntry

        base = FactorBase(
            3**12,
            1,
            11,
            (
                FactorBaseEntry(2, (1,)),
                FactorBaseEntry(3, (0,)),
                FactorBaseEntry(5, (1, 4)),
                FactorBaseEntry(7, (1, 6)),
            ),
        )
        polynomial = Polynomial(3**12, 1, 1, 0)
        expected = reference_positions(polynomial, base, -512, 513, 97)
        for backend in ("list", "array", "bytearray"):
            worker = SieveCollector(
                polynomial,
                base,
                config=SieveConfig(
                    score_policy="fixed",
                    power_plan_bytes=2**20,
                    score_backend=backend,
                    division="bucket",
                    block_width=17,
                    residual_bound=97,
                    max_atoms=4096,
                    max_relations=4096,
                    max_partials=4096,
                    memory_bytes=64 * 2**20,
                ),
                budget=unlimited_budget(),
            )
            result = worker.collect(-512, 513)
            self.assertEqual(result.reason, "complete")
            self.assertEqual(signature(result), expected)
            self.assertGreater(result.stats["plan_replays"], 0)

    def test_independent_full_division_coverage_all_backends_and_tails(self):
        base = build_factor_base(10403, bound=100).factor_base
        for a, b in ((1, 102), (7, 1), (49, 8)):
            polynomial = Polynomial(10403, 1, a, b)
            expected = reference_positions(polynomial, base, -71, 80, 97)
            for policy in ("powers", "fixed"):
                for cap in (0, 256, 2**20):
                    for backend in ("list", "array", "bytearray"):
                        for division in ("bucket", "resieve"):
                            for cutoff in (0, 5):
                                worker = SieveCollector(
                                    polynomial,
                                    base,
                                    config=SieveConfig(
                                        score_policy=policy,
                                        power_plan_bytes=cap,
                                        score_backend=backend,
                                        division=division,
                                        small_prime_cutoff=cutoff,
                                        block_width=17,
                                        residual_bound=97,
                                        max_atoms=4096,
                                        max_partials=4096,
                                        max_relations=4096,
                                        memory_bytes=64 * 2**20,
                                    ),
                                    budget=unlimited_budget(),
                                )
                                result = worker.collect(-71, 80)
                                self.assertEqual(result.reason, "complete")
                                self.assertEqual(signature(result), expected)
                                self.assertLessEqual(
                                    result.stats["plan_bytes"], cap
                                )
                                if cap == 256:
                                    self.assertGreater(
                                        result.stats["plan_refusals"], 0
                                    )

    def test_plan_window_switch_and_polynomial_switch_clear_cached_marks(self):
        base = build_factor_base(10403, bound=100).factor_base
        worker = SieveCollector(
            Polynomial(10403, 1, 1, 102),
            base,
            config=SieveConfig(
                score_policy="fixed",
                power_plan_bytes=2**20,
                block_width=17,
                residual_bound=97,
                memory_bytes=64 * 2**20,
            ),
            budget=unlimited_budget(),
        )
        first = worker.collect(-20, 20)
        self.assertGreater(first.stats["plan_replays"], 0)
        worker.collect(47, 81)
        self.assertEqual(worker._plan_window[:2], (47, 81))
        polynomial = Polynomial(10403, 1, 7, 1)
        roots = tuple(
            polynomial_roots(
                polynomial, base, entry, budget=unlimited_budget()
            )
            for entry in base.entries
        )
        worker.set_polynomial(polynomial, roots)
        self.assertEqual(worker._power_plans, {})
        self.assertIsNone(worker._plan_window)

    def test_work_and_cancellation_preserve_uncommitted_position(self):
        base = build_factor_base(10403, bound=100).factor_base
        polynomial = Polynomial(10403, 1, 7, 1)
        for policy in ("powers", "fixed"):
            budget = unlimited_budget()
            worker = SieveCollector(
                polynomial,
                base,
                config=SieveConfig(
                    score_policy=policy,
                    power_plan_bytes=2**20,
                    residual_bound=97,
                    max_atoms=4096,
                    max_relations=4096,
                    max_partials=4096,
                    memory_bytes=64 * 2**20,
                ),
                budget=budget,
            )
            budget.work_limit = budget.used
            paused = worker.collect(-71, 80)
            self.assertEqual(paused.reason, "work_limit")
            self.assertEqual(paused.next_position, -71)
            self.assertEqual(paused.atoms, ())
            self.assertLessEqual(budget.used, budget.work_limit)
            budget.work_limit = 10**12
            budget.cancelled = lambda: True
            cancelled = worker.collect(paused.next_position, 80)
            self.assertEqual(cancelled.reason, "cancelled")
            self.assertEqual(cancelled.next_position, -71)
            budget.cancelled = None
            complete = worker.collect(cancelled.next_position, 80)
            self.assertEqual(
                signature(complete),
                reference_positions(polynomial, base, -71, 80, 97),
            )

    def test_array_saturation_and_plan_storage_reservation(self):
        base = build_factor_base(10403, bound=100).factor_base
        config = SieveConfig(
            score_backend="array",
            score_policy="fixed",
            memory_bytes=64 * 2**20,
        )
        worker = SieveCollector(
            Polynomial(10403, 1, 1, 102),
            base,
            config=config,
            budget=unlimited_budget(),
        )
        worker._add((0,), 2**40, {})
        self.assertEqual(worker._scores[0], 2**32 - 1)
        self.assertEqual(worker._clip_threshold(2**40), 2**32 - 1)
        with self.assertRaises(MemoryError):
            SieveCollector(
                worker.polynomial,
                base,
                config=replace(
                    config,
                    memory_bytes=worker._workspace,
                    power_plan_bytes=2**20,
                ),
                budget=unlimited_budget(),
            )
        with self.assertRaises(ValueError):
            replace(config, power_plan_bytes=16 * 2**20 + 1)

    def test_checkpoint_rebuilds_disposable_plans_and_preserves_work(self):
        from v2.qs.siqs import SIQSConfig, SIQSJob

        n = 1000003 * 1000033
        config = SIQSConfig(
            base_bound=1000,
            half_width=512,
            batch_width=1,
            memory_bytes=64 * 2**20,
            collector=SieveConfig(
                score_policy="fixed",
                power_plan_bytes=2**20,
                division="bucket",
                residual_bound=1000000,
                max_atoms=8192,
                max_relations=4096,
                max_partials=2048,
            ),
        )
        job = SIQSJob(n, config=config, budget=unlimited_budget())
        paused = job.run(max_blocks=1)
        self.assertIsNone(paused.divisor)
        self.assertTrue(job.engine.collector._power_plans)
        saved = job.checkpoint()
        restored = SIQSJob.from_checkpoint(saved, budget=unlimited_budget())
        self.assertEqual(restored.engine.collector._power_plans, {})
        self.assertEqual(restored.config.collector.power_plan_bytes, 2**20)
        self.assertEqual(
            restored.engine.next_position, job.engine.next_position
        )
        used = restored.budget.used
        result = restored.run()
        self.assertGreater(restored.budget.used, used)
        self.assertIsNotNone(result.divisor)
        self.assertEqual(result.divisor * result.cofactor, n)


if __name__ == "__main__":
    unittest.main()
