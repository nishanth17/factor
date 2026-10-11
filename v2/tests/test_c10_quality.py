"""Quality ordering preserves the flyer population and charged resume."""

import random
import unittest
from dataclasses import replace
from math import prod

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.qs import (
    GrowingDoubleLargeSieveConfig,
    SieveConfig,
    SIQSConfig,
    SIQSJob,
    build_factor_base,
)
from v2.qs.assignment_stream import (
    QUALITY_BATCH_SIZE,
    AssignmentStream,
    quality_workspace,
)


class QualityAssignmentTests(unittest.TestCase):
    def setUp(self):
        self.n = 1000003 * 1000033
        self.base = build_factor_base(self.n, bound=500).factor_base

    def stream(self, policy="quality", count=1000, **options):
        arguments = dict(
            factor_count=3,
            family_count=count,
            pool_size=8,
            seed=29,
            policy=policy,
            budget=Budget(work_limit=10**9),
            memory_bytes=32 * 2**20,
        )
        arguments.update(options)
        return AssignmentStream(self.base, 128, **arguments)

    def test_complete_batches_are_exact_permutations(self):
        control, quality = self.stream("flyer"), self.stream()
        self.assertEqual(len(quality), len(control))
        all_values = [quality[i] for i in range(len(quality))]
        self.assertEqual(len(set(all_values)), len(all_values))
        for start in range(0, len(quality), QUALITY_BATCH_SIZE):
            stop = min(start + QUALITY_BATCH_SIZE, len(quality))
            originals = [control[i] for i in range(start, stop)]
            expected = sorted(
                originals,
                key=lambda values: (
                    abs(prod(values) - quality.target),
                    values,
                ),
            )
            self.assertEqual(all_values[start:stop], expected)
            for values in expected:
                self.assertEqual(len(values), len(set(values)))
                self.assertTrue(all(self.n % p and p != 2 for p in values))

    def test_extension_and_random_access_preserve_prefix(self):
        short, long = self.stream(count=17), self.stream()
        self.assertEqual(short.identity, long.identity)
        indices = list(range(len(short)))
        random.Random(91).shuffle(indices)
        for i in indices:
            self.assertEqual(short[i], long[i])
        with self.assertRaises(IndexError):
            short[len(short)]
        self.assertNotEqual(short.identity, self.stream("flyer").identity)

    def test_one_factor_and_last_partial_batch(self):
        control = self.stream("flyer", factor_count=1, pool_size=7)
        quality = self.stream(factor_count=1, pool_size=7)
        expected = sorted(
            (control[i] for i in range(len(control))),
            key=lambda values: (abs(prod(values) - quality.target), values),
        )
        self.assertEqual([quality[i] for i in range(len(quality))], expected)

    def test_interrupted_batch_is_not_published(self):
        stream = self.stream()
        stream.budget = Budget(work_limit=1)
        with self.assertRaises(BudgetExhaustedError):
            stream[0]
        self.assertEqual(stream._quality_start, -1)
        self.assertEqual(stream._quality_values, ())
        stream.budget = Budget(work_limit=10**9)
        self.assertEqual(stream[0], self.stream()[0])

    def test_reservations_and_cache_rebuild_are_charged(self):
        config = SIQSConfig(assignment_policy="flyer")
        quality = replace(config, assignment_policy="quality")
        self.assertEqual(
            quality.metadata_reserve - config.metadata_reserve,
            quality_workspace(config.factor_count),
        )
        with self.assertRaises(MemoryError):
            self.stream(memory_bytes=self.base.workspace_bytes + 65536)

        first, restored = self.stream(), self.stream()
        before = first.budget.used
        value = first[2]
        first_cost = first.budget.used - before
        before = first.budget.used
        self.assertEqual(first[2], value)
        self.assertEqual(first.budget.used - before, 1)
        before = restored.budget.used
        self.assertEqual(restored[2], value)
        self.assertEqual(restored.budget.used - before, first_cost)

    def test_partial_new_batch_keeps_previous_cache(self):
        stream = self.stream(pool_size=10)
        previous = stream[0]
        saved = stream._quality_values
        preparation = QUALITY_BATCH_SIZE * (
            stream.core_count + stream.target.bit_length() + 5
        )
        one_raw = len(stream.pool) * (stream.core_count + 1)
        stream.budget = Budget(work_limit=preparation + one_raw)
        with self.assertRaises(BudgetExhaustedError):
            stream[QUALITY_BATCH_SIZE]
        self.assertEqual(stream._quality_start, 0)
        self.assertIs(stream._quality_values, saved)

        stream.budget = Budget(work_limit=10**9)
        self.assertEqual(stream[0], previous)
        self.assertEqual(
            stream[QUALITY_BATCH_SIZE],
            self.stream(pool_size=10)[QUALITY_BATCH_SIZE],
        )

    def test_growing_checkpoint_preserves_absent_quotas_and_extension(self):
        collector = GrowingDoubleLargeSieveConfig(
            block_width=32,
            residual_bound=10000,
            large_prime_bound=10000,
            large_product_bound=1000000,
        )
        config = SIQSConfig(
            base_bound=500,
            half_width=128,
            factor_count=3,
            pool_size=8,
            family_count=17,
            assignment_policy="quality",
            max_stalled=100,
            collector=collector,
        )
        self.assertIsNone(config.memory_bytes)
        self.assertIsNone(config.checkpoint_bytes)
        job = SIQSJob(self.n, config=config, budget=Budget(work_limit=10**9))
        self.assertEqual(job.run(max_blocks=1).reason, "paused")
        checkpoint = job.checkpoint()
        restored = SIQSJob.from_checkpoint(
            checkpoint,
            budget=Budget(work_limit=10**9),
            config=replace(config, family_count=28),
            allow_extension=True,
        )
        self.assertGreater(restored.budget.used, job.budget.used)
        self.assertIsNone(restored.config.memory_bytes)
        self.assertIsNone(restored.config.checkpoint_bytes)
        self.assertIsNone(restored.config.collector.max_atoms)
        self.assertEqual(
            restored.assignments.identity, job.assignments.identity
        )
        for i in range(17):
            self.assertEqual(restored.assignments[i], job.assignments[i])
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                checkpoint,
                budget=Budget(work_limit=10**9),
                config=replace(config, assignment_policy="flyer"),
                allow_extension=True,
            )

    def test_job_checkpoint_reconstructs_quality_policy(self):
        config = SIQSConfig(
            base_bound=500,
            half_width=128,
            factor_count=3,
            pool_size=8,
            family_count=20,
            assignment_policy="quality",
            max_stalled=100,
            memory_bytes=64 * 2**20,
            collector=SieveConfig(block_width=32, residual_bound=10000),
        )
        job = SIQSJob(self.n, config=config, budget=Budget(work_limit=10**9))
        result = job.run(max_blocks=1)
        self.assertEqual(result.reason, "paused")
        checkpoint = job.checkpoint()
        spent = job.budget.used
        restored = SIQSJob.from_checkpoint(
            checkpoint, budget=Budget(work_limit=10**9)
        )

        self.assertEqual(restored.config.assignment_policy, "quality")
        self.assertGreater(restored.budget.used, spent)
        self.assertEqual(
            restored.assignments.identity, job.assignments.identity
        )
        self.assertEqual(restored.assignments[0], job.assignments[0])
        result = restored.run()
        self.assertEqual((result.divisor or 1) * result.cofactor, self.n)
        if result.divisor is not None:
            self.assertTrue(1 < result.divisor < self.n)
            self.assertEqual(self.n % result.divisor, 0)


if __name__ == "__main__":
    unittest.main()
