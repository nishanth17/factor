"""Independent arithmetic and boundary controls for the performance audit."""

import ast
import json
import random
import unittest
from collections import Counter
from dataclasses import replace
from math import isqrt, prod
from unittest.mock import patch

from v2.budget import Budget, BudgetExhaustedError
from v2.preprocessing import (
    _POWER_MODULI,
    integer_root,
    power_residue_possible,
)
from v2.qs.linear_algebra import DependencySolver, filter_matrix
from v2.qs.parallel import CollectionPool, ParallelConfig, ParallelSIQSJob
from v2.schedules import SieveContext
from v2.work_budget import PollingBudget


def allowance(work=10**12):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


class PollingTests(unittest.TestCase):
    def test_exact_work_refusal_and_explicit_publication_check(self):
        original = allowance(3)
        budget = PollingBudget(original)
        budget.consume(2)
        with self.assertRaises(BudgetExhaustedError):
            budget.consume(2)
        self.assertEqual(original.used, 2)
        original.work_limit = 10
        budget.consume()
        original.cancelled = lambda: True
        with self.assertRaises(BudgetExhaustedError):
            budget.consume(0)
        self.assertEqual(original.used, 3)
        self.assertEqual(original.reason, "cancelled")

    def test_external_poll_interval_has_a_finite_action_bound(self):
        calls = 0

        def cancelled():
            nonlocal calls
            calls += 1
            return calls == 2

        original = allowance()
        original.cancelled = cancelled
        budget = PollingBudget(original)
        for _ in range(64):
            budget.consume(7)
        with self.assertRaises(BudgetExhaustedError):
            budget.consume(7)
        self.assertEqual(original.used, 64 * 7)
        self.assertEqual(calls, 2)
        for value in (True, -1, 1.5):
            with self.assertRaises((TypeError, ValueError)):
                budget.consume(value)
        with self.assertRaises(ValueError):
            PollingBudget(original, interval=65)


class PreprocessingRepairTests(unittest.TestCase):
    def test_small_prime_screens_never_reject_actual_powers(self):
        generator = random.Random(36102026)
        for exponent, modulus in _POWER_MODULI.items():
            self.assertTrue(
                all(modulus % d for d in range(2, isqrt(modulus) + 1))
            )
            self.assertEqual((modulus - 1) % exponent, 0)
            for base in range(modulus):
                self.assertTrue(
                    power_residue_possible(base**exponent, exponent)
                )
            for _ in range(20):
                base = generator.getrandbits(128)
                value = base**exponent
                self.assertTrue(power_residue_possible(value, exponent))
                for neighbor in (value - 1, value + 1):
                    if not power_residue_possible(neighbor, exponent):
                        root = integer_root(neighbor, exponent)
                        self.assertNotEqual(root**exponent, neighbor)

    def test_cursor_segments_match_independent_trial_division(self):
        context = SieveContext(300, segment_size=16)
        for lo in range(-3, 250):
            for width in (0, 1, 2, 15, 32):
                hi = lo + width
                expected = [
                    n
                    for n in range(max(lo, 2), hi)
                    if all(n % p for p in range(2, isqrt(n) + 1))
                ]
                self.assertEqual(context.prime_segment(lo, hi), expected)
        iterator = context.primes(2, 100)
        self.assertEqual(next(iterator), 2)
        with self.assertRaises(RuntimeError):
            context.prime_segment(2, 10)
        iterator.close()
        self.assertEqual(context.prime_segment(2, 10), [2, 3, 5, 7])

    def test_cursor_segments_preserve_wide_endpoints(self):
        left = 2**32 - 16
        context = SieveContext(left + 64, segment_size=16)
        self.assertEqual(
            context.prime_segment(left, left + 32),
            list(context.primes(left, left + 32)),
        )


class SparseMatrixRepairTests(unittest.TestCase):
    def test_wide_labels_fit_without_changing_lifted_kernels(self):
        high = 1 << 99999
        rows = high, high, 1, 1
        for weight_two in (False, True):
            matrix = filter_matrix(
                rows,
                weight_two=weight_two,
                budget=allowance(),
                memory_bytes=2**20,
            )
            self.assertIs(matrix.original_rows, rows)
            self.assertEqual(matrix.stats["working_columns"], 2)
            for pivot in ("highest", "lowest"):
                masks = DependencySolver(
                    matrix, pivot=pivot, budget=allowance()
                ).run()
                span = {0}
                for mask in masks:
                    span |= {value ^ mask for value in tuple(span)}
                self.assertEqual(span, {0, 3, 12, 15})
        budget = allowance()
        with self.assertRaises(MemoryError):
            filter_matrix(rows, budget=budget, memory_bytes=1)
        self.assertEqual(budget.used, 0)


class PreparationReuseTests(unittest.TestCase):
    def test_cache_requires_identical_base_rows_and_provenance(self):
        from v2.qs.extraction import _prepare_relations, _VerificationCache
        from v2.tests.test_qs_pipeline import ExtractionTests

        fixture = ExtractionTests()
        fixture.setUp()
        cache = _VerificationCache(2**20)

        def prepare(store=None, base=None):
            budget = allowance()
            result = _prepare_relations(
                fixture.relations,
                fixture.base if base is None else base,
                fixture.store if store is None else store,
                budget=budget,
                verification_cache=cache,
            )
            return result, budget.used

        first, first_work = prepare()
        second, second_work = prepare()
        self.assertEqual(first, second)
        self.assertEqual(cache.hits, len(first.relations))
        self.assertLess(second_work, first_work)
        hits = cache.hits
        self.assertEqual(prepare(base=replace(fixture.base))[0], first)
        self.assertEqual(cache.hits, hits)

        combined = fixture.run.combined_relations[0]
        identity = combined.atom_ids[0]
        original = fixture.store[identity]
        altered = dict(fixture.store)
        altered[identity] = replace(original, sign=-original.sign)
        with self.assertRaises(ValueError):
            prepare(store=altered)
        del altered[identity]
        with self.assertRaises(ValueError):
            prepare(store=altered)
        self.assertLessEqual(cache.used, cache.memory_bytes)
        cache.clear()
        self.assertFalse(cache.records)
        self.assertEqual(cache.used, 0)

    def test_failed_checks_and_small_caps_do_not_publish_attestations(
        self,
    ):
        from v2.qs.extraction import _VerificationCache
        from v2.tests.test_qs_pipeline import ExtractionTests

        fixture = ExtractionTests()
        fixture.setUp()
        relation = fixture.run.full_relations[0]
        for cap in (0, 4096):
            cache = _VerificationCache(cap)
            with self.assertRaises(ValueError):
                cache.check(
                    replace(relation, sign=-relation.sign),
                    fixture.base,
                    fixture.store,
                    allowance(),
                    2**20,
                )
            with self.assertRaises(BudgetExhaustedError):
                cache.check(
                    relation, fixture.base, fixture.store, allowance(0), 2**20
                )
            self.assertFalse(cache.records)
            cache.check(
                relation, fixture.base, fixture.store, allowance(), 2**20
            )
            self.assertLessEqual(cache.used, cap)


class CollisionRepairTests(unittest.TestCase):
    def test_frozen_constructor_control_keeps_its_superclass_context(self):
        from v2.benchmarks.performance_audit import apply_revert
        from v2.qs.sss import SSSCollector
        from v2.tests.test_sss import search_collector

        original = SSSCollector.__init__
        try:
            apply_revert("sssf_trees")
            collector = search_collector(mode="sssf", filter_divisor=2)
            self.assertEqual(
                collector.batch.radical, prod(collector.factor_base.primes)
            )
            self.assertTrue(collector.assignment(0))
        finally:
            SSSCollector.__init__ = original

    def test_sssf_two_stage_detection_preserves_admissible_candidates(self):
        from v2.tests.test_sss import search_collector, trial_residual

        for divisor in (1, 2, 10, 1024):
            collector = search_collector(
                mode="sssf", filter_divisor=divisor, filter_bound=1000
            )
            primes = collector.factor_base.primes
            cut = max(1, len(primes) // divisor)
            for index in range(8):
                candidates = collector.assignment(index)
                expected = []
                for position, value in candidates:
                    if not value:
                        expected.append(position)
                        continue
                    first = trial_residual(value, primes[:cut])
                    residual = trial_residual(first, primes[cut:])
                    if first < 1000 and residual <= (
                        collector.config.residual_bound
                    ):
                        expected.append(position)
                stats = dict(
                    filter_rejections=0,
                    generated_candidates=0,
                    tree_rejections=0,
                )
                collector._prepare(index, stats)
                self.assertEqual(collector._assignment, tuple(expected))

    def test_assignment_order_and_work_match_the_immutable_control(self):
        from v2.benchmarks.performance_audit import BASELINE, checked_baseline
        from v2.qs import sss
        from v2.tests.test_sss import search_collector

        source = checked_baseline(BASELINE)["source"]["v2/qs/sss.py"]
        node = next(
            n
            for n in ast.parse(source).body
            if isinstance(n, ast.FunctionDef)
            and n.name == "collision_candidates"
        )
        namespace = dict(vars(sss), Counter=Counter)
        exec(
            compile(
                ast.Module(body=[node], type_ignores=[]), str(BASELINE), "exec"
            ),
            namespace,
        )
        control = namespace["collision_candidates"]
        collector = search_collector()
        for index in range(32):
            collector.budget = allowance()
            current = collector.assignment(index)
            work = collector.budget.used
            collector.budget = allowance()
            with patch.object(sss, "collision_candidates", control):
                expected = collector.assignment(index)
            self.assertEqual(current, expected)
            self.assertEqual(work, collector.budget.used)
        roots = list(collector._roots)
        entry = roots[collector.small_count]
        roots[collector.small_count] = replace(
            entry, roots=entry.roots[:1] * 2
        )
        arguments = (
            collector.polynomial,
            tuple(roots),
            collector.small_count,
            (0, 1, 2),
            collector.coefficients,
        )
        self.assertEqual(
            sss.collision_candidates(*arguments, budget=allowance()),
            control(*arguments, budget=allowance()),
        )


class ChunkedWorkerTests(unittest.TestCase):
    def test_resume_settings_preserve_reserved_geometry(self):
        job = ParallelSIQSJob(
            4001 * 5003, config=self.configuration(), budget=allowance()
        )
        job.run(max_assignments=1, fixed_work=True)
        original = job.config
        for changes in (
            dict(batch_width=64),
            dict(max_batch_atoms=64),
            dict(checkpoint_bytes=original.checkpoint_bytes + 4096),
        ):
            job.config = replace(original, **changes)
            for action in (job.checkpoint, job.run):
                with self.assertRaises(ValueError):
                    action()
        job.config = replace(
            original, assignment_work=20_000_000, poll_interval=1
        )
        self.assertEqual(
            job.run(max_assignments=1, fixed_work=True).reason, "paused"
        )

    def test_sparse_norm_bound_admits_a_large_base_worker_reservation(self):
        config = replace(
            ParallelConfig(),
            base_bound=10000,
            half_width=512,
            factor_count=3,
            batch_width=256,
            family_count=4,
            pool_size=16,
        )
        job = ParallelSIQSJob(
            1000000007 * 1000000009, config=config, budget=allowance()
        )
        job._setup()
        job._memory(4)
        self.assertLess(job.stats["peak_owned_bytes"], config.memory_bytes)
        dense_atom = 4096 + 256 * len(job.base.entries)
        dense_atom += 128 * (
            job.n.bit_length()
            + 2 * config.factor_count * config.base_bound.bit_length()
            + 32
        )
        dense_owned = config.parent_memory_bytes + 4 * (
            config.worker_memory_bytes
            + 3 * 256 * dense_atom
            + job.base.workspace_bytes
        )
        self.assertGreater(dense_owned, config.memory_bytes)
        result = job.run(max_assignments=1, fixed_work=True)
        self.assertEqual(result.reason, "paused")
        from v2.qs.relations import verify_atomic

        for atom in job.engine.collector._atoms.values():
            verify_atomic(
                atom, job.base, residual_bound=config.collector.residual_bound
            )

    def test_paused_solver_releases_pool_budget_adapter(self):
        from v2.qs.linear_algebra import DependencySolver

        job = ParallelSIQSJob(
            4001 * 5003,
            config=replace(ParallelConfig(), family_count=1),
            budget=allowance(),
        )
        original_run = DependencySolver.run

        def pause(solver):
            job.budget.work_limit = job.budget.used
            return original_run(solver)

        with patch.object(DependencySolver, "run", pause):
            result = job.run(fixed_work=True)
        self.assertEqual(result.reason, "work_limit")
        self.assertIsNotNone(job.engine.solver)
        self.assertIs(job.engine.solver.budget, job.budget)
        job.budget.work_limit = 10**12
        self.assertEqual(job.run(fixed_work=True).reason, "factor_found")

    def configuration(self, **options):
        return replace(
            ParallelConfig(),
            family_count=4,
            pool_size=8,
            batch_width=32,
            max_batch_atoms=32,
            **options,
        )

    def test_chunking_removes_whole_family_batch_refusals(self):
        config = self.configuration()
        coarse = ParallelSIQSJob(
            4001 * 5003,
            config=replace(config, batch_width=0),
            budget=allowance(),
        )
        self.assertEqual(coarse.run(fixed_work=True).reason, "batch_limit")
        expected = None
        for mode, workers in (("serial", 1), ("thread", 2), ("process", 2)):
            job = ParallelSIQSJob(coarse.n, config=config, budget=allowance())
            with CollectionPool(mode, workers) as pool:
                result = job.run(pool=pool, fixed_work=True)
            self.assertEqual(result.reason, "factor_found")
            self.assertEqual(result.divisor * result.cofactor, job.n)
            self.assertEqual(
                result.stats["scanned"], len(job.assignments) * 513
            )
            signature = job.engine.collector._atoms, job.budget.used
            if expected is None:
                expected = signature
            else:
                self.assertEqual(signature, expected)

    def test_cross_worker_checkpoint_retains_gray_and_block_progress(self):
        config = self.configuration(factor_count=3)
        job = ParallelSIQSJob(4001 * 5003, config=config, budget=allowance())
        with CollectionPool("thread", 2) as pool:
            self.assertEqual(
                job.run(pool=pool, max_assignments=3, fixed_work=True).reason,
                "paused",
            )
        restored = ParallelSIQSJob.from_checkpoint(
            json.loads(json.dumps(job.checkpoint())), budget=allowance()
        )
        with CollectionPool("process", 4) as pool:
            result = restored.run(pool=pool, fixed_work=True)
        uninterrupted = ParallelSIQSJob(
            job.n, config=config, budget=allowance()
        )
        uninterrupted.run(fixed_work=True)
        self.assertTrue(result.stats["schedule_complete"])
        self.assertEqual(
            restored.engine.collector._atoms,
            uninterrupted.engine.collector._atoms,
        )
        self.assertGreater(restored.budget.used, job.budget.used)

    def test_poll_work_and_rejection_of_future_checkpoint_rows(
        self,
    ):
        from v2.tests.test_qs_parallel import mutate

        config = self.configuration()
        jobs = [
            ParallelSIQSJob(
                4001 * 5003,
                config=replace(config, poll_interval=interval),
                budget=allowance(),
            )
            for interval in (1, 64)
        ]
        for job in jobs:
            job.run(max_assignments=2, fixed_work=True)
        self.assertEqual(jobs[0].budget.used, jobs[1].budget.used)
        self.assertEqual(
            jobs[0].engine.collector._atoms, jobs[1].engine.collector._atoms
        )
        checkpoint = jobs[0].checkpoint()
        for change in (
            lambda p: p.update(version=True),
            lambda p: p.update(next_assignment=0),
        ):
            with self.assertRaises(ValueError):
                ParallelSIQSJob.from_checkpoint(
                    mutate(checkpoint, change), budget=allowance()
                )
        original = jobs[0].budget
        original.work_limit = original.used
        result = jobs[0].run(fixed_work=True)
        self.assertEqual(result.reason, "work_limit")
        self.assertEqual(result.cofactor, jobs[0].n)
        original.work_limit = 10**12
        self.assertEqual(jobs[0].run(fixed_work=True).reason, "factor_found")


if __name__ == "__main__":
    unittest.main()
