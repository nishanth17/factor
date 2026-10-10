"""Restrictive allowance, bounded reuse, and refund/resume regressions."""

import threading
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs.parallel import CollectionPool, ParallelConfig, ParallelSIQSJob


def allowance(work=200000000):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


class FollowupTests(unittest.TestCase):
    def test_resieve_reserves_support_from_actual_window_and_keeps_caps(self):
        from v2.qs import build_factor_base, qs_polynomial
        from v2.qs.sieve_collector import SieveCollector, SieveConfig
        from v2.tests.test_qs_sieve import signature

        base = build_factor_base(
            1000003 * 1000033, bound=10000, budget=allowance()
        ).factor_base
        polynomial = qs_polynomial(base)
        cfg = SieveConfig(
            block_width=4096,
            memory_bytes=64 * 2**20,
            division="resieve",
            score_policy="powers",
        )
        worker = SieveCollector(
            polynomial, base, config=cfg, budget=allowance()
        )
        dense = 4096 * (512 + 128 * len(base.entries))

        self.assertGreater(worker._workspace + dense, cfg.memory_bytes)

        expected = SieveCollector(
            polynomial,
            base,
            config=replace(cfg, division="bucket"),
            budget=allowance(),
        ).collect(-50, 51)

        actual = worker.collect(-50, 51)

        self.assertEqual(signature(actual), signature(expected))
        self.assertLessEqual(actual.workspace_bytes, cfg.memory_bytes)
        self.assertFalse(worker._resieved)
        self.assertEqual(worker._scratch_bytes, 0)
        worker.config = replace(cfg, memory_bytes=worker._workspace + 1)

        refused = worker.collect(100, 151)

        self.assertEqual(refused.reason, "memory_limit")
        self.assertEqual(refused.stats["scanned"], 0)
        worker.config = cfg

        recovered = worker.collect(100, 151)

        expected = SieveCollector(
            polynomial,
            base,
            config=replace(cfg, division="bucket"),
            budget=allowance(),
        ).collect(100, 151)

        self.assertEqual(
            {k: v for k, v in signature(recovered).items() if 100 <= k < 151},
            signature(expected),
        )

    def test_disabled_stage_bounds_do_not_reserve_unused_tables(self):
        cfg = PortfolioConfig(
            pm1_attempts=0, pm1_b2=10**12, ecm_tiers=((10**6, 10**12, 0),)
        )

        self.assertEqual(
            cfg.max_hi, 1 + max(cfg.trial_bound, cfg.max_input_bits)
        )

        run = factorize_bounded(4001 * 5003, config=cfg, budget=allowance())

        self.assertTrue(run.result.complete)
        self.assertEqual(run.result.reconstruct(), 4001 * 5003)

        restored = factorize_bounded(
            4001 * 5003,
            config=cfg,
            checkpoint=run.checkpoint,
            budget=allowance(),
        )

        self.assertEqual(restored.result, run.result)
        with self.assertRaises(MemoryError):
            replace(cfg, ecm_tiers=((10**6, 10**12, 1),))

    def test_local_pools_need_no_multiprocessing_synchronization(self):
        with patch(
            "v2.qs.parallel.multiprocessing.get_context",
            side_effect=AssertionError,
        ):
            for mode, workers in (("serial", 1), ("thread", 2)):
                with CollectionPool(mode, workers) as pool:
                    result = ParallelSIQSJob(
                        4001 * 5003, budget=allowance()
                    ).run(pool=pool)

                    self.assertEqual(
                        result.divisor * result.cofactor, 4001 * 5003
                    )

    def test_local_workers_reuse_a_checked_base_without_revalidating(self):
        from v2.qs.factor_base import FactorBase

        for mode, workers in (("serial", 1), ("thread", 2)):
            job = ParallelSIQSJob(4001 * 5003, budget=allowance())
            job._setup()
            with patch.object(
                FactorBase,
                "__post_init__",
                side_effect=AssertionError("reconstructed trusted base"),
            ):
                with CollectionPool(mode, workers) as pool:
                    result = job.run(pool=pool)

            self.assertEqual(result.divisor * result.cofactor, job.n)

    def test_lease_ceiling_is_not_a_minimum_allowance(self):
        cfg = ParallelConfig(family_count=4, pool_size=8)

        for mode, workers, work in (
            ("serial", 1, 1000000),
            ("thread", 4, 15000000),
        ):
            job = ParallelSIQSJob(
                4001 * 5003, config=cfg, budget=allowance(work)
            )
            with CollectionPool(mode, workers) as pool:
                result = job.run(pool=pool)

            self.assertEqual(result.reason, "factor_found")
            self.assertEqual(result.divisor * result.cofactor, job.n)
            self.assertLessEqual(job.budget.used, work)

    def test_parent_waits_for_refunds_before_declaring_work_exhaustion(self):
        cfg = ParallelConfig(family_count=2, pool_size=8, assignment_work=2)
        job = ParallelSIQSJob(4001 * 5003, config=cfg, budget=allowance())
        released = threading.Event()
        original_memory = job._memory
        original_merge = job._merge

        def memory(workers):
            original_memory(workers)
            job.budget.work_limit = job.budget.used + 3

        def collect(task, state):
            index = task[2]
            if index == 1 and (not released.wait(3)):
                raise AssertionError("parent did not retry pending admission")
            return dict(
                id=index,
                atoms=(),
                divisor=None,
                reason="complete",
                work=1 if index == 0 else 0,
                scanned=0,
                workspace_bytes=0,
                cpu_seconds=0,
                seconds=0,
            )

        def merge(result):
            if job.next_assignment == 0:
                try:
                    job.budget.consume(2)
                except BudgetExhaustedError:
                    released.set()
                    raise

            original_merge(result)

        with CollectionPool("thread", 2) as pool:
            with (
                patch.object(job, "_memory", side_effect=memory),
                patch.object(job, "_merge", side_effect=merge),
                patch("v2.qs.parallel._collect", side_effect=collect),
            ):
                result = job.run(pool=pool, fixed_work=True, max_assignments=2)

            self.assertFalse(pool.lock.locked())

        self.assertTrue(released.is_set())
        self.assertEqual(result.reason, "families_exhausted")
        self.assertEqual(job.next_assignment, 2)
        self.assertFalse(job.pending)
        self.assertEqual(job.budget.used, job.budget.work_limit)

    def test_pending_solver_also_waits_for_running_lease_refunds(self):
        cfg = ParallelConfig(family_count=2, pool_size=8, assignment_work=2)
        job = ParallelSIQSJob(4001 * 5003, config=cfg, budget=allowance())
        released = threading.Event()
        original_memory, original_merge = (job._memory, job._merge)
        failures, finished = ([0], [False])

        def solve(final):
            if finished[0]:
                return None
            job.engine.solver = object()

            try:
                job.budget.consume(2)
            except BudgetExhaustedError:
                failures[0] += 1
                if failures[0] >= 2:
                    released.set()
                raise

            finished[0] = True
            job.engine.solver = None
            return None

        def memory(workers):
            original_memory(workers)
            job.budget.work_limit = job.budget.used + 3
            job.engine._solve = solve

        def collect(task, state):
            index = task[2]
            if index == 1 and (not released.wait(3)):
                raise AssertionError(
                    "parent did not resume the pending solver"
                )
            return dict(
                id=index,
                atoms=(),
                divisor=None,
                reason="complete",
                work=1 if index == 0 else 0,
                scanned=0,
                workspace_bytes=0,
                cpu_seconds=0,
                seconds=0,
            )

        def merge(result):
            if job.next_assignment == 0:
                solve(False)
            original_merge(result)

        with CollectionPool("thread", 2) as pool:
            with (
                patch.object(job, "_memory", side_effect=memory),
                patch.object(job, "_merge", side_effect=merge),
                patch("v2.qs.parallel._collect", side_effect=collect),
            ):
                result = job.run(pool=pool, fixed_work=True, max_assignments=2)

        self.assertGreaterEqual(failures[0], 2)
        self.assertEqual(result.reason, "families_exhausted")
        self.assertEqual(job.next_assignment, 2)
        self.assertEqual(job.budget.used, job.budget.work_limit)

    def test_base_digest_is_cached_per_immutable_instance(self):
        import hashlib
        import json

        from v2.qs.factor_base import build_factor_base
        from v2.qs.families import _identity

        base = build_factor_base(
            1000003 * 1000033, bound=1000, budget=allowance()
        ).factor_base
        expected = hashlib.sha256(
            json.dumps(
                [
                    base.n,
                    base.multiplier,
                    base.bound,
                    [[e.prime, e.square_roots] for e in base.entries],
                ],
                separators=(",", ":"),
            ).encode()
        ).hexdigest()

        self.assertEqual(_identity(base), expected)
        with patch("v2.qs.families.json.dumps", side_effect=AssertionError):
            self.assertEqual(_identity(base), expected)
        changed = replace(base, bound=1001)

        self.assertIsNone(changed._family_identity)
        self.assertNotEqual(_identity(changed), expected)

    def test_worker_counts_its_shared_immutable_base_once(self):
        from v2.qs.families import PolynomialFamily
        from v2.qs.parallel import _collect, _FixedExportCollector
        from v2.qs.sieve_collector import SieveConfig

        cfg = ParallelConfig(
            base_bound=10000,
            half_width=512,
            factor_count=3,
            family_count=4,
            pool_size=16,
            batch_width=1,
            max_batch_atoms=512,
            collector=SieveConfig(
                score_policy="powers", division="bucket", residual_bound=10000
            ),
        )
        job = ParallelSIQSJob(
            1000003 * 1000033, config=cfg, budget=allowance()
        )
        job._setup()
        base = job.base
        family = PolynomialFamily(
            base,
            job._assignment_primes(0),
            budget=allowance(),
            memory_bytes=32 * 2**20,
        )
        step = family._make_step(0, None)
        collector = _FixedExportCollector(
            step.polynomial,
            base,
            budget=allowance(),
            config=replace(
                cfg.collector, memory_bytes=32 * 2**20, max_atoms=512
            ),
            precomputed_roots=step.roots,
        )
        cap = (
            family.workspace_bytes
            + collector._workspace
            - base.workspace_bytes
            + 65536
        )

        self.assertLess(cap, family.workspace_bytes + collector._workspace)
        cfg = replace(cfg, worker_memory_bytes=cap)
        task = (
            (base.n, base.multiplier, base.bound, base.entries),
            job._assignment_primes(0),
            0,
            cfg,
            50000000,
            None,
            True,
        )
        with CollectionPool() as pool:
            result = _collect(task, pool.state)

        self.assertEqual(result["reason"], "complete")
        self.assertLessEqual(result["workspace_bytes"], cap)


if __name__ == "__main__":
    unittest.main()
