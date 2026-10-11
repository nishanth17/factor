"""Allocator-grown DLP retains exact kernels and charged checkpoint state."""

import copy
import json
import unittest
from dataclasses import asdict, replace
from types import SimpleNamespace
from unittest.mock import patch

from v2.execution.budget import BudgetExhaustedError
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import (
    DoubleLargeSieveConfig,
    GrowingDoubleLargeSieveConfig,
    Polynomial,
    SieveCollector,
    SIQSConfig,
    SIQSJob,
    build_factor_base,
    qs_polynomial,
    verify_combined,
)
from v2.qs.large_prime_checkpoint import restore_store
from v2.qs.large_primes import LargePrimeForest, graph_reserve
from v2.qs.storage import RelationPublicationError
from v2.tests.test_qs_dlp import admit, allowance, encoded
from v2.tests.test_qs_pipeline import dense_kernel, span
from v2.tests.test_siqs import mutate


def growing_worker(**options):
    base = build_factor_base(91, bound=3).factor_base
    values = dict(
        block_width=1,
        residual_bound=1000,
        large_prime_bound=1000,
        large_product_bound=1000000,
    )
    values.update(options)
    return SieveCollector(
        qs_polynomial(base),
        base,
        config=GrowingDoubleLargeSieveConfig(**values),
        budget=allowance(),
    )


class GrowingStorageTests(unittest.TestCase):
    def test_absent_quotas_and_legacy_defaults_are_distinct(self):
        growing = growing_worker().config
        for name in (
            "memory_bytes",
            "max_atoms",
            "max_partials",
            "max_relations",
            "split_call_limit",
        ):
            self.assertIsNone(getattr(growing, name))
        self.assertEqual(asdict(growing)["storage_policy"], "grow")
        config = SIQSConfig(collector=growing)
        self.assertIsNone(config.memory_bytes)
        self.assertIsNone(config.checkpoint_bytes)
        self.assertEqual(
            config.metadata_reserve, 65536 + 512 * config.polynomial_limit
        )
        capped = DoubleLargeSieveConfig(
            large_prime_bound=1000, large_product_bound=1000000
        )
        self.assertEqual(capped.max_partials, 256)
        self.assertEqual(capped.max_atoms, 2048)
        self.assertNotIn("storage_policy", asdict(capped))
        self.assertEqual(
            SIQSConfig(collector=capped).memory_bytes, 32 * 1024**2
        )
        self.assertEqual(
            SIQSConfig(collector=capped).checkpoint_bytes, 1024**2
        )
        self.assertEqual(replace(growing, max_atoms=100000).max_atoms, 100000)

    def test_growing_storage_refuses_unsupported_worker_and_sss_adapters(self):
        from v2.qs.parallel import ParallelConfig
        from v2.qs.sss import SSSConfig

        for constructor in (ParallelConfig, SSSConfig):
            with self.assertRaisesRegex(ValueError, "serial SIQS"):
                constructor(collector=growing_worker().config)

    def test_growing_graph_preserves_complete_incidence_kernel(self):
        edges = [(1, 7), (3, 5), (5, 11), (3, 11), (13, 13), (3, 5)]
        forest, cycles = LargePrimeForest(edge_limit=None), []
        for index, (left, right) in enumerate(edges):
            plan = forest.plan(left, right, str(index), None, allowance())
            self.assertFalse(plan.evicted)
            if plan.cycle is not None:
                cycles.append(
                    sum(1 << int(i) for i in plan.cycle) | (1 << index)
                )
            forest.commit(plan)
        self.assertEqual(
            span(cycles), dense_kernel([(1 << a) ^ (1 << b) for a, b in edges])
        )
        self.assertGreater(
            graph_reserve(65537, edge_limit=None), 256 * 1024**2
        )
        with self.assertRaises(MemoryError):
            graph_reserve(65537)

    def test_growth_keeps_partials_and_explicit_fifo_keeps_owned_paths(self):
        for cap in (None, 3):
            collector = growing_worker(max_partials=cap)
            star, _, _ = admit(collector, -7, scalar=41)
            for position, pair in ((-5, (3, 11)), (1, (3, 5)), (-4, (5, 11))):
                admit(collector, position, pair)
            owned = set(collector._combined[0].atom_ids)
            for position, residual in ((5, 67), (8, 233), (23, 499)):
                _, refusal, _ = admit(collector, position, scalar=residual)
                self.assertIsNone(refusal)
            self.assertEqual(star.relation_id in collector._atoms, cap is None)
            self.assertTrue(owned <= set(collector._atoms))
            verify_combined(
                collector._combined[0], collector.factor_base, collector._atoms
            )
            payload = encoded(collector)
            for _ in range(2):
                restored = growing_worker(max_partials=cap)
                before = restored.budget.used
                restore_store(payload, restored, restored.budget)
                self.assertGreater(restored.budget.used, before)
                self.assertEqual(encoded(restored), payload)
                self.assertIsNone(restored._graph.edge_limit)

    def test_explicit_refusals_and_interrupted_publication_preserve_store(
        self,
    ):
        for field, limit, reason in (
            ("max_atoms", 1, "atom_limit"),
            ("max_relations", 0, "relation_limit"),
            ("memory_bytes", None, "memory_limit"),
        ):
            collector = growing_worker()
            admit(collector, 1, (3, 5))
            before = encoded(collector)
            value = collector._workspace if field == "memory_bytes" else limit
            collector.config = replace(collector.config, **{field: value})
            _, refusal, _ = admit(
                collector, 2, (3, 5), polynomial=Polynomial(91, 1, 9, 1, 3)
            )
            self.assertEqual(refusal, reason)
            self.assertEqual(encoded(collector), before)
        collector = growing_worker()
        admit(collector, 1, (3, 5))
        before = encoded(collector)
        with patch(
            "v2.qs.large_prime_collector.combine_relations",
            side_effect=BudgetExhaustedError,
        ):
            with self.assertRaises(BudgetExhaustedError):
                admit(
                    collector, 2, (3, 5), polynomial=Polynomial(91, 1, 9, 1, 3)
                )
        self.assertEqual(encoded(collector), before)
        broken = copy.deepcopy(before)
        broken["atoms"][0][-1] = [3, 9]
        restored = growing_worker()
        with self.assertRaises(ValueError):
            restore_store(broken, restored, restored.budget)
        self.assertFalse(restored._atoms)

    def test_allocator_failure_cannot_publish_a_resumable_store(self):
        collector = growing_worker()
        admit(collector, 1, (3, 5))
        snapshot = encoded(collector)
        commit = collector._graph.commit

        def interrupted(plan):
            commit(plan)
            raise MemoryError("injected allocator failure after graph commit")

        with patch.object(collector._graph, "commit", side_effect=interrupted):
            with self.assertRaises(RelationPublicationError):
                admit(collector, -4, (5, 11))
        with self.assertRaises(RelationPublicationError):
            encoded(collector)
        with self.assertRaises(RelationPublicationError):
            collector.collect(0, 1)
        with self.assertRaises(RelationPublicationError):
            admit(collector, -5, (3, 11))

        restored = growing_worker()
        restore_store(snapshot, restored, restored.budget)
        admit(restored, -4, (5, 11))
        admit(restored, -5, (3, 11))
        verify_combined(
            restored._combined[0], restored.factor_base, restored._atoms
        )

    def test_outer_portfolio_resolves_child_cap_and_resume_keeps_growth(self):
        siqs = SIQSConfig(
            base_bound=200,
            half_width=256,
            collector=GrowingDoubleLargeSieveConfig(
                residual_bound=40000,
                large_prime_bound=20000,
                large_product_bound=5120000,
            ),
        )
        config = PortfolioConfig(
            trial_bound=2,
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=(),
            memory_bytes=64 * 1024**2,
            siqs=siqs,
        )
        self.assertEqual(
            config.siqs.memory_bytes,
            config.memory_bytes - config.workspace_reserve - 8192,
        )
        self.assertIsNone(config.siqs.collector.max_partials)
        n = 4001 * 5003
        first = factorize_bounded(n, config=config, budget=allowance(150000))
        self.assertEqual(first.result.reconstruct(), n)
        resumed = factorize_bounded(
            n, config=config, checkpoint=first.checkpoint, budget=allowance()
        )
        self.assertTrue(resumed.result.complete)
        self.assertEqual(resumed.result.reconstruct(), n)
        self.assertGreater(resumed.work_used, first.work_used)

    def test_uncapped_checkpoint_counts_live_solver_and_preparation(self):
        from v2.qs.checkpoint import _reserve_growing_checkpoint

        collector = growing_worker()
        admit(collector, 1, (3, 5))
        config = SIQSConfig(
            mode="qs",
            base_bound=3,
            half_width=8,
            factor_count=1,
            memory_bytes=4 * 2**20,
            collector=collector.config,
        )
        engine = SimpleNamespace(
            config=replace(collector.config, memory_bytes=4 * 2**20),
            collector=collector,
            prepared=None,
            solver=None,
        )
        job = SimpleNamespace(config=config, engine=engine)
        _reserve_growing_checkpoint(job)

        # Model each independently retained reservation at a solver pause.
        # Refusal must precede any snapshot allocation and leave atoms intact.
        before = encoded(collector)
        for stage in ("prepared", "solver"):
            if stage == "prepared":
                engine.prepared = SimpleNamespace(
                    workspace_bytes=4 * 2**20, shared_workspace_bytes=0
                )
            else:
                engine.solver = SimpleNamespace(
                    matrix=SimpleNamespace(workspace_bytes=4 * 2**20)
                )
            with self.assertRaises(MemoryError):
                _reserve_growing_checkpoint(job)
            self.assertEqual(encoded(collector), before)
            setattr(engine, stage, None)

        # Completed solves release these containers; a historical high-water
        # mark must not become an artificial checkpoint refusal.
        engine._workspace_peak = 4 * 2**20
        _reserve_growing_checkpoint(job)

    def test_dynamic_restore_reserves_input_and_preserves_public_policy(self):
        config = SIQSConfig(
            mode="qs",
            base_bound=3,
            half_width=8,
            factor_count=1,
            memory_bytes=4 * 2**20,
            collector=growing_worker().config,
        )
        job = SIQSJob(91, config=config, budget=allowance())
        job.run(max_blocks=1)
        child_limit = job.engine.config.memory_limit
        snapshot = job.checkpoint()
        for _ in range(2):
            restored = SIQSJob.from_checkpoint(
                snapshot, config=config, budget=allowance()
            )
            self.assertEqual(restored.config, config)
            self.assertIsNone(restored.config.checkpoint_bytes)
            self.assertEqual(restored.engine.config.memory_limit, child_limit)
            self.assertIs(
                restored.engine.collector.config, restored.engine.config
            )
            self.assertEqual(
                restored.stats["capacity"]["metadata_reserve_bytes"],
                config.metadata_reserve,
            )
            snapshot = restored.checkpoint()

        oversized = dict(snapshot, blob=snapshot["blob"] + " " * 2**20)
        with patch("v2.qs.checkpoint.json.loads", side_effect=AssertionError):
            with self.assertRaises(MemoryError):
                SIQSJob.from_checkpoint(
                    oversized, config=config, budget=allowance()
                )

    def test_version_five_full_job_restore_is_charged_and_checked(self):
        config = SIQSConfig(
            mode="qs",
            base_bound=3,
            half_width=8,
            factor_count=1,
            collector=growing_worker().config,
        )
        job = SIQSJob(91, config=config, budget=allowance())
        job.run(max_blocks=1)
        job.budget.seconds = (
            0  # Emergency snapshots remain possible at a stop.
        )
        checkpoint = job.checkpoint()
        self.assertEqual(checkpoint["version"], 5)
        self.assertIsNone(
            json.loads(checkpoint["blob"])["config"]["memory_bytes"]
        )
        restored = SIQSJob.from_checkpoint(checkpoint, budget=allowance())
        self.assertEqual(restored.config, config)
        self.assertGreater(
            restored.budget.used, checkpoint["resources"]["work_used"]
        )
        again = SIQSJob.from_checkpoint(
            restored.checkpoint(), budget=allowance()
        )
        self.assertEqual(again.config, config)
        broken = mutate(
            checkpoint,
            lambda p: p["config"]["collector"].update(storage_policy="other"),
        )
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(broken, budget=allowance())


if __name__ == "__main__":
    unittest.main()
