"""Stable mixed-row selections, exact cache binding and charged recovery."""

import json
import tempfile
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.benchmarks.p38_r3_experiments import history_filter, live_compaction
from v2.budget import BudgetExhaustedError
from v2.qs import linear_algebra
from v2.qs.checkpoint import _restore_store, _store
from v2.qs.extraction import (
    DependencyExtractor,
    _VerificationCache,
    extract_dependency,
    prepare_relations,
)
from v2.qs.linear_algebra import DependencySolver, filter_matrix
from v2.tests import test_qs_pipeline as pipeline_tests
from v2.tests.test_qs import unlimited_budget
from v2.tests.test_qs_pipeline import dense_kernel, span
from v2.tests.test_qs_sieve import collector


class StableRowTests(unittest.TestCase):
    def setUp(self):
        self.fixture = pipeline_tests.ExtractionTests()
        self.fixture.setUp()

    def test_combined_indices_survive_later_full_admission_and_restore(self):
        f = self.fixture
        run = collector(f.polynomial, f.base)
        prefix = ()
        prior_full = 0
        combined_seen = grew_after_combined = False

        for position in range(-128, 129):
            result = run.collect(position, position + 1)

            self.assertEqual(run.matrix_relations[: len(prefix)], prefix)
            if combined_seen and len(run._full) > prior_full:
                grew_after_combined = True
            combined_seen |= bool(run._combined)
            prior_full = len(run._full)
            prefix = run.matrix_relations
            del result

        self.assertTrue(combined_seen and grew_after_combined)
        payload = json.loads(json.dumps(_store(run)))
        restored = collector(f.polynomial, f.base)
        _restore_store(payload, restored, unlimited_budget())

        self.assertEqual(restored.matrix_relations, run.matrix_relations)
        legacy = dict(payload)
        del legacy["row_order"]
        restored = collector(f.polynomial, f.base)
        _restore_store(legacy, restored, unlimited_budget())

        self.assertEqual(
            restored.matrix_relations, tuple(run._full + run._combined)
        )

    def test_bad_row_orders_and_tampered_payloads_are_reverified(self):
        f = self.fixture
        run = collector(f.polynomial, f.base)
        run.collect(-128, 129)
        payload = json.loads(json.dumps(_store(run)))

        for order in (
            [0] * len(payload["row_order"]),
            payload["row_order"][:-1],
            [len(payload["row_order"])] + payload["row_order"][1:],
        ):
            with self.assertRaises(ValueError):
                _restore_store(
                    dict(payload, row_order=order),
                    collector(f.polynomial, f.base),
                    unlimited_budget(),
                )

        payload["atoms"][0][2] *= -1
        with self.assertRaises(ValueError):
            _restore_store(
                payload, collector(f.polynomial, f.base), unlimited_budget()
            )

    def test_store_reset_clears_prefix_and_all_cache_attestations(self):
        f = self.fixture
        run = collector(f.polynomial, f.base)
        run.collect(-128, 129)
        run._preparation_cache = _VerificationCache(2**20)
        run._preparation_cache.tested_dependencies.add(frozenset({1}))
        run.set_polynomial(run.polynomial, run._roots, retain_relations=False)

        self.assertEqual(run.matrix_relations, ())
        self.assertFalse(run._preparation_cache.tested_dependencies)

    def test_frozen_legacy_pending_solver_replays(
        self,
    ):
        from v2.benchmarks.p38_r3 import load_control, modules
        from v2.qs.siqs import SIQSJob

        with tempfile.TemporaryDirectory(prefix="r3-legacy-") as directory:
            control = modules(load_control(directory))
            old = control["qs.siqs"]
            budget = control["budget"].Budget(
                work_limit=10**9, seconds=30, cpu_seconds=30
            )
            job = old.SIQSJob(
                4001 * 5003,
                config=old.SIQSConfig(
                    base_bound=200,
                    half_width=256,
                    collector=old.SieveConfig(residual_bound=500),
                ),
                budget=budget,
            )
            budget.reason = "work_limit"
            with patch.object(
                control["qs.pipeline"].DependencySolver,
                "run",
                side_effect=control["budget"].BudgetExhaustedError(
                    "work_limit"
                ),
            ):
                self.assertEqual(job.run().reason, "work_limit")

            checkpoint = job.checkpoint()

            self.assertEqual(checkpoint["version"], 1)

            restored = SIQSJob.from_checkpoint(
                checkpoint, budget=unlimited_budget()
            )

            result = restored.run()

            self.assertEqual({result.divisor, result.cofactor}, {4001, 5003})

    def test_new_version_requires_the_mixed_order_field(self):
        from v2.qs.siqs import SIQSJob
        from v2.tests.test_siqs import configuration, mutate

        job = SIQSJob(
            4001 * 5003, config=configuration(), budget=unlimited_budget()
        )

        self.assertEqual(job.run(max_blocks=1).reason, "paused")
        checkpoint = job.checkpoint()

        self.assertEqual(checkpoint["version"], 3)
        damaged = mutate(
            checkpoint, lambda payload: payload["store"].pop("row_order")
        )
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(damaged, budget=unlimited_budget())

    def test_sss_and_parallel_new_versions_require_mixed_order(self):
        from v2.qs.parallel import ParallelSIQSJob
        from v2.qs.sss import SSSConfig, SSSJob
        from v2.tests.test_qs_parallel import configuration, mutate
        from v2.tests.test_sss_dispatch import reseal

        job = SSSJob(
            4001 * 4003,
            config=SSSConfig(base_bound=400),
            budget=unlimited_budget(),
        )
        job._setup()
        checkpoint = job.checkpoint()

        self.assertEqual(checkpoint["version"], 3)
        payload = json.loads(checkpoint["blob"])
        del payload["store"]["row_order"]
        with self.assertRaises(ValueError):
            SSSJob.from_checkpoint(
                reseal(checkpoint, payload), budget=unlimited_budget()
            )
        job = ParallelSIQSJob(
            4001 * 5003, config=configuration(), budget=unlimited_budget()
        )
        job._setup()
        checkpoint = job.checkpoint()

        self.assertEqual(json.loads(checkpoint["blob"])["version"], 4)
        damaged = mutate(
            checkpoint, lambda payload: payload["store"].pop("row_order")
        )
        with self.assertRaises(ValueError):
            ParallelSIQSJob.from_checkpoint(damaged, budget=unlimited_budget())

    def test_opt_in_cadence_keeps_resume_and_final_extraction(self):
        from v2.qs.siqs import SIQSJob
        from v2.tests.test_siqs import configuration

        config = configuration(filter_row_growth=32, tested_dependencies=True)
        job = SIQSJob(4001 * 5003, config=config, budget=unlimited_budget())
        job.run(max_blocks=1)
        checkpoint = job.checkpoint()

        restored = SIQSJob.from_checkpoint(
            checkpoint, budget=unlimited_budget()
        )

        self.assertEqual(restored.engine.filter_row_growth, 32)
        self.assertTrue(restored.engine.tested_dependencies)

        result = restored.run()

        self.assertEqual({result.divisor, result.cofactor}, {4001, 5003})
        for growth in (0, True, 4097):
            with self.assertRaises((TypeError, ValueError)):
                configuration(filter_row_growth=growth)
        with self.assertRaises(TypeError):
            configuration(tested_dependencies=1)


class TestedDependencyTests(unittest.TestCase):
    def setUp(self):
        self.fixture = pipeline_tests.ExtractionTests()
        self.fixture.setUp()
        f = self.fixture
        self.prepared = prepare_relations(
            f.relations, f.base, f.store, budget=unlimited_budget()
        )

        dependencies = DependencySolver(
            filter_matrix(self.prepared.rows), budget=unlimited_budget()
        ).run()
        self.mask = next(
            mask
            for mask in dependencies
            if extract_dependency(
                self.prepared, mask, budget=unlimited_budget()
            ).divisor
            is None
        )

    def test_cached_selection_survives_index_shift_and_retains_payload(self):
        cache = _VerificationCache(2**20)
        first = DependencyExtractor(
            self.prepared,
            (self.mask,),
            budget=unlimited_budget(),
            tested_cache=cache,
        )

        self.assertIsNone(first.run())
        self.assertEqual(len(first.trials), 1)
        f = self.fixture
        # Move the last row to the front: same selection, different bit mask.
        rotated = self.prepared.relations[-1:] + self.prepared.relations[:-1]
        prepared = prepare_relations(
            rotated, f.base, f.store, budget=unlimited_budget()
        )
        mask = ((self.mask << 1) & ((1 << len(rotated)) - 1)) | (
            self.mask >> (len(rotated) - 1)
        )
        second = DependencyExtractor(
            prepared, (mask,), budget=unlimited_budget(), tested_cache=cache
        )

        self.assertIsNone(second.run())
        self.assertEqual(second.cache_skips, 1)
        self.assertFalse(second.trials)
        self.assertLessEqual(cache.used, cache.memory_bytes)
        key = cache.dependency_key(prepared, mask, unlimited_budget())
        base, relation, atoms = next(iter(key))
        changed = replace(relation, sign=-relation.sign)
        altered = frozenset(
            (base, changed if row == relation else row, source)
            for base, row, source in key
        )

        self.assertNotEqual(altered, key)
        self.assertNotIn(altered, cache.tested_dependencies)
        if atoms:
            altered_atom = replace(atoms[0], sign=-atoms[0].sign)

            self.assertEqual(altered_atom.relation_id, atoms[0].relation_id)
            self.assertNotEqual(
                (base, relation, (altered_atom,) + atoms[1:]),
                (base, relation, atoms),
            )

    def test_zero_cap_and_refused_trial_publish_no_dependency(self):
        cache = _VerificationCache(0)
        extractor = DependencyExtractor(
            self.prepared,
            (self.mask,),
            budget=unlimited_budget(),
            tested_cache=cache,
        )
        extractor.run()

        self.assertFalse(cache.tested_dependencies)
        cache = _VerificationCache(2**20)
        budget = unlimited_budget()
        budget.work_limit = 0
        extractor = DependencyExtractor(
            self.prepared, (self.mask,), budget=budget, tested_cache=cache
        )
        with self.assertRaises(BudgetExhaustedError):
            extractor.run()

        self.assertEqual(extractor.next_dependency, 0)
        self.assertFalse(cache.tested_dependencies)
        extractor.budget = unlimited_budget()

        self.assertIsNone(extractor.run())
        self.assertEqual(extractor.next_dependency, 1)

    def test_stale_prepared_identity_cannot_hide_changed_row_payload(self):
        cache = _VerificationCache(2**20)
        extractor = DependencyExtractor(
            self.prepared,
            (self.mask,),
            budget=unlimited_budget(),
            tested_cache=cache,
        )
        extractor.run()
        index = (self.mask & -self.mask).bit_length() - 1
        rows = list(self.prepared.relations)
        rows[index] = replace(rows[index], sign=-rows[index].sign)
        damaged = replace(self.prepared, relations=tuple(rows))
        with self.assertRaises(ValueError):
            DependencyExtractor(
                damaged,
                (self.mask,),
                budget=unlimited_budget(),
                tested_cache=cache,
            ).run()

    def test_equal_parity_is_not_an_equal_selection(self):
        f = self.fixture
        altered = dict(f.store)
        combined = next(
            row for row in self.prepared.relations if hasattr(row, "atom_ids")
        )
        atom = altered[combined.atom_ids[0]]
        altered[atom.relation_id] = replace(
            atom, exponents=atom.exponents + ((99991, 2),)
        )
        with self.assertRaises(ValueError):
            prepare_relations(
                (combined,), f.base, altered, budget=unlimited_budget()
            )


class PivotCounterTests(unittest.TestCase):
    def test_incremental_counts_and_refusal_match_actual_pivots(self):
        rows = (7, 11, 13, 14, 3, 5, 6, 0)
        solver = DependencySolver(
            filter_matrix(rows), budget=unlimited_budget()
        )
        while solver.next_row < len(solver.matrix.rows):
            solver.step()
            actual = sum(row.bit_count() for row, _ in solver.pivots.values())

            self.assertEqual(solver.pivot_nonzeros, actual)
            self.assertEqual(solver.peak_nonzeros, actual)

        self.assertEqual(span(solver.run()), dense_kernel(rows))
        refused = DependencySolver(
            filter_matrix(rows), budget=unlimited_budget()
        )
        refused.budget.work_limit = 0
        with self.assertRaises(BudgetExhaustedError):
            refused.step()

        self.assertEqual(refused.pivot_nonzeros, 0)
        refused.budget = unlimited_budget()

        self.assertEqual(span(refused.run()), dense_kernel(rows))


class MatrixChallengerTests(unittest.TestCase):
    def test_histories_disjoint_batches_and_live_maps_preserve_exact_kernel(
        self,
    ):
        import random

        generator = random.Random(38303)
        fixtures = [
            (0, 3, 5, 6, 1, 9),
            (3, 3, 5, 5),
            ((1 << 200) | (1 << 300),) * 3,
        ] + [
            tuple(generator.randrange(64) for _ in range(8)) for _ in range(20)
        ]

        for rows in fixtures:
            expected = dense_kernel(rows)

            for size, use_history in (
                (1, True),
                (32, True),
                (1, False),
                (32, False),
            ):
                matrix = history_filter(
                    rows,
                    linear_algebra,
                    unlimited_budget(),
                    8 * 2**20,
                    batch_size=size,
                    use_history=use_history,
                )
                for pivot in ("highest", "lowest"):
                    solver = DependencySolver(
                        matrix, pivot=pivot, budget=unlimited_budget()
                    )

                    self.assertEqual(span(solver.run()), expected)

            matrix = filter_matrix(rows, weight_two=True)
            compact = live_compaction(matrix, unlimited_budget(), 8 * 2**20)

            self.assertEqual(span(DependencySolver(compact).run()), expected)
            inverse = compact.stats["inverse_columns"]

            for remapped, mask in zip(compact.rows, compact.masks):
                original = 0
                for index, row in enumerate(rows):
                    if mask & (1 << index):
                        original ^= row

                self.assertEqual(
                    sum(
                        1 << column
                        for index, column in enumerate(inverse)
                        if remapped & (1 << index)
                    ),
                    original,
                )

    def test_challenger_storage_and_history_bounds_refuse_before_work(self):
        budget = unlimited_budget()
        with self.assertRaises(MemoryError):
            history_filter((3, 5, 6), linear_algebra, budget, 1)
        self.assertEqual(budget.used, 0)
        with self.assertRaises(ValueError):
            history_filter((1 << 4096,), linear_algebra, budget, 2**30)
        self.assertEqual(budget.used, 0)


if __name__ == "__main__":
    unittest.main()
