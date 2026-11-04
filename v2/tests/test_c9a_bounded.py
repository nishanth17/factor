"""Bounded staging, exact promotion and independent restore regressions."""

import copy
import json
import tempfile
import unittest
from collections import Counter
from dataclasses import asdict, replace
from unittest.mock import patch

from v2 import qs
from v2.benchmarks.qs.c9a import extended
from v2.benchmarks.qs.c9a.configuration import configuration, decode
from v2.execution.budget import BudgetExhaustedError
from v2.qs.extended_checkpoint import pack_store, restore_store
from v2.qs.indexed_relations import FactorRecord
from v2.qs.residual_relations import ResidualAtom
from v2.qs.yield_memory import integer_reserve
from v2.tests.test_c9a_extended import atom_at, policy
from v2.tests.test_qs_dlp import allowance


def bounded(**changes):
    values = dict(asdict(policy()), extension_policy="extended-v2")
    values.update(changes)
    return qs.BoundedYieldSieveConfig(**values)


def worker(**changes):
    base = qs.build_factor_base(91, bound=3).factor_base
    return qs.SieveCollector(
        qs.qs_polynomial(base),
        base,
        config=bounded(**changes),
        budget=allowance(),
    )


class BoundedYieldTests(unittest.TestCase):
    def test_shared_slabs_and_contexts_remain_packed_after_restore(self):
        engine = worker(recompute_prime_limit=0)
        first = atom_at(engine, 1)
        second = atom_at(engine, -21)
        second = ResidualAtom(
            replace(second.polynomial),
            second.position,
            second.sign,
            tuple(second.exponents),
            primes=second.residual_primes,
        )
        self.assertIsNot(first.polynomial, second.polynomial)

        engine._admit(first, Counter())
        first_slab_bytes = engine._slab_bytes
        engine._admit(second, Counter())
        self.assertEqual(engine._slab_bytes, first_slab_bytes)
        self.assertEqual(engine._slab_bytes, engine._arena.allocated_bytes)
        self.assertEqual(len(engine._arena.polynomials), 1)
        saved = json.loads(json.dumps(pack_store(engine)))
        clone = worker(recompute_prime_limit=0)
        restore_store(saved, clone, clone.budget)

        self.assertEqual(pack_store(clone), saved)
        self.assertEqual(clone._workspace, engine._workspace)
        self.assertEqual(clone._slab_bytes, engine._slab_bytes)
        context = next(iter(clone._atoms.values())).polynomial
        for atom in clone._atoms.values():
            self.assertIs(atom.polynomial, context)
            self.assertIsInstance(atom.exponents, FactorRecord)
        self.assertIsInstance(clone._rows[0].exponents, FactorRecord)
        qs.verify_combined(
            clone._rows[0], clone.factor_base, clone._atoms, budget=allowance()
        )

    def test_slab_allocation_is_charged_after_unpublished_budget_refusal(self):
        engine = worker(recompute_prime_limit=0)
        atom = atom_at(engine, 1)
        original = engine._arena.plan
        workspace = engine._workspace

        def stop_after_packing(*args):
            original(*args)
            raise BudgetExhaustedError()

        with patch.object(
            engine._arena, "plan", side_effect=stop_after_packing
        ):
            with self.assertRaises(BudgetExhaustedError):
                engine._admit(atom, Counter())
        self.assertFalse(engine._atoms)
        self.assertGreater(engine._slab_bytes, 0)
        self.assertEqual(engine._workspace, workspace + engine._slab_bytes)

        slabs = engine._slab_bytes
        engine._admit(atom, Counter())
        self.assertEqual(len(engine._atoms), 1)
        self.assertEqual(engine._slab_bytes, slabs)

    def test_memory_refusal_precedes_any_slab_allocation(self):
        engine = worker()
        engine.config = replace(
            engine.config,
            memory_bytes=(engine._workspace + engine._scratch_bytes + 1),
        )
        workspace = engine._workspace

        self.assertEqual(
            engine._admit(atom_at(engine, 1), Counter()), "memory_limit"
        )
        self.assertFalse(engine._atoms)
        self.assertEqual(engine._arena.allocated_bytes, 0)
        self.assertEqual(engine._workspace, workspace)

    def test_conversion_rows_are_packed_and_restore_ownership_matches(self):
        engine = worker()
        engine._admit(atom_at(engine, 1, opaque=True), Counter())
        engine._admit(atom_at(engine, -21), Counter())
        self.assertTrue(engine._conversion_events)
        self.assertIsInstance(engine._rows[0].exponents, FactorRecord)
        saved = json.loads(json.dumps(pack_store(engine)))
        clone = worker()
        restore_store(saved, clone, clone.budget)
        self.assertEqual(pack_store(clone), saved)
        self.assertEqual(clone._workspace, engine._workspace)
        self.assertEqual(clone._slab_bytes, engine._slab_bytes)

    def test_integer_ownership_remains_finite_and_scales_with_digit_storage(
        self,
    ):
        small = integer_reserve(1)
        large = integer_reserve(1 << 100000)
        self.assertGreater(large, 100 * small)
        self.assertLess(large, 100000)

    def test_reset_releases_slab_and_context_ownership(self):
        engine = worker(recompute_prime_limit=0)
        initial = engine._workspace
        engine._admit(atom_at(engine, 1), Counter())
        engine._admit(atom_at(engine, -21), Counter())
        self.assertGreater(engine._slab_bytes, 0)

        engine.set_polynomial(
            engine.polynomial, engine._roots, retain_relations=False
        )
        self.assertFalse(engine._atoms)
        self.assertFalse(engine._arena.polynomials)
        self.assertEqual(engine._slab_bytes, 0)
        self.assertEqual(engine._workspace, initial)

    def test_context_interning_checks_coefficients_after_identity_match(self):
        engine = worker()
        engine._admit(atom_at(engine, 1), Counter())
        polynomial = engine.polynomial
        forged = qs.Polynomial(
            polynomial.n,
            polynomial.multiplier,
            polynomial.a,
            polynomial.b + polynomial.a,
            polynomial.square_coefficient,
        )
        object.__setattr__(forged, "identity", polynomial.identity)

        with self.assertRaises(ValueError):
            engine._arena.polynomial(forged)
        self.assertEqual(len(engine._arena.polynomials), 1)

    def test_cold_restore_preserves_existing_file_extents(self):
        with tempfile.TemporaryDirectory() as directory:
            engine = worker(spill_directory=directory, spill_bytes=2**20)
            engine._admit(atom_at(engine, 1, opaque=True), Counter())
            engine._admit(atom_at(engine, -21), Counter())
            saved = json.loads(json.dumps(pack_store(engine)))
            used = engine._arena.spill.used
            clone = worker(spill_directory=directory, spill_bytes=2**20)

            restore_store(saved, clone, clone.budget)
            self.assertEqual(pack_store(clone), saved)
            self.assertEqual(clone._arena.spill.used, used)
            qs.verify_combined(
                clone._rows[0],
                clone.factor_base,
                clone._atoms,
                budget=allowance(),
            )
            clone._arena.spill.close()
            engine._arena.spill.close()

    def test_fifo_staging_does_not_fill_incidence_and_pair_promotion_is_exact(
        self,
    ):
        engine, stats = worker(opaque_cache_entries=2), Counter()
        for x in (1, 2, 3):
            engine._admit(atom_at(engine, x, opaque=True), stats)
        self.assertFalse(engine._atoms)
        self.assertFalse(engine._graph.pivots)
        self.assertEqual(list(engine._opaque_staging.atoms), [53, 39])
        self.assertEqual(stats["opaque_staging_evictions"], 1)
        engine._admit(atom_at(engine, -22, opaque=True), stats)
        self.assertEqual(len(engine._rows), 1)
        self.assertEqual(len(engine._atoms), 2)
        self.assertEqual(list(engine._opaque_staging.atoms), [39])
        qs.verify_combined(
            engine._rows[0],
            engine.factor_base,
            engine._atoms,
            budget=allowance(),
        )

    def test_proven_factorization_promotes_a_staged_opaque(self):
        engine = worker()
        engine._admit(atom_at(engine, 1, opaque=True), Counter())
        engine._admit(atom_at(engine, -21), Counter())
        self.assertFalse(engine._opaque_staging.atoms)
        self.assertEqual(engine._conversion_events[0][1:], [15, [3, 5]])
        self.assertEqual(len(engine._rows), 1)
        qs.verify_combined(
            engine._rows[0],
            engine.factor_base,
            engine._atoms,
            budget=allowance(),
        )

    def test_byte_cap_and_zero_capacity_refuse_only_unclosed_history(self):
        for changes in (
            dict(opaque_cache_entries=0),
            dict(opaque_cache_bytes=1),
        ):
            engine, stats = worker(**changes), Counter()
            engine._admit(atom_at(engine, 1, opaque=True), stats)
            self.assertFalse(engine._opaque_staging.atoms)
            self.assertFalse(engine._atoms)
            self.assertEqual(stats["opaque_staging_rejections"], 1)

    def test_staging_roundtrip_and_arithmetic_corruption_rejection(self):
        engine = worker()
        for x in (1, 2):
            engine._admit(atom_at(engine, x, opaque=True), Counter())
        saved = json.loads(json.dumps(pack_store(engine)))
        clone = worker()
        restore_store(saved, clone, clone.budget)
        self.assertEqual(pack_store(clone), saved)
        clone._admit(atom_at(clone, -21), Counter())
        self.assertTrue(clone._rows)
        bad = copy.deepcopy(saved)
        bad["opaque_staging"][0][-1] = 21
        with self.assertRaises(ValueError):
            clone = worker()
            restore_store(bad, clone, clone.budget)

    def test_refusal_after_first_promotion_preserves_checked_substep(self):
        from v2.qs.extended_collector import admit_retained

        engine = worker()
        first = atom_at(engine, 1, opaque=True)
        incoming = atom_at(engine, -21, opaque=True)
        engine._admit(first, Counter())

        def refuse_second(owner, atom, stats):
            if atom.relation_id == incoming.relation_id:
                raise BudgetExhaustedError()
            return admit_retained(owner, atom, stats)

        with patch(
            "v2.qs.extended_collector.admit_retained",
            side_effect=refuse_second,
        ):
            with self.assertRaises(BudgetExhaustedError):
                engine._admit(incoming, Counter())
        self.assertIn(first.relation_id, engine._atoms)
        self.assertFalse(engine._opaque_staging.atoms)
        saved = pack_store(engine)
        # A staged duplicate of the retained opaque must fail independently
        # of the original record's identity, including compact token keys.
        bad = copy.deepcopy(saved)
        bad["opaque_staging"] = [
            [
                [
                    incoming.polynomial.a,
                    incoming.polynomial.b,
                    incoming.polynomial.square_coefficient,
                ],
                incoming.position,
                incoming.sign,
                saved["atoms"][0][3],
                [],
                incoming.opaque,
            ]
        ]
        with self.assertRaisesRegex(ValueError, "retained incidence"):
            clone = worker()
            restore_store(bad, clone, clone.budget)
        clone = worker()
        restore_store(saved, clone, clone.budget)
        clone._admit(incoming, Counter())
        self.assertEqual(len(clone._rows), 1)
        qs.verify_combined(
            clone._rows[0], clone.factor_base, clone._atoms, budget=allowance()
        )

    def test_cache_mutation_cannot_reuse_an_old_proof(self):
        engine = worker()
        first = atom_at(engine, 1, opaque=True)
        engine._admit(first, Counter())
        object.__setattr__(first, "opaque", 21)
        with self.assertRaises(ValueError):
            engine._admit(atom_at(engine, -21, opaque=True), Counter())

    def test_schema14_roundtrip_with_staging_and_coupled_policies(self):
        cfg = qs.SIQSConfig(
            base_bound=100,
            half_width=128,
            max_half_width=128,
            memory_bytes=128 * 2**20,
            multiplier=3,
            polynomial_policy="parity2",
            multiplier_policy="parity2_wide",
            assignment_policy="windowed_pair",
            root_update_policy="batch",
            collector=bounded(),
        )
        job = qs.SIQSJob(10007 * 10009, config=cfg, budget=allowance())
        job.run(max_blocks=1)
        saved = job.checkpoint()
        self.assertEqual(saved["version"], 14)
        clone = qs.SIQSJob.from_checkpoint(
            saved, config=cfg, budget=allowance()
        )
        self.assertEqual(clone.divisor, job.divisor)
        self.assertEqual(
            pack_store(clone.engine.collector),
            pack_store(job.engine.collector),
        )
        self.assertEqual(
            clone.run(max_blocks=1).divisor, job.run(max_blocks=1).divisor
        )
        with self.assertRaises(ValueError):
            qs.SIQSJob.from_checkpoint(
                saved,
                config=replace(
                    cfg,
                    collector=replace(cfg.collector, opaque_cache_entries=1),
                ),
                budget=allowance(),
            )

        pending = qs.SIQSJob(
            1000000007 * 1000000087, config=cfg, budget=allowance()
        )
        self.assertEqual(pending.run(max_blocks=1).reason, "paused")
        resumed = qs.SIQSJob.from_checkpoint(
            pending.checkpoint(), config=cfg, budget=allowance()
        )
        a, b = pending.run(max_blocks=1), resumed.run(max_blocks=1)
        self.assertEqual(
            (a.divisor, a.next_position), (b.divisor, b.next_position)
        )
        self.assertEqual(
            pack_store(pending.engine.collector),
            pack_store(resumed.engine.collector),
        )

    def test_bounded_profile_preserves_incumbent_reporting_and_three_lp(self):
        incumbent = configuration(60, "C0")
        candidate = extended.configuration(
            60,
            admission_base_bound=103814,
            profile="bounded-report",
            incumbent=incumbent,
        )
        self.assertIsInstance(candidate.collector, qs.BoundedYieldSieveConfig)
        report = incumbent.collector.candidate_bound or max(
            incumbent.collector.residual_bound,
            incumbent.collector.large_product_bound,
        )
        self.assertEqual(candidate.collector.candidate_bound, report)
        self.assertEqual(candidate.collector.residual_arity, 3)
        self.assertEqual(decode(asdict(candidate), qs), candidate)

    def test_sparse_batch_reservation_preserves_verified_recovery(self):
        base = qs.build_factor_base(10007 * 10009, bound=1000).factor_base
        old = replace(policy(), opaque_pairs=False)
        new = replace(bounded(), opaque_pairs=False)
        stores = []
        sizes = []
        for cfg in (old, new):
            engine = qs.SieveCollector(
                qs.qs_polynomial(base), base, config=cfg, budget=allowance()
            )
            result = engine.collect(1, 200)
            self.assertEqual(result.reason, "complete")
            stores.append(
                [
                    (a.relation_id, tuple(a.exponents), a.residual_product)
                    for a in engine._atoms.values()
                ]
            )
            sizes.append(engine._batch_reserve)
            for atom in engine._atoms.values():
                qs.verify_atomic(
                    atom,
                    base,
                    residual_bound=cfg.residual_bound,
                    large_prime_bound=cfg.large_prime_bound,
                    large_product_bound=cfg.large_product_bound,
                    residual_arity=3,
                    budget=allowance(),
                )
        self.assertEqual(stores[0], stores[1])
        self.assertLess(sizes[1], sizes[0])


if __name__ == "__main__":
    unittest.main()
