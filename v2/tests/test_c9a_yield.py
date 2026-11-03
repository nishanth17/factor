"""Independent admission, cycle-space and arithmetic/resume oracles."""

import copy
import json
import random
import unittest
from collections import Counter, defaultdict
from dataclasses import FrozenInstanceError, asdict, replace
from math import isqrt
from unittest.mock import patch

from v2.benchmarks.qs.c9a.control import control_modules
from v2.execution.budget import BudgetExhaustedError
from v2.qs import (
    AtomicRelation,
    GrowingDoubleLargeSieveConfig,
    Polynomial,
    QSJob,
    SieveCollector,
    SIQSConfig,
    SIQSJob,
    YieldSieveConfig,
    build_factor_base,
    qs_polynomial,
    verify_atomic,
)
from v2.qs.large_prime_checkpoint import pack_store, restore_store
from v2.qs.large_primes import LargePrimeForest
from v2.qs.packed_relations import PackedAtomicRelation, pack_relation
from v2.qs.short_cycles import ShortCycleIndex
from v2.qs.storage import RelationPublicationError
from v2.tests.test_qs_dlp import admit
from v2.tests.test_qs_pipeline import dense_kernel, span
from v2.tests.test_siqs import allowance, mutate


def policy(**changes):
    values = dict(
        large_prime_bound=1000,
        large_product_bound=1000000,
        residual_bound=1000,
        block_width=32,
        division="full",
        score_policy="powers",
        memory_bytes=64 * 2**20,
    )
    values.update(changes)
    return YieldSieveConfig(**values)


def collector(**changes):
    base = build_factor_base(91, bound=3).factor_base
    return SieveCollector(
        qs_polynomial(base),
        base,
        config=policy(**changes),
        budget=allowance(10**13),
    )


def store_json(worker):
    return json.loads(json.dumps(pack_store(worker)))


class AdmissionTests(unittest.TestCase):
    def test_endpoint_product_rejection_precedes_split_but_keeps_slp(self):
        from v2.qs.large_prime_collector import divide_large_residual

        target = collector(small_prime_bound=2, large_prime_bound=5,
                           large_product_bound=25, residual_bound=41)
        counts = Counter()
        # x=1 gives30=2*15. The composite15 cannot have p<=2 and q<=5,
        # even though it lies below the independent scalar allowance41.
        self.assertEqual(divide_large_residual(target, 1, 30, 15,
                                              ((2, 1),), counts), (None, None))
        self.assertEqual(target._split_calls, 0)
        self.assertEqual(counts["endpoint_product_rejections"], 1)
        # x=-7 gives-82=-2*41; a scalar prime41 remains admissible.
        atom, divisor = divide_large_residual(target, -7, -82, 41,
                                              ((2, 1),), counts)
        self.assertIsNone(divisor)
        self.assertEqual(atom.residual, 41)

    def test_a_and_multiplier_exceptions_keep_exact_norm_factors(self):
        base = build_factor_base(91, multiplier=3, bound=3).factor_base
        polynomial = Polynomial(91, 3, 2, 17)
        worker = SieveCollector(
            polynomial,
            base,
            config=policy(small_prime_bound=3),
            budget=allowance(),
        )
        # x=2 gives U=21 and U²-3N=168=2³*3*7. A=2 is
        # already in the full norm; multiplier prime 3 is an admissible LP.
        atom = AtomicRelation(polynomial, 2, 1, ((2, 3),), large_primes=(3, 7))
        self.assertTrue(
            verify_atomic(
                atom,
                base,
                residual_bound=1000,
                large_prime_bound=1000,
                large_product_bound=1000000,
                small_prime_bound=3,
            )
        )
        self.assertIsNone(worker._admit(atom, Counter()))
        self.assertIn(atom.relation_id, worker._atoms)

    def test_independent_trial_division_domain(self):
        base = build_factor_base(1000003 * 1000033, bound=50).factor_base
        poly = qs_polynomial(base)
        config = policy(
            large_prime_bound=211,
            large_product_bound=211**2,
            residual_bound=211,
            small_prime_bound=101,
        )
        expected = set()
        for position in range(-240, 241):
            remaining = abs(poly.value(position))
            if not remaining:
                continue
            for prime in base.primes:
                while remaining % prime == 0:
                    remaining //= prime
            factors = []
            value = remaining
            for divisor in range(2, isqrt(value) + 1):
                while value % divisor == 0:
                    factors.append(divisor)
                    value //= divisor
            if value > 1:
                factors.append(value)
            if (
                not factors
                or len(factors) == 1
                and factors[0] <= 211
                or len(factors) == 2
                and factors[0] <= 101
                and factors[1] <= 211
            ):
                expected.add(position)
        worker = SieveCollector(
            poly, base, config=config, budget=allowance(10**13)
        )
        result = worker.collect(-240, 241)
        self.assertEqual(result.reason, "complete")
        self.assertEqual({a.position for a in result.atoms}, expected)

    def test_equality_square_order_certainty_and_restore_policy(self):
        worker = collector(small_prime_bound=3)
        pair = admit(worker, 1, (3, 5))[0]
        square = admit(worker, 0, (3, 3))[0]
        for atom in (pair, square):
            self.assertTrue(
                verify_atomic(
                    atom,
                    worker.factor_base,
                    residual_bound=1000,
                    large_prime_bound=5,
                    large_product_bound=15,
                    small_prime_bound=3,
                )
            )
            with self.assertRaises(ValueError):
                verify_atomic(
                    atom,
                    worker.factor_base,
                    residual_bound=1000,
                    large_prime_bound=5,
                    large_product_bound=15,
                    small_prime_bound=2,
                )
        with self.assertRaises(ValueError):
            replace(pair, large_primes=(5, 3))
        with self.assertRaises(ValueError):
            verify_atomic(
                replace(square, large_primes=(9, 9)),
                worker.factor_base,
                large_prime_bound=9,
                large_product_bound=81,
                small_prime_bound=9,
            )
        encoded = store_json(worker)
        target = collector(small_prime_bound=2)
        with self.assertRaises(ValueError):
            restore_store(encoded, target, allowance())
        self.assertFalse(target._atoms)

    def test_split_interruption_charges_attempt_without_publication(self):
        worker = collector(small_prime_bound=3)
        stats = defaultdict(int)
        from v2.qs.large_prime_collector import divide_large_residual

        worker.budget.work_limit = worker.budget.used + 30
        with self.assertRaises(BudgetExhaustedError):
            divide_large_residual(worker, 1, 30, 15, ((2, 1),), stats)
        self.assertEqual(worker._split_calls, 1)
        self.assertFalse(worker._atoms)


class GraphAndRecordTests(unittest.TestCase):
    def test_interrupted_shortcut_publication_requires_prior_checkpoint(self):
        target = collector(cycle_policy="short-v1", record_policy="packed-v1")
        saved = store_json(target)
        with patch.object(
            target._short_cycles, "add", side_effect=MemoryError
        ):
            with self.assertRaises(RelationPublicationError):
                admit(target, -5, (3, 11))
        with self.assertRaises(RelationPublicationError):
            pack_store(target)
        restored = collector(
            cycle_policy="short-v1", record_policy="packed-v1"
        )
        restore_store(saved, restored, allowance())
        self.assertEqual(store_json(restored), saved)

    def test_shortcut_cycles_preserve_independent_incidence_kernel(self):
        generator = random.Random(901711)
        for _ in range(30):
            forest, index, edges, cycles = (
                LargePrimeForest(),
                ShortCycleIndex(40),
                [],
                [],
            )
            for i in range(10):
                a, b = generator.randrange(1, 7), generator.randrange(1, 7)
                plan = forest.plan(a, b, str(i), 100, allowance())
                if plan.cycle is not None:
                    path, visits = index.shorter_path(
                        a, b, plan.cycle, allowance()
                    )
                    self.assertLessEqual(visits, 40)
                    self.assertLessEqual(len(path), len(plan.cycle))
                    plan = replace(plan, cycle=path)
                    cycles.append(sum(1 << int(j) for j in path) | (1 << i))
                forest.commit(plan)
                index.add(str(i), a, b)
                edges.append((1 << a) ^ (1 << b))
            self.assertEqual(span(cycles), dense_kernel(edges))
            self.assertEqual(len(span(cycles)), 1 << len(cycles))

    def test_prior_closing_edge_shortcut_and_packed_roundtrip(self):
        worker = collector(cycle_policy="short-v1", record_policy="packed-v1")
        for position, pair in (
            (-5, (3, 11)),
            (-4, (5, 11)),
            (1, (3, 5)),
            (-21, (3, 5)),
            (0, (3, 3)),
        ):
            admit(worker, position, pair)
        self.assertEqual(
            [len(r.atom_ids) for r in worker._combined], [3, 2, 1]
        )
        self.assertIsInstance(
            next(iter(worker._atoms.values())), PackedAtomicRelation
        )
        encoded = store_json(worker)
        self.assertEqual(encoded["record_policy"], "packed-v1")
        target = collector(cycle_policy="short-v1", record_policy="packed-v1")
        restore_store(encoded, target, allowance(10**13))
        self.assertEqual(store_json(target), encoded)
        self.assertEqual(
            target._short_cycles.neighbors, worker._short_cycles.neighbors
        )
        bad = copy.deepcopy(encoded)
        bad["atoms"][0][3]["u32be64"] += "="
        empty = collector(cycle_policy="short-v1", record_policy="packed-v1")
        with self.assertRaises(ValueError):
            restore_store(bad, empty, allowance())
        self.assertFalse(empty._atoms)

    def test_packed_payload_aliasing_and_eviction_lifetimes(self):
        worker = collector(
            record_policy="packed-v1", cycle_policy="short-v1", max_partials=1
        )
        atom = admit(worker, -7, scalar=41)[0]
        packed = worker._atoms[atom.relation_id]
        self.assertNotIn("exponents", packed.__dict__)
        self.assertEqual(asdict(packed), asdict(atom))
        with self.assertRaises(FrozenInstanceError):
            packed._packed_exponents = b""
        admit(worker, 1, (3, 5))
        self.assertNotIn(atom.relation_id, worker._atoms)
        self.assertNotIn(atom.relation_id, worker._short_cycles.edges)
        self.assertEqual(packed.exponents, atom.exponents)
        self.assertTrue(
            verify_atomic(packed, worker.factor_base, residual_bound=1000)
        )
        self.assertIs(pack_relation(packed, allowance()), packed)

    def test_search_cap_and_budget_leave_index_unchanged(self):
        index = ShortCycleIndex(1)
        index.add("a", 3, 5)
        index.add("b", 5, 7)
        before = copy.deepcopy(index.__dict__)
        self.assertEqual(
            index.shorter_path(3, 7, ("x", "y", "z"), allowance())[0],
            ("x", "y", "z"),
        )
        self.assertEqual(index.__dict__, before)
        with self.assertRaises(BudgetExhaustedError):
            index.shorter_path(3, 7, None, allowance(0))
        self.assertEqual(index.__dict__, before)


def parity_witness(**changes):
    base = build_factor_base(10403, bound=40).factor_base
    poly = qs_polynomial(base)
    job = QSJob(
        poly,
        base,
        200,
        201,
        config=policy(**changes),
        budget=allowance(),
        collection_continues=True,
        row_excess=32,
        filter_row_growth=32,
    )
    for position, exponents in (
        (37, ((2, 1), (7, 3), (13, 1))),
        (149, ((2, 1), (7, 1), (13, 1), (17, 2))),
    ):
        atom = AtomicRelation(poly, position, 1, exponents)
        verify_atomic(atom, base)
        job.collector._admit(atom, Counter())
    return job


class PreparationTests(unittest.TestCase):
    def test_target_does_not_delay_an_immediate_square(self):
        job = parity_witness(
            preparation_policy="target-v1", duplicate_policy="pair-v1"
        )
        square = AtomicRelation(job.collector.polynomial, 0, 1, ())
        verify_atomic(square, job.collector.factor_base)
        job.collector._rows = [square]
        job.collector._full = [square]
        job.collector._atoms = {square.relation_id: square}
        self.assertIn(job._solve(False), (101, 103))

    def test_trivial_equal_parity_pair_remains_available(self):
        job = parity_witness(duplicate_policy="pair-v1")
        first, last = job.collector._rows
        mirror = AtomicRelation(first.polynomial, -241, 1, first.exponents)
        verify_atomic(mirror, job.collector.factor_base)
        job.collector._admit(mirror, Counter())
        job.collector._rows = [first, mirror, last]
        self.assertEqual(job._solve(False), 101)
        self.assertEqual(job.stats["pair_trivial"], 1)
        self.assertEqual(job._pair_index.checked_count, 3)
        self.assertEqual(len(job.collector._rows), 3)

    def test_target_deferral_and_equal_parity_extracts_without_discarding(
        self,
    ):
        job = parity_witness(preparation_policy="target-v1")
        self.assertIsNone(job._solve(False))
        self.assertEqual(job.stats["filter_calls"], 0)
        self.assertEqual(job._solve(True), 101)
        job = parity_witness(
            preparation_policy="target-v1",
            duplicate_policy="pair-v1",
            record_policy="packed-v1",
        )
        self.assertEqual(job._solve(False), 101)
        self.assertEqual(job.stats["filter_calls"], 0)
        self.assertEqual(job.stats["pair_dependencies"], 1)
        self.assertEqual(len(job.collector._rows), 2)

    def test_pair_refusal_keeps_pending_cursor_and_charges_resume(self):
        job = parity_witness(duplicate_policy="pair-v1")
        job.budget.work_limit = job.budget.used + 30
        with self.assertRaises(BudgetExhaustedError):
            job._solve(False)
        self.assertLess(job._pair_index.checked_count, 2)
        used = job.budget.used
        job.budget.work_limit = 10**9
        self.assertEqual(job._solve(False), 101)
        self.assertGreater(job.budget.used, used)

    def test_target_still_extracts_at_each_storage_cap(self):
        for reason in ("atom_limit", "relation_limit", "memory_limit"):
            job = parity_witness(preparation_policy="target-v1")
            job.storage_reason = reason
            self.assertEqual(job.run().divisor, 101)
            self.assertTrue(job.storage_solve_done)


class CompatibilityTests(unittest.TestCase):
    def test_disabled_path_matches_immutable_current_mainline(self):
        import v2.qs as current
        from v2.benchmarks.qs.c9a.configuration import decode

        control, budgets = control_modules()
        config = SIQSConfig(
            base_bound=100,
            half_width=32,
            max_half_width=32,
            collector=GrowingDoubleLargeSieveConfig(
                residual_bound=1000,
                large_prime_bound=1000,
                large_product_bound=1000000,
                memory_bytes=64 * 2**20,
                block_width=16,
            ),
            memory_bytes=64 * 2**20,
        )
        outcomes = []
        for module in (current, control):
            budget = (
                allowance(10**9)
                if module is current
                else budgets.Budget(
                    work_limit=10**9, seconds=None, cpu_seconds=None
                )
            )
            job = module.SIQSJob(
                4001 * 5003,
                seed=7,
                config=decode(asdict(config), module),
                budget=budget,
            )
            result = job.run(max_blocks=3)
            store = job.engine.collector
            outcomes.append(
                (
                    result.divisor,
                    result.reason,
                    budget.used,
                    [asdict(a) for a in store._atoms.values()],
                    [asdict(r) for r in store._rows],
                    result.stats["engine"]["collector_counts"],
                )
            )
        self.assertEqual(*outcomes)

    def test_policy_bound_resume_and_corrupt_or_legacy_schema(self):
        config = SIQSConfig(
            base_bound=100,
            half_width=32,
            max_half_width=32,
            memory_bytes=64 * 2**20,
            solve_policy="collection",
            collector=policy(
                small_prime_bound=101,
                cycle_policy="short-v1",
                duplicate_policy="pair-v1",
                record_policy="packed-v1",
                preparation_policy="target-v1",
            ),
        )
        job = SIQSJob(4001 * 5003, config=config, budget=allowance())
        job.run(max_blocks=1)
        saved = job.checkpoint()
        self.assertEqual(saved["version"], 11)
        restored = SIQSJob.from_checkpoint(
            saved, budget=allowance(), config=config
        )
        self.assertEqual(
            store_json(restored.engine.collector),
            store_json(job.engine.collector),
        )
        self.assertGreater(
            restored.budget.used, saved["resources"]["work_used"]
        )
        for key in (
            "small_prime_bound",
            "preparation_policy",
            "cycle_policy",
            "duplicate_policy",
            "record_policy",
        ):
            bad = mutate(saved, lambda p: p["config"]["collector"].pop(key))
            with self.assertRaises(ValueError):
                SIQSJob.from_checkpoint(bad, budget=allowance())
        bad = mutate(saved, lambda p: p.update(version=7))
        bad["version"] = 7
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(bad, budget=allowance())
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                saved,
                budget=allowance(),
                allow_extension=True,
                config=replace(
                    config,
                    collector=replace(config.collector, small_prime_bound=103),
                ),
            )


if __name__ == "__main__":
    unittest.main()
