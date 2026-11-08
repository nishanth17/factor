"""Independent exact oracles for combined useful-yield mechanisms."""

import copy
import json
import random
import unittest
from collections import Counter
from dataclasses import FrozenInstanceError, replace
from math import gcd
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.common import utils
from v2.execution.budget import BudgetExhaustedError
from v2.qs import (
    AtomicRelation,
    CombinedYieldSieveConfig,
    FactorBase,
    QSJob,
    SieveCollector,
    SIQSConfig,
    SIQSJob,
    build_factor_base,
    qs_polynomial,
    verify_atomic,
    verify_combined,
)
from v2.qs.cofactor_policy import (
    KnownEndpoints,
    batch_divisors,
    split_cofactor,
    squfof,
)
from v2.qs.extraction import joint_exponent_product
from v2.qs.indexed_relations import (
    FactorRecord,
    RecordArena,
    decode_factors,
    encode_factors,
)
from v2.qs.large_prime_checkpoint import pack_store, restore_store
from v2.qs.large_primes import LargePrimeForest
from v2.qs.linear_algebra import (
    DependencySolver,
    filter_matrix,
    verify_dependency,
)
from v2.qs.short_cycles import ShortCycleIndex
from v2.qs.storage import RelationPublicationError
from v2.qs.yield_support import CheckedStore, UsefulTarget, residual_floor
from v2.tests.test_qs_dlp import admit, allowance
from v2.tests.test_qs_pipeline import dense_kernel, span
from v2.tests.test_siqs import mutate


def policy(**changes):
    values = dict(
        large_prime_bound=1000,
        large_product_bound=1000000,
        residual_bound=1000,
        block_width=32,
        division="full",
        score_policy="powers",
        memory_bytes=128 * 2**20,
    )
    values.update(changes)
    return adapt_resource_configuration(CombinedYieldSieveConfig, **values)


def collector(**changes):
    base = build_factor_base(91, bound=3).factor_base
    return SieveCollector(
        qs_polynomial(base), base, config=policy(**changes), budget=allowance()
    )


def encoded(worker):
    return json.loads(json.dumps(pack_store(worker)))


def pair_job(**changes):
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


class CofactorTests(unittest.TestCase):
    def test_known_divisor_skips_whole_cofactor_primality(self):
        from v2.qs.yield_collector import divide_yield_residual

        worker = collector(batch_policy="off")
        prime = worker._checked.prime

        def endpoint_only(value, budget):
            self.assertNotEqual(value, 15)
            return prime(value, budget)

        with patch.object(worker._checked, "prime", side_effect=endpoint_only):
            atom, divisor = divide_yield_residual(
                worker, 1, 30, 15, ((2, 1),), Counter(), hint=3
            )
        self.assertIsNone(divisor)
        self.assertEqual(atom.large_primes, (3, 5))

    def test_shared_batch_against_pairwise_divisor_oracle(self):
        generator = random.Random(21179)
        primes = (3, 5, 7, 11, 13, 17, 19)
        for _ in range(80):
            values = tuple(
                generator.choice(primes) * generator.choice(primes)
                for _ in range(8)
            )
            hints = batch_divisors(values, allowance())
            for value, hint in zip(values, hints):
                possible = any(
                    1 < gcd(value, other) < value for other in values
                )
                if possible:
                    self.assertIsNotNone(hint)
                if hint is not None:
                    self.assertTrue(utils.valid_divisor(hint, value))
        self.assertEqual(batch_divisors((15, 15), allowance()), (None, None))
        self.assertEqual(batch_divisors((15, 21, 35), allowance()), (3, 3, 5))

    def test_known_endpoint_product_saturation_and_prime_squares(self):
        known = KnownEndpoints(3)
        known.add((3, 5, 7, 11), allowance())
        self.assertEqual(known.primes, [3, 5, 7])
        self.assertEqual(known.divisor(15, allowance()), 3)
        self.assertEqual(known.divisor(49, allowance()), 7)
        self.assertIsNone(known.divisor(11 * 13, allowance()))
        with patch("v2.qs.cofactor_policy.gcd", side_effect=AssertionError):
            self.assertIsNone(known.divisor(7, allowance()))

    def test_square_form_recurrence_and_finite_failure(self):
        numbers = (15, 21, 35, 77, 91, 143, 221, 437, 1009 * 1013, 22117019)
        successes = 0
        for number in numbers:
            divisor = squfof(number, allowance())
            if divisor is not None:
                self.assertTrue(utils.valid_divisor(divisor, number))
                successes += 1
        self.assertGreaterEqual(successes, 5)
        self.assertEqual(squfof(101**2, allowance()), 101)
        self.assertIsNone(squfof(101, allowance(), steps=1))
        with self.assertRaises(BudgetExhaustedError):
            squfof(22117019, allowance(0))

    def test_resolved_splitters_prove_both_endpoints(self):
        worker = collector()
        for method in ("rho-v2", "portfolio-v2", "ecm-v2"):
            for residual in (3 * 5, 3**2, 97 * 101, 3 * 5 * 7):
                pair, reason, _ = split_cofactor(
                    residual, 1000, allowance(), worker._checked, method
                )
                if pair:
                    self.assertEqual(pair[0] * pair[1], residual)
                    self.assertTrue(
                        all(
                            utils.classify_prime(p) is utils.Primality.PROVEN
                            for p in pair
                        )
                    )
                elif residual == 105:
                    self.assertIn(reason, ("not_two_primes", "split_failed"))
                else:
                    self.fail((method, residual, reason))

    def test_domain_floor_handles_incomplete_custom_and_multiplier_bases(self):
        for multiplier in (1, 3):
            base = build_factor_base(
                10403, multiplier=multiplier, bound=40
            ).factor_base
            self.assertEqual(residual_floor(base, allowance()), 40)
            for omitted in base.entries[1:]:
                incomplete = FactorBase(
                    base.n,
                    base.multiplier,
                    base.bound,
                    tuple(e for e in base.entries if e is not omitted),
                )
                self.assertEqual(
                    residual_floor(incomplete, allowance()),
                    omitted.prime,
                )
        with self.assertRaises(TypeError):
            residual_floor(base, allowance(), 0)

    def test_hole_rejection_does_not_start_split_or_primality(self):
        from v2.qs.yield_collector import divide_yield_residual

        base = build_factor_base(10403, bound=40).factor_base
        worker = SieveCollector(
            qs_polynomial(base),
            base,
            config=policy(residual_bound=100),
            budget=allowance(),
        )
        stats = Counter()
        with patch.object(
            worker._checked,
            "prime",
            side_effect=AssertionError("unnecessary prime test"),
        ):
            self.assertEqual(
                divide_yield_residual(worker, 1, 105, 105, (), stats),
                (None, None),
            )
        self.assertEqual(worker._split_calls, 0)
        self.assertEqual(stats["domain_hole_rejections"], 1)


class CompactStoreTests(unittest.TestCase):
    def test_codec_boundaries_large_exponents_and_corruption(self):
        base = build_factor_base(10403, bound=100).factor_base
        pairs = tuple(
            (entry.prime, (1, 2, 127, 128, 65536)[i % 5])
            for i, entry in enumerate(base.entries)
        )
        data = encode_factors(pairs, base)
        self.assertEqual(
            tuple(decode_factors(data, base.primes, len(pairs))), pairs
        )
        arena = RecordArena(base)
        record = arena.factors(pairs, allowance())
        self.assertIsInstance(record, FactorRecord)
        self.assertEqual(len(record), len(pairs))
        self.assertEqual(tuple(record), pairs)
        with self.assertRaises(FrozenInstanceError):
            record._count = 1
        for bad in (
            b"\x00",
            b"\x80\x00",
            b"\x03\x01",
            data[:-1],
            data + b"\x02",
        ):
            with self.assertRaises(ValueError):
                tuple(decode_factors(bad, base.primes, len(pairs)))

    def test_slabs_preserve_snapshots_and_integer_cycle_provenance(self):
        worker = collector(duplicate_policy="pair-v1")
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
        for row in worker._combined:
            self.assertIsInstance(row._handles, bytes)
            self.assertEqual(len(row._handles), 8 * len(row.atom_ids))
            self.assertTrue(
                verify_combined(row, worker.factor_base, worker._atoms)
            )
        snapshot = worker._combined[0]
        original = (snapshot.atom_ids, tuple(snapshot.exponents))
        worker.set_polynomial(
            worker.polynomial, worker._roots, retain_relations=False
        )
        self.assertEqual(
            (snapshot.atom_ids, tuple(snapshot.exponents)), original
        )

    def test_unowned_eviction_releases_handles_without_aliasing(self):
        worker = collector(max_partials=1)
        original = admit(worker, -7, scalar=41)[0]
        packed = worker._atoms[original.relation_id]
        admit(worker, 1, (3, 5))
        self.assertNotIn(original.relation_id, worker._arena.slots)
        self.assertEqual(tuple(packed.exponents), original.exponents)
        self.assertTrue(
            verify_atomic(packed, worker.factor_base, residual_bound=1000)
        )

    def test_roundtrip_policy_binding_and_early_decode_refusal(self):
        worker = collector(duplicate_policy="pair-v1")
        for x, pair in ((-5, (3, 11)), (-4, (5, 11)), (1, (3, 5))):
            admit(worker, x, pair)
        payload = encoded(worker)
        restored = collector(duplicate_policy="pair-v1")
        restore_store(payload, restored, allowance())
        self.assertEqual(encoded(restored), payload)
        for change in (
            lambda p: p.update(factor_base_identity="foreign"),
            lambda p: p["atoms"][0][3].update(count=10**9),
            lambda p: p["atoms"][0][3].update({"index-v2": "AA=="}),
            lambda p: p["combined"][0][0].append(1000000),
            lambda p: p.update(cofactor_endpoints=[15]),
        ):
            bad, empty = (
                copy.deepcopy(payload),
                collector(duplicate_policy="pair-v1"),
            )
            change(bad)
            with self.assertRaises((ValueError, TypeError)):
                restore_store(bad, empty, allowance())
            self.assertFalse(empty._atoms)
        empty = collector(duplicate_policy="pair-v1")
        empty._workspace = empty.config.memory_limit - 1
        with patch(
            "v2.qs.indexed_relations.import_factors",
            side_effect=AssertionError("decoded before preflight"),
        ):
            with self.assertRaises(MemoryError):
                restore_store(payload, empty, allowance())

    def test_failed_auxiliary_publication_poisoned_not_resumable(self):
        worker = collector()
        with patch.object(
            worker._arena, "commit", side_effect=MemoryError("injected")
        ):
            with self.assertRaises(RelationPublicationError):
                admit(worker, 1, (3, 5))
        with self.assertRaises(RelationPublicationError):
            worker.collect(2, 3)


class ReuseAndDuplicateTests(unittest.TestCase):
    def test_mutated_payload_and_foreign_target_cannot_reuse_proof(self):
        worker = collector()
        source = admit(worker, 1, (3, 5))[0]
        atom = worker._atoms[source.relation_id]
        worker._checked.atom(atom, allowance())
        with patch(
            "v2.qs.yield_support._verify_atomic",
            side_effect=AssertionError("repeated exact verification"),
        ):
            worker._checked.atom(atom, allowance())
        object.__setattr__(atom, "sign", -atom.sign)
        with self.assertRaises(ValueError):
            worker._checked.atom(atom, allowance())
        foreign = build_factor_base(143, bound=3).factor_base
        with self.assertRaises(ValueError):
            CheckedStore(foreign, worker.config).atom(source, allowance())

    def test_mutated_packed_bytes_cannot_reuse_proof(self):
        worker = collector()
        source = admit(worker, 1, (3, 5))[0]
        atom = worker._atoms[source.relation_id]
        record = atom.exponents
        record._block[record._offset] = 0
        with self.assertRaises(ValueError):
            worker._checked.atom(atom, allowance())

    def test_canonical_mirror_elision_preserves_distinct_square_roots(self):
        job = pair_job()
        first, second = job.collector._rows
        mirror = AtomicRelation(
            first.polynomial, -241, 1, tuple(first.exponents)
        )
        stats = Counter()
        job.collector._admit(mirror, stats)
        self.assertEqual(stats["congruence_duplicates"], 1)
        self.assertEqual(len(job.collector._rows), 2)
        self.assertEqual(job._solve(False), 101)
        self.assertIs(job.collector._rows[1], second)
        self.assertEqual(job.stats["filter_calls"], 0)

    def test_zero_and_storage_cap_override_target(self):
        for reason in ("atom_limit", "relation_limit", "memory_limit"):
            job = pair_job(duplicate_policy="identity")
            self.assertIsNone(job._solve(False))
            job.storage_reason = reason
            self.assertEqual(job.run().divisor, 101)
        job = pair_job(duplicate_policy="identity")
        square = AtomicRelation(job.collector.polynomial, 0, 1, ())
        job.collector._admit(square, Counter())
        self.assertIn(job._solve(False), (101, 103))

    def test_active_columns_not_width_and_feedback_is_bounded(self):
        target = UsefulTarget()
        target.columns = 1 | (1 << 5000)
        self.assertEqual(target.target(32), 34)
        target.filtered(
            100,
            dict(output_rows=30, output_columns=50, zero_dependencies=0),
            32,
            4,
            1000,
        )
        self.assertTrue(104 <= target.next_target <= 125)

    def test_pair_budget_refusal_keeps_pending_cursor(self):
        job = pair_job()
        job.budget.work_limit = job.budget.used + 3
        with self.assertRaises(BudgetExhaustedError):
            job._solve(False)
        self.assertLess(job._pair_index.checked_count, 2)
        job.budget.work_limit = 10**13
        self.assertEqual(job._solve(False), 101)


class GraphAndMatrixTests(unittest.TestCase):
    def test_direct_cycle_does_not_walk_forest(self):
        worker = collector(duplicate_policy="identity")
        admit(worker, -5, (3, 11))
        with patch.object(worker._graph, "plan", side_effect=AssertionError):
            admit(worker, -15, (3, 11))
        self.assertTrue(
            verify_combined(
                worker._combined[0], worker.factor_base, worker._atoms
            )
        )

    def test_direct_and_two_edge_fallback_short_circuit(self):
        index = ShortCycleIndex(2, strategy="short-v2")
        for i in range(100):
            index.add(str(i), 1, i + 3, weight=(i,))
        index.add("direct", 1, 2, weight=(1,))
        self.assertEqual(
            index.shorter_path(1, 2, ("a", "b"), allowance()), (("direct",), 1)
        )
        index.remove("direct")
        self.assertEqual(
            index.shorter_path(1, 2, ("a", "b"), allowance()), (("a", "b"), 1)
        )

    def test_bidirectional_cycle_space_against_dense_incidence_kernel(self):
        rng = random.Random(601231)
        for _ in range(80):
            forest, index, edges, cycles = (
                LargePrimeForest(),
                ShortCycleIndex(30, strategy="short-v2"),
                [],
                [],
            )
            for i in range(12):
                a, b = rng.randrange(1, 8), rng.randrange(1, 8)
                plan = forest.plan(a, b, str(i), 100, allowance())
                if plan.cycle is not None:
                    path, visits = index.shorter_path(
                        a, b, plan.cycle, allowance()
                    )
                    self.assertLessEqual(visits, 30)
                    self.assertEqual(len(set(path)), len(path))
                    self.assertLessEqual(len(path), len(plan.cycle))
                    plan = replace(plan, cycle=path)
                    cycles.append(sum(1 << int(j) for j in path) | (1 << i))
                forest.commit(plan)
                index.add(str(i), a, b, weight=(rng.randrange(20),))
                edges.append((1 << a) ^ (1 << b))
            self.assertEqual(span(cycles), dense_kernel(edges))
            self.assertEqual(len(span(cycles)), 1 << len(cycles))

    def test_fill_order_preserves_dense_kernel_and_light_view_lifts(self):
        rng = random.Random(78071)
        for _ in range(60):
            rows = tuple(rng.randrange(64) for _ in range(10))
            matrix = filter_matrix(
                rows,
                weight_two=True,
                merge_policy="fill-v2",
                budget=allowance(),
            )
            kernels = DependencySolver(matrix, budget=allowance()).run()
            self.assertEqual(span(kernels), dense_kernel(rows))
            view = filter_matrix(
                rows,
                weight_two=True,
                merge_policy="fill-v2",
                light_rows=True,
                surplus=1,
                budget=allowance(),
            )
            for mask in DependencySolver(view, budget=allowance()).run():
                self.assertTrue(verify_dependency(mask, rows))

    def test_joint_product_against_independent_modular_powers(self):
        rng = random.Random(35183)
        for modulus in (91, 10403, 1000003 * 1000033):
            for _ in range(30):
                totals = {
                    p: 2 * rng.randrange(1, 65536)
                    for p in (2, 3, 5, 7, 11, 13)
                }
                expected = 1
                for prime, exponent in totals.items():
                    expected = (
                        expected * pow(prime, exponent // 2, modulus) % modulus
                    )
                self.assertEqual(
                    joint_exponent_product(totals, modulus, allowance()),
                    expected,
                )


class CombinedPipelineTests(unittest.TestCase):
    def test_light_view_retries_full_store_inside_one_solve(self):
        import v2.qs.pipeline as pipeline

        job = pair_job(duplicate_policy="identity", preparation_policy="fixed")
        job.row_excess = 0
        calls = []
        original = pipeline.filter_matrix

        def filter_view(*args, **kwargs):
            calls.append(kwargs["light_rows"])
            matrix = original(*args, **kwargs)
            stats = dict(matrix.stats)
            stats["selected_rows"] = 1 if kwargs["light_rows"] else 2
            return replace(matrix, stats=stats)

        with patch.object(pipeline, "filter_matrix", side_effect=filter_view):
            with patch.object(
                pipeline.DependencyExtractor, "run", return_value=None
            ):
                self.assertIsNone(job._solve(False))
        self.assertEqual(calls, [True, False])
        self.assertEqual(len(job.collector.matrix_relations), 2)
        self.assertFalse(job._retry_full_view)

    def test_batch_and_synchronous_collection_have_identical_exact_atoms(self):
        results = []
        for batch in ("off", "shared-v2"):
            worker = collector(
                batch_policy=batch,
                duplicate_policy="identity",
                split_policy="rho-v2",
            )
            result = worker.collect(-30, 70)
            for atom in result.atoms:
                self.assertTrue(
                    verify_atomic(
                        atom,
                        worker.factor_base,
                        residual_bound=1000,
                        large_prime_bound=1000,
                        large_product_bound=1000000,
                    )
                )
            results.append(
                (
                    result.next_position,
                    result.divisor,
                    [
                        (a.relation_id, tuple(a.exponents), a.residual_primes)
                        for a in result.atoms
                    ],
                )
            )
        self.assertEqual(*results)

    def test_full_policy_roundtrip_and_corrupt_progress(self):
        config = SIQSConfig(
            base_bound=100,
            half_width=256,
            max_half_width=256,
            solve_policy="collection",
            collector=policy(block_width=16),
        )
        job = SIQSJob(
            10007 * 10009, seed=73, config=config, budget=allowance()
        )
        job.run(max_blocks=1)
        checkpoint = json.loads(json.dumps(job.checkpoint()))
        self.assertEqual(checkpoint["version"], 12)
        restored = SIQSJob.from_checkpoint(
            checkpoint, budget=allowance(), config=config
        )
        self.assertEqual(
            restored.engine.next_position, job.engine.next_position
        )
        self.assertEqual(
            encoded(restored.engine.collector), encoded(job.engine.collector)
        )
        for change in (
            lambda p: p["config"]["collector"].update(yield_policy="yield-v1"),
            lambda p: p["engine"]["yield_progress"]["active_target"].update(
                columns="0xffff"
            ),
        ):
            with self.assertRaises((ValueError, KeyError)):
                SIQSJob.from_checkpoint(
                    mutate(checkpoint, change), budget=allowance()
                )
        for candidate in (job, restored):
            result = candidate.run(max_blocks=20)
            if result.divisor is not None:
                self.assertEqual(
                    result.divisor * result.cofactor, 10007 * 10009
                )
            self.assertGreaterEqual(result.stats["workspace_bytes"], 0)


if __name__ == "__main__":
    unittest.main()
