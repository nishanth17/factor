"""Independent multigraph kernels and exact, bounded DLP provenance."""

import copy
import json
import random
import unittest
from collections import defaultdict
from dataclasses import asdict, replace
from math import gcd
from unittest.mock import patch

from v2 import arithmetic
from v2.budget import Budget, BudgetExhaustedError
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import (
    AtomicRelation,
    DoubleLargeSieveConfig,
    Polynomial,
    SieveCollector,
    SieveConfig,
    SIQSConfig,
    SIQSJob,
    build_factor_base,
    combine_relations,
    extract_dependency,
    prepare_relations,
    qs_polynomial,
    verify_atomic,
    verify_combined,
)
from v2.qs.large_prime_checkpoint import pack_store, restore_store
from v2.qs.large_primes import (
    CycleTooLongError,
    LargePrimeForest,
    graph_reserve,
    split_two_primes,
)
from v2.tests.test_qs_pipeline import dense_kernel, span
from v2.tests.test_siqs import mutate


def allowance(work=10**13, **options):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None, **options)


def configuration(**options):
    values = dict(
        block_width=1,
        residual_bound=1000,
        large_prime_bound=1000,
        large_product_bound=1000000,
        max_partials=32,
    )
    values.update(options)
    return DoubleLargeSieveConfig(**values)


def worker(*, backend="python-int", **options):
    base = build_factor_base(91, bound=3, backend=backend).factor_base
    return SieveCollector(
        qs_polynomial(base),
        base,
        config=configuration(**options),
        budget=allowance(),
    )


def encoded(collector):
    return json.loads(
        json.dumps(pack_store(collector), default=arithmetic.json_integer)
    )


def admit(collector, position, pair=(), *, polynomial=None, scalar=1):
    """Independent trial division of the complete integer norm by base 2."""
    polynomial = polynomial or collector.polynomial
    value = polynomial.a * polynomial.value(position)
    remaining = abs(value) // polynomial.square_coefficient**2
    exponent = 0
    while remaining % 2 == 0:
        remaining //= 2
        exponent += 1
    expected = pair[0] * pair[1] if pair else scalar
    if remaining != expected:
        raise AssertionError((remaining, expected))
    atom = AtomicRelation(
        polynomial,
        position,
        -1 if value < 0 else 1,
        ((2, exponent),) if exponent else (),
        scalar,
        large_primes=pair,
    )
    verify_atomic(
        atom,
        collector.factor_base,
        residual_bound=1000,
        large_prime_bound=1000,
        large_product_bound=1000000,
        budget=collector.budget,
    )
    stats = defaultdict(int)
    refusal = collector._admit(atom, stats)
    return atom, refusal, stats


class GraphTests(unittest.TestCase):
    def test_complete_spans_against_boolean_incidence_oracle(self):
        generator = random.Random(59203)
        fixtures = [
            [(1, 7), (3, 5), (5, 11), (3, 11), (13, 13), (3, 5)],
            [(2, 3), (5, 7), (2, 3), (5, 7), (11, 11)],
        ]
        fixtures += [
            [
                (generator.randrange(1, 9), generator.randrange(1, 9))
                for _ in range(10)
            ]
            for _ in range(40)
        ]
        for edges in fixtures:
            forest, cycles = LargePrimeForest(), []
            for i, (a, b) in enumerate(edges):
                plan = forest.plan(a, b, str(i), 64, allowance())
                if plan.cycle is not None:
                    cycles.append(
                        sum(1 << int(j) for j in plan.cycle) | (1 << i)
                    )
                forest.commit(plan)
            expected = dense_kernel([(1 << a) ^ (1 << b) for a, b in edges])
            self.assertEqual(span(cycles), expected)
            self.assertEqual(len(expected), 1 << len(cycles))

    def test_eviction_preserves_complete_retained_span_and_ownership(self):
        generator = random.Random(1289)
        forest, retained, cycles = LargePrimeForest(), {}, []
        evictions = 0
        for i in range(35):
            a, b = generator.randrange(1, 22), generator.randrange(1, 22)
            plan = forest.plan(a, b, str(i), 3, allowance())
            for identity in plan.evicted:
                self.assertFalse(any(identity in cycle for cycle in cycles))
                del retained[identity]
            evictions += len(plan.evicted)
            retained[str(i)] = (a, b)
            if plan.cycle is not None:
                cycles.append(plan.cycle + (str(i),))
            forest.commit(plan)
            self.assertLessEqual(len(forest.unowned), 3)
            identities = list(retained)
            rows = [(1 << a) ^ (1 << b) for a, b in retained.values()]
            masks = [
                sum(1 << identities.index(j) for j in cycle)
                for cycle in cycles
            ]
            self.assertEqual(span(masks), dense_kernel(rows))
        self.assertGreater(evictions, 0)

    def test_planning_is_transactional_and_long_paths_are_explicit(self):
        forest = LargePrimeForest()
        for i in range(256):
            forest.commit(forest.plan(i + 1, i + 2, str(i), 300, allowance()))
        saved = copy.deepcopy(forest.__dict__)
        with self.assertRaises(CycleTooLongError):
            forest.plan(1, 257, "closing", 300, allowance())
        self.assertEqual(forest.__dict__, saved)
        for work in (0, 1, 3, 20):
            with self.assertRaises(BudgetExhaustedError):
                forest.plan(400, 401, "new", 256, allowance(work))
            self.assertEqual(forest.__dict__, saved)
        for args in (
            (0, 0, "bad"),
            (10**12 + 1, 2, "bad"),
            (1, 1, ""),
            (1, 1, "0"),
        ):
            with self.assertRaises(ValueError):
                forest.plan(*args, 256, allowance())
        self.assertEqual(graph_reserve(256), 32768 + 4096 * 256)


class RelationTests(unittest.TestCase):
    def test_disconnected_triangle_and_nontrivial_repeated_prime_loop(self):
        collector = worker()
        admit(collector, -7, scalar=41)
        atoms = [
            admit(collector, x, pair)[0]
            for x, pair in ((-5, (3, 11)), (1, (3, 5)), (-4, (5, 11)))
        ]
        triangle = collector._combined[0]
        self.assertEqual(
            set(triangle.atom_ids), {a.relation_id for a in atoms}
        )
        self.assertEqual(triangle.square_correction, 165 % 91)
        self.assertEqual(triangle.exponents, ((2, 2),))
        self.assertEqual(len(collector.partial_ids), 1)
        admit(collector, 0, (3, 3))
        loop = collector._combined[-1]
        prepared = prepare_relations(
            (loop,), collector.factor_base, collector._atoms
        )
        result = extract_dependency(prepared, 1, budget=allowance())
        self.assertEqual(result.x**2 % 91, result.y**2 % 91)
        self.assertEqual(
            {gcd(result.x - result.y, 91), gcd(result.x + result.y, 91)},
            {7, 13},
        )
        self.assertEqual(result.divisor * (91 // result.divisor), 91)

    def test_parallel_edges_and_known_square_shared_with_large_prime(self):
        for backend in ("python-int", "gmpy2-mpz"):
            try:
                collector = worker(backend=backend)
            except arithmetic.BackendUnavailableError:
                continue
            admit(collector, -7, scalar=41)
            atom, _, _ = admit(collector, 1, (3, 5))
            polynomial = Polynomial(collector.factor_base.n, 1, 9, 1, 3)
            admit(collector, 2, (3, 5), polynomial=polynomial)
            row = collector._combined[0]
            self.assertEqual(row.square_correction, 45)
            self.assertEqual(row.exponents, ((2, 2),))
            result = extract_dependency(
                prepare_relations(
                    (row,), collector.factor_base, collector._atoms
                ),
                1,
                budget=allowance(),
            )
            self.assertEqual(
                {gcd(result.x - result.y, 91), gcd(result.x + result.y, 91)},
                {7, 13},
            )
            with self.assertRaises(ValueError):
                prepare_relations(
                    (atom,), collector.factor_base, collector._atoms
                )
            with self.assertRaises(ValueError):
                verify_combined(
                    replace(row, square_correction=44),
                    collector.factor_base,
                    collector._atoms,
                )

    def test_certainty_shapes_and_exact_reconstruction(self):
        collector = worker()
        with self.assertRaises(ValueError):
            AtomicRelation(
                collector.polynomial, 1, 1, ((2, 1),), 3, large_primes=(3, 5)
            )
        for pair in ((5, 3), (3,), (3, 5, 7), (1, 9)):
            with self.assertRaises((ValueError, TypeError)):
                AtomicRelation(
                    collector.polynomial, 1, 1, ((2, 1),), large_primes=pair
                )
        valid = AtomicRelation(
            collector.polynomial, 1, 1, ((2, 1),), large_primes=(3, 5)
        )
        with self.assertRaises(ValueError):
            verify_atomic(valid, collector.factor_base)
        # 19^2 - 91 = 2 * 9 * 15: exact product alone is insufficient.
        composite = AtomicRelation(
            collector.polynomial, 9, 1, ((2, 1),), large_primes=(9, 15)
        )
        with self.assertRaises(ValueError):
            verify_atomic(
                composite,
                collector.factor_base,
                large_prime_bound=100,
                large_product_bound=10000,
            )
        with self.assertRaises(ValueError):
            combine_relations((valid,), collector.factor_base)

    def test_candidate_coverage_against_independent_trial_division(self):
        base = build_factor_base(1009 * 1013, bound=40).factor_base
        polynomial = qs_polynomial(base)
        primes = {p for p in range(2, 201) if all(p % d for d in range(2, p))}
        expected, slp = set(), set()
        for position in range(-64, 65):
            remaining = abs(polynomial.value(position))
            if not remaining:
                continue
            for prime in base.primes:
                while remaining % prime == 0:
                    remaining //= prime
            if remaining == 1 or remaining in primes and remaining <= 97:
                slp.add(position)
                expected.add(position)
            elif any(
                remaining % p == 0 and remaining // p in primes for p in primes
            ):
                expected.add(position)

        for division in ("full", "roots", "bucket", "resieve"):
            for narrow in (False, True):
                collector = SieveCollector(
                    polynomial,
                    base,
                    budget=allowance(),
                    config=configuration(
                        block_width=32,
                        residual_bound=97,
                        large_prime_bound=200,
                        large_product_bound=40000,
                        candidate_bound=97 if narrow else 0,
                        max_partials=512,
                        division=division,
                    ),
                )
                result = collector.collect(-64, 65)
                self.assertEqual(result.reason, "complete")
                positions = {atom.position for atom in result.atoms}
                if narrow:
                    self.assertTrue(slp <= positions <= expected)
                else:
                    self.assertEqual(positions, expected)

    def test_split_bounds_certification_quota_and_cancelled_attempt(self):
        for residual, expected in ((9, (3, 3)), (15, (3, 5)), (77, (7, 11))):
            pair, _, _ = split_two_primes(residual, 100, allowance())
            self.assertEqual(pair, expected)
        for residual, bound in ((13, 100), (27, 100), (77, 10)):
            self.assertFalse(split_two_primes(residual, bound, allowance())[0])
        collector = worker(split_call_limit=1)
        collector.collect(0, 1)
        result = collector.collect(1, 2)
        self.assertEqual(collector._split_calls, 1)
        self.assertEqual(result.stats["split_quota_rejections"], 1)
        collector = worker()
        with patch(
            "v2.qs.large_prime_collector.split_two_primes",
            side_effect=BudgetExhaustedError,
        ):
            collector.budget.reason = "cancelled"
            result = collector.collect(0, 1)
        self.assertEqual(result.next_position, 0)
        self.assertEqual(collector._split_calls, 1)
        self.assertFalse(collector._atoms)


class StoreTests(unittest.TestCase):
    def mixed(self):
        collector = worker()
        admit(collector, -7, scalar=41)
        admit(collector, 0, (3, 3))
        polynomial = Polynomial(91, 1, 9, 1, 3)
        admit(collector, 1, polynomial=polynomial)
        admit(collector, 1, (3, 5))
        admit(collector, 2, (3, 5), polynomial=polynomial)
        return collector

    def test_mixed_order_and_owned_provenance_roundtrip(self):
        collector = self.mixed()
        payload = encoded(collector)
        self.assertEqual(payload["row_order"], [1, 0, 2])
        restored = worker()
        before = restored.budget.used
        restore_store(payload, restored, restored.budget)
        self.assertGreater(restored.budget.used, before)
        self.assertEqual(encoded(restored), payload)
        self.assertEqual(restored._rows, collector._rows)
        self.assertEqual(restored._graph.unowned, collector._graph.unowned)
        self.assertGreater(restored._scratch_peak_bytes, 0)
        self.assertLessEqual(
            restored._workspace + restored._scratch_peak_bytes,
            restored.config.memory_bytes,
        )

    def test_eviction_keeps_triangle_atoms_and_duplicate_is_ignored(self):
        collector = worker(max_partials=3)
        star, _, _ = admit(collector, -7, scalar=41)
        for x, pair in ((-5, (3, 11)), (1, (3, 5)), (-4, (5, 11))):
            admit(collector, x, pair)
        owned = set(collector._combined[0].atom_ids)
        for x, residual in ((5, 67), (8, 233), (23, 499)):
            admit(collector, x, scalar=residual)
        self.assertNotIn(star.relation_id, collector._atoms)
        self.assertTrue(owned <= set(collector._atoms))
        payload = encoded(collector)
        _, refusal, stats = admit(collector, -5, (3, 11))
        self.assertIsNone(refusal)
        self.assertEqual(stats["duplicates"], 1)
        self.assertEqual(encoded(collector), payload)
        restored = worker(max_partials=3)
        restore_store(payload, restored, restored.budget)
        self.assertEqual(encoded(restored), payload)

    def test_caps_and_fault_injection_preserve_store(self):
        for field, cap, reason in (
            ("max_atoms", 1, "atom_limit"),
            ("max_relations", 0, "relation_limit"),
            ("memory_bytes", 0, "memory_limit"),
        ):
            collector = worker()
            admit(collector, 1, (3, 5))
            saved = encoded(collector)
            cap = collector._workspace if field == "memory_bytes" else cap
            collector.config = replace(collector.config, **{field: cap})
            _, refusal, _ = admit(
                collector, 2, (3, 5), polynomial=Polynomial(91, 1, 9, 1, 3)
            )
            self.assertEqual(refusal, reason)
            self.assertEqual(encoded(collector), saved)
        collector = worker()
        admit(collector, 1, (3, 5))
        saved = encoded(collector)
        with patch(
            "v2.qs.large_prime_collector.combine_relations",
            side_effect=BudgetExhaustedError,
        ):
            with self.assertRaises(BudgetExhaustedError):
                admit(
                    collector, 2, (3, 5), polynomial=Polynomial(91, 1, 9, 1, 3)
                )
        self.assertEqual(encoded(collector), saved)

    def test_corruption_and_cancelled_restore_publish_nothing(self):
        payload = encoded(self.mixed())
        changes = [
            lambda p: p["row_order"].__setitem__(0, 0),
            lambda p: p["forest"].append(p["forest"][0]),
            lambda p: p["combined"][0].__setitem__(4, 2),
            lambda p: p["atoms"][1].__setitem__(5, [3, 5]),
            lambda p: p.__setitem__("split_calls", 1048577),
            lambda p: p["combined"][1][0].append(p["combined"][1][0][0]),
            lambda p: p["forest"].pop(),
        ]
        for change in changes:
            bad = copy.deepcopy(payload)
            change(bad)
            restored = worker()
            with self.assertRaises((ValueError, TypeError)):
                restore_store(bad, restored, restored.budget)
            self.assertFalse(restored._atoms)
            self.assertFalse(restored._graph.edges)
        for work in (0, 2, 15, 100):
            restored = worker()
            with self.assertRaises(BudgetExhaustedError):
                restore_store(payload, restored, allowance(work))
            self.assertFalse(restored._atoms)


class JobTests(unittest.TestCase):
    def test_dlp_checkpoint_charges_resume_and_preserves_slp_version(self):
        config = SIQSConfig(
            mode="qs",
            base_bound=3,
            half_width=1,
            batch_width=1,
            collector=configuration(),
        )
        job = SIQSJob(91, seed=7, config=config, budget=allowance())
        self.assertEqual(job.run(max_blocks=1).reason, "paused")
        checkpoint = json.loads(json.dumps(job.checkpoint()))
        self.assertEqual(checkpoint["version"], 4)
        restored = SIQSJob.from_checkpoint(checkpoint, budget=allowance())
        self.assertGreater(
            restored.budget.used, checkpoint["resources"]["work_used"]
        )
        self.assertGreaterEqual(
            restored.budget.prior_wall, checkpoint["resources"]["wall_used"]
        )
        self.assertGreaterEqual(
            restored.budget.prior_cpu, checkpoint["resources"]["cpu_used"]
        )
        result = restored.run()
        self.assertEqual({result.divisor, result.cofactor}, {7, 13})
        self.assertEqual(result.divisor * result.cofactor, 91)
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                mutate(
                    checkpoint,
                    lambda p: p["store"].__setitem__("split_calls", 131073),
                ),
                budget=allowance(),
            )
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                checkpoint,
                config=replace(
                    config,
                    collector=replace(
                        config.collector, large_product_bound=999999
                    ),
                ),
                allow_extension=True,
                budget=allowance(),
            )
        slp = SIQSJob(
            91,
            config=replace(config, collector=SieveConfig(block_width=1)),
            budget=allowance(),
        )
        slp.run(max_blocks=1)
        self.assertEqual(slp.checkpoint()["version"], 3)
        self.assertNotIn("large_primes", asdict(SieveConfig()))
        self.assertNotIn("large_prime_bound", asdict(SieveConfig()))

    def test_portfolio_resume_retains_nested_dlp_and_all_cofactors(self):
        siqs = SIQSConfig(
            base_bound=200,
            half_width=256,
            collector=configuration(
                block_width=256,
                residual_bound=40000,
                large_prime_bound=20000,
                large_product_bound=5120000,
                max_atoms=4096,
                max_relations=4096,
                max_partials=4096,
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
        n = 4001 * 5003
        first = factorize_bounded(n, config=config, budget=allowance(150000))
        self.assertEqual(first.result.reconstruct(), n)
        self.assertFalse(first.result.complete)
        resumed = factorize_bounded(
            n,
            config=config,
            checkpoint=first.checkpoint,
            budget=allowance(),
        )
        self.assertTrue(resumed.result.complete)
        self.assertEqual(resumed.result.reconstruct(), n)
        self.assertGreater(resumed.work_used, first.work_used)

    def test_unsupported_collectors_reject_explicitly(self):
        from v2.qs.parallel import _ExportCollector

        collector = worker()
        with self.assertRaises(ValueError):
            _ExportCollector(
                collector.polynomial,
                collector.factor_base,
                config=configuration(),
                budget=allowance(),
            )


if __name__ == "__main__":
    unittest.main()
