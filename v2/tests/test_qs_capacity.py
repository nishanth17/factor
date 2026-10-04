"""Reachability, exact external squares and bounded continuation contracts."""

import itertools
import json
import unittest
from dataclasses import replace
from math import prod

from v2.budget import Budget, BudgetExhaustedError
from v2.qs.assignment_stream import AssignmentStream, unrank_combination
from v2.qs.capacity import capacity_report, validate_extension
from v2.qs.external_square import CoefficientExhaustedError, external_square
from v2.qs.factor_base import build_factor_base
from v2.qs.families import PolynomialFamily
from v2.qs.linear_algebra import filter_matrix, verify_dependency
from v2.qs.polynomial import Polynomial, polynomial_roots
from v2.qs.reference_collector import collect_block
from v2.qs.relations import verify_atomic, verify_combined
from v2.qs.sieve_collector import SieveCollector, SieveConfig
from v2.qs.siqs import SIQSConfig, SIQSJob


def budget(work=10**12):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


class CapacityTests(unittest.TestCase):
    def setUp(self):
        self.base = build_factor_base(
            4001 * 5003, bound=300, budget=budget()
        ).factor_base

    def stream(self, policy="nearest", count=100, **kwargs):
        options = dict(
            factor_count=3,
            family_count=count,
            pool_size=8,
            seed=17,
            policy=policy,
            budget=budget(),
            memory_bytes=16 * 2**20,
        )
        options.update(kwargs)
        return AssignmentStream(self.base, 256, **options)

    def test_combinadic_oracle_and_unique_prefix_extension(self):
        for size in range(1, 9):
            for count in range(1, size + 1):
                expected = tuple(itertools.combinations(range(size), count))
                self.assertEqual(
                    tuple(
                        unrank_combination(size, count, i)
                        for i in range(len(expected))
                    ),
                    expected,
                )
        for policy in ("nearest", "flyer"):
            short = self.stream(policy, 5)
            longer = self.stream(policy, 1000)
            self.assertEqual(short.identity, longer.identity)
            self.assertEqual(short.workspace_bytes, longer.workspace_bytes)
            self.assertEqual(
                [short[i] for i in range(5)], [longer[i] for i in range(5)]
            )
            values = [longer[i] for i in range(len(longer))]
            self.assertEqual(len(values), len(set(values)))
            self.assertEqual(len(values), len(set(map(prod, values))))
            for values in values:
                self.assertEqual(values, tuple(sorted(set(values))))
                self.assertTrue(all(p in self.base.primes for p in values))

    def test_flyer_is_exact_best_product_in_its_disjoint_domain(self):
        stream = self.stream("flyer")
        for i in range(len(stream)):
            values = stream[i]
            core = prod(p for p in values if p in stream.pool)
            flyer = next(p for p in values if p not in stream.pool)
            self.assertEqual(
                flyer,
                min(
                    stream.flyers,
                    key=lambda p: (abs(core * p - stream.target), p),
                ),
            )

    def test_larger_families_stream_gray_roots_without_history(self):
        primes = tuple(p for p in self.base.primes if p != 2)[:12]
        family = PolynomialFamily(self.base, primes, budget=budget())
        self.assertEqual(family.count, 2048)
        for _ in range(5):
            step = family.next()
            expected = tuple(
                polynomial_roots(
                    step.polynomial, self.base, e, budget=budget()
                )
                for e in self.base.entries
            )
            self.assertEqual(step.roots, expected)
        restored = PolynomialFamily.from_checkpoint(
            self.base, family.checkpoint(), budget=family.budget
        )
        self.assertEqual(restored.next().polynomial, family.next().polynomial)
        self.assertLess(len(json.dumps(family.checkpoint())), 4096)

    def test_search_quota_has_constant_metadata_and_separate_gray_limit(self):
        small = SIQSConfig(
            assignment_policy="nearest",
            family_count=3,
            factor_count=12,
            pool_size=32,
            polynomials_per_family=16,
        )
        large = replace(small, family_count=1000000)
        self.assertEqual(small.metadata_reserve, large.metadata_reserve)
        self.assertEqual(large.polynomial_limit, 16000000)
        report = capacity_report(self.base, large)
        self.assertEqual(report["base_cardinality"], len(self.base.entries))
        self.assertEqual(report["a_maximum"], prod(self.base.primes[-12:]))

    def test_external_square_sieve_matches_independent_full_division(self):
        polynomial, cursor, certainty = external_square(
            self.base, 8, budget=budget(), cursor=307
        )
        self.assertGreater(polynomial.square_coefficient, self.base.bound)
        self.assertEqual(cursor, polynomial.square_coefficient + 4)
        self.assertEqual(certainty, "proven_prime")
        reference = collect_block(
            polynomial,
            self.base,
            -64,
            65,
            residual_bound=10000,
            budget=budget(),
            memory_bytes=64 * 2**20,
            max_relations=1024,
        )
        config = SieveConfig(
            residual_bound=10000,
            max_atoms=4096,
            max_partials=4096,
            max_relations=4096,
            memory_bytes=64 * 2**20,
            score_policy="conservative",
        )
        collector = SieveCollector(
            polynomial, self.base, config=config, budget=budget()
        )
        result = collector.collect(-64, 65)
        expected = {a.position: a for a in reference.relations}
        self.assertTrue(expected)
        self.assertEqual({a.position: a for a in result.atoms}, expected)
        for atom in result.atoms:
            self.assertTrue(
                verify_atomic(
                    atom, self.base, residual_bound=10000, budget=budget()
                )
            )
            self.assertEqual(
                atom.u**2 - self.base.n_prime,
                atom.sign
                * prod(p**e for p, e in atom.exponents)
                * atom.residual
                * atom.square_correction**2,
            )
        for combined in result.combined_relations:
            self.assertTrue(
                verify_combined(
                    combined,
                    self.base,
                    collector._atoms,
                    budget=budget(),
                    memory_bytes=64 * 2**20,
                )
            )
        with self.assertRaises(ValueError):
            Polynomial(self.base.n, 1, polynomial.a, polynomial.b, 2)
        with self.assertRaises(BudgetExhaustedError):
            external_square(self.base, 8, budget=budget(0))
        with self.assertRaises(CoefficientExhaustedError):
            external_square(
                self.base, 8, budget=budget(), cursor=319, trials=1
            )

    def test_external_mpqs_factor_and_checkpoint_store(self):
        cfg = SIQSConfig(
            mode="mpqs",
            external_coefficients=True,
            base_bound=200,
            half_width=256,
            family_count=1000,
            max_stalled=1000,
            collector=SieveConfig(
                residual_bound=10000,
                max_atoms=4096,
                max_relations=2048,
                max_partials=512,
            ),
        )
        job = SIQSJob(self.base.n, config=cfg, budget=budget())
        job.run(max_blocks=1)
        restored = SIQSJob.from_checkpoint(job.checkpoint(), budget=budget())
        self.assertEqual(restored.coefficient_cursor, job.coefficient_cursor)
        self.assertEqual(
            restored.engine.collector._atoms, job.engine.collector._atoms
        )
        result = restored.run()
        self.assertEqual(result.reason, "factor_found")
        self.assertEqual(result.divisor * result.cofactor, self.base.n)

    def test_extend_exhausted_stream_retains_store_and_spent_resources(self):
        cfg = SIQSConfig(
            base_bound=100,
            half_width=8,
            assignment_policy="nearest",
            factor_count=2,
            pool_size=8,
            family_count=1,
            polynomials_per_family=1,
            max_stalled=1000,
            max_trivial=100000,
        )
        job = SIQSJob(1000003, config=cfg, budget=budget())
        self.assertEqual(job.run().reason, "families_exhausted")
        checkpoint = job.checkpoint()
        extended = replace(cfg, family_count=4)
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                checkpoint, budget=budget(), config=extended
            )
        resumed = SIQSJob.from_checkpoint(
            checkpoint, budget=budget(), config=extended, allow_extension=True
        )
        self.assertEqual(resumed.family_index, job.family_index)
        self.assertEqual(
            resumed.engine.collector._atoms, job.engine.collector._atoms
        )
        self.assertGreater(resumed.budget.used, job.budget.used)
        self.assertEqual(resumed.run().reason, "families_exhausted")
        full = SIQSJob(job.n, config=extended, budget=budget())
        self.assertEqual(full.run().cofactor, job.n)
        self.assertEqual(
            resumed.engine.collector._atoms, full.engine.collector._atoms
        )
        self.assertEqual(resumed.stats["polynomials"], 4)
        self.assertEqual(resumed.seen, set())
        for invalid in (
            replace(extended, base_bound=110),
            replace(extended, pool_size=9),
            replace(extended, polynomials_per_family=2),
            replace(extended, checkpoint_bytes=2**21),
            replace(extended, family_count=1),
        ):
            with self.assertRaises(ValueError):
                validate_extension(extended, invalid)

    def test_extension_reserves_growth_of_verified_row_cache(self):
        original = SIQSConfig(
            assignment_policy="nearest",
            collector=SieveConfig(max_relations=64),
        )
        enlarged = replace(
            original,
            collector=replace(original.collector, max_relations=128),
        )
        with self.assertRaisesRegex(ValueError, "checkpoint/cache"):
            validate_extension(original, enlarged)
        enlarged = replace(
            enlarged, memory_bytes=original.memory_bytes + 64 * 2048
        )
        validate_extension(original, enlarged)
        report = capacity_report(self.base, enlarged)
        self.assertEqual(report["preparation_cache_reserve_bytes"], 128 * 2048)
        saturated = replace(
            enlarged,
            collector=replace(enlarged.collector, max_relations=2048),
            memory_bytes=enlarged.memory_bytes + 2 * 2**20,
        )
        validate_extension(enlarged, saturated)
        validate_extension(
            saturated,
            replace(
                saturated,
                collector=replace(saturated.collector, max_relations=4096),
            ),
        )

    def test_sparse_incidence_work_preserves_verified_dependencies(self):
        rows = tuple(1 << (i % 500) for i in range(2000))
        allowance = budget(2_000_000)
        filtered = filter_matrix(
            rows, budget=allowance, memory_bytes=32 * 2**20, weight_two=True
        )
        self.assertLess(allowance.used, 2_000_000)
        for mask in filtered.zero_dependencies:
            self.assertTrue(verify_dependency(mask, rows))
        self.assertEqual(filtered.original_rows, rows)
        wide = (1 << 99999,) * 32
        compact = filter_matrix(
            wide,
            budget=budget(2_000_000),
            memory_bytes=8 * 2**20,
            weight_two=True,
        )
        self.assertEqual(compact.stats["working_columns"], 1)
        for mask in compact.zero_dependencies:
            self.assertTrue(verify_dependency(mask, wide))

    def test_99_digit_reachability_and_probable_coefficient(
        self,
    ):
        n = int(
            "333125016877336815730234617089075553677225251597802911492140000"
            "164835979295625386188653593061163021"
        )
        base = build_factor_base(
            n, bound=10000, budget=budget(), memory_bytes=32 * 2**20
        ).factor_base
        cfg = SIQSConfig(
            base_bound=10000,
            half_width=8192,
            assignment_policy="flyer",
            factor_count=13,
            pool_size=52,
            polynomials_per_family=16,
        )
        report = capacity_report(base, cfg)
        self.assertTrue(report["target_in_product_envelope"])
        polynomial, _, certainty = external_square(base, 8192, budget=budget())
        self.assertEqual(certainty, "probable_prime")
        self.assertGreater(polynomial.square_coefficient, 2**64)
        target = report["a_target"]
        self.assertLess(abs(polynomial.a - target) * 1000, target)
        for x in (-100, 0, 731):
            self.assertEqual(
                polynomial.u_value(x) ** 2 - n,
                polynomial.a * polynomial.value(x),
            )

    def test_checkpoint_size_independent_of_unvisited_search(self):
        cfg = SIQSConfig(
            assignment_policy="nearest",
            base_bound=200,
            family_count=100,
            half_width=256,
        )
        lengths = []
        for count in (100, 1000000):
            job = SIQSJob(
                self.base.n,
                config=replace(cfg, family_count=count),
                budget=budget(),
            )
            job.run(max_blocks=1)
            lengths.append(len(job.checkpoint()["blob"]))
            self.assertEqual(job.seen, set())
        self.assertLess(abs(lengths[0] - lengths[1]), 128)

    def test_external_candidate_allowance_can_extend_after_refusal(self):
        cfg = SIQSConfig(
            mode="mpqs",
            external_coefficients=True,
            base_bound=100,
            half_width=256,
            coefficient_trials=1,
        )
        job = SIQSJob(1000003 * 1000033, config=cfg, budget=budget())
        self.assertEqual(job.run().reason, "coefficient_limit")
        used = job.budget.used
        resumed = SIQSJob.from_checkpoint(
            job.checkpoint(),
            budget=budget(),
            allow_extension=True,
            config=replace(cfg, coefficient_trials=4096),
        )
        self.assertGreaterEqual(resumed.budget.used, used)
        self.assertIn(
            resumed.run(max_blocks=1).reason, ("paused", "factor_found")
        )
        self.assertGreater(resumed.coefficient_cursor, 0)
