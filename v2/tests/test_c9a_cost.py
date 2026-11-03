"""Exact early exits, cache lifetimes and bounded cost regressions."""

import copy
import unittest
from collections import Counter, OrderedDict
from dataclasses import asdict
from types import SimpleNamespace
from unittest.mock import patch

from v2 import qs
from v2.benchmarks.qs.c9a import extended
from v2.common import utils
from v2.execution.budget import BudgetExhaustedError
from v2.qs.cofactor_policy import CofactorBatch, split_cofactor, split_residual
from v2.qs.extended_checkpoint import pack_store, restore_store
from v2.qs.extended_collector import atom_reserve, recover
from v2.qs.indexed_relations import IndexedFactors
from v2.qs.linear_algebra import (
    DependencySolver,
    filter_matrix,
    verify_dependency,
)
from v2.qs.storage import RelationPublicationError
from v2.qs.structured_matrix import panel_reduce
from v2.qs.yield_memory import reserve_batch
from v2.tests.test_c9a_bounded import bounded, worker
from v2.tests.test_c9a_extended import atom_at
from v2.tests.test_qs_dlp import allowance
from v2.tests.test_qs_pipeline import dense_kernel, span


class CountingPivots(dict):
    scans = 0

    def __iter__(self):
        self.scans += 1
        return super().__iter__()


class ChunkAllowance:
    def __init__(self, limit=10000):
        self.limit, self.used, self.polls = limit, 0, []

    def consume(self, amount=1):
        if self.used + amount > self.limit:
            raise BudgetExhaustedError()
        self.polls.append(amount)
        self.used += amount


def lookahead_worker(limit=10000):
    owner = SimpleNamespace(
        budget=ChunkAllowance(limit),
        config=SimpleNamespace(
            smooth_batch_size=0, division="bucket", large_product_bound=100
        ),
        _scores=[int(i in (10, 70)) for i in range(200)],
        _piecewise_passes=lambda offset: True,
    )
    owner._divide_double = (
        lambda worker, x, value, residual, factors, stats, hint=None: (
            x,
            None,
        )
    )
    owner._divide = lambda x, offset, stats: owner._divide_double(
        owner, x, x + 1, 1, (), stats
    )
    return owner


class CostRegressionTests(unittest.TestCase):
    def test_proven_prime_hole_rejects_before_hints_or_primality(self):
        base = qs.build_factor_base(1009 * 1013, bound=100).factor_base
        engine = qs.SieveCollector(
            qs.qs_polynomial(base),
            base,
            config=bounded(residual_bound=10),
            budget=allowance(),
        )
        found = False
        for position in range(1, 20):
            value = engine.polynomial.value(position)
            residual, factors = abs(value), []
            for prime in base.primes:
                exponent = 0
                while residual % prime == 0:
                    residual //= prime
                    exponent += 1
                if exponent:
                    factors.append((prime, exponent))
            if not 10 < residual < engine._residual_floor**2:
                continue
            with (
                patch.object(
                    engine._checked, "prime", side_effect=AssertionError
                ),
                patch.object(
                    engine._endpoints, "divisor", side_effect=AssertionError
                ),
            ):
                self.assertEqual(
                    recover(
                        engine,
                        position,
                        value,
                        residual,
                        tuple(factors),
                        Counter(),
                    ),
                    (None, None),
                )
            found = True
            break
        self.assertTrue(found)

    def test_two_prime_floor_keeps_the_last_useful_rho_attempt(self):
        checked = SimpleNamespace(
            prime=lambda p, budget: utils.classify_prime(p)
            is utils.Primality.PROVEN
        )
        with (
            patch("v2.qs.cofactor_policy.squfof", return_value=None),
            patch(
                "v2.qs.cofactor_policy.factorize_rho",
                side_effect=(None, None, 11),
            ) as rho,
        ):
            factors, _, _ = split_residual(
                11 * 17,
                100,
                allowance(),
                checked,
                "rho-first-v3",
                arity=2,
                dispatch=False,
            )
        self.assertEqual(factors, (11, 17))
        self.assertEqual(
            [call.kwargs["seed"] for call in rho.call_args_list], [7, 11, 23]
        )

    def test_unmatched_staging_skips_unused_congruence_and_keeps_fifo(self):
        engine = worker(
            opaque_cache_entries=1, duplicate_policy="congruence-v2"
        )
        original = engine._opaque_staging.atoms
        with patch.object(
            engine._checked,
            "congruence",
            side_effect=AssertionError("unmatched staging needs no exact key"),
        ):
            engine._admit(atom_at(engine, 1, opaque=True), Counter())
            engine._admit(atom_at(engine, 2, opaque=True), Counter())
        self.assertIs(engine._opaque_staging.atoms, original)
        self.assertEqual(list(original), [53])

    def test_staging_eviction_refusal_preserves_every_original_entry(self):
        engine = worker(opaque_cache_entries=1)
        engine._admit(atom_at(engine, 1, opaque=True), Counter())
        staging = engine._opaque_staging
        prior, used = dict(staging.atoms), staging.used
        atom = atom_at(engine, 2, opaque=True)
        with self.assertRaises(BudgetExhaustedError):
            staging.stage(
                atom,
                ChunkAllowance(len(atom.exponents) + 1),
                Counter(),
                engine,
            )
        self.assertEqual(dict(staging.atoms), prior)
        self.assertEqual(staging.used, used)
        self.assertFalse(engine._publication_failed)

    def test_failed_staging_publication_cannot_be_reused_or_serialized(self):
        class FailingMap(OrderedDict):
            def __setitem__(self, key, value):
                super().__setitem__(key, value)
                raise MemoryError("injected after insertion")

        engine = worker()
        engine._opaque_staging.atoms = FailingMap()
        with self.assertRaises(RelationPublicationError):
            engine._admit(atom_at(engine, 1, opaque=True), Counter())
        self.assertTrue(engine._publication_failed)
        with self.assertRaises(RelationPublicationError):
            pack_store(engine)
        with self.assertRaises(RelationPublicationError):
            engine._admit(atom_at(engine, 2, opaque=True), Counter())

    def test_proven_floor_avoids_impossible_three_prime_recovery(self):
        engine = worker(batch_policy="off", split_policy="rho-first-v3")
        engine._residual_floor = 11
        engine._checked = SimpleNamespace(
            primes={},
            prime=lambda value, budget: False,
            atom=lambda atom, budget: None,
            unit=lambda value, budget: 1,
        )
        for residual, primes, arity in (
            (11 * 13, (11, 13), 2),
            (11**3, (11, 11, 11), 3),
        ):
            with patch(
                "v2.qs.cofactor_policy.split_residual",
                return_value=(primes, "double", 0),
            ) as split:
                atom, divisor = recover(
                    engine, 1, residual, residual, (), Counter()
                )
            self.assertEqual(split.call_args.kwargs["arity"], arity)
            self.assertEqual(atom.residual_product, residual)
            self.assertIsNone(divisor)

    def test_narrower_split_memo_restores_and_enforces_its_arity(self):
        engine = worker()
        key = (187, 1000, engine.config.split_policy, 2, True)
        engine._split_memo.remember(key, ((11, 17), "double"))
        saved = pack_store(engine)
        clone = worker()
        restore_store(saved, clone, clone.budget)
        self.assertEqual(pack_store(clone), saved)

        forged = copy.deepcopy(saved)
        forged["memo"][0][0][0] = 5 * 11 * 17
        forged["memo"][0][1] = [5, 11, 17]
        with self.assertRaises(ValueError):
            clone = worker()
            restore_store(forged, clone, clone.budget)

    def test_disabled_recomputation_does_not_enumerate_indexed_factors(self):
        class CountedFactors(IndexedFactors):
            scans = 0

            def __iter__(self):
                type(self).scans += 1
                return super().__iter__()

        engine = worker(recompute_prime_limit=0)
        source = atom_at(engine, 1)
        factors = CountedFactors(
            tuple(
                (engine.factor_base._columns[p], e)
                for p, e in source.exponents
            ),
            engine.factor_base.primes,
        )
        atom = qs.AtomicRelation(engine.polynomial, 1, 1, factors)
        CountedFactors.scans = 0
        packed, _, _ = engine._arena.plan(atom, None, engine.budget)
        self.assertEqual(CountedFactors.scans, 0)
        self.assertEqual(tuple(packed.exponents), tuple(factors))

    def test_shared_context_does_not_recompute_private_ownership(self):
        engine = worker()
        atom = atom_at(engine, 1)
        engine._admit(atom, Counter())
        with patch(
            "v2.qs.yield_memory.polynomial_reserve",
            side_effect=AssertionError("shared context is already owned"),
        ):
            self.assertGreater(
                atom_reserve(
                    atom,
                    engine.config,
                    engine.factor_base,
                    arena=engine._arena,
                ),
                0,
            )

    def test_batch_bound_reuses_only_successfully_owned_high_water(self):
        engine = worker()
        with patch.object(engine, "_bounds", return_value=(1, 2**20)):
            reserve_batch(engine, 0, 1)
        used, owned = engine.budget.used, engine._batch_reserve
        with patch.object(engine, "_bounds", return_value=(1, 2**10)):
            reserve_batch(engine, 0, 1)
        self.assertEqual(
            (engine.budget.used, engine._batch_reserve), (used, owned)
        )

        bits = engine._batch_bits
        original = engine.config
        engine.config = SimpleNamespace(batch_size=64, memory_limit=0)
        with patch.object(engine, "_bounds", return_value=(1, 2**40)):
            with self.assertRaises(MemoryError):
                reserve_batch(engine, 0, 1)
        self.assertEqual(engine._batch_bits, bits)
        self.assertEqual(engine._batch_reserve, owned)
        engine.config = original
        with patch.object(engine, "_bounds", return_value=(1, 2**40)):
            reserve_batch(engine, 0, 1)
        self.assertGreater(engine._batch_bits, bits)
        self.assertGreater(engine._batch_reserve, owned)

    def test_cost_profile_changes_only_the_measured_expensive_choices(self):
        old = extended.configuration(
            60, admission_base_bound=103814, profile="bounded-report"
        )
        new = extended.configuration(
            60, admission_base_bound=103814, profile="cost-report"
        )
        expected = asdict(old)
        expected["collector"].update(
            solver_panel_bits=0,
            split_policy="rho-first-v3",
            batch_policy="off",
        )
        self.assertEqual(asdict(new), expected)

    def test_rho_first_short_circuits_squfof_and_keeps_its_fallback(self):
        checked = SimpleNamespace(
            prime=lambda p, budget: (
                utils.classify_prime(p) is utils.Primality.PROVEN
            )
        )
        events = []

        def rho(*args, **kwargs):
            events.append("rho")
            return None

        def squfof(*args, **kwargs):
            events.append("squfof")
            return 1009

        with (
            patch("v2.qs.cofactor_policy.factorize_rho", side_effect=rho),
            patch("v2.qs.cofactor_policy.squfof", side_effect=squfof),
        ):
            pair, _, _ = split_cofactor(
                1009 * 1000003, 2000000, allowance(), checked, "rho-first-v3"
            )
        self.assertEqual(pair, (1009, 1000003))
        self.assertEqual(events, ["rho", "rho", "squfof"])

        with (
            patch("v2.qs.cofactor_policy.factorize_rho", return_value=17),
            patch(
                "v2.qs.cofactor_policy.squfof",
                side_effect=AssertionError("unneeded fallback"),
            ),
        ):
            pair, _, _ = split_cofactor(
                17 * 19, 100, allowance(), checked, "rho-first-v3"
            )
        self.assertEqual(pair, (17, 19))

    def test_three_prime_recovery_reuses_a_paid_proper_split(self):
        checked = SimpleNamespace(
            prime=lambda p, budget: (
                utils.classify_prime(p) is utils.Primality.PROVEN
            )
        )
        inputs = []

        def rho(number, **kwargs):
            inputs.append(number)
            return {1353: 3, 451: 11}[number]

        with patch("v2.qs.cofactor_policy.factorize_rho", side_effect=rho):
            primes, _, _ = split_residual(
                1353, 100, allowance(), checked, "rho-first-v3", dispatch=False
            )
        self.assertEqual(primes, (3, 11, 41))
        self.assertEqual(inputs, [1353, 451])

    def test_panel_presence_avoids_scanning_or_sorting_pivots(self):
        pivots = CountingPivots(
            {1: (1, 1), 2: (3, 2), 8: (10, 4), 32: (34, 8), 128: (130, 16)}
        )
        solver = SimpleNamespace(
            pivots=pivots,
            _panel_pivots=171,
            pivot="highest",
            budget=allowance(),
            matrix=SimpleNamespace(
                stats={"panel_bits": 4, "working_columns": 8},
                original_rows=tuple(range(8)),
            ),
        )
        expected = panel_reduce(solver, 171, 7)
        table = solver._panel_cache
        for _ in range(100):
            self.assertEqual(panel_reduce(solver, 171, 7), expected)
        self.assertEqual(pivots.scans, 0)

        pivots[512] = (512, 32)
        solver._panel_pivots |= 512
        self.assertEqual(panel_reduce(solver, 171, 7), expected)
        self.assertEqual(pivots.scans, 0)
        self.assertIs(solver._panel_cache, table)

    def test_noncontiguous_panels_match_independent_dense_kernels(self):
        rows = (3, 9, 10, 129, 136, 130, 4097, 4104, 4098, 4224)
        expected = span(dense_kernel(rows))
        for pivot in ("highest", "lowest"):
            matrix = filter_matrix(rows, panel_bits=3, budget=allowance())
            solver = DependencySolver(matrix, pivot=pivot, budget=allowance())
            dependencies = solver.run()
            self.assertEqual(solver._panel_pivots, sum(solver.pivots))
            self.assertEqual(span(dependencies), expected)
            for mask in dependencies:
                verify_dependency(mask, rows)

    def test_lookahead_prepays_a_chunk_without_reordering_candidates(self):
        owner, batch = lookahead_worker(), CofactorBatch(2)
        self.assertEqual(
            batch.recover(owner, 10, 10, 200, 1, Counter()), (10, None)
        )
        self.assertEqual(owner.budget.polls, [64])
        self.assertEqual(
            batch.recover(owner, 70, 70, 200, 1, Counter()), (70, None)
        )
        self.assertEqual(owner.budget.polls, [64])

    def test_refused_lookahead_publishes_no_records_or_divider_change(self):
        owner, batch = lookahead_worker(limit=63), CofactorBatch(2)
        divider = owner._divide_double
        with self.assertRaises(BudgetExhaustedError):
            batch.recover(owner, 10, 10, 200, 1, Counter())
        self.assertFalse(batch.records)
        self.assertIs(owner._divide_double, divider)
        self.assertEqual(owner.budget.used, 0)


if __name__ == "__main__":
    unittest.main()
