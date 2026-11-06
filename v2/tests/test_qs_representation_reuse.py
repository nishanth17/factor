"""Mutation, arithmetic and publication oracles for retained reuse."""

import unittest
from collections import Counter
from dataclasses import replace
from unittest.mock import patch

from v2.execution.budget import BudgetExhaustedError
from v2.execution.work_budget import PollingBudget
from v2.qs import Polynomial, SieveCollector
from v2.qs.cold_factors import RecomputedFactors
from v2.qs.factor_witness import factor_snapshot
from v2.qs.indexed_relations import FactorRecord, IndexedFactors
from v2.qs.polynomial import PolynomialRoots
from v2.qs.power_sieve import HenselInverses, prime_power_roots
from v2.qs.yield_support import payload_snapshot
from v2.tests.test_c9a_bounded import worker
from v2.tests.test_c9a_extended import atom_at
from v2.tests.test_qs_dlp import allowance


class FactorWitnessTests(unittest.TestCase):
    def test_packed_and_lazy_witnesses_do_not_decode_and_detect_slab_damage(
        self,
    ):
        engine = worker(recompute_prime_limit=0)
        for position in (1, -21):
            engine._admit(atom_at(engine, position), Counter())
        atom, row = next(iter(engine._atoms.values())), engine._rows[0]
        with patch.object(
            FactorRecord, "__iter__", side_effect=AssertionError
        ):
            for value in (atom, row):
                snapshot = payload_snapshot(value)
                engine._checked._known(value, engine.budget)
                self.assertEqual(snapshot, payload_snapshot(value))
        factors = atom.exponents
        factors._block[factors._offset] ^= 1
        with self.assertRaises(ValueError):
            engine._checked.atom(atom, engine.budget)

    def test_recomputed_witness_binds_its_own_provenance_without_iteration(
        self,
    ):
        engine = worker(recompute_prime_limit=13)
        source = atom_at(engine, 1)
        engine._admit(source, Counter())
        atom = engine._atoms[source.relation_id]
        self.assertIs(type(atom.exponents), RecomputedFactors)
        with patch.object(
            RecomputedFactors, "__iter__", side_effect=AssertionError
        ):
            engine._checked.atom(atom, engine.budget)
            previous = factor_snapshot(atom.exponents)
            polynomial = replace(atom.polynomial, b=atom.polynomial.b + 1)
            object.__setattr__(atom.exponents, "polynomial", polynomial)
            self.assertNotEqual(previous, factor_snapshot(atom.exponents))
            with self.assertRaises(ValueError):
                engine._checked.atom(atom, engine.budget)

    def test_indexed_input_is_deeply_immutable_and_bypass_is_detected(self):
        values = [[1, 2]]
        factors = IndexedFactors(values, (2,))
        snapshot = factor_snapshot(factors)
        values[0][1] = 5
        self.assertEqual(tuple(factors), ((2, 2),))
        self.assertEqual(snapshot, factor_snapshot(factors))
        object.__setattr__(factors, "indices", ([1, 2],))
        with self.assertRaises(ValueError):
            factor_snapshot(factors)

    def test_raw_proof_retirement_preserves_packed_canonical_representative(
        self,
    ):
        engine = worker(
            recompute_prime_limit=0, duplicate_policy="congruence-v2"
        )
        raw = atom_at(engine, 1)
        engine._checked.atom(raw, engine.budget)
        self.assertIn(id(raw), engine._checked.records)
        engine._admit(raw, Counter())
        packed = engine._atoms[raw.relation_id]
        self.assertNotIn(id(raw), engine._checked.records)
        self.assertIn(id(packed), engine._checked.records)
        key, divisor = engine._checked.congruence(packed, engine.budget)
        self.assertIsNone(divisor)
        self.assertIsInstance(key[2], bytes)
        self.assertIs(engine._checked.canonical[key], packed)
        self.assertTrue(
            engine._checked.duplicate(key, engine._atoms, engine.budget)
        )
        self.assertEqual(
            engine._checked.used,
            sum(value[3] for value in engine._checked.records.values())
            + 512 * (len(engine._checked.primes) + len(engine._checked.units)),
        )


class HenselReuseTests(unittest.TestCase):
    def test_unfunded_helpers_fall_back_to_identical_scores_and_wide_masks(
        self,
    ):
        cached = worker(score_policy="powers", division="bucket")
        required = cached._workspace - cached._reuse_workspace
        plain = SieveCollector(
            cached.polynomial,
            cached.factor_base,
            config=replace(cached.config, memory_bytes=required),
            budget=allowance(),
        )
        self.assertIsNone(plain._shared_lift_inverses)
        self.assertEqual(plain._hit_masks, ())
        self.assertLessEqual(plain._workspace, required)
        for engine in (plain, cached):
            engine._sieve_direct(1, 17, Counter())
            self.assertEqual(engine._recovery_mask(1024), 1 << 1024)
        self.assertEqual(plain._scores, cached._scores)
        self.assertEqual(plain._hits, cached._hits)

    def test_cross_polynomial_reuse_matches_independent_modular_enumeration(
        self,
    ):
        shared = HenselInverses(8)
        inversions = 0
        from v2.qs import power_sieve

        actual_pow = power_sieve.pow

        def counted(*arguments):
            nonlocal inversions
            inversions += 1
            return actual_pow(*arguments)

        with patch.object(power_sieve, "pow", side_effect=counted):
            for a, b in ((1, 0), (1, 1), (7, 1), (49, 8)):
                polynomial = Polynomial(10403, 1, a, b)
                for prime in (2, 3, 5, 7, 11):
                    roots = PolynomialRoots(
                        prime,
                        tuple(
                            x
                            for x in range(prime)
                            if polynomial.value(x) % prime == 0
                        ),
                        False,
                    )
                    marks = list(
                        prime_power_roots(
                            polynomial,
                            roots,
                            512,
                            1,
                            allowance(),
                            inverses={},
                            shared_inverses=shared,
                        )
                    )
                    for modulus, residues, _ in marks:
                        expected = {
                            x
                            for x in range(modulus)
                            if polynomial.value(x) % modulus == 0
                        }
                        self.assertEqual(set(residues), expected)
        self.assertLessEqual(inversions, 4)
        self.assertLessEqual(len(shared.values), 5)

    def test_singular_parity_multiplier_and_cache_refusal_remain_exact(self):
        for n, multiplier in ((9, 1), (101, 3), (10403, 1)):
            polynomial = Polynomial(n, multiplier, 1, 0)
            for prime in (2, 3, 5, 7):
                roots = PolynomialRoots(
                    prime,
                    tuple(
                        x
                        for x in range(prime)
                        if polynomial.value(x) % prime == 0
                    ),
                    False,
                )
                options = dict(lo=-15, hi=18, inverses={})
                expected = list(
                    prime_power_roots(
                        polynomial,
                        roots,
                        4096,
                        1,
                        allowance(),
                        **options,
                    )
                )
                for cap in (0, 2):
                    shared = HenselInverses(cap)
                    observed = list(
                        prime_power_roots(
                            polynomial,
                            roots,
                            4096,
                            1,
                            allowance(),
                            lo=-15,
                            hi=18,
                            inverses={},
                            shared_inverses=shared,
                        )
                    )
                    self.assertEqual(observed, expected)
                    self.assertLessEqual(len(shared.values), cap)


class RejectSpanTests(unittest.TestCase):
    def test_spans_preserve_exact_refusal_cursor_and_reporting(self):
        for limit in (0, 1, 9, 18, 73, 200):
            outputs = []
            for interval in (None, 4, 64):
                engine = worker()
                ledger = allowance(limit)
                engine.budget = (
                    ledger
                    if interval is None
                    else PollingBudget(ledger, interval=interval)
                )

                def sieve(lo, hi, stats):
                    for i in range(hi - lo):
                        engine._scores[i] = int(lo + i in (1, 5, 15, 70))
                    stats["blocks"] += 1
                    return 1

                def divide(position, offset, stats):
                    engine.budget.consume()
                    return None, None

                with (
                    patch.object(engine, "_sieve", side_effect=sieve),
                    patch.object(engine, "_divide", side_effect=divide),
                ):
                    result = engine.collect(0, 100)
                outputs.append(
                    (result, ledger.used, engine._reporting.state())
                )
            self.assertEqual(outputs[0], outputs[1])
            self.assertEqual(outputs[0], outputs[2])

    def test_span_cancellation_and_overdraw_are_atomic(self):
        calls = []
        ledger = allowance(65)
        ledger.cancelled = lambda: calls.append(ledger.used) or len(calls) > 1
        budget = PollingBudget(ledger)
        budget.consume_span(64)
        with self.assertRaises(BudgetExhaustedError):
            budget.consume_span(1)
        self.assertEqual((ledger.used, calls), (64, [0, 64]))
        ledger.cancelled = None
        with self.assertRaises(BudgetExhaustedError):
            budget.consume_span(2)
        self.assertEqual(ledger.used, 64)
        with self.assertRaises(ValueError):
            budget.consume_span(65)
