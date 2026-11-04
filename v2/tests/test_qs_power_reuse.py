"""Independent modular-root, traversal and bounded-transfer oracles."""

import unittest

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.qs import Polynomial, build_factor_base
from v2.qs.polynomial import polynomial_roots
from v2.qs.power_reuse import RootTransfers
from v2.qs.power_sieve import prime_power_roots
from v2.qs.sparse_queue import build
from v2.tests.test_qs import unlimited_budget
from v2.tests.test_qs_loop_reuse import enumerated_queue, queue_worker


class RootTransferTests(unittest.TestCase):
    def test_complete_shifted_pairs_match_independent_lifting_and_enumeration(
        self,
    ):
        base = build_factor_base(10403, bound=43).factor_base
        for a, values in (
            (7, (-50, -8, -1, 1, 6, 8, 15, 29)),
            (49, (8, -8, 41, 57, -90)),
        ):
            cache = RootTransfers(2**20)
            for b in values:
                polynomial = Polynomial(10403, 1, a, b)
                cache.bind(polynomial)
                for entry in base.entries:
                    roots = polynomial_roots(polynomial, base, entry)
                    prime = entry.prime
                    maximum = prime**3
                    reference_budget = unlimited_budget()
                    expected = list(
                        prime_power_roots(
                            polynomial, roots, maximum, 1, reference_budget
                        )
                    )
                    budget = unlimited_budget()
                    extra = cache.extra_work
                    actual = list(
                        prime_power_roots(
                            polynomial,
                            roots,
                            maximum,
                            1,
                            budget,
                            root_transfers=cache,
                        )
                    )
                    self.assertEqual(actual, expected)
                    self.assertEqual(
                        budget.used,
                        reference_budget.used + cache.extra_work - extra,
                    )
                    for modulus, residues, _ in actual:
                        if modulus > 2000:
                            continue
                        enumerated = {
                            x
                            for x in range(modulus)
                            if polynomial.value(x) % modulus == 0
                        }
                        self.assertEqual(set(residues), enumerated)
                self.assertLessEqual(cache.used, cache.byte_limit)
            self.assertGreater(cache.builds, 0)
            self.assertGreater(cache.replays, cache.builds)

    def test_window_pruning_does_not_cache_incomplete_pairs(self):
        base = build_factor_base(10403, bound=43).factor_base
        cache = RootTransfers(2**20)
        for b in (1, -1, 8, -8, 15):
            polynomial = Polynomial(10403, 1, 7, b)
            cache.bind(polynomial)
            for entry in base.entries:
                roots = polynomial_roots(polynomial, base, entry)
                for lo, hi in ((0, 1), (-9, 10), (-73, 40)):
                    arguments = polynomial, roots, entry.prime**4, 1
                    expected = list(
                        prime_power_roots(
                            *arguments, unlimited_budget(), lo=lo, hi=hi
                        )
                    )
                    actual = list(
                        prime_power_roots(
                            *arguments,
                            unlimited_budget(),
                            lo=lo,
                            hi=hi,
                            root_transfers=cache,
                        )
                    )
                    self.assertEqual(actual, expected)
            self.assertTrue(
                all(len(row[1]) == 2 for row in cache.values.values())
            )

    def test_target_a_and_square_provenance_changes_retire_old_pairs(self):
        cache = RootTransfers(2**20)
        for polynomial in (
            Polynomial(10403, 1, 49, 8),
            Polynomial(10403, 1, 49, 8, 7),
            Polynomial(10403, 1, 1, 1),
            Polynomial(1009 * 1013, 1, 1, 1),
            Polynomial(1009 * 1013, 3, 1, 1),
        ):
            cache.bind(polynomial)
            self.assertFalse(cache.values)
            self.assertEqual(cache.used, 0)
            # A unit root pair is unnecessary for exercising retirement.
            cache.values[25] = (polynomial.b, (1, 24), 1, 4096, 1)
            cache.used = 4096
        cache.bind(object())
        self.assertIsNone(cache.polynomial)
        self.assertFalse(cache.values)

    def test_unfunded_cache_and_refused_update_keep_exact_fallback_state(self):
        for invalid in (False, -1, 1.5, 2**40):
            with self.assertRaises((TypeError, ValueError)):
                RootTransfers(invalid)
        base = build_factor_base(10403, bound=43).factor_base
        first = Polynomial(10403, 1, 7, 1)
        entry = next(
            e
            for e in base.entries
            if e.prime > 7 and len(polynomial_roots(first, base, e).roots) == 2
        )
        roots = polynomial_roots(first, base, entry)
        empty = RootTransfers(65536)
        empty.bind(first)
        budget, reference = unlimited_budget(), unlimited_budget()
        actual = list(
            prime_power_roots(
                first, roots, entry.prime**3, 1, budget, root_transfers=empty
            )
        )
        expected = list(
            prime_power_roots(first, roots, entry.prime**3, 1, reference)
        )
        self.assertEqual((actual, budget.used), (expected, reference.used))
        self.assertEqual(
            (empty.values, empty.used, empty.extra_work), ({}, 0, 0)
        )

        cache = RootTransfers(2**20)
        cache.bind(first)
        modulus, residues, _ = expected[1]
        refused = Budget(work_limit=0, seconds=None, cpu_seconds=None)
        cache.remember(first, entry.prime, modulus, residues, refused)
        self.assertEqual(
            (cache.values, cache.used, cache.extra_work), ({}, 0, 0)
        )
        cache.remember(
            first, entry.prime, modulus, residues, unlimited_budget()
        )
        second = Polynomial(10403, 1, 7, -1)
        cache.bind(second)
        initial = polynomial_roots(second, base, entry).roots
        before = dict(cache.values), cache.used
        budget = Budget(work_limit=1, seconds=None, cpu_seconds=None)
        with self.assertRaises(BudgetExhaustedError):
            cache.shifted(second, entry.prime, modulus, initial, 1, budget)
        self.assertEqual((cache.values, cache.used), before)
        self.assertEqual(budget.used, 1)
        self.assertEqual(cache.replays, 0)

    def test_switched_queue_event_order_matches_uncached_enumeration(self):
        worker = queue_worker(7)
        worker._power_transfers = RootTransfers(2**20)
        worker._power_transfers.bind(worker.polynomial)
        lo, hi = -257, 258
        for b in (0, 1, -1, 102):
            polynomial = Polynomial(worker.polynomial.n, 1, 1, b)
            roots = tuple(
                polynomial_roots(polynomial, worker.factor_base, entry)
                for entry in worker.factor_base.entries
            )
            worker.set_polynomial(polynomial, roots)
            maximum = max(abs(polynomial.value(x)) for x in range(lo, hi))
            expected, reserved = enumerated_queue(
                worker, (lo, hi, maximum), unlimited_budget()
            )
            _, queues, actual = build(worker, (lo, hi, maximum), 2**20)
            self.assertEqual([list(queue) for queue in queues], expected)
            self.assertEqual(actual, reserved)


if __name__ == "__main__":
    unittest.main()
