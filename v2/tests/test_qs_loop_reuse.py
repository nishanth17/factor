"""Exact score, event-order and reservation oracles for loop reuse."""

import unittest
from collections import Counter
from dataclasses import replace
from unittest.mock import patch

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.work_budget import PollingBudget
from v2.qs import Polynomial, SieveCollector, SieveConfig, build_factor_base
from v2.qs.families import verify_polynomial_roots
from v2.qs.polynomial import PolynomialRoots, polynomial_roots
from v2.qs.power_sieve import prime_power_roots
from v2.qs.sparse_queue import ENTRY_BYTES, build
from v2.tests.test_qs import unlimited_budget


def queue_worker(width, marking="sparse"):
    base = build_factor_base(1009 * 1013, bound=100).factor_base
    polynomial = Polynomial(base.n, 1, 1, 0)
    worker = SieveCollector(
        polynomial,
        base,
        budget=unlimited_budget(),
        config=SieveConfig(
            division="bucket",
            marking=marking,
            score_policy="powers",
            block_width=width,
            sparse_queue_bytes=2**20,
            memory_bytes=64 * 2**20,
        ),
    )
    # This oracle isolates event grouping charges from optional higher-root
    # cache construction, whose extra inverse work is tested independently.
    worker._power_transfers = None
    return worker


def enumerated_queue(worker, window, budget):
    """Enumerate positions independently of quotient/remainder stepping."""
    lo, hi, maximum = window
    roots, width, split = (
        worker._roots,
        worker.config.block_width,
        worker._sparse_split,
    )
    count = (hi - lo + width - 1) // width
    queues = [[] for _ in range(count)]
    reserved = 65536 + 1024 * count + 256 + 16 * split
    budget.consume(len(roots) + count + 1)
    for index in range(split, len(roots)):
        root = roots[index]
        for modulus, residues, weight in prime_power_roots(
            worker.polynomial,
            root,
            maximum,
            worker._logs[index],
            budget,
            lo=lo,
            hi=hi,
            inverses={},
        ):
            for residue in residues:
                for position in range(
                    lo + (residue - lo) % modulus, hi, modulus
                ):
                    budget.consume(4)
                    block, offset = divmod(position - lo, width)
                    queues[block].extend(
                        (
                            offset,
                            index,
                            (weight << 1) | (modulus == root.prime),
                        )
                    )
                    reserved += ENTRY_BYTES
    budget.consume(count + 1)
    return queues, reserved


class QueueLoopTests(unittest.TestCase):
    def test_poll_mutation_of_shared_limit_keeps_refused_event_exact(self):
        observations = []
        for grouped in (False, True):
            worker = queue_worker(7)
            lo, hi = -63, 66
            maximum = max(
                abs(worker.polynomial.value(x)) for x in range(lo, hi)
            )
            calls = []
            ledger = unlimited_budget()

            def shorten_limit(shared):
                calls.append(shared.used)
                # Change the limit during the first event poll after setup.
                # A group must not reserve multiple events before seeing it.
                if len(calls) == 3:
                    shared.work_limit = shared.used + 4

            budget = PollingBudget(ledger, interval=7, poll=shorten_limit)
            try:
                if grouped:
                    worker.budget = budget
                    build(worker, (lo, hi, maximum), 2**20)
                else:
                    enumerated_queue(worker, (lo, hi, maximum), budget)
            except BudgetExhaustedError:
                refused = True
            else:
                refused = False

            observations.append((ledger.used, ledger.reason, calls, refused))

        self.assertEqual(observations[0], observations[1])
        self.assertTrue(observations[0][-1])

    def test_flat_event_order_and_complete_work_match_position_enumeration(
        self,
    ):
        for width in (4, 7, 16):
            worker = queue_worker(width)
            lo, hi = -257, 258
            maximum = max(
                abs(worker.polynomial.value(x)) for x in range(lo, hi)
            )
            window = lo, hi, maximum
            reference_budget = unlimited_budget()
            expected, reserved = enumerated_queue(
                worker, window, reference_budget
            )
            ledger = unlimited_budget()
            worker.budget = PollingBudget(ledger, interval=7)
            split, queues, actual = build(worker, window, 2**20)
            self.assertEqual(split, worker._sparse_split)
            self.assertEqual([list(queue) for queue in queues], expected)
            self.assertEqual(actual, reserved)
            self.assertEqual(ledger.used, reference_budget.used)

    def test_grouped_refusal_and_cancel_poll_at_the_same_logical_boundaries(
        self,
    ):
        for work_limit in (0, 150, 300, 2000, 100000):
            for cancel_after in (None, 2, 4):
                observations = []
                for grouped in (False, True):
                    worker = queue_worker(7)
                    lo, hi = -63, 66
                    maximum = max(
                        abs(worker.polynomial.value(x)) for x in range(lo, hi)
                    )
                    calls = []
                    ledger = Budget(
                        work_limit=work_limit, seconds=None, cpu_seconds=None
                    )
                    ledger.cancelled = lambda: calls.append(ledger.used) or (
                        cancel_after is not None and len(calls) >= cancel_after
                    )
                    budget = PollingBudget(ledger, interval=7)
                    try:
                        if grouped:
                            worker.budget = budget
                            build(worker, (lo, hi, maximum), 2**20)
                        else:
                            enumerated_queue(worker, (lo, hi, maximum), budget)
                    except BudgetExhaustedError:
                        refused = True
                    else:
                        refused = False
                    observations.append(
                        (ledger.used, ledger.reason, calls, refused)
                    )
                self.assertEqual(observations[0], observations[1])
                if work_limit == 100000 and cancel_after is None:
                    self.assertFalse(observations[0][-1])


class RootAndScoreLoopTests(unittest.TestCase):
    def test_native_root_pairs_are_validated_without_sorting(self):
        base = build_factor_base(10403, bound=100).factor_base
        polynomial = Polynomial(10403, 1, 49, 8)
        roots = tuple(
            polynomial_roots(polynomial, base, entry) for entry in base.entries
        )
        with patch(
            "v2.qs.families.sorted", side_effect=AssertionError, create=True
        ):
            self.assertTrue(verify_polynomial_roots(polynomial, base, roots))
            index = next(
                i for i, item in enumerate(roots) if len(item.roots) == 2
            )
            for pair in (
                roots[index].roots[::-1],
                (roots[index].roots[0],) * 2,
            ):
                damaged = list(roots)
                damaged[index] = replace(roots[index], roots=pair)
                with self.assertRaises(ValueError):
                    verify_polynomial_roots(polynomial, base, tuple(damaged))

    def test_full_period_reuses_roots_and_terminal_level_needs_no_target(self):
        roots = PolynomialRoots(3, (0,), False)

        class TerminalPolynomial:
            @property
            def n_prime(self):
                raise AssertionError("no lift needs target arithmetic")

        marks = list(
            prime_power_roots(
                TerminalPolynomial(),
                roots,
                3,
                1,
                unlimited_budget(),
                lo=-3,
                hi=3,
            )
        )
        self.assertEqual(marks, [(3, (0,), 1)])
        self.assertIs(marks[0][1], roots.roots)

    def test_fused_marks_and_period_overwrite_match_exact_valuations(self):
        base = build_factor_base(10403, bound=100).factor_base
        for a, b in ((1, 102), (7, 1), (49, 8)):
            polynomial = Polynomial(10403, 1, a, b)
            for marking in (
                "sparse",
                "density+presieve",
                "density+presieve_powers",
            ):
                for width in (17, 64):
                    worker = SieveCollector(
                        polynomial,
                        base,
                        budget=unlimited_budget(),
                        config=SieveConfig(
                            block_width=64,
                            score_policy="powers",
                            division="bucket",
                            marking=marking,
                            memory_bytes=64 * 2**20,
                        ),
                    )
                    worker._scores[:] = [-777] * 64
                    worker._hits[:] = [-1] * 64
                    lo, hi = -31, -31 + width
                    worker._sieve(lo, hi, Counter())
                    for offset, position in enumerate(range(lo, hi)):
                        value = abs(polynomial.value(position))
                        score = bits = 0
                        for index, prime in enumerate(base.primes):
                            remaining, exponent = value, 0
                            while remaining % prime == 0:
                                remaining //= prime
                                exponent += 1
                            score += exponent * (prime - 1).bit_length()
                            if exponent:
                                bits |= 1 << index
                        self.assertEqual(worker._scores[offset], score)
                        self.assertEqual(worker._hits[offset], bits)


class GroupReservationTests(unittest.TestCase):
    def test_group_counts_each_action_and_keeps_overdraw_atomic(self):
        ledger = Budget(work_limit=300, seconds=None, cpu_seconds=None)
        polls = []
        ledger.cancelled = lambda: polls.append(ledger.used) or False
        budget = PollingBudget(ledger)
        budget.consume_actions(64, 4)
        self.assertEqual(
            (ledger.used, polls, budget._remaining), (256, [0], 0)
        )
        budget.consume_actions(10, 4)
        self.assertEqual((ledger.used, polls), (296, [0, 256]))
        with self.assertRaises(BudgetExhaustedError):
            budget.consume_actions(2, 4)
        self.assertEqual(ledger.used, 296)
        with self.assertRaises(ValueError):
            budget.consume_actions(65, 4)
