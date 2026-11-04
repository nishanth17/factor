"""Independent full-mask and exact collector/resume oracles for word hits."""

import unittest
from collections import Counter
from dataclasses import replace

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.qs import Polynomial, build_factor_base
from v2.qs.hit_words import WORD_MASK, WordHits
from v2.qs.polynomial import polynomial_roots
from v2.qs.sieve_collector import SieveCollector, SieveConfig
from v2.tests.test_qs import unlimited_budget
from v2.tests.test_qs_sieve import signature


class HitWordTests(unittest.TestCase):
    def test_words_epochs_and_cached_reads_match_independent_big_masks(self):
        budget = unlimited_budget()
        hits = WordHits(181, [0] * 17, budget)
        for initial in (0, 1, 0, 1):
            hits.begin(initial, budget)
            expected = [initial] * 17
            scores = [0] * 17
            for index in (0, 1, 59, 60, 61, 119, 120, 179, 180):
                offsets = range(index % 3, 17, index % 5 + 1)
                hits.add(index, offsets, scores, 7)
                for offset in offsets:
                    expected[offset] |= 1 << index
            self.assertEqual(hits[:], expected)
            hits.seal()
            before = hits.work
            self.assertEqual(hits.mask(5, budget), expected[5])
            self.assertEqual(hits.work - before, len(hits.words) + 1)
            self.assertEqual(hits.mask(5, budget), expected[5])
            self.assertEqual(hits.work - before, len(hits.words) + 1)
            hits.add_one(120, 5)
            expected[5] |= 1 << 120
            self.assertEqual(list(hits), expected)
            self.assertEqual(hits[-1], expected[-1])

    def test_periodic_initializers_tail_and_private_list_view(self):
        budget = unlimited_budget()
        hits = WordHits(131, [0] * 9, budget)
        hits.begin(1, budget)
        hits.add(125, range(9))
        for values in ([1, 3, 5], [1 << 125, 7, 1 << 60], [0, 1, 0]):
            hits.assign_prefix(values, budget)
            self.assertEqual(hits[:3], values)
            self.assertEqual(hits[3:], [0] * 6)
        hits[::2] = [1 << i for i in (1, 60, 90, 120, 130)]
        expected = [2, 1, 1 << 60, 0, 1 << 90, 0, 1 << 120, 0, 1 << 130]
        self.assertEqual(hits, expected)
        hits[:] = [-1] * 9
        self.assertEqual(hits[:], [-1] * 9)
        hits.begin(0, budget)
        self.assertEqual(hits[:], [0] * 9)

    def test_epoch_wrap_refusal_and_candidate_read_refusal(self):
        budget = unlimited_budget()
        hits = WordHits(121, [0] * 7, budget)
        hits.begin(0, budget)
        hits.add(120, range(7))
        hits.epoch = WORD_MASK
        refused = Budget(work_limit=8, seconds=None, cpu_seconds=None)
        with self.assertRaises(BudgetExhaustedError):
            hits.begin(0, refused)
        self.assertEqual(hits.epoch, WORD_MASK)
        hits.begin(1, budget)
        self.assertEqual(hits[:], [1] * 7)
        hits.seal()
        refused = Budget(work_limit=0, seconds=None, cpu_seconds=None)
        with self.assertRaises(BudgetExhaustedError):
            hits.mask(0, refused)
        self.assertEqual(hits[0], 1)
        hits.cache_epoch = WORD_MASK
        hits.begin(0, budget)
        self.assertEqual(hits[:], [0] * 7)


def make_worker(width=64, kernels=(), sparse=0):
    base = build_factor_base(1000003 * 1000033, bound=2000).factor_base
    return SieveCollector(
        Polynomial(base.n, 1, 1, 0),
        base,
        budget=unlimited_budget(),
        config=SieveConfig(
            division="bucket",
            score_policy="powers",
            block_width=width,
            marking="+".join(kernels) if kernels else "sparse",
            sparse_queue_bytes=sparse,
            memory_bytes=64 * 2**20,
            residual_bound=2003,
        ),
    )


class HitCollectorTests(unittest.TestCase):
    def test_every_mask_and_score_matches_list_and_modular_enumeration(self):
        for width in (17, 64):
            for kernels in (
                (),
                ("presieve", "unroll4", "parity_two"),
                ("presieve_powers",),
            ):
                for sparse in (0, 2**20):
                    worker = make_worker(width, kernels, sparse)
                    reference = make_worker(width, kernels, sparse)
                    reference._hits = [0] * width
                    self.assertIsNotNone(worker._word_hits)
                    worker._set_plan_interval(-79, 82)
                    reference._set_plan_interval(-79, 82)
                    for lo, hi in ((-79, -79 + width), (3, 16)):
                        actual_stats, expected_stats = Counter(), Counter()
                        self.assertEqual(
                            worker._sieve(lo, hi, actual_stats),
                            reference._sieve(lo, hi, expected_stats),
                        )
                        self.assertEqual(worker._scores, reference._scores)
                        self.assertEqual(
                            worker._hits[: hi - lo], reference._hits[: hi - lo]
                        )
                        self.assertEqual(actual_stats, expected_stats)
                        for offset, position in enumerate(range(lo, hi)):
                            mask = sum(
                                1 << index
                                for index, entry in enumerate(
                                    worker.factor_base.entries
                                )
                                if worker.polynomial.value(position)
                                % entry.prime
                                == 0
                            )
                            self.assertEqual(worker._hits[offset], mask)

    def test_refused_block_and_resumed_publication_match_uninterrupted(self):
        worker = make_worker(kernels=("presieve", "unroll4", "parity_two"))
        reference = make_worker(kernels=("presieve", "unroll4", "parity_two"))
        reference._hits = [0] * 64
        expected = reference.collect(-79, 82)
        worker.budget = Budget(work_limit=1, seconds=None, cpu_seconds=None)
        refused = worker.collect(-79, 82)
        self.assertEqual(refused.next_position, -79)
        self.assertEqual(refused.atoms, ())
        worker.budget = unlimited_budget()
        actual = worker.collect(refused.next_position, 82)
        self.assertEqual(actual.reason, expected.reason)
        self.assertEqual(signature(actual), signature(expected))
        self.assertEqual(actual.next_position, expected.next_position)

    def test_optional_construction_declines_unfunded_memory(self):
        worker = make_worker()
        count, width = len(worker.factor_base.entries), 64
        self.assertLessEqual(WordHits.reservation(count, width), 32 * 2**20)
        self.assertEqual(len(worker._word_hits.descriptors), count)
        limit = worker._workspace + 512 * count + 32768 - 1
        fallback = SieveCollector(
            worker.polynomial,
            worker.factor_base,
            budget=unlimited_budget(),
            config=replace(worker.config, memory_bytes=limit),
        )
        self.assertIsNone(fallback._word_hits)
        self.assertIsInstance(fallback._hits, list)
        with self.assertRaises(BudgetExhaustedError):
            WordHits(
                count,
                [0] * width,
                Budget(work_limit=0, seconds=None, cpu_seconds=None),
            )

    def test_polynomial_switch_rebuilds_words_and_a_support_exactly(self):
        worker = make_worker(kernels=("presieve", "unroll4", "parity_two"))
        reference = make_worker(kernels=("presieve", "unroll4", "parity_two"))
        reference._hits = [0] * 64
        for b in (1, -1, 3, -3):
            polynomial = Polynomial(worker.factor_base.n, 1, 2, b)
            roots = tuple(
                polynomial_roots(polynomial, worker.factor_base, entry)
                for entry in worker.factor_base.entries
            )
            worker.set_polynomial(polynomial, roots)
            reference.set_polynomial(polynomial, roots)
            actual = worker.collect(-79, 82)
            expected = reference.collect(-79, 82)
            self.assertEqual(actual.reason, expected.reason)
            self.assertEqual(signature(actual), signature(expected))
            self.assertEqual(actual.next_position, expected.next_position)


if __name__ == "__main__":
    unittest.main()
