"""Byte epochs preserve full recovery masks and fund every wraparound clear."""

import unittest
from collections import Counter

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.qs.byte_hit_words import ByteStampWordHits
from v2.qs.hit_words import WORD_MASK, WordHits
from v2.tests.test_qs import unlimited_budget
from v2.tests.test_qs_hit_words import make_worker
from v2.tests.test_qs_sieve import signature


class ByteStampTests(unittest.TestCase):
    def test_epochs_wrap_without_stale_bits_or_changed_cached_reads(self):
        budget = unlimited_budget()
        compact = ByteStampWordHits(131, [0] * 9, budget)
        reference = WordHits(131, [0] * 9, budget)
        self.assertTrue(all(isinstance(s, bytearray) for s in compact.stamps))
        for epoch in range(270):
            expected = [epoch % 2] * 9
            for words in (compact, reference):
                words.begin(epoch % 2, budget)
                for index in (0, 59, 60, 119, 120, 130):
                    words.add(index, range(index % 3, 9, 3))
                words.seal()
            for index in (0, 59, 60, 119, 120, 130):
                for offset in range(index % 3, 9, 3):
                    expected[offset] |= 1 << index
            self.assertEqual(compact[:], expected)
            self.assertEqual(reference[:], expected)
            self.assertEqual(compact.mask(3, budget), expected[3])
            self.assertEqual(compact.mask(3, budget), expected[3])
        self.assertEqual(compact.epoch_resets, 1)
        self.assertLessEqual(compact.epoch, 255)

    def test_wrap_refusal_keeps_epoch_until_all_clears_are_funded(self):
        budget = unlimited_budget()
        words = ByteStampWordHits(121, [0] * 7, budget)
        words.begin(0, budget)
        words.add(120, range(7))
        words.epoch = 255
        refused = Budget(work_limit=8, seconds=None, cpu_seconds=None)
        with self.assertRaises(BudgetExhaustedError):
            words.begin(1, refused)
        self.assertEqual(words.epoch, 255)
        self.assertEqual(words.epoch_resets, 0)
        self.assertEqual(refused.used, 8)
        words.begin(1, budget)
        self.assertEqual(words[:], [1] * 7)
        self.assertEqual(words.epoch_resets, 1)
        self.assertEqual(words.epoch, 1)

        before = words.work
        words.epoch = 255
        with self.assertRaises(ValueError):
            words.begin(WORD_MASK + 1, budget)
        self.assertEqual(words.work, before)
        self.assertEqual(words.epoch, 255)

    def test_periodic_prefix_general_masks_and_negative_oracle_views(self):
        budget = unlimited_budget()
        compact = ByteStampWordHits(131, [0] * 9, budget)
        reference = WordHits(131, [0] * 9, budget)
        for values in ([1, 3, 5], [1 << 125, 7, 1 << 60], [0, 1, 0]):
            for words in (compact, reference):
                words.assign_prefix(values, budget)
            self.assertEqual(compact[:], reference[:])
        for words in (compact, reference):
            words[::2] = [1 << i for i in (1, 60, 90, 120, 130)]
        self.assertEqual(compact[:], reference[:])
        for words in (compact, reference):
            words[:] = [-1] * 9
        self.assertEqual(compact[:], reference[:])
        compact.cache_epoch = WORD_MASK
        compact.begin(0, budget)
        self.assertEqual(compact[:], [0] * 9)

    def test_collector_switch_wrap_and_interruption_match_native_stamps(self):
        kernels = ("presieve", "unroll4", "parity_two")
        for sparse, packed in ((0, False), (2**20, False), (2**20, True)):
            selected = kernels + (("packed_queue",) if packed else ())
            actual = make_worker(64, selected + ("byte_stamps",), sparse)
            expected = make_worker(64, selected, sparse)
            self.assertIsInstance(actual._word_hits, ByteStampWordHits)
            actual._word_hits.epoch = 254
            for lo, hi in ((-79, -15), (3, 16), (-31, 33)):
                a, b = Counter(), Counter()
                self.assertEqual(
                    actual._sieve(lo, hi, a), expected._sieve(lo, hi, b)
                )
                self.assertEqual(actual._scores, expected._scores)
                self.assertEqual(actual._hits[:], expected._hits[:])
                self.assertEqual(a, b)
            self.assertGreater(actual._word_hits.epoch_resets, 0)

            actual.budget = Budget(
                work_limit=1, seconds=None, cpu_seconds=None
            )
            stopped = actual.collect(-79, 82)
            self.assertEqual(stopped.next_position, -79)
            actual.budget = unlimited_budget()
            left, right = actual.collect(-79, 82), expected.collect(-79, 82)
            self.assertEqual(signature(left), signature(right))
            self.assertEqual(left.next_position, right.next_position)
