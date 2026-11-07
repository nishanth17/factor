"""Funded cross-polynomial inverse reuse and exact scarce-memory fallback."""

import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.qs import Polynomial, SieveCollector
from v2.qs.polynomial import polynomial_roots
from v2.tests.test_qs import unlimited_budget
from v2.tests.test_qs_hit_words import make_worker
from v2.tests.test_qs_sieve import signature


class WideInverseTests(unittest.TestCase):
    def test_larger_cache_keeps_scores_masks_and_atomic_relations_exact(self):
        # Exercise the actual old/new capacity branch with a small complete
        # arithmetic oracle, rather than allocating a large training base.
        with patch("v2.qs.sieve_collector.MAX_SHARED_INVERSE_BYTES", 1):
            reference = make_worker(64, sparse=2**20)
            cached = make_worker(64, kernels=("wide_inverses",), sparse=2**20)
        self.assertIsNone(reference._shared_lift_inverses)
        self.assertIsNotNone(cached._shared_lift_inverses)
        reference._power_transfers = cached._power_transfers = None
        for b in (0, 1, -1, 3):
            polynomial = Polynomial(reference.factor_base.n, 1, 1, b)
            roots = tuple(
                polynomial_roots(polynomial, reference.factor_base, entry)
                for entry in reference.factor_base.entries
            )
            for worker in (reference, cached):
                worker.set_polynomial(
                    polynomial, roots, retain_relations=False
                )
            left, right = reference.collect(-79, 82), cached.collect(-79, 82)
            self.assertEqual(signature(left), signature(right))
            self.assertEqual(reference._scores, cached._scores)
            self.assertEqual(reference._hits[:], cached._hits[:])
            self.assertLessEqual(cached._workspace, cached.config.memory_limit)
        self.assertGreater(len(cached._shared_lift_inverses.values), 0)
        self.assertLessEqual(
            len(cached._shared_lift_inverses.values),
            len(cached.factor_base.entries),
        )

    def test_optional_inverse_cache_leaves_the_next_switch_funded(self):
        with patch("v2.qs.sieve_collector.MAX_SHARED_INVERSE_BYTES", 1):
            reference = make_worker(64)
        count = len(reference.factor_base.entries)
        mandatory = reference._workspace - reference._reuse_workspace
        inverse_bytes = 256 + 512 * count
        switch_bytes = 512 * count + 32768
        for extra, enabled in ((-1, False), (1, True)):
            config = replace(
                reference.config,
                marking="wide_inverses",
                memory_bytes=mandatory + inverse_bytes + switch_bytes + extra,
            )
            worker = SieveCollector(
                reference.polynomial,
                reference.factor_base,
                config=config,
                budget=unlimited_budget(),
            )
            self.assertEqual(worker._shared_lift_inverses is not None, enabled)
            roots = tuple(
                polynomial_roots(worker.polynomial, worker.factor_base, entry)
                for entry in worker.factor_base.entries
            )
            worker.set_polynomial(worker.polynomial, roots)
            self.assertLessEqual(
                worker._workspace + switch_bytes, config.memory_limit
            )


if __name__ == "__main__":
    unittest.main()
