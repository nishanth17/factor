"""Funded larger derived-root storage with exact arithmetic fallback."""

import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.qs import Polynomial, SieveCollector
from v2.qs.polynomial import polynomial_roots
from v2.qs.power_reuse import MAX_EXTENDED_TRANSFER_BYTES, RootTransfers
from v2.tests.test_qs import unlimited_budget
from v2.tests.test_qs_hit_words import make_worker
from v2.tests.test_qs_sieve import signature


class WideTransferTests(unittest.TestCase):
    def test_capacity_requires_an_explicit_finite_larger_grant(self):
        with self.assertRaises(ValueError):
            RootTransfers(33 * 2**20)
        for invalid in (False, -1, 1.5, MAX_EXTENDED_TRANSFER_BYTES + 1):
            with self.assertRaises((TypeError, ValueError)):
                RootTransfers(0, capacity_limit=invalid)
        cache = RootTransfers(
            MAX_EXTENDED_TRANSFER_BYTES,
            capacity_limit=MAX_EXTENDED_TRANSFER_BYTES,
        )
        self.assertEqual(cache.byte_limit, MAX_EXTENDED_TRANSFER_BYTES - 65536)
        self.assertFalse(cache.values)

    def test_larger_cache_preserves_scores_masks_atoms_and_switches(self):
        # Shrink only the legacy admission ceiling; use the production
        # larger-capacity branch with a complete small arithmetic oracle.
        with patch("v2.qs.sieve_collector.MAX_TRANSFER_BYTES", 65536):
            reference = make_worker(64, sparse=2**20)
            cached = make_worker(64, kernels=("wide_transfers",), sparse=2**20)
        self.assertEqual(reference._power_transfers.byte_limit, 0)
        self.assertGreater(cached._power_transfers.byte_limit, 0)
        for a, b in ((1, 0), (1, 1), (1, -1), (1, 3), (2, 1)):
            polynomial = Polynomial(reference.factor_base.n, 1, a, b)
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
        self.assertGreater(cached._power_transfers.builds, 0)
        self.assertGreater(cached._power_transfers.replays, 0)
        self.assertLessEqual(
            cached._power_transfers.used, cached._power_transfers.byte_limit
        )

    def test_larger_cache_never_spends_the_next_polynomial_switch(self):
        template = make_worker(64)
        options = replace(template.config, division="roots", marking="sparse")
        with (
            patch("v2.qs.sieve_collector.MAX_TRANSFER_BYTES", 65536),
            patch("v2.qs.sieve_collector.MAX_SHARED_INVERSE_BYTES", 1),
        ):
            reference = SieveCollector(
                template.polynomial,
                template.factor_base,
                config=options,
                budget=unlimited_budget(),
            )
            count = len(template.factor_base.entries)
            mandatory = reference._workspace - reference._reuse_workspace
            desired = min(MAX_EXTENDED_TRANSFER_BYTES, 65536 + 8192 * count)
            switch = 512 * count + 32768
            for extra, enabled in ((-1, False), (1, True)):
                config = replace(
                    options,
                    marking="wide_transfers",
                    memory_bytes=mandatory + desired + switch + extra,
                )
                worker = SieveCollector(
                    template.polynomial,
                    template.factor_base,
                    config=config,
                    budget=unlimited_budget(),
                )
                self.assertEqual(worker._power_transfers is not None, enabled)
                self.assertLessEqual(
                    worker._workspace + switch, config.memory_limit
                )


if __name__ == "__main__":
    unittest.main()
