"""One-class cache replay against independent Hensel and queue oracles."""

import unittest

from v2.benchmarks.infrastructure.resource_control import Budget
from v2.benchmarks.infrastructure.resource_control import (
    historical_replace as replace,
)
from v2.execution.budget import BudgetExhaustedError
from v2.qs import Polynomial, SIQSJob, build_factor_base
from v2.qs.polynomial import polynomial_roots
from v2.qs.power_reuse import RootTransfers
from v2.qs.power_sieve import prime_power_roots
from v2.qs.sparse_queue import build
from v2.tests.test_qs import unlimited_budget
from v2.tests.test_qs_hit_words import make_worker
from v2.tests.test_qs_loop_reuse import enumerated_queue, queue_worker
from v2.tests.test_qs_sieve import signature


class PrunedTransferTests(unittest.TestCase):
    def test_pruned_replay_requires_an_explicit_boolean(self):
        self.assertFalse(RootTransfers(2**20).pruned_replay)
        for value in (0, 1, None, "yes"):
            with self.assertRaises(ValueError):
                RootTransfers(2**20, pruned_replay=value)

    def test_pruned_classes_preserve_every_root_order_and_paid_action(self):
        base = build_factor_base(10403, bound=43).factor_base
        for a, values in (
            (7, (1, -1, 8, -8, 15, 29)),
            (49, (8, -8, 41, 57, -90)),
        ):
            cache = RootTransfers(2**20, pruned_replay=True)
            first = Polynomial(10403, 1, a, values[0])
            cache.bind(first)
            # Seed complete pairs independently of any observation window.
            for entry in base.entries:
                roots = polynomial_roots(first, base, entry)
                list(
                    prime_power_roots(
                        first,
                        roots,
                        entry.prime**4,
                        1,
                        unlimited_budget(),
                        root_transfers=cache,
                    )
                )
            for b in values:
                polynomial = Polynomial(10403, 1, a, b)
                cache.bind(polynomial)
                for entry in base.entries:
                    roots = polynomial_roots(polynomial, base, entry)
                    for lo, hi in ((0, 1), (-9, 10), (-73, 40), (31, 87)):
                        args = polynomial, roots, entry.prime**4, 1
                        direct, replay = unlimited_budget(), unlimited_budget()
                        extra = cache.extra_work
                        expected = list(
                            prime_power_roots(*args, direct, lo=lo, hi=hi)
                        )
                        actual = list(
                            prime_power_roots(
                                *args,
                                replay,
                                lo=lo,
                                hi=hi,
                                root_transfers=cache,
                            )
                        )
                        self.assertEqual(actual, expected)
                        self.assertEqual(
                            replay.used,
                            direct.used + cache.extra_work - extra,
                        )
                self.assertTrue(
                    all(len(row[1]) == 2 for row in cache.values.values())
                )
                self.assertLessEqual(cache.used, cache.byte_limit)
            self.assertGreater(cache.pruned_replays, 0)

    def test_interrupted_pair_replay_never_publishes_partial_cache_state(self):
        base = build_factor_base(10403, bound=43).factor_base
        cache = RootTransfers(2**20, pruned_replay=True)
        first = Polynomial(10403, 1, 7, 1)
        cache.bind(first)
        entry = next(
            e
            for e in base.entries
            if e.prime > 7 and len(polynomial_roots(first, base, e).roots) == 2
        )
        roots = polynomial_roots(first, base, entry)
        list(
            prime_power_roots(
                first,
                roots,
                entry.prime**3,
                1,
                unlimited_budget(),
                root_transfers=cache,
            )
        )
        second = Polynomial(10403, 1, 7, -1)
        cache.bind(second)
        initial = polynomial_roots(second, base, entry).roots
        before = dict(cache.values), cache.used, cache.extra_work
        budget = Budget(work_limit=1, seconds=None, cpu_seconds=None)
        with self.assertRaises(BudgetExhaustedError):
            cache.shifted_class(
                second,
                entry.prime,
                entry.prime**2,
                initial,
                initial[0],
                1,
                budget,
            )
        self.assertEqual((cache.values, cache.used, cache.extra_work), before)
        self.assertEqual((budget.used, cache.pruned_replays), (1, 0))

    def test_collector_masks_atoms_and_sparse_event_order_remain_exact(self):
        reference = make_worker(64, sparse=2**20)
        cached = make_worker(64, kernels=("pruned_transfers",), sparse=2**20)
        self.assertTrue(cached._power_transfers.pruned_replay)
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
            self.assertEqual(
                signature(reference.collect(-79, 82)),
                signature(cached.collect(-79, 82)),
            )
            self.assertEqual(reference._scores, cached._scores)
            self.assertEqual(reference._hits[:], cached._hits[:])

        worker = queue_worker(7)
        worker._power_transfers = RootTransfers(2**20, pruned_replay=True)
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

    def test_resume_rebuilds_derived_pairs_without_changing_relations(self):
        from v2.tests.test_current_qs_workers import (
            N,
            allowance,
            configuration,
        )

        config = configuration(extended=True)
        config = replace(
            config,
            checkpoint_bytes=2**20,
            collector=replace(
                config.collector, marking="wide_transfers+pruned_transfers"
            ),
        )
        job = SIQSJob(N, seed=29, config=config, budget=allowance())
        job.run(max_blocks=1)
        restored = SIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )
        self.assertTrue(
            restored.engine.collector._power_transfers.pruned_replay
        )
        expected, actual = job.run(max_blocks=4), restored.run(max_blocks=4)
        self.assertEqual(
            (actual.reason, actual.divisor),
            (expected.reason, expected.divisor),
        )

        def atoms(owner):
            return sorted(
                (a.relation_id, a.sign, tuple(a.exponents), a.residual_product)
                for a in owner.engine.collector._atoms.values()
            )

        self.assertEqual(atoms(restored), atoms(job))
        self.assertEqual(
            actual.stats["engine"]["collected_positions"],
            expected.stats["engine"]["collected_positions"],
        )


if __name__ == "__main__":
    unittest.main()
