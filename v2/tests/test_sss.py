"""Independent smoothness/CRT oracles and interrupted SSS factor recovery."""

import contextlib
import io
import random
import unittest
from dataclasses import replace
from math import prod
from unittest.mock import patch

from v2 import utils
from v2.budget import Budget
from v2.qs import (
    build_factor_base,
    qs_polynomial,
    verify_atomic,
    verify_combined,
)
from v2.qs.smooth_batch import SmoothBatch, product_tree
from v2.qs.sss import SSSCollector, SSSConfig, SSSJob


def allowance():
    return Budget(work_limit=10**10, seconds=None, cpu_seconds=None)


def trial_residual(value, primes):
    """Use scalar division without a product, remainder or gcd algorithm."""
    for prime in primes:
        while value % prime == 0:
            value //= prime
    return value


def search_collector(n=1009 * 1013, **changes):
    config = replace(SSSConfig(base_bound=400), **changes)
    base = build_factor_base(n, bound=config.base_bound).factor_base
    return SSSCollector(
        qs_polynomial(base), base, search_config=config, budget=allowance()
    )


class SmoothBatchTests(unittest.TestCase):
    def test_random_high_powers_and_duplicate_leaves(self):
        primes = (2, 3, 5, 7, 11)
        generator = random.Random(351)
        values = (1, 2**257 * 3**81, 2**257 * 13, 17**32, 45, 45)
        values += tuple(generator.randrange(1, 10**12) for _ in range(128))
        detector = SmoothBatch(primes, budget=allowance())
        self.assertEqual(
            detector.residuals(values),
            tuple(trial_residual(v, primes) for v in values),
        )
        self.assertEqual(detector.residuals(()), ())
        self.assertEqual(
            SmoothBatch((), budget=allowance()).residuals(values), values
        )

    def test_odd_tree_and_leaf_moduli(self):
        for values in ((), (1,), (3, 4, 5, 7, 9)):
            tree = product_tree(values, budget=allowance())
            if values:
                self.assertEqual(tree[-1], (prod(values),))
            else:
                self.assertEqual(tree, ())
        detector = SmoothBatch((2, 3), budget=allowance())
        self.assertEqual(detector.residuals((1, 2, 3, 6, 12, 36)), (1,) * 6)

    def test_tree_limits_refuse_before_work_and_invalid_bases(self):
        for options in (
            dict(max_bits=3),
            dict(max_nodes=2),
            dict(memory_bytes=1),
        ):
            budget = allowance()
            with self.assertRaises(MemoryError):
                product_tree((2, 3, 5), budget=budget, **options)
            self.assertEqual(budget.used, 0)
        for primes in ((2, 4), (3, 2), (2, 2), (100003,)):
            with self.assertRaises(ValueError):
                SmoothBatch(primes, budget=allowance())
        for values in ((0,), (-1,), (True,)):
            with self.assertRaises((TypeError, ValueError)):
                product_tree(values, budget=allowance())


class CollisionTests(unittest.TestCase):
    def test_crt_and_signed_collisions_against_enumeration(self):
        collector = search_collector()
        selected = (0, 1, 2)
        roots = collector._roots
        small_count = collector.small_count
        remaining_roots = roots[small_count:]
        modulus = prod(roots[i].prime for i in selected)
        choices = {i: roots[i].roots[0] for i in selected}
        expected = set()
        for changed in selected:
            if len(roots[changed].roots) != 2:
                continue
            choices[changed] = roots[changed].roots[1]
            position = next(
                x
                for x in range(modulus)
                if all(x % roots[i].prime == choices[i] for i in selected)
            )
            for dropped in (None,) + selected:
                if dropped == changed:
                    continue
                divisor = 1 if dropped is None else roots[dropped].prime
                step = modulus // divisor
                for shift in range(-roots[-1].prime, roots[-1].prime):
                    argument = position + shift * step
                    count = sum(
                        -entry.prime <= shift < entry.prime
                        and collector.polynomial.value(argument) % entry.prime
                        == 0
                        for entry in remaining_roots
                    )
                    if count >= 3:
                        expected.add(
                            (
                                argument,
                                abs(collector.polynomial.value(argument))
                                // step,
                            )
                        )
        from v2.qs.sss import collision_candidates

        actual = collision_candidates(
            collector.polynomial,
            roots,
            collector.small_count,
            selected,
            collector.coefficients,
            budget=allowance(),
        )
        self.assertEqual(set(actual), expected)
        self.assertTrue(any(position < 0 for position, _ in actual))

    def test_seed_repeatability_and_assignment_cap(self):
        first, second = search_collector(), search_collector()
        for index in range(8):
            self.assertEqual(first.assignment(index), second.assignment(index))
        capped = search_collector(max_candidates=1)
        with self.assertRaises(MemoryError):
            capped.assignment(0)
        self.assertEqual(capped._atoms, {})


class SSSJobTests(unittest.TestCase):
    def test_complete_balanced_inputs_and_quiet_library(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            for mode in ("sss", "sssf"):
                for n in (1009 * 1013, 4001 * 4003, 10007 * 10009):
                    job = SSSJob(
                        n,
                        config=SSSConfig(mode=mode, base_bound=400),
                        budget=allowance(),
                    )
                    result = job.run()
                    self.assertTrue(utils.valid_divisor(result.divisor, n))
                    self.assertEqual(result.divisor * result.cofactor, n)
                    self.assertEqual(result.reason, "factor_found")
                    again = job.run()
                    self.assertEqual(again.divisor, result.divisor)
                    self.assertEqual(
                        again.stats["work_used"], result.stats["work_used"]
                    )
        self.assertEqual(output.getvalue(), "")

    def test_candidate_refusal_resume_retains_the_exact_store(self):
        reference = search_collector(n=4001 * 4003)
        expected = reference.collect(0, 3)
        for cut in (1, 20, 100, 300, 600):
            collector = search_collector(n=4001 * 4003)
            polls = [0]

            def cancel():
                polls[0] += 1
                return polls[0] >= cut

            collector.budget.cancelled = cancel
            refused = collector.collect(0, 3)
            self.assertEqual(refused.reason, "cancelled")
            spent = collector.budget.used
            collector.budget.cancelled = None
            resumed = collector.collect(refused.next_position, 3)
            self.assertEqual(resumed.atoms, expected.atoms)
            self.assertEqual(
                resumed.combined_relations, expected.combined_relations
            )
            self.assertGreaterEqual(collector.budget.used, spent)
            for atom in resumed.atoms:
                verify_atomic(
                    atom,
                    collector.factor_base,
                    residual_bound=10000,
                    budget=allowance(),
                )
            for relation in resumed.combined_relations:
                verify_combined(
                    relation,
                    collector.factor_base,
                    collector._atoms,
                    budget=allowance(),
                )

    def test_pause_and_budget_extension_match_uninterrupted_factor(self):
        config = SSSConfig(base_bound=400)
        baseline = SSSJob(4001 * 4003, config=config, budget=allowance()).run()
        budget = allowance()
        job = SSSJob(4001 * 4003, config=config, budget=budget)
        first = job.run(batch_limit=1)
        self.assertEqual(first.reason, "paused")
        spent = budget.used
        budget.work_limit = spent
        refused = job.run()
        self.assertEqual(refused.reason, "work_limit")
        budget.work_limit = 10**10
        self.assertEqual(job.run().divisor, baseline.divisor)
        self.assertGreater(budget.used, spent)
        job.budget = allowance()
        with self.assertRaises(ValueError):
            job.run()

    def test_store_caps_final_solve_and_filter_loss(self):
        config = SSSConfig(
            base_bound=400,
            collector=replace(
                SSSConfig().collector,
                max_relations=0,
            ),
        )
        job = SSSJob(4001 * 4003, config=config, budget=allowance())
        result = job.run()
        self.assertEqual(result.reason, "relation_limit")
        self.assertIsNone(result.divisor)
        self.assertEqual(
            job.run().stats["work_used"], result.stats["work_used"]
        )
        filtered = search_collector(mode="sssf", filter_bound=1)
        result = filtered.collect(0, 3)
        self.assertGreater(result.stats["filter_rejections"], 0)
        self.assertLessEqual(
            result.workspace_bytes, filtered.config.memory_bytes
        )
        for atom in result.atoms:
            verify_atomic(atom, filtered.factor_base, residual_bound=10000)

    def test_exhaustion_setup_failure_and_validation(self):
        job = SSSJob(
            1009 * 1013,
            config=SSSConfig(base_bound=200, search_rounds=1),
            budget=allowance(),
        )
        result = job.run()
        self.assertEqual(result.reason, "search_exhausted")
        self.assertEqual(result.cofactor, job.n)
        self.assertEqual(
            job.run().stats["work_used"], result.stats["work_used"]
        )
        self.assertEqual(SSSJob(10403).run().divisor, 101)
        self.assertEqual(
            SSSJob(10403, budget=Budget(work_limit=0)).run().reason,
            "work_limit",
        )
        self.assertEqual(
            SSSJob(10403, config=SSSConfig(memory_bytes=1)).run().reason,
            "memory_limit",
        )
        refused = SSSJob(1009 * 1013, config=SSSConfig(memory_bytes=100000))
        result = refused.run()
        self.assertEqual(result.reason, "memory_limit")
        self.assertEqual(
            refused.run().stats["work_used"], result.stats["work_used"]
        )
        for changes in (
            dict(mode="other"),
            dict(search_rounds=0),
            dict(selection_size=17),
            dict(filter_bound=-1),
        ):
            with self.assertRaises(ValueError):
                SSSConfig(**changes)
        job = SSSJob(4001 * 4003, budget=allowance())
        with patch(
            "v2.qs.sss.collision_candidates", side_effect=ArithmeticError
        ):
            with self.assertRaises(ArithmeticError):
                job.run()


if __name__ == "__main__":
    unittest.main()
