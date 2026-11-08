"""Independent valuations, scalar identities and interrupted marking."""

import json
import unittest
from collections import defaultdict
from dataclasses import replace
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.execution.budget import BudgetExhaustedError
from v2.qs import Polynomial, SIQSConfig, SIQSJob, build_factor_base
from v2.qs.factor_base import FactorBase, FactorBaseEntry
from v2.qs.polynomial import polynomial_roots
from v2.qs.sieve_collector import SieveCollector, SieveConfig
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_qs import reference_positions, unlimited_budget
from v2.tests.test_qs_sieve import signature


def config(**kwargs):
    values = dict(
        block_width=17,
        division="bucket",
        score_policy="powers",
        residual_bound=97,
        max_atoms=4096,
        max_partials=4096,
        max_relations=4096,
        memory_bytes=64 * 2**20,
    )
    values.update(kwargs)
    return adapt_resource_configuration(SieveConfig, **values)


def worker(poly, base, marking, **kwargs):
    return SieveCollector(
        poly,
        base,
        config=config(marking=marking, **kwargs),
        budget=unlimited_budget(),
    )


class DensityKernelTests(unittest.TestCase):
    def test_scores_hits_and_order_against_exhaustive_division(self):
        # k=9 gives a singular multiplier prime; repeated A support and 2
        # exercise specialized loops without nonsingularity assumptions.
        for k, a, b in (
            (1, 1, 102),
            (1, 7, 1),
            (1, 49, 8),
            (9, 1, 306),
            (9, 3, 0),
        ):
            base = build_factor_base(
                10403, multiplier=k, bound=100
            ).factor_base
            poly = Polynomial(10403, k, a, b)
            for cap in (0, 256, 2**20):
                for width in (1, 7, 17, 64):
                    for policy in ("powers", "fixed"):
                        left = worker(
                            poly,
                            base,
                            "sparse",
                            power_plan_bytes=cap,
                            block_width=width,
                            score_policy=policy,
                        )
                        right = worker(
                            poly,
                            base,
                            "density",
                            power_plan_bytes=cap,
                            block_width=width,
                            score_policy=policy,
                        )
                        expected = reference_positions(poly, base, -71, 80, 97)
                        candidates = []
                        for collector in (left, right):
                            seen = []
                            divide = collector._divide

                            def record(position, offset, stats):
                                seen.append(position)
                                return divide(position, offset, stats)

                            with patch.object(collector, "_divide", record):
                                result = collector.collect(-71, 80)
                            self.assertEqual(result.reason, "complete")
                            self.assertEqual(signature(result), expected)
                            candidates.append(
                                (seen, result, collector.budget.used)
                            )
                        self.assertEqual(candidates[0], candidates[1])

            # The oracle visits every x and divides F(x), independently of
            # roots, lifted residues, hit ranges and marking implementation.
            right = worker(poly, base, "density", block_width=151)
            stats = defaultdict(int, threshold_min=2**63)
            right._sieve(-71, 80, stats)
            for offset, x in enumerate(range(-71, 80)):
                value = abs(a * x * x + 2 * b * x + (b * b - k * 10403) // a)
                self.assertEqual(
                    (a * x + b) ** 2 - k * 10403, a * poly.value(x)
                )
                bits = sum(
                    1 << i for i, p in enumerate(base.primes) if value % p == 0
                )
                self.assertEqual(right._hits[offset], bits)
                if value:
                    score = 0
                    for p, log in zip(base.primes, right._logs):
                        while value % p == 0:
                            score += log
                            value //= p
                    self.assertEqual(right._scores[offset], score)

    def test_backends_policies_and_division_fall_back_exactly(self):
        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 49, 8)
        for backend in ("list", "array", "bytearray"):
            for policy in (
                "adaptive",
                "conservative",
                "candidate",
                "powers",
                "fixed",
            ):
                for division in ("full", "roots", "bucket", "resieve"):
                    pair = [
                        worker(
                            poly,
                            base,
                            m,
                            score_backend=backend,
                            score_policy=policy,
                            division=division,
                            small_prime_cutoff=5,
                        ).collect(-31, 38)
                        for m in ("sparse", "density")
                    ]
                    self.assertEqual(pair[0], pair[1])

        # A clipped byte threshold must still retain a value with 300 twos.
        b = 2**298 - 1
        poly = Polynomial(b * b, 1, 1, b)
        base = build_factor_base(b * b, bound=3).factor_base
        for marking in ("sparse", "density"):
            result = worker(
                poly, base, marking, score_backend="bytearray"
            ).collect(2, 3)
            self.assertEqual(result.atoms[0].exponents, ((2, 300),))

    def test_singular_lift_cap_and_zero_value(self):
        base = FactorBase(
            3**12,
            1,
            11,
            (
                FactorBaseEntry(2, (1,)),
                FactorBaseEntry(3, (0,)),
                FactorBaseEntry(5, (1, 4)),
                FactorBaseEntry(7, (1, 6)),
            ),
        )
        poly = Polynomial(3**12, 1, 1, 0)
        for cap in (0, 256, 2**20):
            pair = [
                worker(poly, base, m, power_plan_bytes=cap).collect(-740, 741)
                for m in ("sparse", "density")
            ]
            self.assertEqual(pair[0], pair[1])
            # The square target yields a proper divisor at its zero norm.
            # Exact collection stops there; compare the committed prefix,
            # not positions beyond a successfully discovered factor.
            self.assertEqual(pair[1].reason, "factor_found")
            self.assertEqual(base.n % pair[1].divisor, 0)
            self.assertEqual(
                signature(pair[1]),
                reference_positions(
                    poly, base, -740, pair[1].next_position, 97
                ),
            )
            self.assertNotIn(-729, signature(pair[1]))
            self.assertNotIn(729, signature(pair[1]))

    def test_exact_reservation_sequence_and_interrupted_prefix(self):
        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 49, 8)
        for cap in (0, 256, 2**20):
            traces = []
            for marking in ("sparse", "density"):
                collector = worker(poly, base, marking, power_plan_bytes=cap)
                trace, consume = [], collector.budget.consume

                def record(amount=1):
                    trace.append(amount)
                    return consume(amount)

                with patch.object(collector.budget, "consume", record):
                    collector.collect(-31, 38)
                traces.append(trace)
            self.assertEqual(traces[0], traces[1])
            cuts = (0, 1, 20, len(trace) // 3, len(trace) // 2, len(trace) - 1)
            for cut in cuts:
                pair = []
                for marking in ("sparse", "density"):
                    collector = worker(
                        poly, base, marking, power_plan_bytes=cap
                    )
                    collector.budget.work_limit = collector.budget.used + sum(
                        trace[:cut]
                    )
                    stopped = collector.collect(-31, 38)
                    collector.budget = unlimited_budget()
                    resumed = collector.collect(stopped.next_position, 38)
                    pair.append((stopped, resumed))
                self.assertEqual(pair[0], pair[1])
                self.assertEqual(
                    signature(pair[1][1]),
                    reference_positions(poly, base, -31, 38, 97),
                )

    def test_cancelled_mark_keeps_first_uncommitted_position(self):
        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 49, 8)
        pair = []
        for marking in ("sparse", "density"):
            collector = worker(poly, base, marking)
            calls = 0

            def cancel():
                nonlocal calls
                calls += 1
                return calls == 23

            collector.budget = Budget(work_limit=10**9, cancelled=cancel)
            stopped = collector.collect(-31, 38)
            self.assertEqual(stopped.next_position, -31)
            self.assertFalse(stopped.atoms)
            pair.append(stopped)
        self.assertEqual(pair[0], pair[1])

    def test_switch_ownership_validation_and_checkpoint_rebuild(self):
        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 1, 102)
        first = worker(poly, base, "density", power_plan_bytes=2**20)
        second = worker(poly, base, "density", power_plan_bytes=2**20)
        first.collect(-31, 38)
        self.assertFalse(second._power_plans)
        old = (first.polynomial, first._roots, dict(first._power_plans))
        next_poly = Polynomial(10403, 1, 7, 1)
        roots = tuple(
            polynomial_roots(next_poly, base, e) for e in base.entries
        )
        original_config = first.config
        first.config = replace(first.config)
        with self.assertRaises(MemoryError):
            first.set_polynomial(next_poly, roots)
        self.assertEqual(
            (first.polynomial, first._roots, first._power_plans), old
        )
        first.config = original_config

        first.budget.work_limit = first.budget.used
        with self.assertRaises(BudgetExhaustedError):
            first.set_polynomial(next_poly, roots)
        self.assertEqual(
            (first.polynomial, first._roots, first._power_plans), old
        )
        first.budget = unlimited_budget()
        with self.assertRaises(ValueError):
            first.set_polynomial(next_poly, roots[:-1])
        first.set_polynomial(next_poly, roots)
        self.assertFalse(first._power_plans)
        self.assertIsNone(first._plan_window)

        with self.assertRaises(BudgetExhaustedError):
            SieveCollector(
                poly,
                base,
                config=config(marking="density"),
                budget=Budget(work_limit=0),
            )
        rebuilt = worker(poly, base, "density").collect(-31, 38)
        self.assertEqual(
            signature(rebuilt), reference_positions(poly, base, -31, 38, 97)
        )

        cfg = SIQSConfig(
            base_bound=1000,
            half_width=512,
            batch_width=1,
            collector=config(marking="density", power_plan_bytes=2**20),
        )
        job = SIQSJob(
            1000003 * 1000033, seed=7, config=cfg, budget=unlimited_budget()
        )
        first = job.run(max_blocks=1)
        self.assertIsNone(first.divisor)
        self.assertTrue(job.engine.collector._power_plans)
        snapshot = json.loads(json.dumps(job.checkpoint()))
        for _ in range(2):
            budget = unlimited_budget()
            resumed = SIQSJob.from_checkpoint(
                snapshot, config=cfg, budget=budget
            )
            self.assertFalse(resumed.engine.collector._power_plans)
            result = resumed.run(max_blocks=1)
            self.assertGreater(budget.used, snapshot["resources"]["work_used"])
            self.assertEqual((result.divisor or 1) * result.cofactor, job.n)
            self.assertEqual(resumed.config, cfg)


if __name__ == "__main__":
    unittest.main()
