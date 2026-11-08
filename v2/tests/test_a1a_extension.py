"""Independent extended-period scores, exact root writes and charged resume."""

import itertools
import json
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.qs import CountedSIQSConfig, Polynomial, SIQSJob, build_factor_base
from v2.qs.collector_kernels import add_unrolled, combine, workspace
from v2.qs.factor_base import FactorBase, FactorBaseEntry
from v2.qs.polynomial import polynomial_roots
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_a1a_density import config, worker
from v2.tests.test_a1a_integration import SELECTED, scored
from v2.tests.test_qs import reference_positions, unlimited_budget
from v2.tests.test_qs_sieve import signature

POWER = "density+presieve_powers+scan_fused"
UNROLL = "density+presieve+unroll4+scan_fused"
BOTH = "density+presieve_powers+unroll4+scan_fused"
ARMS = (POWER, UNROLL, BOTH)


class ExtendedKernelTests(unittest.TestCase):
    def test_unrolled_writes_against_position_predicate(self):
        for width, step, bit in itertools.product(
            (1, 2, 3, 4, 7, 17, 64), (1, 2, 3, 4, 7, 17, 64), (0, 1 << 71)
        ):
            for offset in range(step):
                scores, hits = [3] * width, [2] * width
                add_unrolled(scores, hits, offset, width, step, 7, bit)
                for position in range(width):
                    selected = (
                        position >= offset and (position - offset) % step == 0
                    )
                    self.assertEqual(scores[position], 10 if selected else 3)
                    self.assertEqual(
                        hits[position], (2 | bit) if selected else 2
                    )

    def test_exact_valuations_bits_and_candidate_order(self):
        for k, a, b in ((1, 1, 102), (1, 7, 1), (1, 49, 8), (9, 3, 0)):
            base = build_factor_base(
                10403, multiplier=k, bound=100
            ).factor_base
            poly = Polynomial(10403, k, a, b)
            for width in (1, 7, 17, 64, 151):
                for policy in ("powers", "fixed"):
                    for arm in ARMS:
                        collector = worker(
                            poly,
                            base,
                            arm,
                            block_width=width,
                            score_policy=policy,
                        )
                        end = min(-71 + width, 80)
                        scored(collector, -71, end)
                        for offset, x in enumerate(range(-71, end)):
                            value = abs(
                                a * x * x
                                + 2 * b * x
                                + (b * b - k * 10403) // a
                            )
                            score, bits = 0, 0
                            for index, (p, weight) in enumerate(
                                zip(base.primes, collector._logs)
                            ):
                                if value % p == 0:
                                    bits |= 1 << index
                                while value and value % p == 0:
                                    score += weight
                                    value //= p
                            self.assertEqual(collector._scores[offset], score)
                            self.assertEqual(collector._hits[offset], bits)

            expected = reference_positions(poly, base, -71, 80, 97)
            for cap in (0, 256, 2**20):
                traces = []
                for arm in (SELECTED, *ARMS):
                    collector = worker(
                        poly, base, arm, block_width=64, power_plan_bytes=cap
                    )
                    seen, divide = [], collector._divide

                    def record(position, offset, stats):
                        seen.append(position)
                        return divide(position, offset, stats)

                    with patch.object(collector, "_divide", record):
                        result = collector.collect(-71, 80)
                    self.assertEqual(signature(result), expected)
                    traces.append(seen)
                self.assertTrue(all(trace == traces[0] for trace in traces))

    def test_bound_changes_rotations_backends_cutoffs_and_policies(self):
        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 49, 8)
        for backend, policy, cutoff in itertools.product(
            ("list", "array", "bytearray"),
            ("powers", "fixed", "adaptive", "conservative"),
            (0, 5, 17),
        ):
            options = dict(
                block_width=64,
                score_backend=backend,
                score_policy=policy,
                small_prime_cutoff=cutoff,
            )
            for lo, hi in ((-64, -1), (-7, 18), (18, 19), (-31, -14)):
                expected = scored(
                    worker(poly, base, "sparse", **options), lo, hi
                )
                for arm in ARMS:
                    self.assertEqual(
                        scored(worker(poly, base, arm, **options), lo, hi),
                        expected,
                    )

        collector = worker(poly, base, POWER, block_width=64)
        for maximum in (2, 64, 7, 1, 31):
            expected = worker(poly, base, "sparse", block_width=64)
            with (
                patch.object(collector, "_bounds", return_value=(1, maximum)),
                patch.object(expected, "_bounds", return_value=(1, maximum)),
            ):
                self.assertEqual(
                    scored(collector, -9, 10), scored(expected, -9, 10)
                )
            self.assertTrue(
                all(
                    modulus <= maximum
                    for _, modulus in collector._presieve.levels
                )
            )

    def test_selected_powers_exist_and_period_is_bounded(self):
        base = build_factor_base(10403, multiplier=3, bound=100).factor_base
        poly = Polynomial(10403, 3, 1, 177)
        collector = worker(poly, base, POWER, block_width=4096)
        scored(collector, -11, 13)
        table = collector._presieve
        self.assertLessEqual(len(table.scores), 4096)
        self.assertTrue(
            any(
                modulus > base.primes[index] for index, modulus in table.levels
            )
        )
        self.assertEqual(table.builds, 1)
        scored(collector, -11, 13)
        self.assertEqual((table.builds, table.replays), (1, 2))

    def test_universal_and_singular_fallback_levels_are_not_replaced(self):
        cases = []
        for prime in (3, 71):
            n = prime**12
            base = FactorBase(
                n,
                1,
                prime + 1,
                (FactorBaseEntry(2, (1,)), FactorBaseEntry(prime, (0,))),
            )
            poly = Polynomial(n, 1, prime, 0)
            cases.append((poly, base, -513, 512))
        n = 2**24 + 1
        base = FactorBase(n, 1, 3, (FactorBaseEntry(2, (1,)),))
        cases.append((Polynomial(n, 1, 2, 1), base, -513, 512))
        b = 2**298 - 1
        poly = Polynomial(b * b, 1, 1, b)
        cases.append(
            (poly, build_factor_base(b * b, bound=3).factor_base, 2, 3)
        )
        for poly, base, lo, hi in cases:
            for cutoff in (0, 3):
                for bounds in ((lo, hi), (-1, 1)):
                    expected = scored(
                        worker(
                            poly,
                            base,
                            "sparse",
                            block_width=4096,
                            small_prime_cutoff=cutoff,
                        ),
                        *bounds,
                    )
                    for arm in ARMS:
                        collector = worker(
                            poly,
                            base,
                            arm,
                            block_width=4096,
                            small_prime_cutoff=cutoff,
                        )
                        self.assertEqual(scored(collector, *bounds), expected)

    def test_unrolling_preserves_every_reservation_and_refused_prefix(self):
        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 49, 8)
        traces = []
        for arm in (SELECTED, UNROLL):
            collector = worker(poly, base, arm, block_width=64)
            trace, consume = [], collector.budget.consume

            def record(amount=1):
                trace.append(amount)
                return consume(amount)

            with patch.object(collector.budget, "consume", record):
                collector.collect(-31, 38)
            traces.append(trace)
        self.assertEqual(traces[0], traces[1])
        for cut in (0, 1, 20, len(trace) // 2, len(trace) - 1):
            outcomes = []
            for arm in (SELECTED, UNROLL):
                collector = worker(poly, base, arm, block_width=64)
                collector.budget.work_limit = collector.budget.used + sum(
                    trace[:cut]
                )
                stopped = collector.collect(-31, 38)
                charged = collector.budget.used
                collector.budget = unlimited_budget()
                resumed = collector.collect(stopped.next_position, 38)
                outcomes.append((stopped, resumed, charged))
            self.assertEqual(outcomes[0], outcomes[1])

    def test_power_period_interruption_is_atomic_and_replay_is_charged(self):
        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 49, 8)
        for cut in (1, 2, 5, 11, 25, 50, 100):
            collector = worker(poly, base, BOTH, block_width=64)
            calls = 0

            def cancel():
                nonlocal calls
                calls += 1
                return calls == cut

            collector.budget = Budget(work_limit=10**9, cancelled=cancel)
            stopped = collector.collect(-31, 38)
            used = collector.budget.used
            collector.budget = unlimited_budget()
            resumed = collector.collect(stopped.next_position, 38)
            combined = signature(stopped) | signature(resumed)
            self.assertEqual(
                combined, reference_positions(poly, base, -31, 38, 97)
            )
            if cut > 1:
                self.assertGreater(used, 0)
            self.assertGreater(collector.budget.used, 0)

        collector = worker(poly, base, BOTH, block_width=64)
        scored(collector, -31, 33)
        table = collector._presieve
        old = table.key, table.scores, table.hits, table.levels
        next_poly = Polynomial(10403, 1, 7, 1)
        collector.set_polynomial(
            next_poly,
            tuple(
                polynomial_roots(next_poly, base, entry)
                for entry in base.entries
            ),
        )
        collector.budget = Budget(work_limit=64)
        stopped = collector.collect(-31, 33)
        self.assertEqual(stopped.next_position, -31)
        self.assertEqual(
            (table.key, table.scores, table.hits, table.levels), old
        )

    def test_simultaneous_memory_and_mutually_exclusive_periods(self):
        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 49, 8)
        scalar = worker(poly, base, "sparse", block_width=64)
        reserve = workspace(BOTH, 64, len(base.entries))
        selected = worker(poly, base, BOTH, block_width=64)
        self.assertEqual(selected._workspace, scalar._workspace + reserve)
        uncapped = worker(
            poly,
            base,
            BOTH,
            block_width=64,
            # Optional reuse caches may decline their allocation. Test
            # the required presieve floor independently of those caches.
            memory_bytes=(
                scalar._workspace - scalar._reuse_workspace + reserve - 1
            ),
        )

        self.assertEqual(uncapped._workspace, selected._workspace)
        with self.assertRaises(ValueError):
            config(marking="presieve+presieve_powers")
        with self.assertRaises(ValueError):
            combine(("presieve", "presieve_powers"))


class ExtendedResumeTests(unittest.TestCase):
    def test_counted_checkpoint_restores_saved_arm_and_empty_derived_caches(
        self,
    ):
        for arm in ARMS:
            cfg = CountedSIQSConfig(
                base_count=30,
                base_bound=500,
                factor_count=3,
                batch_width=1,
                collector=config(marking=arm, power_plan_bytes=2**20),
            )
            job = SIQSJob(
                1000003 * 1000033, config=cfg, budget=unlimited_budget()
            )
            first = job.run(max_blocks=1)
            self.assertIsNone(first.divisor)
            saved = json.loads(json.dumps(job.checkpoint()))
            states = []
            for _ in range(2):
                budget = unlimited_budget()
                restored = SIQSJob.from_checkpoint(saved, budget=budget)
                self.assertEqual(restored.config, cfg)
                self.assertIsNone(restored.engine.collector._presieve.key)
                self.assertFalse(restored.engine.collector._power_plans)
                self.assertGreater(
                    budget.used, saved["resources"]["work_used"]
                )
                result = restored.run(max_blocks=2)
                self.assertEqual(
                    (result.divisor or 1) * result.cofactor, job.n
                )
                states.append(
                    (
                        result.divisor,
                        result.reason,
                        restored.engine.next_position,
                    )
                )
            self.assertEqual(states[0], states[1])
            with self.assertRaises(ValueError):
                SIQSJob.from_checkpoint(
                    saved,
                    budget=unlimited_budget(),
                    config=replace(
                        cfg, collector=replace(cfg.collector, marking=SELECTED)
                    ),
                    allow_extension=True,
                )
