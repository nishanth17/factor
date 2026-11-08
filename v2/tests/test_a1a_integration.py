"""Selected kernels, counted policies and charged derived-state rebuilding."""

import itertools
import json
import unittest
from collections import defaultdict
from dataclasses import asdict, replace
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import (
    CountedSIQSConfig,
    Polynomial,
    SIQSJob,
    build_factor_base,
    default_siqs_config,
)
from v2.qs.collector_kernels import combine, workspace
from v2.qs.count_policy import endpoint_for_count
from v2.qs.factor_base import FactorBase, FactorBaseEntry
from v2.qs.polynomial import polynomial_roots
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_a1a_density import config, worker
from v2.tests.test_qs import reference_positions, unlimited_budget
from v2.tests.test_qs_sieve import signature

SELECTED = "density+presieve+scan_fused"
ARMS = tuple(
    combine(parts)
    for size in range(1, 4)
    for parts in itertools.combinations(
        ("density", "presieve", "scan_fused"), size
    )
)


def scored(collector, lo, hi):
    stats = defaultdict(int, threshold_min=2**63)
    threshold = collector._sieve(lo, hi, stats)
    return (
        threshold,
        list(collector._scores[: hi - lo]),
        list(collector._hits[: hi - lo]),
        stats["root_hits"],
        stats["power_hits"],
    )


class SelectedKernelTests(unittest.TestCase):
    def test_composition_against_direct_values_and_candidate_order(self):
        for k, a, b in ((1, 1, 102), (1, 7, 1), (1, 49, 8), (9, 3, 0)):
            base = build_factor_base(
                10403, multiplier=k, bound=100
            ).factor_base
            poly = Polynomial(10403, k, a, b)
            expected = reference_positions(poly, base, -71, 80, 97)

            for cap in (0, 256, 2**20):
                for width in (1, 17, 64):
                    traces = []
                    for arm in ("sparse", *ARMS):
                        collector = worker(
                            poly,
                            base,
                            arm,
                            power_plan_bytes=cap,
                            block_width=width,
                        )
                        visited, divide = [], collector._divide

                        def record(position, offset, stats):
                            visited.append(position)
                            return divide(position, offset, stats)

                        with patch.object(collector, "_divide", record):
                            result = collector.collect(-71, 80)

                        self.assertEqual(result.reason, "complete")
                        self.assertEqual(signature(result), expected, arm)
                        traces.append(visited)
                    self.assertTrue(all(t == traces[0] for t in traces))

            for policy in ("powers", "fixed"):
                collector = worker(
                    poly,
                    base,
                    SELECTED,
                    score_policy=policy,
                    block_width=151,
                )
                scored(collector, -71, 80)
                self.assertGreater(collector._presieve.builds, 0)
                for offset, x in enumerate(range(-71, 80)):
                    value = abs(
                        a * x * x + 2 * b * x + (b * b - k * 10403) // a
                    )
                    bits = sum(
                        1 << i
                        for i, p in enumerate(base.primes)
                        if value % p == 0
                    )
                    self.assertEqual(collector._hits[offset], bits)
                    score = 0
                    for p, log in zip(base.primes, collector._logs):
                        while value and value % p == 0:
                            score += log
                            value //= p
                    self.assertEqual(collector._scores[offset], score)

    def test_score_backends_rotations_tails_cutoffs_and_low_norms(self):
        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 49, 8)
        for backend, policy, cutoff, bounds in itertools.product(
            ("list", "array", "bytearray"),
            ("adaptive", "conservative", "powers", "fixed"),
            (0, 5, 17),
            ((-64, -1), (-7, 18), (18, 19), (-31, -14)),
        ):
            options = dict(
                score_backend=backend,
                score_policy=policy,
                small_prime_cutoff=cutoff,
                block_width=64,
            )
            expected = scored(worker(poly, base, "sparse", **options), *bounds)
            self.assertEqual(
                scored(worker(poly, base, SELECTED, **options), *bounds),
                expected,
            )

        # Period selection must omit primes above the current norm bound;
        # zero/one-hit decisions use the actual tail, not configured width.
        collector = worker(poly, base, SELECTED, block_width=64)
        with patch.object(collector, "_bounds", return_value=(1, 2)):
            scored(collector, -1, 0)
        self.assertTrue(
            all(base.primes[i] <= 2 for i in collector._presieve.indices)
        )
        options = dict(block_width=17, small_prime_cutoff=17)
        self.assertEqual(
            scored(worker(poly, base, SELECTED, **options), -31, -14),
            scored(worker(poly, base, "sparse", **options), -31, -14),
        )

    def test_presieve_rebuild_cancellation_and_simultaneous_memory(self):
        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 49, 8)
        collector = worker(poly, base, SELECTED, block_width=64)
        collector.collect(-31, 38)
        table = collector._presieve
        old = table.key, table.scores, table.hits, table.indices
        next_poly = Polynomial(10403, 1, 7, 1)
        collector.set_polynomial(
            next_poly,
            tuple(polynomial_roots(next_poly, base, e) for e in base.entries),
            retain_relations=False,
        )
        calls = 0

        def cancel():
            nonlocal calls
            calls += 1
            return calls == 3

        collector.budget = Budget(work_limit=10**9, cancelled=cancel)
        stopped = collector.collect(-31, 38)
        self.assertEqual(stopped.next_position, -31)
        self.assertEqual(
            (table.key, table.scores, table.hits, table.indices), old
        )
        self.assertGreater(collector.budget.used, 0)

        collector.budget = unlimited_budget()
        resumed = collector.collect(stopped.next_position, 38)
        self.assertEqual(
            signature(resumed),
            reference_positions(next_poly, base, -31, 38, 97),
        )
        self.assertNotEqual(table.key, old[0])
        fresh = worker(next_poly, base, SELECTED, block_width=64)
        self.assertIsNot(fresh._presieve, table)
        self.assertIsNone(fresh._presieve.key)
        reserve = workspace(SELECTED, 64, len(base.entries))
        scalar = worker(poly, base, "sparse", block_width=64)
        self.assertEqual(fresh._workspace - scalar._workspace, reserve)
        uncapped = worker(
            poly,
            base,
            SELECTED,
            block_width=64,
            # Optional reuse caches may decline their allocation. Test
            # the required presieve floor independently of those caches.
            memory_bytes=(
                scalar._workspace - scalar._reuse_workspace + reserve - 1
            ),
        )
        self.assertEqual(uncapped._workspace, fresh._workspace)

    def test_selector_rejects_unported_or_noncanonical_support(self):
        for marking in (
            "inverse",
            "paired",
            "trusted",
            "presieve+density",
            "density+density",
            "scan_fused+",
            "",
        ):
            with self.assertRaises(ValueError):
                config(marking=marking)

    def test_large_universal_prime_retains_scalar_fallback_weight(self):
        base = FactorBase(
            71**4,
            1,
            100,
            (FactorBaseEntry(2, (1,)), FactorBaseEntry(71, (0,))),
        )
        poly = Polynomial(base.n, 1, 71, 0)
        selected = worker(poly, base, SELECTED, block_width=256)
        self.assertTrue(selected._roots[1].all_positions)
        self.assertEqual(
            scored(selected, -9, 10),
            scored(worker(poly, base, "sparse", block_width=256), -9, 10),
        )
        self.assertNotIn(1, selected._presieve.indices)
        self.assertTrue(all(bits & 2 for bits in selected._hits[:19]))

    def test_composed_saturation_and_interrupted_period_rotation(self):
        b = 2**298 - 1
        poly = Polynomial(b * b, 1, 1, b)
        base = build_factor_base(b * b, bound=3).factor_base
        for backend in ("list", "array", "bytearray"):
            for policy in ("powers", "fixed"):
                result = worker(
                    poly,
                    base,
                    SELECTED,
                    score_backend=backend,
                    score_policy=policy,
                ).collect(2, 3)
                self.assertEqual(result.atoms[0].exponents, ((2, 300),))

        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 49, 8)
        collector = worker(poly, base, SELECTED, block_width=64)
        scored(collector, -31, 33)
        table = collector._presieve
        old = table.key, table.scores, table.hits
        calls = 0

        def cancel_rotation():
            nonlocal calls
            calls += 1
            return calls == 2

        collector.budget = Budget(work_limit=10**9, cancelled=cancel_rotation)
        stopped = collector.collect(-9, 8)
        self.assertEqual(stopped.next_position, -9)
        self.assertEqual(collector.budget.used, 17)
        self.assertEqual((table.key, table.scores, table.hits), old)
        collector.budget = unlimited_budget()
        complete = collector.collect(stopped.next_position, 8)
        self.assertEqual(
            signature(complete), reference_positions(poly, base, -9, 8, 97)
        )


class CountedIntegrationTests(unittest.TestCase):
    def test_counted_active_resume_rebuilds_without_serialized_trust(self):
        cfg = CountedSIQSConfig(
            base_count=30,
            base_bound=500,
            factor_count=3,
            batch_width=1,
            collector=config(marking=SELECTED, power_plan_bytes=2**20),
        )
        job = SIQSJob(
            1000003 * 1000033, seed=7, config=cfg, budget=unlimited_budget()
        )
        first = job.run(max_blocks=1)
        self.assertIsNone(first.divisor)
        self.assertGreater(job.engine.collector._presieve.replays, 0)
        checkpoint = json.loads(json.dumps(job.checkpoint()))
        self.assertNotIn(
            "_kernels", json.loads(checkpoint["blob"])["config"]["collector"]
        )
        states = []
        for _ in range(2):
            budget = unlimited_budget()
            restored = SIQSJob.from_checkpoint(checkpoint, budget=budget)
            self.assertEqual(restored.config, cfg)
            self.assertEqual(restored.base, job.base)
            self.assertIsNone(restored.engine.collector._presieve.key)
            self.assertFalse(restored.engine.collector._power_plans)
            self.assertGreater(
                budget.used, checkpoint["resources"]["work_used"]
            )
            result = restored.run(max_blocks=2)
            self.assertEqual((result.divisor or 1) * result.cofactor, job.n)
            states.append(
                (result.divisor, result.reason, restored.engine.next_position)
            )
        self.assertEqual(states[0], states[1])

        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                checkpoint,
                budget=unlimited_budget(),
                config=replace(
                    cfg, collector=replace(cfg.collector, marking="sparse")
                ),
                allow_extension=True,
            )

    def test_count_resolution_exact_arithmetic_and_refused_rebuild(self):
        n = 1009 * 1013
        primes = [p for p in range(2, 100) if all(p % d for d in range(2, p))]
        for multiplier in (1, 3, 9, 15):
            eligible = [
                p
                for p in primes
                if any(x * x % p == n * multiplier % p for x in range(p))
            ]
            for count in (1, 3, len(eligible)):
                bound, divisor = endpoint_for_count(
                    n, multiplier, count, unlimited_budget(), ceiling=100
                )
                self.assertEqual(
                    (bound, divisor), (eligible[count - 1] + 1, None)
                )

        cfg = CountedSIQSConfig(
            base_count=10, base_bound=500, collector=config(marking=SELECTED)
        )
        job = SIQSJob(n, config=cfg, budget=Budget(work_limit=0))
        result = job.run()
        self.assertEqual(result.reason, "cancelled")
        checkpoint = json.loads(json.dumps(job.checkpoint()))
        resumed = SIQSJob.from_checkpoint(
            checkpoint, budget=unlimited_budget()
        )
        self.assertEqual(resumed.config, cfg)
        complete = resumed.run()
        self.assertEqual((complete.divisor or 1) * complete.cofactor, n)

    def test_scoped_defaults_explicit_configs_and_saved_jobs(self):
        for digits in (40, *range(50, 61)):
            for memory in (None, 128 * 2**20, 2**30):
                cfg = adapt_resource_configuration(
                    default_siqs_config,
                    10 ** (digits - 1),
                    recipe="calibrated"
                    if memory in (None, 2**30)
                    else "reference",
                )
                selected = 50 <= digits <= 55 and memory in (None, 2**30)
                self.assertEqual(
                    cfg.collector.marking,
                    SELECTED + "+parity_two" if selected else "parity_two",
                )
                self.assertNotIn("_kernels", asdict(cfg.collector))
        n = 10**54 + 123
        explicit = default_siqs_config(n)
        self.assertIs(SIQSJob(n, config=explicit).config, explicit)
        self.assertEqual(PortfolioConfig().for_input(n).siqs, explicit)
        first = factorize_bounded(n, budget=Budget(work_limit=0))
        saved = first.checkpoint["payload"]["config"]
        resumed = factorize_bounded(
            n, checkpoint=first.checkpoint, budget=Budget(work_limit=0)
        )
        self.assertEqual(resumed.checkpoint["payload"]["config"], saved)
        self.assertEqual(resumed.result.reconstruct(), n)
