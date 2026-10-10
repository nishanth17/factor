"""Complete offline graph space and exact lifting checks for C1 diagnostics."""

import random
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.benchmarks.qs.c1.c1_feasibility import kernel, xor_selected
from v2.benchmarks.qs.c1.c1_followup import (
    FollowupAudit,
    campaign_decision,
    check_cycle,
    compact_report,
    fundamental_cycles,
    load_followup,
    refine_prefixes,
    validate_records,
)
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.qs import SIQSConfig, build_factor_base, qs_polynomial
from v2.qs.sieve_collector import SieveCollector, SieveConfig
from v2.tests.test_qs import unlimited_budget


class CompleteOfflineGraphTests(unittest.TestCase):
    def test_disconnected_loops_parallel_edges_and_isolated_slp(self):
        pairs = (
            (1, 101),
            (103, 107),
            (107, 109),
            (109, 103),
            (113, 113),
            (127, 131),
            (127, 131),
            (1, 1),
        )
        cycles, stats = fundamental_cycles(pairs, unlimited_budget())
        expected = {
            mask for mask in range(1 << len(pairs)) if check_cycle(pairs, mask)
        }
        actual = {
            xor_selected(cycles, mask) for mask in range(1 << len(cycles))
        }

        self.assertEqual(expected, actual)
        self.assertEqual(stats["cycle_rank"], 4)
        self.assertEqual(stats["components"], 4)
        self.assertEqual(stats["cycles_disconnected_from_slp"], 3)
        self.assertEqual(stats["edges_outside_cycles"], 1)
        self.assertIn(0b1110, actual)

    def test_random_multigraphs_against_independent_dense_elimination(self):
        rng = random.Random(548128)
        for _ in range(100):
            pairs = tuple(
                (rng.randrange(10), rng.randrange(10)) for _ in range(14)
            )
            rows = tuple((1 << a) ^ (1 << b) for a, b in pairs)
            oracle = kernel(rows, unlimited_budget())
            cycles, stats = fundamental_cycles(pairs, unlimited_budget())
            self.assertEqual(len(cycles), len(oracle))
            self.assertEqual(stats["cycle_rank"], len(oracle))
            self.assertTrue(
                all(xor_selected(rows, mask) == 0 for mask in cycles)
            )
            # Independence plus equal dimension establishes complete span.
            self.assertEqual(kernel(cycles, unlimited_budget()), ())

    def test_long_disconnected_path_is_iterative_and_bounded(self):
        pairs = tuple((i, i + 1) for i in range(1000, 3000)) + ((1000, 3000),)
        cycles, _ = fundamental_cycles(pairs, unlimited_budget())
        self.assertEqual(cycles, ((1 << len(pairs)) - 1,))
        with self.assertRaises(MemoryError):
            fundamental_cycles(pairs, unlimited_budget(), memory_bytes=1024)
        with self.assertRaises(BudgetExhaustedError):
            fundamental_cycles(pairs, Budget(work_limit=1))
        with self.assertRaises(BudgetExhaustedError):
            fundamental_cycles(pairs, Budget(cancelled=lambda: True))

    def test_filter_discovery_lifts_original_square_corrections(self):
        records = [
            dict(u=10, sign=1, square=1, residual=23, lp=(23,), exponents=()),
            dict(u=13, sign=1, square=2, residual=23, lp=(23,), exponents=()),
        ]
        report = compact_report(records, 77)
        self.assertEqual(report["proper_divisors"], [7, 11])
        self.assertEqual(report["dependencies"], 1)
        self.assertEqual(report["nontrivial_dependencies"], 1)
        self.assertEqual(report["trivial_dependencies"], 0)
        self.assertEqual(report["post_filter"]["zero_dependencies"], 1)
        self.assertEqual(report["remaining"], [])
        self.assertTrue(report["complete"])
        self.assertEqual(report["factors"], [7, 11])
        with self.assertRaises(ValueError):
            validate_records(
                [dict(records[1], square=1)], 77, unlimited_budget()
            )
        with self.assertRaises(ValueError):
            validate_records(
                [dict(records[1], exponents=((2, 10**100),))],
                77,
                unlimited_budget(),
            )
        records[0].update(polynomial="p", position=0)
        records[1].update(polynomial="p", position=0)
        with self.assertRaises(ValueError):
            validate_records(records, 77, unlimited_budget())

    def test_first_factor_cost_audit_keeps_full_arithmetic_validation(self):
        records = [
            dict(
                u=9, sign=1, square=1, residual=1, lp=(), exponents=((2, 2),)
            ),
            dict(u=10, sign=1, square=1, residual=23, lp=(23,), exponents=()),
            dict(u=13, sign=1, square=2, residual=23, lp=(23,), exponents=()),
        ]
        exhaustive = compact_report(records, 77)
        first = compact_report(records, 77, first_factor=True)

        self.assertEqual(exhaustive["dependencies"], 2)
        self.assertEqual(exhaustive["unextracted_dependencies"], 0)
        self.assertEqual(first["algebraic_dependencies"], 2)
        self.assertEqual(first["dependencies"], 1)
        self.assertEqual(first["unextracted_dependencies"], 1)
        self.assertEqual(first["factors"], [7, 11])
        self.assertEqual(first["remaining"], [])
        # A later atom is still checked even when the first row factors n.
        records[-1]["square"] = 1
        with self.assertRaises(ValueError):
            compact_report(records, 77, first_factor=True)

    def test_censoring_is_explicit_not_zero_yield(self):
        row = dict(u=10, sign=1, square=1, residual=23, lp=(23,), exponents=())
        report = compact_report([row], 77, memory_bytes=1)
        self.assertIn("censored", report)
        self.assertNotIn("dependencies", report)
        report = compact_report([row], 77, budget=Budget(work_limit=0))
        self.assertIn("censored", report)


class FollowupCensusTests(unittest.TestCase):
    def test_retained_prefix_refinement_verifies_both_sides_of_bracket(self):
        config = SIQSConfig(
            collector=SieveConfig(
                division="bucket",
                score_policy="powers",
            )
        )
        audit = FollowupAudit(config, 7, 30)
        audit.records = [
            dict(
                kind="slp",
                block=1,
                u=10,
                sign=1,
                square=1,
                residual=23,
                lp=(23,),
                exponents=(),
            ),
            dict(
                kind="slp",
                block=3,
                u=13,
                sign=1,
                square=2,
                residual=23,
                lp=(23,),
                exponents=(),
            ),
        ]
        snapshots = [
            dict(
                blocks=4,
                costs={},
                reports={
                    policy: compact_report(audit.records, 77)
                    for policy in ("slp", "64", "128")
                },
            )
        ]

        refine_prefixes(audit, snapshots, 77)

        self.assertEqual([p["blocks"] for p in snapshots], [4, 2, 3])
        self.assertFalse(snapshots[1]["reports"]["128"]["complete"])
        self.assertEqual(snapshots[1]["reports"]["128"]["remaining"], [77])
        self.assertEqual(snapshots[2]["reports"]["128"]["factors"], [7, 11])
        self.assertLessEqual(len(snapshots), 9)

    def test_investment_gate_requires_two_inputs_and_slp_shortfall(self):
        witness = dict(blocks=2048, slp_has_factor=False, charged_cpu=1.0)
        first = dict(name="50-0", decision={"64": [witness], "128": []})
        self.assertEqual(campaign_decision([first]), {})
        second = dict(name="50-1", decision={"64": [witness], "128": []})
        self.assertEqual(campaign_decision([first, second]), {"50": ["64"]})
        first["decision"]["64"] = [dict(witness, slp_has_factor=True)]
        second["decision"]["64"] = [dict(witness, slp_has_factor=True)]
        self.assertEqual(campaign_decision([first, second]), {})

    def test_nested_bounds_and_global_sampling_preserve_slp(self):
        config = SIQSConfig(
            base_bound=100,
            half_width=128,
            collector=SieveConfig(
                block_width=32,
                division="bucket",
                score_policy="powers",
                residual_bound=500,
                max_partials=4,
                memory_bytes=256 * 2**20,
            ),
        )
        base = build_factor_base(10403, bound=100).factor_base
        polynomial = qs_polynomial(base)
        plain = SieveCollector(
            polynomial,
            base,
            config=config.collector,
            budget=unlimited_budget(),
        )
        audited = SieveCollector(
            polynomial,
            base,
            config=config.collector,
            budget=unlimited_budget(),
        )
        audit = FollowupAudit(config, 7, 30)
        original = SieveCollector._sieve

        def observed(collector, lo, hi, stats):
            threshold = original(collector, lo, hi, stats)
            audit.observe(collector, lo, hi, threshold)
            return threshold

        expected = plain.collect(-128, 129)
        with patch.object(SieveCollector, "_sieve", observed):
            actual = audited.collect(-128, 129)
        self.assertEqual(expected, actual)
        self.assertEqual(plain.budget.used, audited.budget.used)
        self.assertEqual(audit.product_limit, 128 * 100**2)
        self.assertEqual(audit.inner_limit, 64 * 100**2)
        self.assertTrue(
            all(r["residual"] <= audit.product_limit for r in audit.records)
        )

    def test_geometric_sampling_is_independent_of_block_partition(self):
        config = SIQSConfig(
            collector=SieveConfig(division="bucket", score_policy="powers")
        )

        def positions(widths):
            audit = FollowupAudit(config, 7, 30)
            result, start = [], 0
            for width in widths:
                while audit.next_sample < width:
                    result.append(start + audit.next_sample)
                    audit.next_sample += 1 + audit.sample_gap()
                audit.next_sample -= width
                start += width
            return result

        self.assertEqual(
            positions([4096] * 4), positions([1, 4095, 8000, 4288])
        )
        self.assertGreater(len(positions([16384])), 20)

    def test_versioned_followup_inputs_are_distinct_training_cases(self):
        control, fixtures, configs = load_followup()
        self.assertEqual(len(fixtures), 6)
        self.assertEqual(len({f["n"] for f in fixtures}), 6)
        self.assertEqual(control["product_multipliers"], [64, 128])
        self.assertTrue(all(f["split"] == "training" for f in fixtures))
        self.assertEqual(configs[50].base_bound, 30000)
        self.assertEqual(replace(configs[50]).factor_count, 6)


if __name__ == "__main__":
    unittest.main()
