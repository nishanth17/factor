"""Independent checks for diagnostic sampling, certificates and incidence."""

import random
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.benchmarks.qs.c1.c1_feasibility import (
    ResidualAudit,
    incidence_report,
    kernel,
    load_control,
    split_residual,
    validate_control,
    verify_square,
    xor_selected,
)
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.qs import SIQSConfig, build_factor_base, qs_polynomial
from v2.qs.sieve_collector import SieveCollector, SieveConfig
from v2.tests.test_qs import unlimited_budget


class DiagnosticOracleTests(unittest.TestCase):
    def test_complete_cycle_space_including_disconnected_multigraph(self):
        # A separate triangle, a self-loop, and parallel edges all contribute.
        rows = (1, 1, 2 ^ 4, 4 ^ 8, 8 ^ 2, 16 ^ 16, 32 ^ 64, 32 ^ 64)
        expected = {
            mask
            for mask in range(1 << len(rows))
            if xor_selected(rows, mask) == 0
        }
        basis = kernel(rows, unlimited_budget())
        actual = {xor_selected(basis, mask) for mask in range(1 << len(basis))}

        self.assertEqual(actual, expected)
        self.assertEqual(len(basis), 4)
        self.assertTrue(any(mask & 0b11100 == 0b11100 for mask in actual))

    def test_random_incidence_agrees_with_exhaustive_subset_oracle(self):
        rng = random.Random(54)
        for _ in range(50):
            rows = tuple(rng.randrange(64) for _ in range(9))
            expected = {
                mask for mask in range(512) if xor_selected(rows, mask) == 0
            }
            basis = kernel(rows, unlimited_budget())
            actual = {
                xor_selected(basis, mask) for mask in range(1 << len(basis))
            }

            self.assertEqual(actual, expected)

    def test_square_corrections_and_both_gcd_signs(self):
        records = [
            dict(u=10, sign=1, square=1, residual=23, lp=(23,), exponents=()),
            dict(u=13, sign=1, square=2, residual=23, lp=(23,), exponents=()),
        ]
        self.assertEqual(verify_square(records, 3, 77), {7, 11})
        report = incidence_report(records, 77)
        self.assertEqual(report["independent_lp_constraints"], 1)
        self.assertEqual(report["dependencies"], 1)
        self.assertEqual(report["proper_divisors"], [7, 11])

        for key, value in (("u", 14), ("square", 1), ("residual", 21)):
            damaged = [records[0], dict(records[1], **{key: value})]
            with self.assertRaises(AssertionError):
                verify_square(damaged, 3, 77)

    def test_split_certainty_bounds_squares_and_finite_failure(self):
        for residual, limit, expected in (
            (101 * 103, 103, "dlp"),
            (101**2, 101, "dlp"),
            (101 * 103, 101, "endpoint_bound"),
            (101, 200, "prime"),
            (101**3, 20000, "not_two_primes"),
        ):
            endpoints, reason, evaluations = split_residual(
                residual, limit, unlimited_budget()
            )
            self.assertEqual(reason, expected)
            self.assertLessEqual(evaluations, 4096)
            if reason == "dlp":
                self.assertEqual(endpoints[0] * endpoints[1], residual)

        with self.assertRaises(BudgetExhaustedError):
            split_residual(101 * 103, 200, Budget(work_limit=0))
        with self.assertRaises(BudgetExhaustedError):
            split_residual(101 * 103, 200, Budget(cancelled=lambda: True))
        with patch(
            "v2.benchmarks.qs.c1.c1_feasibility.factorize_rho",
            return_value=None,
        ):
            self.assertEqual(
                split_residual(101 * 103, 200, unlimited_budget())[1],
                "split_failed",
            )


class DiagnosticSamplingTests(unittest.TestCase):
    def config(self):
        return SIQSConfig(
            base_bound=100,
            half_width=128,
            collector=SieveConfig(
                score_policy="powers",
                division="bucket",
                block_width=32,
                residual_bound=500,
                max_partials=4,
                memory_bytes=256 * 2**20,
            ),
        )

    def test_observer_preserves_slp_rows_eviction_and_work(self):
        config = self.config()
        base = build_factor_base(10403, bound=100).factor_base
        polynomial = qs_polynomial(base)
        reference = SieveCollector(
            polynomial,
            base,
            config=config.collector,
            budget=unlimited_budget(),
        )
        observed = SieveCollector(
            polynomial,
            base,
            config=config.collector,
            budget=unlimited_budget(),
        )
        audit = ResidualAudit(config, 7)
        original = SieveCollector._sieve

        def observer(collector, lo, hi, stats):
            threshold = original(collector, lo, hi, stats)
            audit.observe(collector, lo, hi, threshold)
            return threshold

        result = reference.collect(-128, 129)
        with patch.object(SieveCollector, "_sieve", observer):
            checked = observed.collect(-128, 129)

        self.assertEqual(result, checked)
        self.assertEqual(reference.budget.used, observed.budget.used)
        self.assertIsNone(audit.stop)
        self.assertGreater(audit.counts["uniform_positions"], 0)

    def test_reservoir_is_bounded_and_includes_late_rejections(self):
        audit = ResidualAudit(self.config(), 29)
        for value in range(10000):
            audit.sample("rejected", value)

        self.assertEqual(len(audit.samples["rejected"]), 128)
        self.assertEqual(audit.sample_seen["rejected"], 10000)
        self.assertTrue(
            any(v["residual"] > 9000 for v in audit.samples["rejected"])
        )

    def test_split_cap_and_budget_stop_are_not_empty_success(self):
        audit = ResidualAudit(self.config(), 7)
        audit.counts["split_attempts"] = 8192
        self.assertEqual(audit.classify(101 * 103)[1], "split_unexamined")
        self.assertEqual(audit.stop, "split_attempt_limit")

        base = build_factor_base(10403, bound=100).factor_base
        run = SieveCollector(
            qs_polynomial(base),
            base,
            config=self.config().collector,
            budget=unlimited_budget(),
        )
        audit = ResidualAudit(self.config(), 7)
        audit.budget = Budget(work_limit=0)
        audit.observe(run, 0, 1, 0)
        self.assertEqual(audit.stop, "diagnostic_work_limit")

    def test_control_certainty_uses_strict_a10_range(self):
        for factor, expected in (
            (2**64 + 13, "proven_prime"),
            (3317044064679887385961981, "probable_prime"),
        ):
            fixture = dict(n=factor, factors=[factor])
            row = dict(
                factors=[factor],
                remaining=[],
                complete=True,
                divisor=None,
                certainty=[expected],
                work=0,
                stats={},
            )
            validate_control(row, fixture)
            row["certainty"] = [
                "probable_prime"
                if expected == "proven_prime"
                else "proven_prime"
            ]
            with self.assertRaises(AssertionError):
                validate_control(row, fixture)

    def test_frozen_loaders_need_only_versioned_files(self):
        control, fixtures, configs = load_control()
        self.assertEqual(control["seeds"], [7, 29])
        self.assertEqual([f["digits"] for f in fixtures], [30, 40, 60])
        self.assertEqual(configs[30].half_width, 8192)
        self.assertEqual(configs[40].half_width, 65536)
        self.assertEqual(configs[60].base_bound, 100000)
        self.assertEqual(replace(configs[30]).collector.threshold_extra, 0)


if __name__ == "__main__":
    unittest.main()
