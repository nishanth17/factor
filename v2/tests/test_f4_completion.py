"""Completion-study controls must preserve owners and censored outcomes."""

import unittest
from pathlib import Path

from v2.benchmarks.ecm.b12_f4.completion import (
    configuration,
    point_comparison,
    validate_known_factorization,
)
from v2.benchmarks.ecm.b12_f4.selection_child import (
    execution_options,
    selection_options,
    warmup_validation,
)


class CompletionStudyTests(unittest.TestCase):
    def test_serial_control_cannot_accidentally_adopt_current_defaults(self):
        self.assertEqual(
            execution_options(dict(arm="S0"), None), {"ecm_executor": None}
        )
        owner = object()
        self.assertIs(
            execution_options(dict(arm="POLICY"), owner)["ecm_executor"], owner
        )

    def test_legacy_master_signature_receives_no_candidate_keyword(self):
        self.assertEqual(
            execution_options(dict(arm="S0"), None, supports_executor=False),
            {},
        )

    def test_cold_auto_and_retained_warm_are_explicit(self):
        request = dict(arm="AUTO", mode="cold", ecm_memory_bytes=8 * 2**30)
        self.assertEqual(
            execution_options(request, None),
            dict(ecm_executor="auto", ecm_memory_bytes=8 * 2**30),
        )
        with self.assertRaises(ValueError):
            execution_options(dict(request, mode="warm"), None)

    def test_retained_context_owns_the_funded_parent_without_changing_legacy(
        self,
    ):
        self.assertEqual(
            selection_options(dict(seed=7)), dict(seed=7, prefer_warm=False)
        )
        values = selection_options(
            dict(
                seed=7,
                prefer_warm=True,
                parent_memory=2 * 2**30 + 16 * 2**20,
                ecm_memory_bytes=8 * 2**30,
            )
        )
        self.assertEqual(values["parent_memory"], 2 * 2**30 + 16 * 2**20)
        self.assertEqual(values["memory_bytes"], 8 * 2**30)
        self.assertTrue(values["prefer_warm"])

    def test_resource_variants_keep_accepted_sixty_digit_policies(self):
        from v2.portfolio import PortfolioConfig

        n = 10**59 + 123
        root = Path(__file__).resolve().parents[2]
        original = PortfolioConfig().for_input(n)
        funded = configuration(n, root)
        self.assertEqual(original.ecm_tiers, funded.ecm_tiers)
        self.assertEqual(funded.siqs.solve_policy, original.siqs.solve_policy)
        self.assertEqual(funded.siqs.collector, original.siqs.collector)
        self.assertEqual(funded.siqs.base_count, original.siqs.base_count)
        self.assertEqual(funded.siqs.checkpoint_bytes, 64 * 2**20)
        self.assertGreater(
            funded.siqs.memory_bytes, original.siqs.memory_bytes
        )

    def test_stock_configuration_preserves_actual_default_resource_policy(
        self,
    ):
        from dataclasses import asdict

        from v2.portfolio import PortfolioConfig

        n = 10**59 + 123
        root = Path(__file__).resolve().parents[2]
        actual = configuration(n, root, funded=False)
        self.assertEqual(
            asdict(actual), asdict(PortfolioConfig().for_input(n))
        )
        self.assertEqual(actual.siqs.memory_bytes, 2**30)
        self.assertEqual(actual.siqs.checkpoint_bytes, 2**20)
        self.assertIsNone(configuration(10**69 + 123, root, funded=False).siqs)

    def test_seventy_digit_fallback_is_an_explicit_common_configuration(self):
        from v2.portfolio import PortfolioConfig

        n = 10**69 + 123
        root = Path(__file__).resolve().parents[2]
        original = PortfolioConfig().for_input(n)
        funded = configuration(n, root)
        self.assertIsNone(original.siqs)
        self.assertEqual(original.ecm_tiers, funded.ecm_tiers)
        self.assertEqual(funded.siqs.base_count, 10500)
        self.assertEqual(funded.siqs.multiplier_policy, "parity2_wide")
        self.assertLessEqual(
            funded.siqs.memory_bytes + funded.workspace_reserve + 8192,
            funded.memory_bytes,
        )

    def test_warmup_oracle_detects_false_prime_despite_exact_reconstruction(
        self,
    ):
        row = dict(
            n=221,
            endpoint="complete",
            reconstruct=221,
            factorization_complete=True,
            complete=True,
            factors=[[13, 1, "proven_prime"], [17, 1, "proven_prime"]],
            remaining=[],
            reason="complete",
            seconds=1,
            work=100,
            events=object(),
        )
        record = warmup_validation(row)
        self.assertNotIn("events", record)
        validate_known_factorization(record, {13: 1, 17: 1})

        # Reconstruction alone accepts a false composite-as-prime result.
        # The independently certified oracle must still see and reject it.
        false_prime = dict(record, factors=[[221, 1, "probable_prime"]])
        with self.assertRaisesRegex(ValueError, "independent factor oracle"):
            validate_known_factorization(false_prime, {13: 1, 17: 1})
        partial = dict(
            record,
            factors=[[13, 1, "proven_prime"]],
            remaining=[17],
            factorization_complete=False,
            complete=False,
            reason="wall_limit",
        )
        validate_known_factorization(partial, {13: 1, 17: 1})

    def test_incomplete_fast_return_is_not_a_speedup(self):
        complete = dict(complete=True, seconds=100)
        incomplete = dict(complete=False, seconds=0.1)
        self.assertIsNone(point_comparison(complete, incomplete)["speedup"])
        self.assertEqual(
            point_comparison(complete, dict(complete=True, seconds=50))[
                "speedup"
            ],
            2,
        )
        with self.assertRaises(ValueError):
            point_comparison(
                complete, dict(complete=True, seconds=float("nan"))
            )
