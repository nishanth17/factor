"""Completion-study controls must preserve owners and censored outcomes."""

import unittest
from pathlib import Path

from v2.benchmarks.ecm.b12_f4.completion import configuration, point_comparison
from v2.benchmarks.ecm.b12_f4.selection_child import execution_options


class CompletionStudyTests(unittest.TestCase):
    def test_serial_control_cannot_accidentally_adopt_current_defaults(self):
        self.assertEqual(execution_options(dict(arm="S0"), None),
                         {"ecm_executor": None})
        owner = object()
        self.assertIs(execution_options(dict(arm="POLICY"), owner)
                      ["ecm_executor"], owner)

    def test_cold_auto_and_retained_warm_are_explicit(self):
        request = dict(arm="AUTO", mode="cold", ecm_memory_bytes=8 * 2**30)
        self.assertEqual(execution_options(request, None),
                         dict(ecm_executor="auto", ecm_memory_bytes=8 * 2**30))
        with self.assertRaises(ValueError):
            execution_options(dict(request, mode="warm"), None)

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
        self.assertGreater(funded.siqs.memory_bytes, original.siqs.memory_bytes)

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
        self.assertLessEqual(funded.siqs.memory_bytes + funded.workspace_reserve
                             + 8192, funded.memory_bytes)

    def test_incomplete_fast_return_is_not_a_speedup(self):
        complete = dict(complete=True, seconds=100)
        incomplete = dict(complete=False, seconds=.1)
        self.assertIsNone(point_comparison(complete, incomplete)["speedup"])
        self.assertEqual(point_comparison(complete,
                         dict(complete=True, seconds=50))["speedup"], 2)
        with self.assertRaises(ValueError):
            point_comparison(complete, dict(complete=True, seconds=float("nan")))
