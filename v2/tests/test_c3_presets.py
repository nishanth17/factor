"""Source schedule composition and preservation of protected SIQS bounds."""

import unittest
from dataclasses import replace

from v2.execution.allocation import ECMAllocation
from v2.execution.budget import Budget
from v2.execution.ecm_presets import ECM_PRESETS, with_ecm_preset
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import SIQSConfig


def base_config(**changes):
    config = PortfolioConfig(
        trial_bound=5,
        rho_attempts=0,
        pm1_attempts=0,
        max_input_bits=128,
        segment_size=16,
        memory_bytes=64 * 2**20,
        siqs=SIQSConfig(
            base_bound=100, half_width=64, memory_bytes=32 * 2**20
        ),
        allocation=ECMAllocation(
            "pretest",
            pretest_work=0,
            fallback_work=1000,
            fallback_seconds=1,
            fallback_cpu_seconds=1,
        ),
    )
    return replace(config, **changes)


class PresetTests(unittest.TestCase):
    def test_each_preset_preserves_siqs_and_cumulative_allocation(self):
        base = base_config()
        for name in ECM_PRESETS:
            config = with_ecm_preset(base, name)
            self.assertIs(config.siqs, base.siqs)
            self.assertIs(config.allocation, base.allocation)
            self.assertEqual(config.memory_bytes, base.memory_bytes)
            self.assertEqual(config.ecm_tiers, ECM_PRESETS[name])
            self.assertLessEqual(
                config.siqs.memory_bytes + config.workspace_reserve + 8192,
                config.memory_bytes,
            )
        self.assertEqual(base.ecm_tiers, ((2000, 147396, 32),))

    def test_source_tables_have_exact_integer_bounds_and_escalation(self):
        self.assertEqual(sum(t[2] for t in ECM_PRESETS["alpertron25"]), 415)
        self.assertEqual(
            sum(t[2] for t in ECM_PRESETS["yamaquasi_ecm64"]), 140
        )
        self.assertEqual(ECM_PRESETS["gmp_ecm20"], ((11000, 1900000, 74),))
        for tiers in ECM_PRESETS.values():
            self.assertTrue(
                all(type(x) is int for tier in tiers for x in tier)
            )
            self.assertTrue(
                all(2 <= b1 <= b2 and c > 0 for b1, b2, c in tiers)
            )
        with self.assertRaises(TypeError):
            ECM_PRESETS["unbounded"] = ()

    def test_missing_protection_or_unknown_preset_is_rejected(self):
        for config in (
            base_config(allocation=None),
            base_config(allocation=None, siqs=None),
        ):
            with self.assertRaisesRegex(ValueError, "protected allocation"):
                with_ecm_preset(config, "alpertron20")
        for field in (
            "fallback_work",
            "fallback_seconds",
            "fallback_cpu_seconds",
        ):
            base = base_config()
            base = replace(
                base, allocation=replace(base.allocation, **{field: 0})
            )
            with self.assertRaisesRegex(ValueError, "positive SIQS"):
                with_ecm_preset(base, "alpertron20")
        with self.assertRaisesRegex(ValueError, "unknown"):
            with_ecm_preset(base_config(), "guess")

    def test_unfunded_larger_workspace_cannot_shrink_siqs_storage(self):
        base = base_config()
        tight = replace(
            base,
            memory_bytes=base.siqs.memory_bytes
            + base.workspace_reserve
            + 8192,
        )
        with self.assertRaisesRegex(MemoryError, "coexistence"):
            with_ecm_preset(tight, "gmp_ecm25")
        self.assertEqual(tight.siqs.memory_bytes, 32 * 2**20)

    def test_large_campaign_obeys_pretest_floor_and_restores_exact_tiers(self):
        base = base_config()
        base = replace(
            base, allocation=replace(base.allocation, fallback_work=500_000)
        )
        config = with_ecm_preset(base, "alpertron25")
        first = factorize_bounded(
            1009 * 1013,
            config=config,
            budget=Budget(work_limit=1000, seconds=5, cpu_seconds=5),
        )
        self.assertEqual(first.reason, "insufficient_fallback_work")
        self.assertFalse(any(e["stage"] == "ecm" for e in first.events))
        resumed = factorize_bounded(
            1009 * 1013,
            checkpoint=first.checkpoint,
            budget=Budget(work_limit=2_000_000, seconds=5, cpu_seconds=5),
        )
        self.assertTrue(resumed.result.complete)
        self.assertEqual(resumed.result.reconstruct(), 1009 * 1013)
        saved = resumed.checkpoint["payload"]["config"]
        self.assertEqual(
            tuple(tuple(t) for t in saved["ecm_tiers"]), config.ecm_tiers
        )
        self.assertEqual(
            saved["siqs"]["memory_bytes"], config.siqs.memory_bytes
        )
