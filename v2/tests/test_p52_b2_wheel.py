"""Freeze independent follow-up inputs and isolate all benchmark controls."""

import unittest
from unittest.mock import patch

from v2 import ecm
from v2.benchmarks import p52_b2, p52_b2_wheel


class WheelProtocolTests(unittest.TestCase):
    def test_fresh_proofs_frozen_controls_and_admitted_grid(self):
        protocol, corpus = p52_b2_wheel.inputs()
        _, previous = p52_b2.inputs()
        numbers = {fixture["n"] for fixture in corpus["fixtures"]}
        self.assertEqual(len(numbers), len(corpus["fixtures"]))
        self.assertFalse(
            numbers & {fixture["n"] for fixture in previous["fixtures"]}
        )
        self.assertEqual(protocol["sampling"], [[3, 9], [5, 31], [8, 63]])
        for case in protocol["cases"]:
            for arm in protocol["training_arms"]:
                engine = p52_b2_wheel.engine_for(arm)
                config = engine.PortfolioConfig(
                    **p52_b2_wheel.options(protocol, case, arm)
                )
                self.assertLess(config.workspace_reserve, config.memory_bytes)
        for arm in ("streamed", "programs", "legacy"):
            engine = p52_b2_wheel.engine_for(arm)
            self.assertNotIn(
                "ecm_pair_wheel", engine.PortfolioConfig.__dataclass_fields__
            )
            scalar = engine.advance_job.__globals__["ecm"].scalar_multiply
            expected = scalar(97, 3, 1, 1009, 2)
            with patch.object(
                ecm,
                "scalar_multiply",
                side_effect=AssertionError("active engine"),
            ):
                self.assertEqual(scalar(97, 3, 1, 1009, 2), expected)
