"""Owned SLP loading, fixed arm controls and complete outcome validation."""

import tempfile
import unittest
from dataclasses import asdict

from v2.benchmarks.c1_implementation import (
    ARMS,
    arm_config,
    decode_config,
    load_slp,
    run_one,
    unstable,
)
from v2.qs import SieveConfig, SIQSConfig


class ImplementationHarnessTests(unittest.TestCase):
    def test_owned_control_and_current_dlp_reconstruct_same_input(self):
        fixture = dict(
            id="tiny", kind="balanced", digits=2, n=91, factors=[7, 13]
        )
        config = SIQSConfig(
            mode="qs",
            base_bound=5,
            half_width=1,
            batch_width=1,
            collector=SieveConfig(block_width=1, residual_bound=25),
        )
        with tempfile.TemporaryDirectory(prefix="c1-control-test-") as root:
            baseline = load_slp(root)
            for arm in ("slp", "dlp"):
                row = run_one(
                    fixture, 7, arm_config(config, arm), arm, 3, baseline
                )
                self.assertEqual(row["factors"], [7, 13])
                self.assertEqual(row["remaining"], [])
                self.assertTrue(row["complete"])
                self.assertGreater(row["seconds"], 0)

    def test_fixed_controls_and_stability(self):
        config = SIQSConfig(
            base_bound=200, collector=SieveConfig(residual_bound=40000)
        )
        for arm in ARMS:
            candidate = arm_config(config, arm)
            self.assertEqual(decode_config(asdict(candidate)), candidate)
            self.assertEqual(candidate.memory_bytes, config.memory_bytes)
            self.assertEqual(
                candidate.collector.max_atoms, config.collector.max_atoms
            )
        self.assertEqual(arm_config(config, "dlp_half").base_bound, 100)
        self.assertEqual(
            arm_config(config, "graph_slp").collector.large_product_bound, 4
        )
        self.assertFalse(unstable([dict(seconds=1) for _ in range(9)]))
        self.assertTrue(
            unstable([dict(seconds=t) for t in (1, 1, 1, 2, 2, 2, 3, 3, 3)])
        )


if __name__ == "__main__":
    unittest.main()
