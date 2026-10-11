"""Default routing, finite admission and numerical checkpoint identity."""

import json
import subprocess
import sys
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path

from v2.execution.allocation import ECMAllocation
from v2.execution.budget import Budget
from v2.execution.ecm_defaults import COMPACT64, ESCALATING140
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import SIQSConfig

LEGACY32 = ((2000, 147396, 32),)


def allowance(work=0):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


class DefaultTests(unittest.TestCase):
    def test_exact_observable_size_boundaries_and_sign(self):
        for n, tiers in (
            (10**29 - 1, LEGACY32),
            (10**29, COMPACT64),
            (10**30 - 1, COMPACT64),
            (10**30, LEGACY32),
            (10**39 - 1, LEGACY32),
            (10**39, ESCALATING140),
            (10**40 - 1, ESCALATING140),
            (10**40, LEGACY32),
        ):
            for signed in (n, -n):
                self.assertEqual(
                    PortfolioConfig().for_input(signed).ecm_tiers, tiers
                )

    def test_explicit_empty_custom_and_old_default_tiers_take_precedence(self):
        for tiers in ((), LEGACY32, ((50, 1000, 7),)):
            config = PortfolioConfig(ecm_tiers=tiers)
            self.assertIs(config.for_input(10**39), config)
            self.assertEqual(config.ecm_tiers, tiers)

    def test_resolved_plan_is_frozen_for_children_and_serialization(self):
        config = PortfolioConfig().for_input(10**39)
        self.assertIs(config.for_input(10**29), config)
        values = asdict(config)
        self.assertNotIn("_automatic_ecm_tiers", values)
        restored = PortfolioConfig(**values)
        self.assertEqual(restored.ecm_tiers, ESCALATING140)
        self.assertIs(restored.for_input(10**29), restored)

    def test_cli_headroom_preserves_all_siqs_settings_and_allocation(self):
        siqs = SIQSConfig(memory_bytes=64 * 2**20)
        policy = ECMAllocation(
            "campaign",
            fallback_work=1000,
            fallback_seconds=1,
            fallback_cpu_seconds=1,
        )
        original = PortfolioConfig(
            memory_bytes=80 * 2**20, siqs=siqs, allocation=policy
        )
        selected = original.for_input(10**39)
        self.assertEqual(selected.ecm_tiers, ESCALATING140)
        self.assertIs(selected.siqs, siqs)
        self.assertIs(selected.allocation, policy)
        for field in asdict(original):
            if field != "ecm_tiers":
                self.assertEqual(
                    getattr(selected, field), getattr(original, field)
                )
        self.assertLessEqual(
            siqs.memory_bytes + selected.workspace_reserve + 8192,
            selected.memory_bytes,
        )

    def test_tight_memory_retains_old_tiers_without_shrinking_siqs(self):
        siqs = SIQSConfig(memory_bytes=64 * 2**20)
        base = PortfolioConfig(memory_bytes=80 * 2**20, siqs=siqs)
        tight = PortfolioConfig(
            memory_bytes=siqs.memory_bytes + base.workspace_reserve + 8192,
            siqs=siqs,
            ecm_chain_mode=base.ecm_chain_mode,
            ecm_program_bytes=base.ecm_program_bytes,
            ecm_chain_bytes=base.ecm_chain_bytes,
        )
        selected = tight.for_input(10**39)
        self.assertEqual(selected.ecm_tiers, LEGACY32)
        self.assertEqual(selected.memory_bytes, tight.memory_bytes)
        self.assertIs(selected.siqs, siqs)

    def test_new_implicit_checkpoint_restores_exact_numerical_schedule(self):
        for n, tiers in ((10**29 + 1, COMPACT64), (10**39 + 1, ESCALATING140)):
            first = factorize_bounded(n, budget=allowance())
            self.assertEqual(first.work_used, 0)
            self.assertEqual(first.result.reconstruct(), n)
            resumed = factorize_bounded(
                n, checkpoint=first.checkpoint, budget=allowance(1000)
            )
            self.assertEqual(resumed.result.reconstruct(), n)
            saved = resumed.checkpoint["payload"]["config"]
            self.assertEqual(tuple(map(tuple, saved["ecm_tiers"])), tiers)
            self.assertLessEqual(resumed.work_used, 1000)
            with self.assertRaisesRegex(ValueError, "incompatible"):
                factorize_bounded(
                    n,
                    config=PortfolioConfig(ecm_tiers=LEGACY32),
                    checkpoint=first.checkpoint,
                    budget=allowance(1000),
                )

    def test_old_fixed32_checkpoint_does_not_adopt_new_default(self):
        n = 10**39 + 1
        first = factorize_bounded(
            n, config=PortfolioConfig(ecm_tiers=LEGACY32), budget=allowance()
        )
        for config in (None, PortfolioConfig()):
            resumed = factorize_bounded(
                n,
                config=config,
                checkpoint=first.checkpoint,
                budget=allowance(1000),
            )
            saved = resumed.checkpoint["payload"]["config"]
            self.assertEqual(tuple(map(tuple, saved["ecm_tiers"])), LEGACY32)
            self.assertEqual(resumed.result.reconstruct(), n)

    def test_protected_fallback_refuses_without_spending_reserve(self):
        n = 10**39 + 1
        config = PortfolioConfig(
            memory_bytes=80 * 2**20,
            siqs=SIQSConfig(memory_bytes=64 * 2**20),
            allocation=ECMAllocation("campaign", fallback_work=1_000_000),
        )
        # Fund mandatory checks, but leave the fallback floor unfunded.
        run = factorize_bounded(n, config=config, budget=allowance(100_000))
        self.assertEqual(run.reason, "insufficient_fallback_work")
        self.assertEqual(run.result.reconstruct(), n)
        self.assertFalse(any(e["stage"] == "ecm" for e in run.events))
        saved = run.checkpoint["payload"]["config"]
        self.assertEqual(tuple(map(tuple, saved["ecm_tiers"])), ESCALATING140)
        self.assertEqual(
            saved["siqs"]["memory_bytes"], config.siqs.memory_bytes
        )

    def test_cli_omitted_tiers_resume_and_explicit_curves(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "state.json"
            command = [sys.executable, "-B", "-m", "v2.factor"]
            for flags, expected in (
                ([], ESCALATING140),
                (["--ecm-curves", "32"], LEGACY32),
            ):
                first = subprocess.run(
                    command
                    + [
                        str(10**39 + 1),
                        "--bounded",
                        "--work-limit",
                        "0",
                        "--checkpoint",
                        str(path),
                    ]
                    + flags,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(first.returncode, 1, first.stderr)
                self.assertTrue(path.exists(), first.stderr)
                saved = json.loads(path.read_text())["payload"]["config"]
                self.assertEqual(
                    tuple(map(tuple, saved["ecm_tiers"])), expected
                )
                resumed = subprocess.run(
                    command
                    + [
                        "--resume",
                        str(path),
                        "--work-limit",
                        "1000",
                        "--checkpoint",
                        str(path),
                    ],
                    capture_output=True,
                    text=True,
                )
                self.assertIn(resumed.returncode, (0, 1), resumed.stderr)
                self.assertNotIn("incompatible", resumed.stderr)
                saved = json.loads(path.read_text())["payload"]["config"]
                self.assertEqual(
                    tuple(map(tuple, saved["ecm_tiers"])), expected
                )

    def test_direct_siqs_request_keeps_empty_ecm_schedule(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "state.json"
            run = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    "-m",
                    "v2.factor",
                    str(10**39 + 1),
                    "--method",
                    "siqs",
                    "--work-limit",
                    "0",
                    "--checkpoint",
                    str(path),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(run.returncode, 1, run.stderr)
            saved = json.loads(path.read_text())["payload"]["config"]
            self.assertEqual(saved["ecm_tiers"], [])
