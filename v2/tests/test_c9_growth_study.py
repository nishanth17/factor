"""C9.2 causal controls and charged interruption contracts."""

import json
import unittest
from dataclasses import asdict
from pathlib import Path
from tempfile import TemporaryDirectory

from v2.benchmarks.qs.c9 import c9_growth_confirmation as confirmation
from v2.benchmarks.qs.c9 import c9_growth_study as study


class GrowthProtocolTests(unittest.TestCase):
    def test_growth_preserves_frozen_candidate_and_polynomial_domains(self):
        original = json.loads(
            (study.INPUTS / "controls/c9_training.json").read_text()
        )["configurations"]
        configs = {a: study.configuration(original, a) for a in study.ARMS}
        fixed, grown = (asdict(configs[a]) for a in ("dlp", "grow"))
        for name in (
            "max_atoms",
            "max_partials",
            "max_relations",
            "split_call_limit",
            "memory_bytes",
        ):
            self.assertIsNone(grown["collector"][name])
            grown["collector"][name] = fixed["collector"][name]
        grown["collector"].pop("storage_policy")
        self.assertEqual(fixed, grown)
        aligned = asdict(configs["grow_aligned"])
        self.assertEqual(aligned["collector"]["residual_bound"], 3000000)
        aligned["collector"]["residual_bound"] = 900000000
        self.assertEqual(aligned, asdict(configs["grow"]))
        for config in configs.values():
            self.assertEqual(config.memory_bytes, study.MEMORY)
            self.assertEqual(study.decode(asdict(config)), config)

    def test_validator_uses_declared_memory_and_rejects_bad_arithmetic(self):
        fixture = dict(n=21, factors=[3, 7])
        row = dict(
            factors=[3, 7],
            remaining=[],
            complete=True,
            divisor=3,
            certainty=["proven_prime", "proven_prime"],
            work=100,
            stats={"workspace_bytes": 700 * 2**20},
        )
        study.validate_result(row, fixture, 1024 * 2**20)
        with self.assertRaises(AssertionError):
            study.validate_result(row, fixture, 256 * 2**20)
        row["divisor"] = 1
        with self.assertRaises(AssertionError):
            study.validate_result(row, fixture, 1024 * 2**20)
        row["divisor"] = 3
        row["certainty"][0] = "probable_prime"
        with self.assertRaises(AssertionError):
            study.validate_result(row, fixture, 1024 * 2**20)
        row.update(
            factors=[],
            remaining=[21],
            complete=False,
            divisor=None,
            certainty=[],
        )
        study.validate_result(row, fixture, 1024 * 2**20)

    def test_expense_restart_keeps_failed_reservation(self):
        with TemporaryDirectory() as directory:
            output = Path(directory)
            selected = {
                "spent_before_confirmation": {"wall": 14000, "cpu": 14000}
            }
            expense = confirmation.Expenses(
                output, selected, {"generation_seconds": 1}
            )

            def interrupted():
                raise KeyboardInterrupt()

            with self.assertRaises(KeyboardInterrupt):
                expense.call("interrupted", 300, interrupted)
            loaded = confirmation.Expenses(
                output, selected, {"generation_seconds": 1}
            )
            self.assertGreaterEqual(loaded.used(), 14301)
            with self.assertRaisesRegex(RuntimeError, "aggregate"):
                loaded.call("forbidden", 100, lambda: None)
            self.assertEqual(len(loaded.data["calls"]), 1)
            self.assertEqual(loaded.data["calls"][0]["status"], "failed")

    def test_no_reveal_or_reselection_during_decode(self):
        self.assertNotEqual(confirmation.GENERATION_SEED, 202610100909)
        self.assertEqual(confirmation.CAPS[50], 360)
        self.assertEqual(confirmation.SEEDS, (7, 29, 47))
        self.assertEqual(study.ARMS, ("slp", "dlp", "grow", "grow_aligned"))


if __name__ == "__main__":
    unittest.main()
