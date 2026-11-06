"""Bundled constants preserve the retained mathematical records."""

import importlib
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class BundledDataTests(unittest.TestCase):
    def test_every_program_matches_immutable_oracle(self):
        for name, filename in (
            ("chain_data", "c6_fast_records.json"),
            ("cf_data", "b3_cf_records.json"),
        ):
            module = importlib.import_module("v2.ecm." + name)
            original = json.loads(
                (ROOT / "benchmarks/inputs/controls" / filename).read_text()
            )
            for family, rows in original["families"].items():
                expected = tuple(
                    (
                        int(scalar),
                        bytes.fromhex(row["code"]),
                        bytes.fromhex(row["masks"]),
                        row["output"],
                        row["slots"],
                    )
                    for scalar, row in rows.items()
                )
                self.assertEqual(getattr(module, family.upper()), expected)

    def test_old_imports_resolve_to_one_implementation(self):
        for name in (
            "auto",
            "defaults",
            "presets",
            "selection",
            "stage_two",
            "worker",
            "worker_state",
            "workers",
        ):
            self.assertIs(
                importlib.import_module("v2.execution.ecm_" + name),
                importlib.import_module("v2.ecm." + name),
            )

    def test_policy_returns_owned_values(self):
        from v2.ecm.worker_policy import policy_data

        first = policy_data()
        expected = first["rows"][0]["digits"]
        first["rows"][0]["digits"] = 1
        self.assertEqual(policy_data()["rows"][0]["digits"], expected)
