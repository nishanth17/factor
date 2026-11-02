"""Independent lower-band design and literal-policy checks, without timing."""

import copy
import importlib
import json
import sys
import tempfile
import unittest
from dataclasses import asdict, replace

from v2.benchmarks.qs.c10a01 import control_v4, lower_plan_v4 as plan
from v2.benchmarks.suites.build_phase_two_corpus import verify_certificates
from v2.qs import CountedSIQSConfig, GrowingDoubleLargeSieveConfig, SIQSConfig
from v2.qs.collector_kernels import workspace


class LowerDesignTests(unittest.TestCase):
    def test_budgets_inputs_and_unchanged_parent(self):
        data = plan.build_plan()
        self.assertTrue(plan.verify(data))
        self.assertEqual(len(data["calls"]), 20)
        self.assertEqual(
            sum(c["reservation_seconds"] for c in data["calls"]), 2860
        )
        self.assertEqual(
            6 * 125 + 6 * 305 + 8 * 35 + 500 + 120 + 120 + 100, 3700
        )
        self.assertEqual(38600 + 3700, 42300)
        self.assertEqual(3124 + 3400 + 31600 + 3700, 41824)
        self.assertEqual(
            data["work"], sum(p["work"] for p in data["phases"].values())
        )
        self.assertEqual(len(data["inputs"]), 7)
        self.assertEqual(len({r["n"] for r in data["inputs"]}), 7)
        self.assertEqual(
            [r["digits"] for r in data["pending_generation"]], [52, 56]
        )
        for digits in (50, 55):
            cells = [
                c
                for c in data["calls"]
                if c["kind"] == "complete" and c["digits"] == digits
            ]
            self.assertEqual(len(cells), 6)
            self.assertEqual(len({c["fixture"] for c in cells}), 3)
        for name in data["corpora"]:
            corpus = json.loads((plan.CORPORA / name).read_text())
            verify_certificates(corpus["certificates"])
        bad = copy.deepcopy(data)
        bad["calls"].pop()
        with self.assertRaises(ValueError):
            plan.verify(bad)
        self.assertFalse(data["execution_ready"])

    def test_current_default_policies_and_explicit_capacity_differences(self):
        with tempfile.TemporaryDirectory() as directory:
            module, _, _ = control_v4.load_control(directory)
            try:
                defaults = importlib.import_module(
                    control_v4.PACKAGE + ".qs.defaults"
                )
                for digits in (50, 51, 52, 55, 56):
                    values = plan.configuration(digits, "control")
                    cfg = defaults.default_siqs_config(10 ** (digits - 1) + 1)
                    collector = replace(
                        cfg.collector,
                        max_relations=12000,
                        max_atoms=65536,
                        max_partials=32768,
                        split_call_limit=1048576,
                    )
                    expected = replace(
                        cfg,
                        collector=collector,
                        memory_bytes=values["memory_bytes"],
                    )
                    self.assertEqual(asdict(expected), values)
                    typ = (
                        module.CountedSIQSConfig
                        if "base_count" in values
                        else module.SIQSConfig
                    )
                    values["collector"] = module.GrowingDoubleLargeSieveConfig(
                        **values["collector"]
                    )
                    typ(**values)
            finally:
                for name in tuple(sys.modules):
                    if name == control_v4.PACKAGE or name.startswith(
                        control_v4.PACKAGE + "."
                    ):
                        del sys.modules[name]

    def test_nominee_retains_lower_marking_and_extra_scratch_bound(self):
        for digits in (50, 51, 52, 55, 56):
            values = plan.configuration(digits, "nominee")
            values["collector"] = GrowingDoubleLargeSieveConfig(
                **values["collector"]
            )
            typ = CountedSIQSConfig if "base_count" in values else SIQSConfig
            cfg = typ(**values)
            self.assertEqual(cfg.multiplier, 0)
            self.assertIn("parity_two", cfg.collector.marking)
            if digits <= 55:
                self.assertIn("presieve", cfg.collector.marking)
                self.assertIn("scan_fused", cfg.collector.marking)
            self.assertLessEqual(
                workspace(cfg.collector.marking, 4096, 5000), 32 * 2**20
            )
        with self.assertRaises(ValueError):
            plan.configuration(54, "control")
