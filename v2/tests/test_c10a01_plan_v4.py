"""Design accounting and isolated-control checks without experimental calls."""

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
    recording_budget,
)
from v2.benchmarks.infrastructure.resource_control import (
    historical_asdict as asdict,
)
from v2.benchmarks.qs.c10a01 import control_v4, plan_v4
from v2.benchmarks.qs.c10a01.configuration_v4 import configuration
from v2.qs import CountedSIQSConfig, GrowingDoubleLargeSieveConfig, SIQSConfig


class ExperimentDesignTests(unittest.TestCase):
    def test_independent_call_and_budget_census(self):
        plan = plan_v4.build_plan()
        self.assertTrue(plan_v4.verify_plan(plan))
        rows = plan["calibration"]["calls"]
        phases = {}
        for row in rows:
            phases.setdefault(row["phase"], []).append(row)
        self.assertEqual(
            {k: len(v) for k, v in phases.items()},
            dict(
                primary=16,
                roots=6,
                population=4,
                ranking=4,
                forced=4,
                observer=4,
            ),
        )
        self.assertEqual(sum(r["seconds"] + 5 for r in rows), 1960)
        self.assertEqual(1960 + 40 * 25 + 180 + 260, 3400)
        self.assertEqual(512 * 131073, 67109376)
        self.assertEqual(512 * 33, 16896)
        self.assertEqual(32 * 33, 1056)
        self.assertEqual(
            (1800 + 5) * 2 + (2400 + 5) * 2 + 600 + 120 + 60, 9200
        )
        self.assertEqual((3600 + 5) * 2 + 300 + 60 + 30, 7600)
        self.assertEqual((7200 + 5) * 2 + 300 + 60 + 30, 14800)
        self.assertEqual(3600 + 3400 + 9200 + 7600 + 14800, 38600)
        self.assertEqual(3124 + 3400 + 9200 + 7600 + 14800, 38124)
        for stage in (plan["calibration"], *plan["pilots"].values()):
            limits = stage["limits"]
            for field in ("wall", "cpu", "work"):
                self.assertEqual(
                    limits["total"][field],
                    sum(p[field] for p in limits["phases"].values()),
                )
        for stage in plan["pilots"].values():
            for digits in {r["digits"] for r in stage["calls"]}:
                pair = [r for r in stage["calls"] if r["digits"] == digits]
                self.assertEqual(len(pair), 2)
                self.assertEqual(pair[0]["work"], pair[1]["work"])
                self.assertEqual(
                    pair[0]["configuration"]["memory_bytes"],
                    pair[1]["configuration"]["memory_bytes"],
                )
        self.assertLess(95**59 * 20, 100**59)
        self.assertGreaterEqual(95**58 * 20, 100**58)

    def test_literal_challengers_and_nomination_menu_construct(self):
        plan = plan_v4.build_plan()
        rows = plan["calibration"]["calls"] + [
            row for stage in plan["pilots"].values() for row in stage["calls"]
        ]
        for row in rows:
            if row["arm"] == "control":
                continue
            values = copy.deepcopy(row["configuration"])
            values["collector"] = adapt_resource_configuration(
                GrowingDoubleLargeSieveConfig, **values["collector"]
            )
            self.assertEqual(
                asdict(
                    adapt_resource_configuration(CountedSIQSConfig, **values)
                ),
                row["configuration"],
            )
        for digits in (60, 62, 65, 70):
            for roots in ("standard", "cached", "batch"):
                for assignments in ("paired_flyer", "windowed_pair"):
                    for ranking in ("parity2_wide", "parity2_norm"):
                        cfg = configuration(
                            digits,
                            "coupled",
                            root_update_policy=roots,
                            assignment_policy=assignments,
                            multiplier_policy=ranking,
                        )
                        self.assertEqual(cfg.root_update_policy, roots)
        with self.assertRaises(ValueError):
            plan_v4.config_values(60, "misspelled")

    def test_drift_omission_and_extra_cells_are_rejected(self):
        plan = plan_v4.build_plan()
        for change in (
            lambda p: p["pilots"].pop("C"),
            lambda p: p["calibration"]["calls"].append(
                copy.deepcopy(p["calibration"]["calls"][0])
            ),
            lambda p: p["calibration"]["calls"][0].update(seconds=31),
            lambda p: p["calibration"]["limits"]["total"].update(wall=3500),
            lambda p: p.update(control_commit="changed"),
        ):
            bad = copy.deepcopy(plan)
            change(bad)
            with self.assertRaises(ValueError):
                plan_v4.verify_plan(bad)
        self.assertEqual(plan["execution_status"], "design_only")
        self.assertTrue(plan["execution_prerequisites"])

    def test_root_ablation_and_forced_form_deltas_are_literal(self):
        rows = plan_v4.build_plan()["calibration"]["calls"]
        for digits in (62, 70):
            roots = [
                r
                for r in rows
                if r["phase"] == "roots" and r["digits"] == digits
            ]
            configs = [dict(r["configuration"]) for r in roots]
            self.assertEqual(
                {c.pop("root_update_policy") for c in configs},
                {"standard", "cached", "batch"},
            )
            self.assertEqual(configs[0], configs[1])
            self.assertEqual(configs[1], configs[2])
        forms = [
            dict(r["configuration"])
            for r in rows
            if r["id"].startswith("forced-60-residue5")
        ]
        self.assertEqual(
            {c.pop("polynomial_policy") for c in forms},
            {"standard", "parity2"},
        )
        self.assertEqual(forms[0], forms[1])
        with patch.object(
            control_v4,
            "configuration",
            return_value=configuration(60, "coupled"),
        ):
            with self.assertRaises(ValueError):
                control_v4.configuration_values(60)


class OwnedControlTests(unittest.TestCase):
    def test_owned_types_literal_configs_and_small_proper_split(self):
        with tempfile.TemporaryDirectory() as folder:
            module, budgets, _ = control_v4.load_control(Path(folder))
            try:
                self.assertIsNot(module.SIQSConfig, SIQSConfig)
                for digits in (60, 62, 65, 70):
                    cfg = control_v4.control_configuration(digits, module)
                    self.assertEqual(
                        asdict(cfg), control_v4.configuration_values(digits)
                    )
                # Ordinary correctness QA, not a timed experiment or pilot.
                result = module.SIQSJob(
                    10403,
                    config=module.SIQSConfig(base_bound=200),
                    budget=recording_budget(
                        budgets.Budget,
                        work_limit=1000000,
                        seconds=None,
                        cpu_seconds=None,
                    ),
                ).run()
                self.assertEqual(result.reason, "factor_found")
                self.assertEqual(result.divisor * result.cofactor, 10403)
                self.assertTrue(1 < result.divisor < 10403)
            finally:
                for name in tuple(sys.modules):
                    if name == control_v4.PACKAGE or name.startswith(
                        control_v4.PACKAGE + "."
                    ):
                        del sys.modules[name]

    def test_corrupt_owned_source_is_rejected_before_loading(self):
        baseline = control_v4.checked_baseline()
        baseline["source"]["v2/__init__.py"] += "# changed\n"
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "control.json"
            path.write_text(json.dumps(baseline))
            with patch.object(control_v4, "BASELINE", path):
                with self.assertRaises(ValueError):
                    control_v4.checked_baseline()
