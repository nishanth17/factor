"""Production A7a routing, literal Y0 recipes and saved execution identity."""

import json
import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.qs.a7a.cell import configuration
from v2.execution.budget import Budget
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs.defaults import default_siqs_config
from v2.qs.production_policy import default_qs_workers, qs_recipe
from v2.qs.siqs import SIQSConfig
from v2.tests.test_current_qs_workers import configuration as small_config
from v2.tests.test_current_qs_workers import settings as small_workers
from v2.tests.test_siqs_cli import cli


class ProductionQSTests(unittest.TestCase):
    def test_automatic_pool_executes_and_drains_on_a_real_small_factorization(
        self,
    ):
        import multiprocessing

        n = 97 * 4001 * 5003
        config = PortfolioConfig(
            trial_bound=100,
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=(),
            siqs=small_config(),
        )
        before = {p.pid for p in multiprocessing.active_children()}
        with patch(
            "v2.qs.production_policy.default_qs_workers",
            return_value=small_workers("process", 2),
        ):
            run = factorize_bounded(
                n,
                seed=29,
                config=config,
                ecm_executor=None,
                budget=Budget(seconds=15, cpu_seconds=60),
            )
        self.assertEqual(run.result.reconstruct(), n)
        self.assertTrue(run.result.complete)
        self.assertTrue(any(e["stage"] == "siqs" for e in run.events))
        self.assertEqual(
            {p.pid for p in multiprocessing.active_children()}, before
        )

    def test_y0_anchors_match_frozen_experiment_policy(self):
        for digits in (40, 50, 55, 60, 65, 70):
            n = 10 ** (digits - 1) + 123
            self.assertEqual(
                asdict(default_siqs_config(n)), asdict(configuration(n, "Y0"))
            )

    def test_nearest_worker_rows_signs_ties_and_extrapolation(self):
        bands = (
            (1, 39, 0, 0),
            (40, 45, 4, 4),
            (46, 52, 4, 8),
            (53, 57, 4, 16),
            (58, 62, 8, 16),
            (63, 67, 14, 64),
            (68, 100, 8, 64),
        )
        for low, high, count, packet in bands:
            for digits in range(low, high + 1):
                for sign in (-1, 1):
                    recipe = qs_recipe(sign * 10 ** (digits - 1))
                    self.assertEqual(
                        (recipe.workers, recipe.packet), (count, packet)
                    )
                    self.assertEqual(recipe.extrapolated, digits > 70)
        self.assertEqual(qs_recipe(10**44).anchor, 40)

    def test_worker_bundle_preserves_selected_features_and_capacities(self):
        for digits, count, packet in (
            (40, 4, 4),
            (50, 4, 8),
            (55, 4, 16),
            (60, 8, 16),
            (65, 14, 64),
            (70, 8, 64),
        ):
            n = 10 ** (digits - 1) + 123
            selected = default_qs_workers(n, default_siqs_config(n))
            self.assertEqual(
                (selected.workers, selected.polynomials_per_task),
                (count, packet),
            )
            self.assertEqual(selected.max_batch_atoms, 2048)
            self.assertEqual(selected.state_cache_bytes, 64 * 2**20)
            self.assertEqual(selected.packet_bytes, 4 * 2**20)
            self.assertEqual(selected.lookahead, 4)
            self.assertEqual(selected.verification, "worker")
            self.assertEqual(selected.stop_backend, "semaphore")
            for flag in (
                "lease_budget",
                "reuse_state",
                "packed_transport",
                "private_stop",
                "coefficient_cache",
            ):
                self.assertTrue(getattr(selected, flag))
            for flag in (
                "observe",
                "whole_family",
                "arrival_order",
                "cpu_snapshot",
                "family_affinity",
                "grouped_lease_charging",
            ):
                self.assertFalse(getattr(selected, flag))

    def test_unsupported_explicit_engines_remain_serial(self):
        n = 10**49 + 123
        base = SIQSConfig()
        for config in (
            None,
            replace(base, mode="mpqs"),
            replace(base, growth_steps=1),
            replace(base, shared_relations=False),
            replace(base, diverse=False),
        ):
            self.assertIsNone(default_qs_workers(n, config))
        self.assertIsNone(default_qs_workers(10403, base))

    def test_zero_budget_serial_override_and_saved_recipe(self):
        n = 10**64 + 123
        config = PortfolioConfig(siqs=default_siqs_config(n))
        # No child pool may be prepared before the relation stage is reached.
        with patch(
            "v2.qs.parallel.CollectionPool", side_effect=AssertionError
        ):
            automatic = factorize_bounded(
                n, config=config, budget=Budget(seconds=0), ecm_executor=None
            )
            serial = factorize_bounded(
                n,
                config=config,
                budget=Budget(seconds=0),
                ecm_executor=None,
                qs_workers=None,
            )
        state = automatic.checkpoint["payload"]["state"]
        self.assertEqual(state["qs_workers"]["workers"], 14)
        self.assertEqual(
            state["qs_routing"],
            dict(policy="a7a-Y0", anchor=65, extrapolated=False),
        )
        self.assertNotIn("qs_workers", serial.checkpoint["payload"]["state"])
        with patch(
            "v2.qs.production_policy.default_qs_workers",
            side_effect=AssertionError("resume reselected workers"),
        ):
            resumed = factorize_bounded(
                n,
                checkpoint=automatic.checkpoint,
                budget=Budget(seconds=0),
                ecm_executor=None,
            )
            legacy = factorize_bounded(
                n,
                checkpoint=serial.checkpoint,
                budget=Budget(seconds=0),
                ecm_executor=None,
            )
        self.assertEqual(
            resumed.checkpoint["payload"]["state"]["qs_workers"],
            state["qs_workers"],
        )
        self.assertNotIn("qs_workers", legacy.checkpoint["payload"]["state"])
        with self.assertRaisesRegex(ValueError, "worker settings changed"):
            factorize_bounded(
                n,
                checkpoint=automatic.checkpoint,
                budget=Budget(seconds=0),
                ecm_executor=None,
                qs_workers=None,
            )
        large = 10**74 + 123
        extrapolated = factorize_bounded(
            large,
            config=PortfolioConfig(siqs=default_siqs_config(large)),
            budget=Budget(seconds=0),
            ecm_executor=None,
        )
        self.assertTrue(
            extrapolated.checkpoint["payload"]["state"]["qs_routing"][
                "extrapolated"
            ]
        )

    def test_cli_auto_and_serial_checkpoint_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            for flag in ("auto", "serial"):
                path = Path(directory) / (flag + ".json")
                result = cli(
                    10**69 + 123,
                    "--method",
                    "siqs",
                    "--seconds",
                    "0",
                    "--qs-workers",
                    flag,
                    "--checkpoint",
                    path,
                )
                self.assertEqual(result.returncode, 1, result.stderr)
                checkpoint = json.loads(path.read_text())
                worker = checkpoint["payload"]["state"].get("qs_workers")
                self.assertEqual(
                    None if worker is None else worker["workers"],
                    8 if flag == "auto" else None,
                )
                resumed = cli(
                    "--resume", path, "--seconds", "0", "--checkpoint", path
                )
                self.assertEqual(resumed.returncode, 1, resumed.stderr)
                self.assertEqual(
                    json.loads(path.read_text())["payload"]["state"].get(
                        "qs_workers"
                    ),
                    worker,
                )


if __name__ == "__main__":
    unittest.main()
