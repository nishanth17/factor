"""Wide-cache admission and independent retained v1 assignment vectors."""

import json
import unittest
from dataclasses import replace
from pathlib import Path

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.execution.ecm_stage_two import StageTwoSchedule
from v2.execution.ecm_workers import arithmetic_options
from v2.execution.portfolio_config import PortfolioConfig
from v2.execution.worker_contract import (
    Manifest,
    Search,
    WorkerSettings,
    digest,
)

VECTORS = (
    Path(__file__).resolve().parents[1]
    / "benchmarks/inputs/controls/b12_f4_v1_assignment_vectors.json"
)


class CacheCapacityTests(unittest.TestCase):
    def test_resource_caps_preserve_retained_mathematical_assignments(self):
        fixture = json.loads(VECTORS.read_text())
        self.assertEqual(fixture["search_version"], "b12-f4-v1")
        for vector in fixture["vectors"]:
            config = adapt_resource_configuration(
                PortfolioConfig, **vector["config"]
            )
            for mebibytes in (16, 64):
                settings = WorkerSettings(
                    decoded_cache_bytes=mebibytes * 2**20,
                    stage_two_cache_bytes=mebibytes * 2**20,
                )
                manifest = Manifest(
                    vector["n"],
                    vector["seed"],
                    vector["lineage"],
                    arithmetic_options(config, settings),
                )
                actual = [
                    [
                        position,
                        manifest.assignment(position)["id"],
                        manifest.assignment(position)["sigma"],
                    ]
                    for position, _, _ in vector["expected"]
                ]
                self.assertEqual(actual, vector["expected"])

    def test_wide_cache_still_requires_funded_simultaneous_ownership(self):
        settings = WorkerSettings(
            workers=4,
            decoded_cache_bytes=64 * 2**20,
            stage_two_cache_bytes=64 * 2**20,
        )
        self.assertGreaterEqual(settings.owned_bytes, 4 * (128 * 2**20))
        with self.assertRaises(TypeError):
            replace(settings, memory_bytes=512 * 2**20)

    def test_larger_capacity_is_finite_and_rejects_oversized_requests(self):
        StageTwoSchedule(1024, 64 * 2**20)
        for keyword in ("decoded_cache_bytes", "stage_two_cache_bytes"):
            with self.assertRaises(ValueError):
                adapt_resource_configuration(
                    WorkerSettings, **{keyword: 64 * 2**20 + 1}
                )
        with self.assertRaises(ValueError):
            StageTwoSchedule(1024, 64 * 2**20 + 1)

    def test_previous_resource_checkpoint_cannot_reset_charges(self):
        settings = WorkerSettings()
        config = PortfolioConfig(
            ecm_tiers=((50, 500, 4),), siqs=None, sss=None
        )
        manifest = Manifest(
            1000003, 7919, "root", arithmetic_options(config, settings)
        )
        checkpoint = Search(manifest, settings, fixed_work=True).snapshot()
        checkpoint["payload"]["version"] = "b12-f4-v6"
        checkpoint["sha256"] = digest(checkpoint["payload"])
        with self.assertRaisesRegex(ValueError, "unsupported"):
            Search.restore(checkpoint, manifest, settings, fixed_work=True)


if __name__ == "__main__":
    unittest.main()
