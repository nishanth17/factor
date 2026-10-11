"""C3 allocation checkpoints retain C9 growing storage and finite ownership."""

import unittest
from dataclasses import replace

from v2.execution.allocation import ECMAllocation
from v2.portfolio import factorize_bounded
from v2.qs import GrowingDoubleLargeSieveConfig, SIQSConfig
from v2.tests.test_c3_allocation import configuration, ledger


class GrowingAllocationIntegrationTests(unittest.TestCase):
    def configuration(self):
        return configuration(
            siqs=SIQSConfig(
                base_bound=100,
                half_width=64,
                collector=GrowingDoubleLargeSieveConfig(
                    residual_bound=10000,
                    large_prime_bound=10000,
                    large_product_bound=1000000,
                ),
            ),
            allocation=ECMAllocation(
                "pretest", pretest_work=0, fallback_work=10000
            ),
        )

    def check_implicit_restore(self, work, expected_reason):
        config = self.configuration()
        number = 1009 * 1013
        first = factorize_bounded(number, config=config, budget=ledger(work))
        self.assertEqual(first.reason, expected_reason)
        self.assertEqual(first.result.reconstruct(), number)
        payload = first.checkpoint["payload"]
        self.assertEqual(payload["version"], 12)
        saved = payload["config"]["siqs"]
        self.assertEqual(saved["collector"]["storage_policy"], "grow")
        self.assertEqual(
            saved["memory_bytes"],
            config.memory_bytes - config.workspace_reserve - 8192,
        )
        self.assertGreater(saved["memory_bytes"], 0)
        for field in ("max_atoms", "max_partials", "max_relations"):
            self.assertIsNone(saved["collector"][field])

        for _ in range(2):
            restored = factorize_bounded(
                number, checkpoint=first.checkpoint, budget=ledger(10**13)
            )
            self.assertTrue(restored.result.complete)
            self.assertEqual(restored.result.reconstruct(), number)
            self.assertGreater(restored.work_used, first.work_used)
            self.assertGreaterEqual(restored.wall_seconds, first.wall_seconds)
            self.assertGreaterEqual(restored.cpu_seconds, first.cpu_seconds)
            self.assertEqual(
                restored.checkpoint["payload"]["config"]["siqs"], saved
            )
            self.assertEqual(
                sum(event["stage"] == "handoff" for event in restored.events),
                1,
            )
        return first, config

    def test_unadmitted_growing_fallback_restores_implicit_configuration(self):
        self.check_implicit_restore(1000, "insufficient_fallback_work")

    def test_live_growing_fallback_restores_and_rejects_changed_policy(self):
        first, config = self.check_implicit_restore(30000, "work_limit")
        current = first.checkpoint["payload"]["state"]["current"]
        self.assertTrue(current["fallback_admitted"])
        self.assertIn("siqs_checkpoint", current)

        changed = replace(
            config,
            siqs=replace(
                config.siqs,
                collector=replace(config.siqs.collector, max_partials=256),
            ),
        )
        with self.assertRaisesRegex(ValueError, "incompatible checkpoint"):
            factorize_bounded(
                1009 * 1013,
                config=changed,
                checkpoint=first.checkpoint,
                budget=ledger(),
            )
