"""Independent queue events, refusal boundaries and derived-state resume."""

import json
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.execution.budget import BudgetExhaustedError
from v2.execution.work_budget import PollingBudget
from v2.qs import Polynomial, SieveCollector, SIQSJob, build_factor_base
from v2.qs.sieve_collector import SieveConfig
from v2.qs.sparse_queue import build
from v2.tests.test_qs import reference_positions, unlimited_budget
from v2.tests.test_qs_loop_reuse import enumerated_queue, queue_worker
from v2.tests.test_qs_sieve import signature
from v2.tests.test_siqs import allowance, configuration


class PackedQueueTests(unittest.TestCase):
    def test_explicit_large_memory_grant_retains_a_finite_ceiling(self):
        from v2.qs.current_parallel import QSWorkerConfig

        self.assertEqual(
            QSWorkerConfig(memory_bytes=22 * 2**30).memory_bytes,
            22 * 2**30,
        )
        with self.assertRaises(ValueError):
            QSWorkerConfig(memory_bytes=22 * 2**30 + 1)

    def test_packed_words_match_independent_events_and_paid_work(self):
        for width in (4, 7, 16):
            worker = queue_worker(width)
            worker.config = replace(worker.config, marking="packed_queue")
            lo, hi = -257, 258
            maximum = max(
                abs(worker.polynomial.value(x)) for x in range(lo, hi)
            )
            reference_budget = unlimited_budget()
            expected, reserved = enumerated_queue(
                worker, (lo, hi, maximum), reference_budget
            )
            ledger = unlimited_budget()
            worker.budget = PollingBudget(ledger, interval=7)
            split, queues, actual = build(worker, (lo, hi, maximum), 2**20)
            offset_bits = (width - 1).bit_length()
            index_bits = (len(worker._roots) - 1).bit_length()
            decoded = []
            for queue in queues:
                row = []
                for word in queue:
                    row.extend(
                        (
                            word & ((1 << offset_bits) - 1),
                            (word >> offset_bits) & ((1 << index_bits) - 1),
                            word >> (offset_bits + index_bits),
                        )
                    )
                decoded.append(row)
                self.assertEqual(queue.itemsize, 8)
            self.assertEqual(decoded, expected)
            self.assertEqual(split, worker._sparse_split)
            # Preserve the older conservative reservation and exact refusal
            # envelope while the physical queue drops from12 to8 bytes/hit.
            self.assertEqual(actual, reserved)
            self.assertEqual(ledger.used, reference_budget.used)

    def test_signed_singular_windows_and_refused_queue_keep_coverage(self):
        base = build_factor_base(10403, bound=100).factor_base
        for a, b in ((1, 102), (7, 1), (49, 8)):
            polynomial = Polynomial(10403, 1, a, b)
            expected = reference_positions(polynomial, base, -71, 80, 97)
            for cap in (1, 2**20):
                config = SieveConfig(
                    memory_bytes=64 * 2**20,
                    score_policy="powers",
                    division="bucket",
                    block_width=17,
                    residual_bound=97,
                    sparse_queue_bytes=cap,
                    marking="packed_queue",
                    max_atoms=4096,
                    max_relations=4096,
                    max_partials=4096,
                )
                collector = SieveCollector(
                    polynomial, base, config=config, budget=unlimited_budget()
                )
                result = collector.collect(-71, 80)
                self.assertEqual(signature(result), expected)
                self.assertLessEqual(collector._workspace, config.memory_limit)

    def test_interrupted_construction_releases_memory_before_replay(self):
        worker = queue_worker(7)
        worker.config = replace(worker.config, marking="packed_queue")
        before = worker._workspace
        with patch(
            "v2.qs.sparse_queue.build",
            side_effect=BudgetExhaustedError("work"),
        ):
            result = worker.collect(-63, 66)
        self.assertEqual(result.next_position, -63)
        self.assertEqual(worker._workspace, before)
        self.assertIsNone(worker._sparse_queue_key)
        replay = worker.collect(-63, 66)
        self.assertEqual(replay.reason, "complete")

    def test_checkpoint_rebuilds_derived_queue_and_pins_layout(self):
        config = configuration(
            sign_policy="swap_high",
            collector=SieveConfig(
                memory_bytes=64 * 2**20,
                score_policy="powers",
                division="bucket",
                sparse_queue_bytes=2**20,
                marking="packed_queue",
                residual_bound=500,
                max_atoms=4096,
                max_relations=4096,
                max_partials=4096,
            ),
        )
        job = SIQSJob(4001 * 5003, config=config, seed=73, budget=allowance())
        job.run(max_blocks=1)
        saved = json.loads(json.dumps(job.checkpoint()))
        restored = SIQSJob.from_checkpoint(saved, budget=allowance())
        self.assertEqual(restored.config, config)
        left, right = job.run(), restored.run()
        self.assertEqual(
            (left.divisor, left.cofactor), (right.divisor, right.cofactor)
        )
        self.assertEqual(left.divisor * left.cofactor, job.n)
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                saved,
                budget=allowance(),
                config=replace(
                    config,
                    collector=replace(config.collector, marking="sparse"),
                ),
            )


if __name__ == "__main__":
    unittest.main()
