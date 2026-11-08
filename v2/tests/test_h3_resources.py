"""Deterministic resource ownership; no elapsed-time assertions."""

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.resources import (
    ProcessSample, RemovedResourceLimitError, atomic_save,
)


class Clock:
    def __init__(self):
        self.value = 0.0

    def __call__(self):
        return self.value


class Provider:
    def __init__(self):
        self.rows = [ProcessSample(os.getpid(), 1, 100, 0)]

    def physical_memory(self):
        return 1000

    def sample(self, roots):
        return list(self.rows)


class ResourceTests(unittest.TestCase):
    def setUp(self):
        self.wall, self.cpu, self.provider = Clock(), Clock(), Provider()

    def budget(self, **kwargs):
        return Budget(
            wall_clock=self.wall, cpu_clock=self.cpu,
            provider=self.provider, **kwargs,
        )

    def test_defaults_and_work_counter(self):
        budget = self.budget()
        self.assertIsNone(budget.seconds)
        self.assertIsNone(budget.cpu_seconds)
        self.assertEqual(budget.rss_bytes, 800)
        budget.consume(10**100)
        self.assertEqual(budget.used, 10**100)

    def test_removed_work_argument(self):
        with self.assertRaisesRegex(RemovedResourceLimitError, "work_limit"):
            self.budget(work_limit=10**20)

    def test_precedence_and_payload(self):
        budget = self.budget(seconds=1, cpu_seconds=1, rss_bytes=100)
        self.wall.value = self.cpu.value = 2
        budget.cancelled = lambda: True
        with self.assertRaisesRegex(BudgetExhaustedError, "cancelled"):
            budget.consume()
        budget.cancelled = None
        for name, limit in (("wall", 1), ("cpu", 1), ("rss", 100)):
            with self.assertRaisesRegex(BudgetExhaustedError, name + "_limit"):
                budget.publication()
            self.assertEqual(budget.stop["limit"], name)
            self.assertEqual(budget.stop["limit_value"], limit)
            self.assertEqual(budget.stop["scope"], "factorization")
            self.assertGreaterEqual(budget.stop["observed"], limit)
            setattr(budget, {"wall": "seconds", "cpu": "cpu_seconds",
                             "rss": "rss_bytes"}[name], None)

    def test_sampled_current_tree_rss(self):
        budget = self.budget(rss_bytes=500)
        self.provider.rows.append(ProcessSample(123, 2, 450, 1))
        with self.assertRaisesRegex(BudgetExhaustedError, "rss_limit"):
            budget.publication()
        self.assertEqual(budget.stop["observed"], 550)
        self.provider.rows.pop()
        budget.publication()
        self.assertEqual(budget.observed_rss, 100)

    def test_live_ack_exit_no_double_charge(self):
        budget = self.budget()
        self.provider.rows.append(ProcessSample(123, 2, 200, 3))
        budget.publication()
        self.assertEqual(budget.cpu_used, 3)
        budget.observe_child((123, 2), 3)
        budget.observe_child((123, 2), 5)
        budget.observe_child((123, 2), 5)
        self.assertEqual(budget.cpu_used, 5)
        self.provider.rows.pop()
        budget.publication()
        self.assertEqual(budget.cpu_used, 5)
        budget.observe_child((123, 3), 2)
        self.assertEqual(budget.cpu_used, 7)

    def test_close_start_excludes_idle(self):
        budget = self.budget()
        self.wall.value, self.cpu.value = 4, 2
        budget.close()
        self.wall.value, self.cpu.value = 104, 52
        self.assertEqual((budget.wall_used, budget.cpu_used), (4, 2))
        budget.start()
        self.wall.value, self.cpu.value = 107, 53
        budget.close().close()
        self.assertEqual((budget.wall_used, budget.cpu_used), (7, 3))

    def test_saved_history_and_retained_owner(self):
        budget = self.budget(cpu_seconds=60)
        budget.restore_charges(work=10, wall=20, cpu=40)
        self.cpu.value = 15
        self.assertEqual(budget.cpu_used, 55)
        budget.restore_charges(work=10, wall=20, cpu=40)
        self.assertEqual(budget.cpu_used, 55)
        fresh = self.budget(cpu_seconds=60)
        fresh.restore_charges(work=10, wall=20, cpu=40)
        self.assertEqual(fresh.cpu_used, 40)

    def test_atomic_optional_oversized_and_interrupted_save(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.json"
            old = {"published": "old"}
            atomic_save(path, old)
            large = {"published": "x" * (2**20 + 1)}
            self.assertIsNone(atomic_save(path, large))
            self.assertEqual(json.loads(path.read_text()), large)
            atomic_save(path, old)
            event = atomic_save(path, large, checkpoint_bytes=2**20)
            self.assertEqual(event["reason"], "checkpoint_size_skip")
            self.assertGreater(event["requested"], event["limit_value"])
            self.assertEqual(json.loads(path.read_text()), old)
            with patch("os.replace", side_effect=OSError("injected")):
                with self.assertRaisesRegex(OSError, "injected"):
                    atomic_save(path, large)
            self.assertEqual(json.loads(path.read_text()), old)
            self.assertEqual(list(Path(directory).iterdir()), [path])
