"""Deterministic resource ownership; no elapsed-time assertions."""

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.resources import (
    ProcessSample,
    RemovedResourceLimitError,
    atomic_save,
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
            wall_clock=self.wall,
            cpu_clock=self.cpu,
            provider=self.provider,
            **kwargs,
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
            setattr(
                budget,
                {"wall": "seconds", "cpu": "cpu_seconds", "rss": "rss_bytes"}[
                    name
                ],
                None,
            )

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


class StorageMigrationTests(unittest.TestCase):
    def test_internal_limit_arguments_are_named_errors(self):
        from v2.portfolio import PortfolioConfig
        from v2.qs import SieveConfig, SIQSConfig
        from v2.qs.extended_policy import BoundedYieldSieveConfig

        for factory, argument in (
            (SieveConfig, "max_atoms"),
            (SieveConfig, "max_relations"),
            (SieveConfig, "memory_bytes"),
            (SIQSConfig, "checkpoint_bytes"),
            (PortfolioConfig, "memory_bytes"),
            (BoundedYieldSieveConfig, "spill_bytes"),
        ):
            with self.subTest(factory=factory, argument=argument):
                with self.assertRaisesRegex(
                    RemovedResourceLimitError, argument
                ):
                    factory(**{argument: None})

    def test_graph_65536_is_not_an_admission_ceiling(self):
        from v2.qs.large_primes import LargePrimeForest, graph_reserve

        class BoundaryDict(dict):
            def __len__(self):
                return 65536

        graph = LargePrimeForest()
        graph.edges = BoundaryDict()
        plan = graph.plan_link(101, 103, "boundary", Budget(rss_bytes=None))
        graph.commit_link(plan)
        self.assertIn("boundary", graph.edges)
        self.assertGreater(graph_reserve(65537), graph_reserve(65536))

    def test_matrix_estimate_is_diagnostic(self):
        from v2.qs.linear_algebra import filter_matrix

        with patch(
            "v2.qs.linear_algebra.matrix_workspace", return_value=10**30
        ):
            matrix = filter_matrix((3, 3), budget=Budget(rss_bytes=None))
        self.assertEqual(matrix.original_rows, (3, 3))
        self.assertEqual(matrix.workspace_bytes, 10**30)

    def test_old_checkpoint_refused_and_saved_resume(self):
        from v2.execution.resources import UnsupportedResourceVersionError
        from v2.qs.siqs import SIQSConfig, SIQSJob

        job = SIQSJob(
            1009 * 1013, config=SIQSConfig(), budget=Budget(rss_bytes=None)
        )
        job.run(max_blocks=1)
        snapshot = job.checkpoint()
        old = dict(snapshot)
        del old["resource_version"]
        with self.assertRaises(UnsupportedResourceVersionError):
            SIQSJob.from_checkpoint(old, budget=Budget(rss_bytes=None))
        restored = SIQSJob.from_checkpoint(
            snapshot, budget=Budget(rss_bytes=None)
        )
        self.assertEqual(
            restored.engine.next_position, job.engine.next_position
        )
        self.assertGreaterEqual(restored.budget.used, job.budget.used)

    def test_stage_share_hands_off_with_unfunded_reserve(self):
        import contextlib
        import io

        from v2.execution.allocation import ECMAllocation
        from v2.portfolio import PortfolioConfig, factorize_bounded
        from v2.qs.siqs import SIQSConfig

        output = io.StringIO()
        config = PortfolioConfig(
            allocation=ECMAllocation(
                "pretest", pretest_work=0, fallback_work=10**100
            ),
            siqs=SIQSConfig(mode="qs", base_bound=100),
            ecm_chain_mode="off",
        )
        with (
            contextlib.redirect_stdout(output),
            contextlib.redirect_stderr(output),
        ):
            run = factorize_bounded(
                1009 * 1013,
                config=config,
                budget=Budget(rss_bytes=None),
                ecm_executor=None,
            )
        self.assertEqual(output.getvalue(), "")
        self.assertEqual(run.result.reconstruct(), 1009 * 1013)
        self.assertNotIn(
            run.reason, ("pretest_exhausted", "insufficient_fallback_work")
        )


class ChildReaperTests(unittest.TestCase):
    def test_live_ack_exit_and_crash_use_one_total(self):
        from types import SimpleNamespace

        provider = Provider()
        provider.rows.append(ProcessSample(123, 10, 50, 2))
        budget = Budget(
            provider=provider, wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        process = SimpleNamespace(
            pid=123, _popen=SimpleNamespace(pid=123, returncode=None)
        )
        budget.track_process(process)
        budget.acknowledge(123, 3)
        usage = SimpleNamespace(ru_utime=4, ru_stime=1)
        with patch("os.wait4", return_value=(123, 9, usage)) as wait:
            self.assertEqual(process._popen.poll(), -9)
            self.assertEqual(process._popen.poll(), -9)
            wait.assert_called_once_with(123, os.WNOHANG)
        # The initial retained CPU was a baseline, so only 5 - 2 is new.
        self.assertEqual(budget.cpu_used, 3)
        budget.acknowledge(123, 5)
        self.assertEqual(budget.cpu_used, 3)

    def test_retained_process_rebinds_and_excludes_idle_cpu(self):
        from types import SimpleNamespace

        provider = Provider()
        provider.rows.append(ProcessSample(123, 10, 50, 2))
        first = Budget(
            provider=provider, wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        process = SimpleNamespace(
            pid=123, _popen=SimpleNamespace(pid=123, returncode=None)
        )
        first.track_process(process)
        first.close()
        provider.rows[-1] = ProcessSample(123, 10, 50, 7)
        second = Budget(
            provider=provider, wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        second.track_process(process)
        usage = SimpleNamespace(ru_utime=8, ru_stime=1)
        with patch("os.wait4", return_value=(123, 0, usage)):
            self.assertEqual(process._popen.poll(0), 0)
        self.assertEqual(first.cpu_used, 0)
        self.assertEqual(second.cpu_used, 2)
