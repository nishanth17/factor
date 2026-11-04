"""Explicit oversubscription QA; no automatic worker policy change."""

import json
import multiprocessing
import unittest
from unittest.mock import patch

from v2.benchmarks.ecm.b12_f4.profile_runner import signature
from v2.execution.ecm_workers import ECMExecutor
from v2.execution.worker_contract import WorkerSettings
from v2.tests.test_ecm_workers import allowance, configuration


class SixteenWorkerTests(unittest.TestCase):
    def test_spawn_durable_restart_and_exact_serial_manifest(self):
        before = {child.pid for child in multiprocessing.active_children()}
        n = 1009 * 1013
        config = configuration(ecm_tiers=((20, 200, 32),))
        settings = WorkerSettings(
            workers=16,
            allow_wide_workers=True,
            startup_seconds=5,
            assignment_cpu=2,
            memory_bytes=4 * 2**30,
            checkpoint_bytes=4 * 2**20,
        )
        with ECMExecutor(seed=7919, force_workers=False) as serial:
            reference = serial.run(
                n, config=config, budget=allowance(), fixed_work=True
            )
        with ECMExecutor(
            seed=7919, settings=settings, force_workers=True
        ) as parallel:
            run = parallel.run(
                n, config=config, budget=allowance(), fixed_work=True
            )
            self.assertEqual(
                run.stats["termination_reason"], "curves_exhausted"
            )
            self.assertEqual(len(parallel.slots), 16)
            self.assertEqual(
                signature(run.checkpoint, n, 32),
                signature(reference.checkpoint, n, 32),
            )
        checkpoint = json.loads(json.dumps(run.checkpoint))
        with ECMExecutor(
            seed=7919, settings=settings, force_workers=True
        ) as resumed:
            replay = resumed.run(
                n,
                config=config,
                budget=allowance(),
                checkpoint=checkpoint,
                fixed_work=True,
            )
            self.assertEqual(
                signature(replay.checkpoint, n, 32),
                signature(reference.checkpoint, n, 32),
            )
            self.assertEqual(replay.stats["worker_dispatches"], 0)
            self.assertGreaterEqual(
                replay.checkpoint["budget"]["cpu"],
                checkpoint["budget"]["cpu"],
            )
        self.assertEqual(
            {child.pid for child in multiprocessing.active_children()}, before
        )


class FixedWorkTerminationTests(unittest.TestCase):
    def test_a_factor_does_not_hide_an_incomplete_fixed_work_search(self):
        calls = []

        def report(engine, task, **kwargs):
            calls.append(task["assignment"]["id"])
            first = len(calls) == 1
            return dict(
                id=task["assignment"]["id"],
                token=task["token"],
                attempt=task["attempt"],
                stage=task["stage"],
                factor=5 if first else None,
                reason="complete" if first else "cpu_limit",
                work=1,
            )

        with patch("v2.execution.ecm_worker.CurveEngine.execute", report):
            with ECMExecutor(force_workers=False) as owner:
                run = owner.run(
                    35,
                    config=configuration(),
                    budget=allowance(),
                    fixed_work=True,
                )
        self.assertEqual(run.divisor, 5)
        self.assertEqual(run.reason, "factor_found")
        self.assertEqual(run.stats["termination_reason"], "cpu_limit")
        self.assertEqual(run.stats["completed_after"], 1)
        self.assertEqual(len(calls), 2)
