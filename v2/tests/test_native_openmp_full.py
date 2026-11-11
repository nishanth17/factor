"""Actual whole-curve OpenMP admission; optional provided runtime required."""

import json
import os
import subprocess
import sys
import threading
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.ecm.native_curve_team import OpenMPCurveTeam
from v2.ecm.native_packet import SOURCE
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.ecm_native import NativeECMExecutor
from v2.execution.ecm_openmp import OpenMPExecutor
from v2.execution.portfolio_config import (
    _execution_context,
    close_execution_context,
)
from v2.execution.stage_jobs import new_job
from v2.execution.worker_contract import WorkerSettings
from v2.tests.test_h3_resources import Clock, Provider
from v2.tests.test_native_correctness import (
    FIXTURES,
    deterministic_budget,
    selected,
)
from v2.tests.test_native_full import drive, signature


@unittest.skipUnless(
    (SOURCE / "build/openmp/manifest.json").exists(),
    "provided whole-curve OpenMP build required",
)
class WholeOpenMPTests(unittest.TestCase):
    def environment(self, workers=4, schedule="dynamic", chunk=1, case=None):
        case = case or dict(n=2**199 - 189, seed=7, b1=31, b2=233)
        config = selected("c-gmp-full", case["b1"], case["b2"])
        budget = deterministic_budget()
        context, programs = _execution_context({}, config, budget)
        source = programs or context
        self.addCleanup(close_execution_context, source)
        team = OpenMPCurveTeam(
            workers=workers, schedule=schedule, chunk_size=chunk
        )
        self.addCleanup(team.close)
        team._grow(team.capacity, case["n"], budget)
        job = new_job("ecm", case["n"], case["seed"], case["b1"], case["b2"])
        return config, budget, source, team, job

    def finish(self, config, budget, source, team, job):
        for _ in range(10000):
            team.advance_jobs([(0, job)], budget, source, config)
            if job["done"]:
                return
        self.fail("finite OpenMP curve did not finish")

    def test_every_canonical_boundary_matches_independent_python(self):
        for case in json.loads(FIXTURES.read_text())["cases"]:
            expected, trace, _ = drive(case, "python")
            for workers in (1, 2, 4):
                for schedule in ("static", "dynamic"):
                    for chunk in (1, 3):
                        with self.subTest(
                            case=case["id"],
                            workers=workers,
                            schedule=schedule,
                            chunk=chunk,
                        ):
                            config, budget, source, team, job = (
                                self.environment(
                                    workers, schedule, chunk, case
                                )
                            )
                            team.atomic_actions = True
                            actual_trace = []
                            for _ in range(10000):
                                team.advance_jobs(
                                    [(0, job)], budget, source, config
                                )
                                info = team.bindings[0]._read()
                                if info.boundary or info.done:
                                    actual_trace.append(signature(job))
                                if job["done"]:
                                    break
                            self.assertEqual(job, expected)
                            self.assertEqual(actual_trace, trace)
                            self.assertEqual(
                                team.lib.b4_omp_participants(team.owner),
                                (1 << workers) - 1,
                            )
                            self.assertEqual(team.actual[0], workers)
                            team.close()
                            # libomp retains its own workers; these counters
                            # never misrepresent them as owner-joined threads.
                            self.assertEqual(team.lifecycle[:2], (0, 0))
                            self.assertGreater(team.lifecycle[2], 0)

    def test_context_alias_refusal_and_constrained_team_do_no_work(self):
        config, budget, source, team, job = self.environment()
        packet = team.packets[0]
        packet.bind_job(job, budget, source, config)
        team.contexts[0] = team.contexts[1] = packet.context
        before = packet._read().setup_ops
        self.assertEqual(
            team._native_step(
                team.owner,
                team.contexts,
                2,
                4,
                0,
                1,
                5,
                128,
                team.statuses,
                team.actual,
            ),
            -1,
        )
        self.assertEqual(packet._read().setup_ops, before)
        script = """
from v2.tests.test_native_openmp_full import WholeOpenMPTests
from v2.ecm.native_team import TeamRefusedError
case = WholeOpenMPTests()
try:
    config, budget, source, team, job = case.environment()
    try:
        team.advance_jobs([(0, job)], budget, source, config)
    except TeamRefusedError:
        assert team.actual[0] == 1
        assert team.bindings[0]._read().setup_ops == 0
        assert not job['done']
    else:
        raise AssertionError('constrained team silently accepted')
finally:
    case.doCleanups()
"""
        env = dict(os.environ, OMP_THREAD_LIMIT="1", OMP_DYNAMIC="FALSE")
        result = subprocess.run(
            [sys.executable, "-B", "-c", script],
            env=env,
            capture_output=True,
            text=True,
            timeout=20,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_cancel_failure_and_close_drain_before_reuse(self):
        config, budget, source, team, job = self.environment()
        team._native_inject(team.owner, -1, 0, 0)
        with self.assertRaises(RuntimeError):
            team.advance_jobs([(0, job)], budget, source, config)
        self.assertFalse(job["done"])
        # Clear the prior failed generation before the diagnostic waiter can
        # inspect its cancellation flag, independently of thread scheduling.
        self.assertEqual(team._native_reset(team.owner), 0)
        team._native_inject(team.owner, -1, -1, 1)
        errors = []

        def run():
            try:
                team.advance_jobs([(0, job)], budget, source, config)
            except BaseException as error:
                errors.append(error)

        caller = threading.Thread(target=run)
        caller.start()
        try:
            team._native_wait(team.owner, 4)
            self.assertEqual(team._native_shutdown(team.owner), -1)
            with self.assertRaises(RuntimeError):
                team.close()
            with self.assertRaises(RuntimeError):
                team.bindings[0].close()
            other = OpenMPCurveTeam(workers=1)
            try:
                with self.assertRaisesRegex(RuntimeError, "another native"):
                    other.advance_jobs([(0, job)], budget, source, config)
            finally:
                other.close()
            team.cancel()
        finally:
            team._native_release(team.owner)
            caller.join(10)
        self.assertFalse(caller.is_alive())
        self.assertEqual(len(errors), 1)
        self.assertIsInstance(errors[0], BudgetExhaustedError)
        self.assertFalse(job["done"])
        self.assertTrue(all(not p.call_lock.locked() for p in team.packets))
        team._native_inject(team.owner, -1, -1, 0)
        self.finish(config, budget, source, team, job)
        team.bindings[0].sync(allow_result=True)
        self.assertEqual(
            job, drive(dict(n=job["n"], seed=7, b1=31, b2=233), "python")[0]
        )
        team.close()
        team.close()
        self.assertEqual(team.packets, [])

    def test_root_limits_after_native_work_precede_publication(self):
        for cause in ("wall", "cpu", "rss"):
            config, _, source, team, job = self.environment()
            wall, cpu, provider = Clock(), Clock(), Provider()
            budget = Budget(
                seconds=10,
                cpu_seconds=10,
                rss_bytes=500,
                wall_clock=wall,
                cpu_clock=cpu,
                provider=provider,
            )
            native_step = team._native_step

            def expire(*args):
                result = native_step(*args)
                if cause == "wall":
                    wall.value = 10
                elif cause == "cpu":
                    cpu.value = 10
                else:
                    provider.rows = [replace(provider.rows[0], rss=500)]
                return result

            with patch.object(team, "_native_step", expire):
                before = budget.used
                with self.assertRaisesRegex(
                    BudgetExhaustedError, cause + "_limit"
                ):
                    team.advance_jobs([(0, job)], budget, source, config)
                self.assertGreater(budget.used, before)
                self.assertGreater(team.bindings[0]._read().setup_ops, 0)
                self.assertFalse(job["done"])

    def test_tails_assignment_history_replay_and_backend_identity(self):
        config = replace(
            selected("c-gmp-full", 31, 233), ecm_tiers=((31, 233, 37),)
        )
        expected = None
        for workers in (1, 2, 4):
            for schedule in ("static", "dynamic"):
                owner = OpenMPExecutor(
                    seed=73,
                    settings=WorkerSettings(workers=workers),
                    schedule=schedule,
                    capacity=32,
                )
                try:
                    result = owner.run(
                        2**199 - 189,
                        config=config,
                        budget=deterministic_budget(),
                        fixed_work=True,
                    )
                    completed = result.checkpoint["payload"]["completed"]
                    self.assertEqual(len(completed), 37)
                    self.assertEqual(owner.team.actual[0], workers)
                    if expected is None:
                        expected = completed
                    self.assertEqual(completed, expected)
                finally:
                    owner.close()
        owner = OpenMPExecutor(seed=73, settings=WorkerSettings(workers=2))
        self.addCleanup(owner.close)
        before = [0]

        def cancel():
            return owner.team is not None and owner.team.slices > before[0]

        paused = owner.run(
            2**199 - 189,
            config=config,
            budget=Budget(provider=Provider(), cancelled=cancel),
            fixed_work=True,
        )
        self.assertEqual(paused.reason, "cancelled")
        checkpoint = json.loads(json.dumps(paused.checkpoint))
        resumed = owner.run(
            2**199 - 189,
            config=config,
            budget=deterministic_budget(),
            fixed_work=True,
            checkpoint=checkpoint,
        )
        self.assertEqual(resumed.checkpoint["payload"]["completed"], expected)
        control = NativeECMExecutor(
            seed=73, settings=WorkerSettings(workers=2)
        )
        try:
            with self.assertRaises(ValueError):
                control.run(
                    2**199 - 189,
                    config=config,
                    budget=deterministic_budget(),
                    checkpoint=checkpoint,
                    fixed_work=True,
                )
        finally:
            control.close()
