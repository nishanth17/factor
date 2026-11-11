"""Whole-ECM C execution, canonical oracles and interruption recovery."""

import copy
import hashlib
import json
import signal
import threading
import unittest
from contextlib import ExitStack
from dataclasses import replace
from unittest.mock import patch

from v2.ecm import core
from v2.ecm.chains import ChainPlan
from v2.ecm.native_curve import NativeCurvePacket
from v2.ecm.native_curve_team import NativeCurveTeam, OpenMPCurveTeam
from v2.ecm.native_packet import SOURCE
from v2.ecm.native_team import TeamRefusedError
from v2.execution import stage_jobs
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.ecm_native import NativeECMExecutor
from v2.execution.portfolio_config import (
    _execution_context,
    close_execution_context,
)
from v2.execution.worker_contract import WorkerSettings
from v2.tests.test_h3_resources import Clock, Provider
from v2.tests.test_native_correctness import (
    FIXTURES,
    deterministic_budget,
    selected,
)


def signature(job):
    return hashlib.sha256(json.dumps(job, sort_keys=True).encode()).hexdigest()


def drive(case, kernel, *, job=None, atomic=False, stop=None, config=None):
    config = config or selected(kernel, case["b1"], case["b2"])
    budget = deterministic_budget()
    context, programs = _execution_context({}, config, budget)
    source = programs or context
    packet = getattr(source, "native_packet", None)
    if isinstance(packet, NativeCurvePacket):
        packet.atomic_actions = atomic
    job = (
        copy.deepcopy(job)
        if job is not None
        else stage_jobs.new_job(
            "ecm", case["n"], case["seed"], case["b1"], case["b2"]
        )
    )
    trace = []
    try:
        for _ in range(100000):
            stage_jobs.advance_job(job, budget, source, config)
            if isinstance(packet, NativeCurvePacket):
                if job["done"] and not atomic:
                    packet.sync(allow_result=True)
                info = packet._read()
                if atomic and (info.boundary or info.done):
                    trace.append(signature(job))
            else:
                trace.append(signature(job))
            if job["done"] or stop is not None and stop(job, packet):
                break
        else:
            raise AssertionError("finite native curve did not finish")
    finally:
        close_execution_context(source)
    stats = None if packet is None else dict(getattr(packet, "statistics", {}))
    return job, trace, stats


@unittest.skipUnless(
    (SOURCE / "build/manifest.json").exists(), "scalar native build required"
)
class WholeScalarTests(unittest.TestCase):
    def test_auto_chain_refill_checkpoint_resumes_in_fresh_owner(self):
        corpus = FIXTURES.parents[1] / "corpora/b4a_production_3p_v1.json"
        cases = json.loads(corpus.read_text())["fixtures"]
        for case in cases:
            if case["population"] != "failed" or case["cluster"] != 0:
                continue
            with self.subTest(case=case["id"]):
                case = dict(case, seed=73, b1=127, b2=4000)
                config = replace(
                    selected("c-gmp-full", 127, 4000), ecm_chain_mode="auto"
                )
                expected, trace, _ = drive(
                    case, "python", config=replace(config, ecm_kernel="python")
                )
                actual, checkpoints, _ = drive(
                    case, "c-gmp-full", atomic=True, config=config
                )
                self.assertEqual(actual, expected)
                self.assertEqual(checkpoints, trace)

                def after_refills(job, packet):
                    return (
                        job["phase"] == "stage_two"
                        and job["cursor"]["left"] > 255
                    )

                paused, _, _ = drive(
                    case,
                    "c-gmp-full",
                    atomic=True,
                    config=config,
                    stop=after_refills,
                )
                self.assertFalse(paused["done"])
                self.assertGreater(paused["cursor"]["left"], 255)
                resumed, _, _ = drive(
                    case, "c-gmp-full", job=paused, atomic=True, config=config
                )
                self.assertEqual(resumed, expected)

    def test_portfolio_checkpoint_resume_after_actual_full_native_deadline(
        self,
    ):
        from v2.portfolio import factorize_bounded
        from v2.tests.test_native_correctness import after_native_slice

        case = json.loads(FIXTURES.read_text())["cases"][0]
        config = selected("c-gmp-full")
        wall, cpu = Clock(), Clock()
        budget = Budget(
            seconds=10, wall_clock=wall, cpu_clock=cpu, provider=Provider()
        )
        advance = NativeCurvePacket.advance_job

        def expire(packet, job, allowance, context, selected_config):
            packet.bind_job(job, allowance, context, selected_config)
            with after_native_slice(
                packet, lambda: setattr(wall, "value", 10)
            ):
                return advance(
                    packet, job, allowance, context, selected_config
                )

        with patch.object(NativeCurvePacket, "advance_job", expire):
            paused = factorize_bounded(
                case["n"],
                seed=1,
                config=config,
                budget=budget,
                ecm_executor=None,
            )
        self.assertEqual(paused.reason, "wall_limit")
        saved = json.loads(json.dumps(paused.checkpoint))
        resumed = factorize_bounded(
            case["n"],
            seed=1,
            config=config,
            budget=deterministic_budget(),
            checkpoint=saved,
            ecm_executor=None,
        )
        control = factorize_bounded(
            case["n"],
            seed=1,
            config=selected("python"),
            budget=deterministic_budget(),
            ecm_executor=None,
        )
        self.assertEqual(resumed.result, control.result)
        self.assertEqual(resumed.result.reconstruct(), case["n"])
        self.assertGreaterEqual(resumed.work_used, paused.work_used)

    def test_restore_preparation_failure_preserves_original_canonical_job(
        self,
    ):
        case = dict(n=2**199 - 189, seed=7, b1=31, b2=233)
        saved, _, _ = drive(
            case, "python", stop=lambda job, packet: bool(job["powers"])
        )
        for error in (
            MemoryError("allocation"),
            BudgetExhaustedError("cancelled"),
        ):
            config = selected("c-gmp-full", 31, 233)
            budget = deterministic_budget()
            context, programs = _execution_context({}, config, budget)
            source = programs or context
            packet = source.native_packet
            job = copy.deepcopy(saved)
            try:
                with patch.object(
                    packet.blocks, "register", side_effect=error
                ):
                    with self.assertRaises(type(error)):
                        packet.bind_job(job, budget, source, config)
            finally:
                close_execution_context(source)
            self.assertEqual(job, saved)

    def test_strict_coroutine_interruption_exports_saved_chunk_input(self):
        import random

        case = dict(n=1009 * 1013 * (10**40 + 121), seed=17, b1=2000, b2=2500)
        config = replace(
            selected("c-gmp-full", 2000, 2500),
            ecm_tiers=((2000, 2500, 8),),
            ecm_chain_mode="auto",
        )
        budget = deterministic_budget()
        context, programs = _execution_context({}, config, budget)
        source = programs or context
        packet = source.native_packet
        packet.actions = 1
        job = stage_jobs.new_job("ecm", case["n"], 17, 2000, 2500)
        job.update(
            phase="stage_one",
            value=[1009, 1],
            a24=2,
            sigma=random.Random(17).randrange(6, 2**63),
        )
        try:
            for _ in range(10000):
                stage_jobs.advance_job(job, budget, source, config)
                if packet._read().strict_state:
                    break
            else:
                self.fail("native strict coroutine was not reached")
        finally:
            close_execution_context(source)
        self.assertFalse(job["done"])
        self.assertEqual(job["value"], [1009, 1])
        self.assertEqual(len(job["powers"]), 16)
        self.assertNotIn("chain_chunks", job)
        resumed, _, _ = drive(case, "c-gmp-full", job=job, config=config)
        # Compare the supplied exceptional point with its independent oracle.
        control_job = stage_jobs.new_job("ecm", case["n"], 17, 2000, 2500)
        control_job.update(
            phase="stage_one", value=[1009, 1], a24=2, sigma=job["sigma"]
        )
        expected, _, _ = drive(
            case,
            "python",
            job=control_job,
            config=replace(config, ecm_kernel="python"),
        )
        self.assertEqual(resumed, expected)

    def environment(self, case):
        config = selected("c-gmp-full", case["b1"], case["b2"])
        budget = deterministic_budget()
        context, programs = _execution_context({}, config, budget)
        source = programs or context
        job = stage_jobs.new_job(
            "ecm", case["n"], case["seed"], case["b1"], case["b2"]
        )
        packet = source.native_packet
        self.addCleanup(close_execution_context, source)
        packet.bind_job(job, budget, source, config)
        return config, budget, source, job, packet

    def test_mid_native_cancellation_and_preparation_failure_resume(self):
        case = dict(n=2**199 - 189, seed=7, b1=31, b2=233)
        expected, _, _ = drive(case, "python")
        for preparation, cancel, error in (
            (0, 5, BudgetExhaustedError),
            (1, 0, MemoryError),
        ):
            config, budget, source, job, packet = self.environment(case)
            packet.lib.b4_curve_faults(packet.context, preparation, cancel)
            with self.assertRaises(error):
                for _ in range(10000):
                    stage_jobs.advance_job(job, budget, source, config)
                    self.assertFalse(job["done"])
            self.assertGreater(packet._read().stage_one_ops, 0)
            close_execution_context(source)
            self.assertFalse(job["done"])
            actual, _, _ = drive(case, "c-gmp-full", job=job)
            self.assertEqual(actual, expected)

    def test_root_stops_after_native_work_before_result_publication(self):
        for cause in ("wall", "cpu", "rss", "cancelled"):
            config = selected("c-gmp-full", 11, 47)
            wall, cpu, provider = Clock(), Clock(), Provider()
            cancelled = [False]
            budget = Budget(
                seconds=10,
                cpu_seconds=10,
                rss_bytes=500,
                wall_clock=wall,
                cpu_clock=cpu,
                provider=provider,
                cancelled=lambda: cancelled[0],
            )
            context, programs = _execution_context({}, config, budget)
            source = programs or context
            packet = source.native_packet
            job = stage_jobs.new_job("ecm", 15, 7, 11, 47)
            packet.bind_job(job, budget, source, config)
            original = packet.lib

            class Library:
                def __getattr__(self, name):
                    return getattr(original, name)

                def b4_step(self, *args):
                    result = original.b4_step(*args)
                    if cause == "wall":
                        wall.value = 10
                    elif cause == "cpu":
                        cpu.value = 10
                    elif cause == "rss":
                        provider.rows = [replace(provider.rows[0], rss=500)]
                    else:
                        cancelled[0] = True
                    return result

            packet.lib = Library()
            try:
                with self.assertRaisesRegex(
                    BudgetExhaustedError,
                    "cancelled" if cause == "cancelled" else cause + "_limit",
                ):
                    stage_jobs.advance_job(job, budget, source, config)
                self.assertFalse(job["done"])
                self.assertGreater(packet._read().setup_ops, 0)
            finally:
                packet.lib = original
                close_execution_context(source)
            self.assertFalse(job["done"])

    def test_native_table_geometry_boundary_rejects_without_fallback(self):
        packet = NativeCurvePacket()
        self.addCleanup(packet.close)
        packet.available(1013, deterministic_budget())
        packet._buffers()
        owners, values = packet._values([1013, 7])
        try:
            for distance, status in ((65536, 0), (65537, -2)):
                self.assertEqual(
                    packet.lib.b4_curve_begin(
                        packet.context,
                        values[0][0],
                        values[1][0],
                        200000,
                        distance**2,
                        16,
                        64,
                        1024,
                    ),
                    status,
                )
            self.assertNotEqual(
                packet.lib.b4_curve_info_get(packet.context, packet.info), 0
            )
        finally:
            for item in reversed(values + owners):
                packet.ffi.release(item)

    def test_all_fixed_fixtures_match_every_canonical_action(self):
        for case in json.loads(FIXTURES.read_text())["cases"]:
            with self.subTest(case=case["id"]):
                expected, trace, _ = drive(case, "python")
                actual, native_trace, stats = drive(
                    case, "c-gmp-full", atomic=True
                )
                self.assertEqual(actual, expected)
                self.assertEqual(native_trace, trace)
                self.assertEqual(actual["factor"], case["factor"])
                self.assertGreater(stats["setup_ops"], 0)

    def test_no_python_ecm_math_or_recovery_can_execute(self):
        case = dict(n=2**199 - 189, seed=7, b1=31, b2=233)
        expected, _, _ = drive(case, "python")
        with ExitStack() as stack:
            for name in (
                "setup_curve",
                "point_double",
                "point_add",
                "scalar_multiply",
            ):
                stack.enter_context(
                    patch.object(core, name, side_effect=AssertionError(name))
                )
            stack.enter_context(
                patch.object(
                    stage_jobs, "gcd", side_effect=AssertionError("Python GCD")
                )
            )
            stack.enter_context(
                patch.object(
                    ChainPlan,
                    "recover",
                    side_effect=AssertionError("Python recovery"),
                )
            )
            actual, _, stats = drive(case, "c-gmp-full")
        self.assertEqual(actual, expected)
        self.assertTrue(
            all(
                stats[name] > 0
                for name in ("setup_ops", "stage_one_ops", "stage_two_ops")
            )
        )

    def test_resume_each_canonical_phase_keeps_durable_progress(self):
        case = dict(n=2**199 - 189, seed=7, b1=31, b2=233)
        expected, _, _ = drive(case, "python")
        phases = set()

        def stop(job, packet):
            return job["phase"] not in phases

        job = None
        for _ in range(12):
            job, _, _ = drive(
                case, "c-gmp-full", job=job, atomic=True, stop=stop
            )
            phases.add(job["phase"])
            if job["done"]:
                break
        self.assertEqual(job, expected)
        self.assertTrue(
            {"stage_one", "baby_steps", "giant_setup", "stage_two"} <= phases
        )

    def test_cancel_after_real_native_work_records_debit_before_stop(self):
        config = selected("c-gmp-full", 31, 233)
        budget = deterministic_budget()
        context, programs = _execution_context({}, config, budget)
        source = programs or context
        packet = source.native_packet
        job = stage_jobs.new_job("ecm", 2**199 - 189, 7, 31, 233)
        packet.bind_job(job, budget, source, config)
        original = packet.lib

        class Library:
            def __getattr__(self, name):
                return getattr(original, name)

            def b4_step(self, *args):
                result = original.b4_step(*args)
                budget.cancelled = lambda: True
                return result

        packet.lib = Library()
        before = budget.used
        try:
            with self.assertRaises(BudgetExhaustedError):
                stage_jobs.advance_job(job, budget, source, config)
            self.assertGreater(budget.used, before)
            self.assertFalse(job["done"])
        finally:
            packet.lib = original
            close_execution_context(source)
        self.assertFalse(job["done"])
        resumed, _, _ = drive(
            dict(n=job["n"], seed=7, b1=31, b2=233), "c-gmp-full", job=job
        )
        control, _, _ = drive(
            dict(n=job["n"], seed=7, b1=31, b2=233), "python"
        )
        self.assertEqual(resumed, control)

    def test_table_limit_is_explicit_and_use_after_close_rejects(self):
        packet = NativeCurvePacket()
        packet.close()
        packet.close()
        with self.assertRaises(RuntimeError):
            packet.bind_job(
                stage_jobs.new_job("ecm", 1013, 7, 11, 47),
                deterministic_budget(),
                None,
                None,
            )

    def test_prac_strict_nonunits_and_exceptional_retry_in_c(
        self,
    ):
        case = dict(n=1009 * 1013 * (10**40 + 121), seed=17, b1=2000, b2=2500)
        base = replace(
            selected("python", 2000, 2500),
            ecm_tiers=((2000, 2500, 8),),
            ecm_chain_mode="auto",
        )
        import random

        for point in ([1009, 1], [0, 1], [1, 0]):
            job = stage_jobs.new_job("ecm", case["n"], 17, 2000, 2500)
            job.update(
                phase="stage_one",
                value=point,
                a24=2,
                sigma=random.Random(17).randrange(6, 2**63),
            )
            expected, trace, _ = drive(case, "python", job=job, config=base)
            native = replace(base, ecm_kernel="c-gmp-full")
            with patch.object(
                ChainPlan,
                "recover",
                side_effect=AssertionError("Python strict recovery"),
            ):
                actual, native_trace, _ = drive(
                    case, "c-gmp-full", job=job, atomic=True, config=native
                )
            self.assertEqual(actual, expected)
            self.assertEqual(native_trace, trace)
            self.assertGreater(actual.get("chain_replays", 0), 0)

    def test_stage_two_empty_ranges_small_distances_and_short_tails(self):
        for b1, b2 in (
            (2, 2),
            (2, 19),
            (3, 19),
            (5, 29),
            (11, 47),
            (31, 31),
            (31, 233),
        ):
            for batch in (1, 3, 64):
                case = dict(n=2**199 - 189, seed=7, b1=b1, b2=b2)
                control = replace(selected("python", b1, b2), gcd_batch=batch)
                expected, trace, _ = drive(case, "python", config=control)
                actual, native_trace, _ = drive(
                    case,
                    "c-gmp-full",
                    atomic=True,
                    config=replace(control, ecm_kernel="c-gmp-full"),
                )
                self.assertEqual(actual, expected)
                self.assertEqual(native_trace, trace)


@unittest.skipUnless(
    (SOURCE / "build/pthread/manifest.json").exists(),
    "pthread native build required",
)
class WholeThreadTests(unittest.TestCase):
    def environment(self, workers=16):
        config = selected("c-gmp-full", 31, 233)
        budget = deterministic_budget()
        context, programs = _execution_context({}, config, budget)
        source = programs or context
        team = NativeCurveTeam(workers=workers)
        self.addCleanup(close_execution_context, source)
        self.addCleanup(team.close)
        team._grow(team.capacity, 2**199 - 189, budget)
        job = stage_jobs.new_job("ecm", 2**199 - 189, 7, 31, 233)
        return config, budget, source, team, job

    def test_extended_openmp_canonical_fixture_actions(self):
        for case in json.loads(FIXTURES.read_text())["cases"]:
            expected, trace, _ = drive(case, "python")
            for count in (8, 10, 12, 14, 16):
                with self.subTest(case=case["id"], workers=count):
                    config = selected("c-gmp-full", case["b1"], case["b2"])
                    budget = deterministic_budget()
                    context, programs = _execution_context({}, config, budget)
                    source = programs or context
                    team = OpenMPCurveTeam(workers=count)
                    team.atomic_actions = True
                    job = stage_jobs.new_job(
                        "ecm", case["n"], case["seed"], case["b1"], case["b2"]
                    )
                    observed = []
                    try:
                        for _ in range(10000):
                            team.advance_jobs(
                                [(0, job)], budget, source, config
                            )
                            info = team.bindings[0]._read()
                            if info.boundary or info.done:
                                observed.append(signature(job))
                            if job["done"]:
                                break
                        else:
                            self.fail("OpenMP fixture did not complete")
                        self.assertEqual(job, expected)
                        self.assertEqual(observed, trace)
                        self.assertEqual(
                            team.worker_activity()["entered_mask"],
                            (1 << count) - 1,
                        )
                    finally:
                        team.close()
                        close_execution_context(source)

    def test_extended_count_bounds_and_actual_workers(self):
        for count in (0, 17, -1, True, 8.0):
            with self.subTest(invalid=count):
                with self.assertRaises(ValueError):
                    NativeCurveTeam(workers=count)
        for count in (8, 10, 12, 14, 16):
            with self.subTest(workers=count):
                config, budget, source, team, job = self.environment(count)
                self.assertEqual(
                    team.lib.b4_pthread_inject(team.owner, 16, -1, 0), -1
                )
                for _ in range(10000):
                    team.advance_jobs([(0, job)], budget, source, config)
                    if job["done"]:
                        break
                else:
                    self.fail("undersupplied team did not complete")
                team.bindings[0].sync(allow_result=True)
                activity = team.worker_activity()
                self.assertEqual(activity["actual"], count)
                self.assertEqual(activity["entered_mask"], (1 << count) - 1)
                ids = [row["identity"] for row in activity["workers"]]
                self.assertEqual(len(set(ids)), count)
                self.assertTrue(all(ids))
                self.assertEqual(
                    sum(row["curves"] for row in activity["workers"]), 1
                )
                self.assertEqual(
                    job,
                    drive(dict(n=job["n"], seed=7, b1=31, b2=233), "python")[
                        0
                    ],
                )
                team.close()
                self.assertEqual(team.lifecycle[:2], (count, count))

    def test_partial_startup_worker_and_metadata_failures_drain(self):
        for create_at, job_at, error in (
            (0, -1, TeamRefusedError),
            (2, -1, TeamRefusedError),
            (15, -1, TeamRefusedError),
            (-1, 0, RuntimeError),
        ):
            config, budget, source, team, job = self.environment()
            team.lib.b4_pthread_inject(team.owner, create_at, job_at, 0)
            with self.assertRaises(error):
                team.advance_jobs([(0, job)], budget, source, config)
            self.assertFalse(job["done"])
            if create_at >= 0:
                self.assertEqual(team.bindings[0]._read().setup_ops, 0)
                self.assertEqual(team.lib.b4_pthread_counter(team.owner, 2), 0)
                self.assertEqual(team.worker_activity()["actual"], 0)
            team.lib.b4_pthread_inject(team.owner, -1, -1, 0)
            for _ in range(10000):
                team.advance_jobs([(0, job)], budget, source, config)
                if job["done"]:
                    break
            else:
                self.fail("team did not recover after drained failure")
            team.bindings[0].sync(allow_result=True)
            self.assertEqual(
                job,
                drive(dict(n=job["n"], seed=7, b1=31, b2=233), "python")[0],
            )
            team.close()
            team.close()
        config, budget, source, team, job = self.environment()
        with patch.object(
            team.blocks,
            "register",
            side_effect=MemoryError("metadata injection"),
        ):
            with self.assertRaises(MemoryError):
                for _ in range(8):
                    team.advance_jobs([(0, job)], budget, source, config)
        self.assertFalse(job["done"])
        team.close()

    def test_signal_cancellation_inside_full_foreign_generation_and_reuse(
        self,
    ):
        config, budget, source, team, job = self.environment()
        team.advance_jobs([(0, job)], budget, source, config)
        self.assertGreater(team.bindings[0]._read().setup_ops, 0)
        team.lib.b4_pthread_inject(team.owner, -1, -1, 1)
        errors = []

        def execute():
            try:
                team.advance_jobs([(0, job)], budget, source, config)
            except BaseException as error:
                errors.append(error)

        caller = threading.Thread(target=execute)
        caller.start()
        try:
            team.lib.b4_pthread_wait(team.owner, 4)
            with self.assertRaises(RuntimeError):
                team.close()
            with self.assertRaises(RuntimeError):
                team.bindings[0].close()
            previous = signal.signal(signal.SIGUSR1, lambda *_: team.cancel())
            try:
                signal.raise_signal(signal.SIGUSR1)
            finally:
                signal.signal(signal.SIGUSR1, previous)
        finally:
            team.lib.b4_pthread_release(team.owner)
            caller.join(10)
        self.assertFalse(caller.is_alive())
        self.assertEqual(len(errors), 1)
        self.assertIsInstance(errors[0], BudgetExhaustedError)
        self.assertFalse(job["done"])
        team.lib.b4_pthread_inject(team.owner, -1, -1, 0)
        for _ in range(10000):
            team.advance_jobs([(0, job)], budget, source, config)
            if job["done"]:
                break
        self.assertTrue(job["done"])
        team.close()
        self.assertEqual(team.lifecycle[:2], (team.workers, team.workers))

    def test_repeated_cancelled_owner_resume_keeps_all_saved_history(self):
        config = replace(
            selected("c-gmp-full", 31, 233), ecm_tiers=((31, 233, 3),)
        )
        owner = NativeECMExecutor(seed=73, settings=WorkerSettings(workers=2))
        self.addCleanup(owner.close)
        checkpoint = None
        for _ in range(2):
            before = 0 if owner.team is None else owner.team.slices
            wall, cpu = Clock(), Clock()

            def cancel():
                if owner.team is not None and owner.team.slices > before:
                    wall.value, cpu.value = 2, 3
                    return True
                return False

            budget = Budget(
                wall_clock=wall,
                cpu_clock=cpu,
                provider=Provider(),
                cancelled=cancel,
            )
            paused = owner.run(
                2**199 - 189,
                config=config,
                budget=budget,
                checkpoint=checkpoint,
                fixed_work=True,
            )
            self.assertEqual(paused.reason, "cancelled")
            self.assertGreater(paused.stats["native_coverage"]["setup_ops"], 0)
            if checkpoint:
                for key in ("work", "wall", "cpu"):
                    self.assertGreater(
                        paused.checkpoint["payload"]["budget"][key],
                        checkpoint["payload"]["budget"][key],
                    )
            checkpoint = json.loads(json.dumps(paused.checkpoint))
        resumed = owner.run(
            2**199 - 189,
            config=config,
            budget=deterministic_budget(),
            checkpoint=checkpoint,
            fixed_work=True,
        )
        self.assertEqual(resumed.reason, "curves_exhausted")
        self.assertEqual(len(resumed.checkpoint["payload"]["completed"]), 3)
        self.assertEqual(resumed.checkpoint["payload"]["runs"], 3)

    def test_full_thread_root_wall_cpu_and_rss_stop_before_publication(self):
        for cause in ("wall", "cpu", "rss"):
            config = selected("c-gmp-full", 31, 233)
            wall, cpu, provider = Clock(), Clock(), Provider()
            budget = Budget(
                seconds=10,
                cpu_seconds=10,
                rss_bytes=500,
                wall_clock=wall,
                cpu_clock=cpu,
                provider=provider,
            )
            context, programs = _execution_context({}, config, budget)
            source = programs or context
            team = NativeCurveTeam(workers=4)
            team._grow(team.capacity, 2**199 - 189, budget)
            job = stage_jobs.new_job("ecm", 2**199 - 189, 7, 31, 233)
            original = team.lib

            class Library:
                def __getattr__(self, name):
                    return getattr(original, name)

                def b4_pthread_step(self, *args):
                    result = original.b4_pthread_step(*args)
                    if cause == "wall":
                        wall.value = 10
                    elif cause == "cpu":
                        cpu.value = 10
                    else:
                        provider.rows = [replace(provider.rows[0], rss=500)]
                    return result

            team.lib = Library()
            before = budget.used
            try:
                with self.assertRaisesRegex(
                    BudgetExhaustedError, cause + "_limit"
                ):
                    team.advance_jobs([(0, job)], budget, source, config)
                self.assertGreater(budget.used, before)
                self.assertGreater(team.bindings[0]._read().setup_ops, 0)
                self.assertFalse(job["done"])
            finally:
                team.lib = original
                team.close()
                close_execution_context(source)

    def test_simultaneous_setup_factors_use_stable_queue_publication(self):
        from v2.ecm.workers import arithmetic_options
        from v2.execution.worker_contract import Manifest

        config = replace(
            selected("c-gmp-full", 11, 47), ecm_tiers=((11, 47, 32),)
        )
        for workers in (1, 2, 4, 6, 8, 10, 12, 14, 16):
            settings = WorkerSettings(
                workers=workers,
                allow_wide_workers=workers not in (1, 2, 4),
            )
            manifest = Manifest(
                15, 73, "root", arithmetic_options(config, settings)
            )
            factors = [
                core.setup_curve(
                    15, manifest.assignment(index)["sigma"]
                ).factor
                for index in range(32)
            ]
            self.assertGreater(
                sum(factor is not None for factor in factors), 1
            )
            expected = next(factor for factor in factors if factor is not None)
            owner = NativeECMExecutor(seed=73, settings=settings, capacity=32)
            try:
                result = owner.run(
                    15, config=config, budget=deterministic_budget()
                )
                self.assertEqual(result.divisor, expected)
                self.assertEqual(result.reason, "factor")
                self.assertGreater(
                    result.stats["native_coverage"]["setup_ops"], 0
                )
            finally:
                owner.close()

    def test_all_fixture_actions_on_all_teams_and_schedules(
        self,
    ):
        for case in json.loads(FIXTURES.read_text())["cases"]:
            expected, trace, _ = drive(case, "python")
            for workers in (1, 2, 4, 6, 8, 10, 12, 14, 16):
                for schedule in ("static", "dynamic"):
                    for chunk in (1, 3):
                        with self.subTest(
                            case=case["id"],
                            workers=workers,
                            schedule=schedule,
                            chunk=chunk,
                        ):
                            config = selected(
                                "c-gmp-full", case["b1"], case["b2"]
                            )
                            budget = deterministic_budget()
                            context, programs = _execution_context(
                                {}, config, budget
                            )
                            source = programs or context
                            team = NativeCurveTeam(
                                workers=workers,
                                schedule=schedule,
                                chunk_size=chunk,
                            )
                            team.atomic_actions = True
                            job = stage_jobs.new_job(
                                "ecm",
                                case["n"],
                                case["seed"],
                                case["b1"],
                                case["b2"],
                            )
                            actual_trace = []
                            try:
                                for _ in range(10000):
                                    team.advance_jobs(
                                        [(0, job)], budget, source, config
                                    )
                                    packet = team.bindings[0]
                                    info = packet._read()
                                    if info.boundary or info.done:
                                        actual_trace.append(signature(job))
                                    if job["done"]:
                                        break
                                else:
                                    self.fail("whole team did not finish")
                                self.assertEqual(job, expected)
                                self.assertEqual(actual_trace, trace)
                                self.assertGreater(team.slices, 0)
                                self.assertGreater(
                                    packet.statistics["setup_ops"], 0
                                )
                            finally:
                                team.close()
                                close_execution_context(source)
                            self.assertEqual(
                                team.lifecycle[:2], (workers, workers)
                            )

    def test_second_tier_logical_cancel_preserves_completed_prefix(self):
        from v2.ecm.workers import ECMExecutor

        config = replace(
            selected("c-gmp-full", 31, 233),
            ecm_tiers=((31, 233, 3), (43, 331, 5)),
        )
        n = 2**199 - 189
        cancelled = [False]
        owner = NativeECMExecutor(
            seed=73, settings=WorkerSettings(workers=4), capacity=4
        )
        advance = NativeCurveTeam.advance_jobs

        def stop_after_second_tier_generation(team, jobs, *args):
            ready = advance(team, jobs, *args)
            if jobs[0][1]["b1"] == 43:
                cancelled[0] = True
            return ready

        try:
            budget = deterministic_budget()
            budget.cancelled = lambda: cancelled[0]
            with patch.object(
                NativeCurveTeam,
                "advance_jobs",
                stop_after_second_tier_generation,
            ):
                paused = owner.run(
                    n, config=config, budget=budget, fixed_work=True
                )
            self.assertEqual(paused.reason, "cancelled")
            prefix = paused.checkpoint["payload"]["completed"]
            self.assertEqual(list(prefix), ["0", "1", "2"])
            self.assertEqual(paused.stats["lost_chunks"], 4)
            self.assertGreater(
                paused.checkpoint["payload"]["budget"]["work"], 0
            )

            resumed = owner.run(
                n,
                config=config,
                budget=deterministic_budget(),
                checkpoint=json.loads(json.dumps(paused.checkpoint)),
                fixed_work=True,
            )
            self.assertEqual(resumed.reason, "curves_exhausted")
            self.assertEqual(
                {
                    k: resumed.checkpoint["payload"]["completed"][k]
                    for k in prefix
                },
                prefix,
            )
            self.assertGreaterEqual(
                resumed.checkpoint["payload"]["budget"]["work"],
                paused.checkpoint["payload"]["budget"]["work"],
            )
        finally:
            owner.close()

        with ECMExecutor(
            seed=73, settings=WorkerSettings(workers=1)
        ) as control:
            expected = control.run(
                n,
                config=replace(config, ecm_kernel="python"),
                budget=deterministic_budget(),
                fixed_work=True,
            )
        actual_rows = resumed.checkpoint["payload"]["completed"]
        expected_rows = expected.checkpoint["search"]["payload"]["completed"]
        self.assertEqual(
            {
                k: {field: row[field] for field in ("id", "sigma", "factor")}
                for k, row in actual_rows.items()
            },
            {
                k: {field: row[field] for field in ("id", "sigma", "factor")}
                for k, row in expected_rows.items()
            },
        )

    def test_full_owner_matches_python_assignments_and_saved_history(self):
        from v2.ecm.workers import ECMExecutor

        config = selected("c-gmp-full", 31, 233)
        config = replace(config, ecm_tiers=((31, 233, 37),))
        control = ECMExecutor(seed=73, settings=WorkerSettings(workers=1))
        try:
            expected = control.run(
                2**199 - 189,
                config=replace(config, ecm_kernel="python"),
                budget=deterministic_budget(),
                fixed_work=True,
            )
        finally:
            control.close()
        for workers in (1, 2, 4, 6, 8, 10, 12, 14, 16):
            for schedule in ("static", "dynamic"):
                owner = NativeECMExecutor(
                    seed=73,
                    settings=WorkerSettings(
                        workers=workers,
                        allow_wide_workers=workers not in (1, 2, 4),
                    ),
                    schedule=schedule,
                )
                try:
                    actual = owner.run(
                        2**199 - 189,
                        config=config,
                        budget=deterministic_budget(),
                        fixed_work=True,
                    )
                    self.assertEqual(
                        {
                            k: {
                                field: row[field]
                                for field in ("id", "sigma", "factor")
                            }
                            for k, row in actual.checkpoint["payload"][
                                "completed"
                            ].items()
                        },
                        {
                            k: {
                                field: row[field]
                                for field in ("id", "sigma", "factor")
                            }
                            for k, row in expected.checkpoint["search"][
                                "payload"
                            ]["completed"].items()
                        },
                    )
                    self.assertTrue(
                        all(
                            value > 0
                            for value in actual.stats[
                                "native_coverage"
                            ].values()
                        )
                    )
                    resumed = owner.run(
                        2**199 - 189,
                        config=config,
                        budget=deterministic_budget(),
                        fixed_work=True,
                        checkpoint=json.loads(json.dumps(actual.checkpoint)),
                    )
                    self.assertEqual(
                        resumed.checkpoint["payload"]["completed"],
                        actual.checkpoint["payload"]["completed"],
                    )
                    self.assertGreaterEqual(
                        resumed.checkpoint["payload"]["budget"]["work"],
                        actual.checkpoint["payload"]["budget"]["work"],
                    )
                finally:
                    owner.close()
