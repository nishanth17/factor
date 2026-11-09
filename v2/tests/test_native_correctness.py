"""Fixed-curve native admission and deterministic H3 interruption checks."""

import copy
import json
import os
import unittest
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from v2.common.arithmetic import BackendUnavailableError
from v2.ecm.native_packet import SOURCE, NativePacket
from v2.ecm.native_team import NativeTeam
from v2.execution import stage_jobs
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.portfolio_config import (
    PortfolioConfig,
    _execution_context,
    close_execution_context,
)
from v2.execution.worker_contract import digest
from v2.portfolio import factorize_bounded
from v2.tests.test_h3_resources import Clock, Provider

FIXTURES = (
    Path(__file__).resolve().parents[1]
    / "benchmarks/inputs/controls/b4a_native_correctness_v1.json"
)


def selected(kernel="c-gmp", b1=11, b2=47):
    return PortfolioConfig(
        trial_bound=2,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=((b1, b2, 1),),
        ecm_chain_mode="off",
        ecm_kernel=kernel,
        siqs=None,
        segment_size=16,
        gcd_batch=4,
    )


def deterministic_budget(**changes):
    return Budget(
        wall_clock=lambda: 0,
        cpu_clock=lambda: 0,
        provider=Provider(),
        **changes,
    )


@contextmanager
def after_native_slice(packet, callback):
    """Inject a stop after real C work and before Python publication."""
    original = packet.lib

    class Library:
        def __getattr__(self, name):
            return getattr(original, name)

        def b4_step(self, *args):
            result = original.b4_step(*args)
            callback()
            return result

    packet.lib = Library()
    try:
        yield
    finally:
        packet.lib = original


def fixed_curve(case, kernel):
    config = selected(kernel, case["b1"], case["b2"])
    budget = deterministic_budget()
    context, programs = _execution_context({}, config, budget)
    job = stage_jobs.new_job(
        "ecm", case["n"], case["seed"], case["b1"], case["b2"]
    )
    trace, saturated, calls = [], [], []
    gcd, finish = stage_jobs.gcd, NativePacket._finish

    def observe_gcd(value, modulus):
        result = gcd(value, modulus)
        if result == modulus:
            saturated.append(job["phase"])
        return result

    def observe_finish(packet, allowance):
        result = finish(packet, allowance)
        calls.append(True)
        return result

    try:
        with (
            patch.object(stage_jobs, "gcd", observe_gcd),
            patch.object(NativePacket, "_finish", observe_finish),
        ):
            for _ in range(256):
                stage_jobs.advance_job(
                    job, budget, programs or context, config
                )
                trace.append(copy.deepcopy(job))
                if job["done"]:
                    return job, trace, saturated, len(calls)
        raise AssertionError("fixed correctness recipe did not finish")
    finally:
        close_execution_context(programs or context)
        budget.close()


class ExplicitBackendTests(unittest.TestCase):
    def test_missing_explicit_team_fails_without_scalar_fallback(self):
        owner = NativeTeam()
        self.addCleanup(owner.close)
        with patch("v2.ecm.native_packet.SOURCE", SOURCE / "missing"):
            with self.assertRaisesRegex(
                BackendUnavailableError, "OpenMP team backend unavailable"
            ):
                owner.execute(
                    [("ladder", 13, (2, 1), 1013, 2)],
                    deterministic_budget(),
                    lambda *_: False,
                )
        self.assertEqual(owner.slices, 0)
        self.assertEqual(owner.packets, [])


@unittest.skipUnless(
    (SOURCE / "build/manifest.json").exists(), "scalar native build required"
)
class NativeCorrectnessTests(unittest.TestCase):
    def test_fixed_full_width_curves_and_nonunit_recovery(self):
        rows = []
        for case in json.loads(FIXTURES.read_text())["cases"]:
            with self.subTest(case=case["id"]):
                expected, trace, saturation, _ = fixed_curve(case, "python")
                actual, native_trace, native_saturation, calls = fixed_curve(
                    case, "c-gmp"
                )
                self.assertEqual(native_trace, trace)
                self.assertEqual(actual, expected)
                self.assertEqual(actual["factor"], case["factor"])
                if "digits" in case:
                    self.assertEqual(len(str(case["n"])), case["digits"])
                    self.assertEqual(case["n"] % case["known_factor"], 0)
                self.assertEqual(native_saturation, saturation)
                self.assertEqual(
                    actual["phase"],
                    "stage_two"
                    if case["stage"] == "not_found"
                    else case["stage"],
                )
                self.assertGreater(calls, 0)
                if case.get("saturated"):
                    self.assertTrue(saturation)
                if actual["factor"] is not None:
                    self.assertTrue(1 < actual["factor"] < case["n"])
                    self.assertEqual(case["n"] % actual["factor"], 0)
                rows.append(
                    dict(
                        case=case["id"],
                        digits=len(str(case["n"])),
                        seed=case["seed"],
                        b1=case["b1"],
                        b2=case["b2"],
                        sigma=actual.get("sigma"),
                        factor=actual["factor"],
                        phase=actual["phase"],
                        native_packets=calls,
                        saturation=saturation,
                        state_sha256=digest(actual),
                        all_boundaries_equal=True,
                    )
                )
        output = os.environ.get("FACTOR_B4A_CORRECTNESS_OUTPUT")
        if output:
            with Path(output).open("a") as stream:
                stream.write(json.dumps(dict(exactness=rows)) + "\n")

    def test_deadline_stop_flag_and_rss_stop_after_real_native_work(self):
        for cause in ("deadline", "flag", "rss"):
            with self.subTest(cause=cause):
                wall, cpu, provider = Clock(), Clock(), Provider()
                stop, slices = [False], []
                budget = Budget(
                    seconds=10,
                    rss_bytes=500,
                    wall_clock=wall,
                    cpu_clock=cpu,
                    provider=provider,
                    cancelled=lambda: stop[0],
                )
                packet = NativePacket()
                self.addCleanup(packet.close)
                packet.available(1013, budget)

                def inject():
                    slices.append(True)
                    if cause == "deadline":
                        wall.value = 10
                    elif cause == "flag":
                        stop[0] = True
                    else:
                        provider.rows = [replace(provider.rows[0], rss=500)]

                with after_native_slice(packet, inject):
                    with self.assertRaisesRegex(
                        BudgetExhaustedError,
                        dict(
                            deadline="wall_limit",
                            flag="cancelled",
                            rss="rss_limit",
                        )[cause],
                    ):
                        packet.ladder(2**4095 + 19, (2, 1), 1013, 2, budget)
                self.assertEqual(len(slices), 1)
                self.assertFalse(packet._lock.locked())
                self.assertEqual(
                    packet.lib.b4_step(packet.context, 5, 128), -1
                )
                self.assertEqual(
                    packet.ladder(13, (2, 1), 1013, 2, deterministic_budget()),
                    stage_jobs.ecm.scalar_multiply(13, 2, 1, 1013, 2),
                )
                packet.close()
                self.assertIsNone(packet.context)
                self.assertEqual(packet.buffers, [])
                self.assertEqual(packet.programs, {})

    def test_deadline_preserves_job_and_resume_matches_uninterrupted(self):
        case = json.loads(FIXTURES.read_text())["cases"][0]
        config = selected()
        wall, cpu = Clock(), Clock()
        budget = Budget(
            seconds=10, wall_clock=wall, cpu_clock=cpu, provider=Provider()
        )
        context, programs = _execution_context({}, config, budget)
        packet = (programs or context).native_packet
        self.addCleanup(close_execution_context, programs or context)
        packet.available(case["n"], budget)
        original = stage_jobs.new_job(
            "ecm", case["n"], case["seed"], case["b1"], case["b2"]
        )
        job, paused = copy.deepcopy(original), False
        saved = None

        def expire():
            wall.value = 10
            cpu.value = 3  # Process CPU includes every native thread once.

        for _ in range(256):
            saved = copy.deepcopy(job)
            try:
                with after_native_slice(packet, expire):
                    stage_jobs.advance_job(
                        job, budget, programs or context, config
                    )
            except BudgetExhaustedError:
                paused = True
                break
        self.assertTrue(paused)
        self.assertEqual(job, saved)
        history = budget.charges()
        self.assertEqual((history["wall"], history["cpu"]), (10, 3))
        used = budget.used
        for _ in range(2):
            budget.restore_charges(work=used, wall=10, cpu=3)
        self.assertEqual(budget.used, used)
        self.assertEqual((budget.wall_used, budget.cpu_used), (10, 3))
        budget.seconds = 20
        for _ in range(256):
            stage_jobs.advance_job(job, budget, programs or context, config)
            if job["done"]:
                break
        expected = fixed_curve(case, "python")[0]
        self.assertTrue(job["done"])
        self.assertEqual(job, expected)
        self.assertEqual(budget.prior_cpu, 0)  # No child/thread surcharge.

    def test_portfolio_checkpoint_resume_after_native_deadline(self):
        case = json.loads(FIXTURES.read_text())["cases"][0]
        config = selected()
        wall, cpu = Clock(), Clock()
        budget = Budget(
            seconds=10, wall_clock=wall, cpu_clock=cpu, provider=Provider()
        )
        finish = NativePacket._finish

        def expire(packet, allowance):
            with after_native_slice(
                packet, lambda: setattr(wall, "value", 10)
            ):
                return finish(packet, allowance)

        with patch.object(NativePacket, "_finish", expire):
            paused = factorize_bounded(
                case["n"], seed=1, config=config, budget=budget
            )
        self.assertEqual(paused.reason, "wall_limit")
        resumed = factorize_bounded(
            case["n"],
            seed=1,
            config=config,
            budget=deterministic_budget(),
            checkpoint=json.loads(json.dumps(paused.checkpoint)),
        )
        expected = factorize_bounded(
            case["n"],
            seed=1,
            config=selected("python"),
            budget=deterministic_budget(),
        )
        self.assertEqual(resumed.result, expected.result)
        self.assertEqual(resumed.result.reconstruct(), case["n"])

    def test_repeated_scalar_close_releases_buffers_without_new_threads(self):
        import ctypes
        import multiprocessing

        from v2.benchmarks.ecm.b4a.ownership import GroupObserver, TaskInfo

        observer = GroupObserver()

        def thread_count():
            info = TaskInfo()
            size = observer.lib.proc_pidinfo(
                os.getpid(), 4, 0, ctypes.byref(info), ctypes.sizeof(info)
            )
            self.assertEqual(size, ctypes.sizeof(info))
            return info.thread_count

        before = thread_count()
        children = [child.pid for child in multiprocessing.active_children()]
        for _ in range(4):
            packet = NativePacket()
            packet.ladder(13, (2, 1), 1013, 2, deterministic_budget())
            packet.close()
            packet.close()
            self.assertIsNone(packet.context)
            self.assertIsNone(packet.ffi)
            self.assertEqual(packet.buffers, [])
            self.assertEqual(packet.output_views, {})
            self.assertEqual(packet.store.used, 0)
            self.assertEqual(thread_count(), before)
        self.assertEqual(
            [child.pid for child in multiprocessing.active_children()],
            children,
        )

    def test_scalar_build_refuses_all_team_sizes_and_schedules(self):
        packet = NativePacket()
        self.addCleanup(packet.close)
        allowance = deterministic_budget()
        packet.available(1013, allowance)
        self.assertEqual(packet.lib.b4_openmp(), 0)
        ffi = packet.ffi
        owners = [
            ffi.new("b4_context *[1]"),
            ffi.new("int[1]"),
            ffi.new("unsigned *"),
        ]
        for owner in owners:
            self.addCleanup(ffi.release, owner)
        contexts, statuses, actual = owners
        for workers in (1, 2, 4):
            for schedule in range(3):
                packet.prepare_ladder(13, (2, 1), 1013, 2, allowance)
                contexts[0] = packet.context
                try:
                    self.assertEqual(
                        packet.lib.b4_team_step(
                            contexts,
                            1,
                            workers,
                            schedule,
                            5,
                            128,
                            statuses,
                            actual,
                        ),
                        -2,
                    )
                    self.assertEqual(actual[0], 0)
                    self.assertEqual(
                        packet.lib.b4_result(packet.context, packet.output), -1
                    )
                finally:
                    packet.end()
