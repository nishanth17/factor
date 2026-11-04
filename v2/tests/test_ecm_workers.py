"""Queued arithmetic/process QA. NEVER part of the pure deterministic lane."""

import copy
import multiprocessing
import os
import time
import unittest
from dataclasses import replace
from math import gcd, prod
from unittest.mock import Mock, patch

from v2.benchmarks.support.prac_oracle import (
    affine_multiply,
    matches,
    twist_point,
)
from v2.ecm import core as ecm
from v2.execution.budget import Budget
from v2.execution.ecm_worker import CurveEngine
from v2.execution.ecm_workers import ECMExecutor, arithmetic_options
from v2.execution.stage_jobs import advance_job, new_job
from v2.execution.worker_contract import (
    LeaseRefusedError,
    Ledger,
    Manifest,
    WorkerSettings,
)
from v2.portfolio import PortfolioConfig, _execution_context, factorize_bounded
from v2.tests.test_ecm_programs import primes_between


def allowance(**changes):
    options = dict(work_limit=200_000_000, seconds=60, cpu_seconds=240)
    options.update(changes)
    return Budget(**options)


def configuration(**changes):
    return replace(
        PortfolioConfig(
            trial_bound=5,
            rho_attempts=0,
            pm1_attempts=0,
            siqs=None,
            ecm_tiers=((20, 200, 4),),
            segment_size=16,
            max_input_bits=512,
            memory_bytes=16 * 2**20,
        ),
        **changes,
    )


def task(plan, position, attempt=1):
    return dict(
        assignment=plan.assignment(position),
        options=plan.options,
        schedule=plan.metadata,
        token=position + 1,
        attempt=attempt,
        work=100_000_000,
        cpu=30,
        deadline=time.monotonic() + 45,
    )


class ArithmeticWorkerTests(unittest.TestCase):
    def test_identical_manifest_serial_one_two_four_workers(self):
        n = 1009 * 1013
        config = configuration()
        expected = None
        before = {p.pid for p in multiprocessing.active_children()}
        for workers, forced in ((1, False), (1, True), (2, True), (4, True)):
            with ECMExecutor(
                seed=7919,
                settings=WorkerSettings(workers=workers),
                force_workers=forced,
            ) as executor:
                result = executor.run(
                    n, config=config, budget=allowance(), fixed_work=True
                )
                records = result.checkpoint["search"]["payload"]["completed"]
                actual = {
                    int(i): (m["id"], m["factor"], m["sigma"], m["phase"])
                    for i, m in records.items()
                }
                self.assertEqual(len(actual), 4)
                for _, divisor, _, _ in actual.values():
                    if divisor is not None:
                        self.assertTrue(1 < divisor < n and n % divisor == 0)
                if expected is not None:
                    self.assertEqual(actual, expected)
                expected = actual
        self.assertEqual(
            {p.pid for p in multiprocessing.active_children()}, before
        )

    def test_stage_one_schedule_against_independent_lcm_ladder(self):
        config = configuration(ecm_tiers=((10, 10, 2),), chunk_size=1)
        plan = Manifest(
            1000003,
            7919,
            "prime-control",
            arithmetic_options(config, WorkerSettings()),
        )
        job = new_job("ecm", plan.n, 6, 10, 10)
        budget = allowance()
        context, programs = _execution_context({}, config, budget)
        advance_job(job, budget, programs or context, config)
        setup = ecm.setup_curve(plan.n, job["sigma"])
        scalar = prod(
            max(p**e for e in range(1, 7) if p**e <= 10)
            for p in primes_between(2, 11)
        )
        point, curve_a, curve_b = twist_point(setup.point, setup.a24, plan.n)
        expected = affine_multiply(scalar, point, plan.n, curve_a, curve_b)
        while job["phase"] == "stage_one" and not job["done"]:
            advance_job(job, budget, programs or context, config)
        self.assertTrue(matches(job["value"], expected, plan.n))
        # Scalar construction is independent of the production prime cursor.
        self.assertEqual(gcd(job["value"][1], plan.n), 1)

    def test_setup_factor_retry_both_stages_and_saturation(self):
        engine = CurveEngine()
        config = configuration(gcd_batch=2)
        plan = Manifest(
            35, 7919, "setup", arithmetic_options(config, WorkerSettings())
        )
        # Constructor fields are read from the actual setup result type.
        template = ecm.setup_curve(1000003, 17)
        with patch(
            "v2.ecm.core.setup_curve", return_value=replace(template, factor=5)
        ):
            result = engine.execute(task(plan, 0), cancelled=lambda: False)
        self.assertEqual(result["factor"], 5)
        with patch(
            "v2.ecm.core.setup_curve",
            return_value=replace(template, factor=None, retry=True),
        ):
            result = engine.execute(task(plan, 1), cancelled=lambda: False)
        self.assertEqual(result["reason"], "complete")
        self.assertIsNone(result["factor"])
        self.assertEqual(engine.preparations, 1)

        # Saturation n=5*7 must recover an individual factor, not publish n.
        job = new_job("ecm", 35, 0, 20, 200)
        job.update(phase="stage_two", terms=[5, 7], product=0)
        context, programs = _execution_context({}, config, allowance())
        advance_job(job, allowance(), programs or context, config)
        self.assertEqual(job["phase"], "term_replay")
        advance_job(job, allowance(), programs or context, config)
        self.assertEqual(job["factor"], 5)

    def test_full_schedule_prac_eligibility_and_once_per_context(self):
        config = configuration(
            ecm_tiers=((2000, 2000, 8),), ecm_chain_mode="auto"
        )
        plan = Manifest(
            10**39 + 7,
            7919,
            "PRAC",
            arithmetic_options(config, WorkerSettings()),
        )
        self.assertEqual(plan.options["ecm_chain_mode"], "reuse")
        self.assertEqual(plan.options["ecm_tiers"][0][2], 8)
        engine = CurveEngine()
        first = engine.execute(task(plan, 0), cancelled=lambda: False)
        second = engine.execute(task(plan, 1), cancelled=lambda: False)
        self.assertEqual(engine.preparations, 1)
        self.assertEqual(
            (first["reason"], second["reason"]), ("complete", "complete")
        )
        # Eligibility does not require every curve to reach a chain chunk.
        changed = Manifest(
            plan.n,
            7919,
            "PRAC",
            dict(plan.options, ecm_tiers=[[2000, 2100, 8]]),
        )
        engine.execute(task(changed, 0), cancelled=lambda: False)
        self.assertEqual(engine.preparations, 2)

    def test_cancellation_and_refusal_before_setup(self):
        with ECMExecutor() as executor:
            cancelled = executor.run(
                1009 * 1013,
                config=configuration(),
                budget=allowance(cancelled=lambda: True),
            )
            self.assertEqual(cancelled.reason, "cancelled")
            self.assertFalse(executor.slots)
        with ECMExecutor(force_workers=True) as executor:
            expired = executor.run(
                1009 * 1013,
                config=configuration(),
                budget=allowance(seconds=0),
            )
            self.assertEqual(expired.reason, "wall_limit")
            self.assertFalse(executor.slots)

    def test_drained_restart_charge_and_count_rejection(self):
        config = configuration()
        with ECMExecutor(seed=7, force_workers=True) as executor:
            partial = executor.run(
                1009 * 1013,
                config=config,
                budget=allowance(cancelled=lambda: True),
            )
        checkpoint = copy.deepcopy(partial.checkpoint)
        with ECMExecutor(
            seed=7, settings=WorkerSettings(workers=2), force_workers=True
        ) as changed:
            with self.assertRaises(ValueError):
                changed.run(
                    1009 * 1013,
                    config=config,
                    budget=allowance(),
                    checkpoint=checkpoint,
                )
        with ECMExecutor(seed=7, force_workers=True) as resumed:
            result = resumed.run(
                1009 * 1013,
                config=config,
                budget=allowance(),
                checkpoint=checkpoint,
            )
            self.assertGreaterEqual(
                result.checkpoint["budget"]["work"],
                checkpoint["budget"]["work"],
            )

    def test_portfolio_reconstruction_preprocessing_and_checkpoint_marker(
        self,
    ):
        config = configuration()
        with ECMExecutor(seed=7) as executor:
            run = factorize_bounded(
                2**8 * 3**3 * 1009,
                config=config,
                seed=29,
                budget=allowance(),
                ecm_executor=executor,
            )
            self.assertFalse(executor.slots)
            self.assertEqual(
                prod(f.value**f.exponent for f in run.result.factors)
                * prod(run.result.remaining),
                2**8 * 3**3 * 1009,
            )
            with self.assertRaises(ValueError):
                factorize_bounded(
                    2**8 * 3**3 * 1009,
                    config=config,
                    budget=allowance(),
                    checkpoint=run.checkpoint,
                )


def _fake_worker(channel, frame_bytes, startup_work, mode):
    from v2.execution.ecm_worker import read_message
    from v2.execution.worker_contract import Decoder, frame

    channel.settimeout(0.01)
    channel.sendall(
        frame(
            dict(
                kind="ready", cpu=time.process_time(), rss=1, work=startup_work
            ),
            frame_bytes,
        )
    )
    decoder = Decoder(frame_bytes)
    deadline = time.monotonic() + 3
    job = None
    while job is None and time.monotonic() < deadline:
        job = read_message(channel, decoder, block=True)
    if job is None:
        os._exit(19)
    packet = job
    job = packet["jobs"][0]
    if mode == "truncated":
        channel.sendall(b"\x00\x00\x00\x50{")
    elif mode == "oversized":
        channel.sendall((frame_bytes + 1).to_bytes(4, "big"))
    elif mode == "ack":
        channel.sendall(
            frame(
                dict(
                    kind="result",
                    packet=packet["packet"],
                    stage=job["stage"],
                    prepared_schedule=packet["schedule"],
                    metadata_cache=[packet["schedule"]],
                    id=job["assignment"]["id"],
                    token=job["token"],
                    attempt=job["attempt"],
                    work=1,
                    cpu=0,
                    rss=1,
                    reason="complete",
                    factor=5,
                    sigma=job["assignment"]["sigma"],
                ),
                frame_bytes,
            )
        )
    os._exit(17)


def crash_before_report(
    channel, frame_bytes, startup_work, settings_data=None
):
    _fake_worker(channel, frame_bytes, startup_work, "crash")


def publish_then_crash(channel, frame_bytes, startup_work, settings_data=None):
    _fake_worker(channel, frame_bytes, startup_work, "ack")


def truncate_then_crash(
    channel, frame_bytes, startup_work, settings_data=None
):
    _fake_worker(channel, frame_bytes, startup_work, "truncated")


def oversize_then_crash(
    channel, frame_bytes, startup_work, settings_data=None
):
    _fake_worker(channel, frame_bytes, startup_work, "oversized")


def exit_during_startup(
    channel, frame_bytes, startup_work, settings_data=None
):
    os._exit(23)


class StartupOwnershipTests(unittest.TestCase):
    """Fake startup; queued here because this module imports arithmetic."""

    def prepare(self, workers=4, startup_seconds=10):
        executor = ECMExecutor(
            settings=WorkerSettings(
                workers=workers, startup_seconds=startup_seconds
            )
        )
        events, processes = [], []

        def process_factory(**options):
            index = len(processes)
            process = Mock()
            process.start.side_effect = lambda: events.append(index)
            process.is_alive.return_value = True
            processes.append(process)
            return process

        context = Mock()
        context.get_start_method.return_value = "spawn"
        context.Process.side_effect = process_factory
        return executor, context, events, processes

    def test_launches_overlap_and_ready_reports_can_arrive_out_of_order(self):
        executor, context, started, _ = self.prepare()
        ledger = Ledger(20_000_000, 100)
        reads = []

        def read(slot):
            self.assertEqual(started, [0, 1, 2, 3])
            reads.append(slot["index"])
            if reads == [0]:
                return None
            return dict(kind="ready", work=1_000_000, cpu=0.125, rss=1)

        with (
            patch(
                "v2.execution.ecm_workers.multiprocessing.get_context",
                return_value=context,
            ),
            patch(
                "v2.execution.ecm_workers.socket.socketpair",
                side_effect=lambda: (Mock(), Mock()),
            ),
            patch.object(executor, "_read", side_effect=read),
        ):
            executor._spawn(ledger, time.monotonic() + 60, lambda: None)

        self.assertEqual(reads, [0, 1, 2, 3, 0])
        self.assertFalse(ledger.active)
        self.assertEqual(ledger.work, 5_000_000)
        self.assertEqual(
            ledger.cpu, 14 + 4 * (0.125 + executor.settings.shutdown_cpu)
        )
        self.assertEqual(
            [s["reported_cpu"] for s in executor.slots], [0.125] * 4
        )
        executor.slots.clear()  # Mock handles only; no children were started.

    def test_constructor_failure_closes_the_unlaunched_socket_pair(self):
        executor, context, _, _ = self.prepare()
        context.Process.side_effect = OSError("injected construction failure")
        parent, child = Mock(), Mock()
        ledger = Ledger(20_000_000, 100)
        with (
            patch(
                "v2.execution.ecm_workers.multiprocessing.get_context",
                return_value=context,
            ),
            patch(
                "v2.execution.ecm_workers.socket.socketpair",
                return_value=(parent, child),
            ),
        ):
            with self.assertRaises(OSError):
                executor._spawn(ledger, time.monotonic() + 60, lambda: None)
        parent.close.assert_called_once()
        child.close.assert_called_once()
        self.assertFalse(executor.slots)
        self.assertEqual(ledger.lost, 1)
        self.assertFalse(ledger.active)

    def test_bad_ready_keeps_every_unreported_startup_lease(self):
        executor, context, _, _ = self.prepare()
        ledger = Ledger(20_000_000, 100)
        with (
            patch(
                "v2.execution.ecm_workers.multiprocessing.get_context",
                return_value=context,
            ),
            patch(
                "v2.execution.ecm_workers.socket.socketpair",
                side_effect=lambda: (Mock(), Mock()),
            ),
            patch.object(executor, "_read", return_value=dict(kind="result")),
        ):
            with self.assertRaisesRegex(ValueError, "ready report"):
                executor._spawn(ledger, time.monotonic() + 60, lambda: None)

        self.assertFalse(ledger.active)
        self.assertEqual(ledger.lost, 4)
        self.assertEqual(
            ledger.cpu, 14 + 4 * (10 + executor.settings.shutdown_cpu)
        )
        self.assertEqual(ledger.work, 5_000_000)
        executor.slots.clear()

    def test_one_ready_report_refunds_only_that_child_before_peer_failure(
        self,
    ):
        executor, context, _, _ = self.prepare()
        ledger = Ledger(20_000_000, 100)
        messages = [
            dict(kind="ready", work=1_000_000, cpu=0.125, rss=1),
            dict(kind="result"),
        ]
        with (
            patch(
                "v2.execution.ecm_workers.multiprocessing.get_context",
                return_value=context,
            ),
            patch(
                "v2.execution.ecm_workers.socket.socketpair",
                side_effect=lambda: (Mock(), Mock()),
            ),
            patch.object(executor, "_read", side_effect=messages),
        ):
            with self.assertRaisesRegex(ValueError, "ready report"):
                executor._spawn(ledger, time.monotonic() + 60, lambda: None)

        self.assertFalse(ledger.active)
        self.assertEqual(ledger.lost, 3)
        self.assertEqual(
            ledger.cpu, 14 + 0.125 + 30 + 4 * executor.settings.shutdown_cpu
        )
        self.assertEqual(ledger.work, 5_000_000)
        executor.slots.clear()

    def test_startup_deadline_is_not_reset_after_launching_other_workers(self):
        executor, context, _, processes = self.prepare(
            workers=2, startup_seconds=1
        )
        ledger = Ledger(20_000_000, 100)
        clock = [0.0]
        factory = context.Process.side_effect

        def slow_start(**options):
            process = factory(**options)

            def advance_clock():
                clock[0] += 0.6

            process.start.side_effect = advance_clock
            return process

        context.Process.side_effect = slow_start
        with (
            patch(
                "v2.execution.ecm_workers.multiprocessing.get_context",
                return_value=context,
            ),
            patch(
                "v2.execution.ecm_workers.socket.socketpair",
                side_effect=lambda: (Mock(), Mock()),
            ),
            patch(
                "v2.execution.ecm_workers.time.monotonic",
                side_effect=lambda: clock[0],
            ),
            patch.object(executor, "_read") as read,
        ):
            with self.assertRaisesRegex(LeaseRefusedError, "startup_timeout"):
                executor._spawn(ledger, 60, lambda: None)

        self.assertEqual(len(processes), 2)
        read.assert_not_called()
        self.assertEqual(ledger.lost, 2)
        self.assertFalse(ledger.active)
        executor.slots.clear()


class ProcessFailureTests(unittest.TestCase):
    def test_peer_startup_exit_closes_all_owned_children(self):
        original = {p.pid for p in multiprocessing.active_children()}
        settings = WorkerSettings(workers=4)
        with patch("v2.execution.ecm_worker.worker_main", exit_during_startup):
            with ECMExecutor(
                settings=settings, force_workers=True
            ) as executor:
                result = executor.run(
                    35, config=configuration(), budget=allowance()
                )
        self.assertEqual(result.reason, "worker_failure")
        self.assertEqual(result.stats["worker_dispatches"], 0)
        self.assertGreaterEqual(result.stats["lost_reports"], 4)
        self.assertGreaterEqual(
            result.stats["charged_child_cpu"],
            4 * (settings.startup_seconds + settings.shutdown_cpu),
        )
        self.assertEqual(
            {p.pid for p in multiprocessing.active_children()}, original
        )

    def test_process_failures_drain_and_conservatively_charge(self):
        original = {p.pid for p in multiprocessing.active_children()}
        for target in (
            crash_before_report,
            truncate_then_crash,
            oversize_then_crash,
            publish_then_crash,
        ):
            with patch("v2.execution.ecm_worker.worker_main", target):
                with ECMExecutor(force_workers=True) as executor:
                    result = executor.run(
                        35, config=configuration(), budget=allowance()
                    )
            self.assertEqual(
                {p.pid for p in multiprocessing.active_children()}, original
            )
            self.assertGreater(result.stats["charged_work"], 0)
            self.assertFalse(
                result.checkpoint["search"]["payload"].get("active")
            )
            if result.divisor is not None:
                self.assertEqual(result.divisor, 5)
            else:
                self.assertEqual(result.reason, "worker_failure")
                self.assertGreater(result.stats["lost_reports"], 0)


class EngagementAndSeedTests(unittest.TestCase):
    def test_ecm_quota_does_not_move_the_opt_in_portfolio_stream(self):
        import random
        from types import SimpleNamespace

        from v2.common import utils
        from v2.portfolio import _advance

        for curves in (1, 9):
            config = configuration(ecm_tiers=((20, 200, curves),))
            generator = random.Random(12345)
            before = generator.getstate()
            state = dict(
                original=35,
                current=dict(n=35, mult=1, stage="ecm", tier=0, job=None),
                events=[],
                dropped_events=0,
            )
            with ECMExecutor() as executor:
                with patch.object(
                    executor,
                    "run",
                    return_value=SimpleNamespace(
                        divisor=None,
                        reason="curves_exhausted",
                        checkpoint={},
                        stats={},
                    ),
                ):
                    _advance(
                        state,
                        config,
                        allowance(),
                        None,
                        generator,
                        {},
                        None,
                        utils.PRIMALITY_POLICY,
                        executor,
                    )
            self.assertEqual(generator.getstate(), before)
            self.assertEqual(state["current"]["tier"], 1)

    def test_engagement_below_threshold_and_unfunded_start_stay_serial(self):
        from v2.execution.worker_contract import EngagementRule

        config = configuration()
        settings = WorkerSettings(workers=4)
        plan = Manifest(
            1009 * 1013, 7, "root", arithmetic_options(config, settings)
        )
        for threshold in (10000, float("nan"), 0.001):
            rule = EngagementRule(plan.schedule, threshold, threshold, (1,))
            with ECMExecutor(seed=7, settings=settings, rule=rule) as executor:
                # Four cold starts cannot borrow from a three-second CPU grant.
                run = executor.run(
                    plan.n, config=config, budget=allowance(cpu_seconds=3)
                )
                self.assertFalse(run.stats["engaged"])
                self.assertFalse(executor.slots)

    def test_cancellation_at_each_existing_logical_stage(self):
        from v2.execution import stage_jobs

        n = 1000000000039 * 1000000000061
        config = configuration(ecm_tiers=((50, 1000, 4),))
        plan = Manifest(
            n, 7, "root", arithmetic_options(config, WorkerSettings())
        )
        original = stage_jobs.advance_job
        for phase in (
            "setup",
            "stage_one",
            "stage_two_setup",
            "baby_steps",
            "giant_setup",
            "stage_two",
        ):
            stopped = [phase == "setup"]

            def advance(job, budget, context, options):
                original(job, budget, context, options)
                if job["phase"] == phase:
                    stopped[0] = True

            known_job = new_job("ecm", n, 17, 50, 1000)
            with (
                patch.object(stage_jobs, "new_job", return_value=known_job),
                patch.object(stage_jobs, "advance_job", side_effect=advance),
            ):
                run = CurveEngine().execute(
                    task(plan, 0), cancelled=lambda: stopped[0]
                )
            self.assertEqual(run["reason"], "cancelled")
            self.assertIsNone(run["factor"])
