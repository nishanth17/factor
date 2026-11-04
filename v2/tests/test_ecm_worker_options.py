"""Queued F4 arithmetic/lifecycle tests; not the pure harness import graph."""

import multiprocessing
import time
import unittest
from unittest.mock import Mock, patch

from v2.execution.budget import BudgetExhaustedError
from v2.execution.ecm_worker import CurveEngine
from v2.execution.ecm_worker_state import ActionBudget, DecodedPrograms
from v2.execution.ecm_workers import ECMExecutor, arithmetic_options
from v2.execution.worker_contract import (
    EngagementRule,
    Ledger,
    Manifest,
    Search,
    WorkerSettings,
)
from v2.tests.test_ecm_workers import allowance, configuration, task


class PollingAndCacheTests(unittest.TestCase):
    def test_every_work_reservation_is_checked_with_deferred_clock_polling(
        self,
    ):
        calls = []
        budget = ActionBudget(
            poll_every=4,
            work_limit=9,
            seconds=60,
            cpu_seconds=60,
            cancelled=lambda: calls.append(1) and False,
        )
        for _ in range(9):
            budget.consume(1)
        self.assertEqual(budget.used, 9)
        self.assertEqual(len(calls), 3)
        with self.assertRaisesRegex(BudgetExhaustedError, "work_limit"):
            budget.consume(1)
        self.assertEqual(budget.used, 9)

    def test_forced_boundary_check_and_bounded_cancellation_precedence(self):
        stop = [False]
        budget = ActionBudget(
            poll_every=4,
            work_limit=100,
            seconds=60,
            cpu_seconds=60,
            cancelled=lambda: stop[0],
        )
        budget.consume(1)
        stop[0] = True
        with self.assertRaisesRegex(BudgetExhaustedError, "cancelled"):
            for _ in range(4):
                budget.consume(1)
        self.assertLessEqual(budget.used, 4)
        budget.used = 100
        with self.assertRaisesRegex(BudgetExhaustedError, "cancelled"):
            budget.consume(1)
        with self.assertRaisesRegex(BudgetExhaustedError, "cancelled"):
            budget.enter_phase("publication")

    def test_decoding_is_immutable_bounded_and_refusal_keeps_fallback(self):
        packed = Mock()
        packed.program_segment.return_value = [2, 3, 5, 7]
        packed.power_values.return_value = [8, 9, 5, 7]
        cache = DecodedPrograms(packed, 4096)
        budget = allowance()
        first = cache.program_segment(2, 11, budget, bound=10)
        second = cache.program_segment(2, 11, budget, bound=10)
        self.assertIs(first, second)
        self.assertEqual(second, (2, 3, 5, 7))
        self.assertEqual(cache.power_values(2, 11, 10, 1, 2), (9, 5))
        packed.program_segment.assert_called_once()
        with self.assertRaises(TypeError):
            second[0] = 17
        refused = DecodedPrograms(packed, 512)
        self.assertEqual(
            refused.program_segment(2, 11, budget, bound=10), [2, 3, 5, 7]
        )
        self.assertFalse(refused.blocks)
        self.assertEqual(refused.refusals, 1)
        for index in range(80):
            cache.program_segment(20 + index, 21 + index, budget)
            self.assertLessEqual(cache.used_bytes, cache.memory_bytes)
            self.assertLessEqual(len(cache.blocks), 64)

    def test_contexts_reuse_across_moduli_and_evict_before_refused_allocation(
        self,
    ):
        settings = WorkerSettings(metadata_slots=2, decoded_cache_bytes=65536)
        engine = CurveEngine(settings)
        config = configuration(ecm_program_bytes=65536)
        first = Manifest(1000003, 7, "a", arithmetic_options(config, settings))
        second = Manifest(10000019, 7, "b", first.options)
        for plan in (first, second, first):
            result = engine.execute(task(plan, 0), cancelled=lambda: False)
            self.assertEqual(result["reason"], "complete")
        self.assertEqual(engine.preparations, 1)
        self.assertGreater(engine.decode_hits, 0)
        changed = Manifest(
            first.n, 7, "a", dict(first.options, ecm_tiers=[[21, 200, 4]])
        )
        engine.execute(task(changed, 0), cancelled=lambda: False)
        self.assertEqual(len(engine.contexts), 2)
        third = Manifest(
            first.n, 7, "a", dict(first.options, ecm_tiers=[[22, 200, 4]])
        )
        with patch(
            "v2.portfolio._execution_context",
            side_effect=MemoryError("injected"),
        ):
            result = engine.execute(task(third, 0), cancelled=lambda: False)
        self.assertEqual(result["reason"], "error")
        self.assertEqual(len(engine.contexts), 1)
        self.assertNotIn(first.metadata, engine.contexts)
        self.assertIsNone(engine.identity)


class PacketAndStageTests(unittest.TestCase):
    def test_packets_reserve_each_job_and_shrink_to_funding_and_frames(self):
        settings = WorkerSettings(
            packet_curves=4,
            assignment_work=10,
            assignment_cpu=2,
            frame_bytes=2048,
        )
        config = configuration()
        plan = Manifest(
            1000003, 7, "root", arithmetic_options(config, settings)
        )
        search = Search(plan, settings, fixed_work=True)
        with ECMExecutor(settings=settings) as owner:
            ledger = Ledger(1_000_000, 5)
            packet, issued = owner._packet(
                plan,
                search,
                ledger,
                time.monotonic() + 60,
                dict(index=0, metadata_cache=[plan.metadata]),
            )
        self.assertEqual(len(issued), 3)
        self.assertEqual([job["cpu"] for job in packet["jobs"]], [2, 2, 1])
        self.assertEqual(len(ledger.active), 3)
        self.assertEqual(sum(cpu for _, _, cpu in ledger.active.values()), 5)
        for token in list(ledger.active):
            ledger.lose(token)
        self.assertEqual(ledger.cpu, 5)
        self.assertEqual(ledger.work, 3 * 10 + 4 * settings.frame_bytes)

    def test_packet_admission_keeps_ready_work_for_other_idle_slots(self):
        settings = WorkerSettings(workers=4, packet_curves=4)
        plan = Manifest(
            1000003, 7, "root", arithmetic_options(configuration(), settings)
        )
        search = Search(plan, settings, fixed_work=True)
        ledger = Ledger(200_000_000, 240)
        with ECMExecutor(settings=settings) as owner:
            packets = []
            for index in range(4):
                packet, issued = owner._packet(
                    plan,
                    search,
                    ledger,
                    time.monotonic() + 60,
                    dict(index=index, metadata_cache=[plan.metadata]),
                    peers=4 - index,
                )
                packets.append(packet)
                self.assertEqual(len(issued), 1)
        self.assertEqual(
            len({p["jobs"][0]["assignment"]["id"] for p in packets}), 4
        )
        self.assertEqual(len(search.active), 4)

    def test_boundary_continues_without_repeating_setup_or_stage_one(self):
        config, settings = (
            configuration(),
            WorkerSettings(stage_policy="staged"),
        )
        plan = Manifest(
            1000003, 7919, "root", arithmetic_options(config, settings)
        )
        engine = CurveEngine(settings)
        first = engine.execute(
            dict(task(plan, 0), stage="stage_one"), cancelled=lambda: False
        )
        self.assertEqual(first["reason"], "boundary")
        self.assertNotIn("baby", first["residue"])
        continuation = dict(
            task(plan, 0), stage="stage_two", residue=first["residue"]
        )
        with patch(
            "v2.ecm.core.setup_curve",
            side_effect=AssertionError("setup repeated"),
        ):
            result = CurveEngine(settings).execute(
                continuation, cancelled=lambda: False
            )
        whole = CurveEngine().execute(task(plan, 0), cancelled=lambda: False)
        self.assertEqual(
            (result["factor"], result["phase"], result["sigma"]),
            (whole["factor"], whole["phase"], whole["sigma"]),
        )
        self.assertEqual(result["boundary_sha256"], first["residue"]["sha256"])

    def test_drained_stage_restart_keeps_boundary_and_prior_charges(self):
        settings = WorkerSettings(stage_policy="staged", residue_capacity=1)
        config = configuration(ecm_tiers=((20, 200, 1),))
        stopped = [False]
        original = Search.receive

        def acknowledge(search, position, token, message, **kwargs):
            original(search, position, token, message, **kwargs)
            if message["reason"] == "boundary":
                stopped[0] = True

        with patch.object(
            Search, "receive", side_effect=acknowledge, autospec=True
        ):
            with ECMExecutor(seed=7, settings=settings) as owner:
                partial = owner.run(
                    1000003,
                    config=config,
                    budget=allowance(cancelled=lambda: stopped[0]),
                )
        self.assertEqual(partial.reason, "cancelled")
        self.assertEqual(
            len(partial.checkpoint["search"]["payload"]["boundaries"]), 1
        )
        with ECMExecutor(seed=7, settings=settings) as owner:
            with patch(
                "v2.ecm.core.setup_curve",
                side_effect=AssertionError("stage one repeated"),
            ):
                final = owner.run(
                    1000003,
                    config=config,
                    budget=allowance(),
                    checkpoint=partial.checkpoint,
                )
        self.assertEqual(final.reason, "curves_exhausted")
        self.assertEqual(
            final.checkpoint["search"]["payload"]["attempts"], [1]
        )
        self.assertEqual(
            final.checkpoint["search"]["payload"]["continuation_attempts"], [1]
        )
        self.assertGreater(
            final.checkpoint["budget"]["work"],
            partial.checkpoint["budget"]["work"],
        )

    def test_same_declared_curves_with_packets_stages_and_worker_counts(self):
        expected = None
        variants = (
            (1, 1, "whole", False),
            (1, 4, "staged", False),
            (1, 4, "whole", True),
            (2, 2, "staged", True),
            (4, 1, "staged", True),
        )
        before = {child.pid for child in multiprocessing.active_children()}
        for workers, batch, stages, forced in variants:
            settings = WorkerSettings(
                workers=workers,
                packet_curves=batch,
                stage_policy=stages,
                residue_capacity=4,
                poll_every=4,
                decoded_cache_bytes=65536,
            )
            with ECMExecutor(
                seed=7919, settings=settings, force_workers=forced
            ) as owner:
                run = owner.run(
                    1000003,
                    config=configuration(),
                    budget=allowance(),
                    fixed_work=True,
                )
            records = run.checkpoint["search"]["payload"]["completed"]
            actual = {
                int(i): (r["id"], r["factor"], r["sigma"], r["phase"])
                for i, r in records.items()
            }
            self.assertEqual(len(actual), 4)
            if expected is not None:
                self.assertEqual(actual, expected)
            expected = actual
            self.assertLessEqual(run.stats["residue_queue_peak"], 4)
        self.assertEqual(
            {child.pid for child in multiprocessing.active_children()}, before
        )

    def test_engagement_threshold_and_warm_state_with_bound_cost_identity(
        self,
    ):
        settings = WorkerSettings()
        config = configuration()
        plan = Manifest(
            1000003, 7, "root", arithmetic_options(config, settings)
        )
        with ECMExecutor(settings=settings) as probe:
            cost_identity = probe.cost_identity
        cold = EngagementRule(
            plan.schedule, 999, 999, (0.1,), execution=cost_identity
        )
        with ECMExecutor(seed=7, settings=settings, rule=cold) as serial:
            run = serial.run(plan.n, config=config, budget=allowance())
            self.assertFalse(run.stats["engaged"])
            self.assertFalse(serial.slots)
            self.assertEqual(run.stats["engagement"]["reason"], "threshold")
        warm = EngagementRule(
            plan.schedule, 0.001, 0.001, (0.1,), execution=cost_identity
        )
        with ECMExecutor(seed=7, settings=settings, rule=warm) as owner:
            first = owner.run(plan.n, config=config, budget=allowance())
            pids = [slot["process"].pid for slot in owner.slots]
            second = owner.run(plan.n, config=config, budget=allowance())
            self.assertTrue(first.stats["engaged"])
            self.assertEqual(first.stats["engagement"]["state"], "cold")
            self.assertEqual(second.stats["engagement"]["state"], "warm")
            self.assertEqual(second.stats["preparations"], 0)
            self.assertEqual(
                [slot["process"].pid for slot in owner.slots], pids
            )
            self.assertGreaterEqual(second.stats["reported_child_cpu"], 0)

    def test_cancellation_before_publication_is_not_hidden_by_poll_interval(
        self,
    ):
        plan = Manifest(
            35,
            7,
            "root",
            arithmetic_options(configuration(), WorkerSettings()),
        )
        stop = [False]

        def finish(job, budget, context, config):
            job.update(done=True, factor=5)
            stop[0] = True

        with patch("v2.execution.stage_jobs.advance_job", side_effect=finish):
            result = CurveEngine(WorkerSettings(poll_every=64)).execute(
                task(plan, 0), cancelled=lambda: stop[0]
            )
        self.assertEqual(result["reason"], "cancelled")
        self.assertIsNone(result["factor"])


class PoolShutdownTests(unittest.TestCase):
    def test_shutdown_deadlines_apply_to_the_pool_and_signals_are_broadcast(
        self,
    ):
        clock, events = [0.0], []

        class Process:
            def __init__(self, pid):
                self.pid, self.alive, self.killing = pid, True, False

            def is_alive(self):
                return self.alive

            def terminate(self):
                events.append(("terminate", self.pid))

            def kill(self):
                events.append(("kill", self.pid))
                self.killing = True

            def join(self, timeout):
                events.append(("join", self.pid))
                clock[0] += timeout
                if self.killing:
                    self.alive = False

            def close(self):
                events.append(("closed", self.pid))

        owner = ECMExecutor(settings=WorkerSettings(workers=4))
        owner.slots = [
            dict(process=Process(i), channel=Mock()) for i in range(4)
        ]
        with (
            patch.object(owner, "_send"),
            patch(
                "v2.execution.ecm_workers.time.monotonic",
                side_effect=lambda: clock[0],
            ),
        ):
            owner._destroy()
        self.assertEqual(clock[0], 3)
        self.assertFalse(owner.slots)
        first_kill = next(i for i, e in enumerate(events) if e[0] == "kill")
        self.assertEqual(
            [e for e in events[:first_kill] if e[0] == "terminate"],
            [("terminate", i) for i in range(4)],
        )


def packet_factor_worker(
    channel, frame_bytes, startup_work, settings_data=None
):
    """Finite fake service for real spawn/healthy reuse QA; no ECM calls."""
    from v2.execution.ecm_worker import read_message
    from v2.execution.worker_contract import Decoder, frame

    channel.settimeout(0.01)
    decoder = Decoder(frame_bytes)
    channel.sendall(
        frame(
            dict(
                kind="ready", cpu=time.process_time(), rss=1, work=startup_work
            ),
            frame_bytes,
        )
    )
    while True:
        packet = read_message(channel, decoder, block=True)
        if packet is None or packet["kind"] == "cancel":
            continue
        if packet["kind"] == "close":
            break
        for index, job in enumerate(packet["jobs"]):
            channel.sendall(
                frame(
                    dict(
                        kind="result",
                        packet=packet["packet"],
                        id=job["assignment"]["id"],
                        token=job["token"],
                        attempt=job["attempt"],
                        stage=job["stage"],
                        sigma=job["assignment"]["sigma"]
                        if index == 0
                        else None,
                        reason="complete" if index == 0 else "cancelled",
                        factor=5 if index == 0 else None,
                        work=1 if index == 0 else 0,
                        cpu=0,
                        rss=1,
                        prepared_schedule=packet["schedule"],
                        metadata_cache=[packet["schedule"]],
                    ),
                    frame_bytes,
                )
            )
    channel.close()


class RetainedPoolTests(unittest.TestCase):
    def test_first_factor_drains_packet_and_retains_healthy_workers(self):
        before = {child.pid for child in multiprocessing.active_children()}
        settings = WorkerSettings(workers=2, packet_curves=4)
        with patch(
            "v2.execution.ecm_worker.worker_main", packet_factor_worker
        ):
            with ECMExecutor(settings=settings, force_workers=True) as owner:
                first = owner.run(
                    35, config=configuration(), budget=allowance()
                )
                pids = [slot["process"].pid for slot in owner.slots]
                second = owner.run(
                    35, config=configuration(), budget=allowance()
                )
                self.assertEqual(
                    [slot["process"].pid for slot in owner.slots], pids
                )
                self.assertEqual((first.divisor, second.divisor), (5, 5))
                self.assertEqual(second.stats["startup_seconds"], 0)
                self.assertTrue(
                    all(not slot["active"] for slot in owner.slots)
                )
        self.assertEqual(
            {child.pid for child in multiprocessing.active_children()}, before
        )


class ParentHeadroomTests(unittest.TestCase):
    def test_live_leases_leave_parent_cpu_to_process_results(self):
        settings = WorkerSettings(workers=4, memory_bytes=2**31)
        owner = ECMExecutor(settings=settings)
        plan = Manifest(
            35,
            7,
            "root",
            arithmetic_options(
                configuration(ecm_tiers=((10, 100, 4),)),
                settings,
            ),
        )
        search = Search(plan, settings, fixed_work=True)
        ledger = Ledger(10**8, 5)
        for slot in range(4):
            packet, issued = owner._packet(
                plan,
                search,
                ledger,
                10**9,
                dict(index=slot, metadata_cache=[]),
                peers=4 - slot,
            )
            self.assertEqual(len(issued), 1)
            self.assertEqual(len(packet["jobs"]), 1)
        self.assertEqual(len(ledger.active), 4)
        ledger.charge(parent_cpu=0.01)
        self.assertLess(ledger.cpu + ledger.parent_cpu, ledger.cpu_limit)
