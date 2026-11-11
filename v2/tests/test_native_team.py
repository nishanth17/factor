"""Exact pthread queues, publication and cumulative resume checks."""

import copy
import os
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import RecordingBudget
from v2.ecm import core
from v2.ecm.native_packet import SOURCE
from v2.ecm.native_team import TEAM_OVERHEAD, NativeTeam
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.ecm_native import NativeECMExecutor
from v2.execution.ecm_workers import ECMExecutor
from v2.execution.portfolio_config import PortfolioConfig
from v2.execution.resources import RemovedResourceLimitError
from v2.execution.worker_contract import WorkerSettings

ENV = dict(
    OMP_WAIT_POLICY="PASSIVE",
    KMP_BLOCKTIME="0",
    OMP_DYNAMIC="FALSE",
    OMP_MAX_ACTIVE_LEVELS="1",
    OMP_THREAD_LIMIT="4",
    KMP_STACKSIZE="4M",
    LIBOMP_NUM_HIDDEN_HELPER_THREADS="0",
)


def budget(**kwargs):
    return RecordingBudget(
        **dict(dict(work_limit=10**8, seconds=8, cpu_seconds=16), **kwargs)
    )


def team(workers=2, schedule="dynamic", capacity=8):
    return NativeTeam(
        workers=workers,
        schedule=schedule,
        capacity=capacity,
    )


def settings(workers=2):
    return WorkerSettings(workers=workers)


def config(kernel="c-gmp"):
    return PortfolioConfig(
        ecm_tiers=((11, 47, 3),),
        ecm_kernel=kernel,
        siqs=None,
        segment_size=16,
        max_input_bits=512,
    )


class NativeTeamAdmissionTests(unittest.TestCase):
    def test_retired_memory_and_invalid_configuration_refuse(self):
        with self.assertRaises(RemovedResourceLimitError):
            NativeTeam(memory_bytes=TEAM_OVERHEAD)
        for options in (
            dict(workers=0),
            dict(workers=17),
            dict(workers=True),
            dict(schedule="auto"),
            dict(capacity=33),
        ):
            with self.assertRaises(ValueError):
                NativeTeam(**options)
        with self.assertRaises(RemovedResourceLimitError):
            NativeECMExecutor(
                settings=replace(settings(), memory_bytes=100 * 2**20),
                capacity=8,
            )

    def test_lazy_owner_preserves_cheap_portfolio_exits(self):
        from v2.portfolio import factorize_bounded

        owner = NativeECMExecutor(settings=settings())
        try:
            with patch.object(
                NativeTeam,
                "execute",
                side_effect=AssertionError("eager native team"),
            ):
                run = factorize_bounded(
                    9, config=config(), budget=budget(), ecm_executor=owner
                )
            self.assertTrue(run.result.complete)
            self.assertEqual(run.result.reconstruct(), 9)
            self.assertIsNone(owner.team)
        finally:
            owner.close()


class NativeObservedResourceTests(unittest.TestCase):
    def test_unlimited_root_and_saved_history_without_loading_native(self):
        selected = replace(config(), ecm_tiers=((11, 47, 0),))
        owner = NativeECMExecutor(settings=settings())
        self.addCleanup(owner.close)
        first_budget = Budget(rss_bytes=None)
        first = owner.run(1013, config=selected, budget=first_budget)
        self.assertEqual(first.reason, "curves_exhausted")
        self.assertIsNone(owner.team)
        restored = Budget(rss_bytes=None)
        second = owner.run(
            1013, config=selected, budget=restored, checkpoint=first.checkpoint
        )
        self.assertEqual(second.reason, "curves_exhausted")
        saved = first.checkpoint["payload"]["budget"]
        self.assertGreaterEqual(restored.used, saved["work"])
        self.assertGreaterEqual(restored.prior_wall, saved["wall"])
        self.assertGreaterEqual(restored.prior_cpu, saved["cpu"])
        self.assertEqual(second.checkpoint["payload"]["runs"], 2)

    def test_root_stop_is_structured_and_precedes_native_loading(self):
        owner = NativeECMExecutor(settings=settings())
        self.addCleanup(owner.close)
        run = owner.run(1013, config=config(), budget=Budget(seconds=0))
        self.assertEqual(run.reason, "wall_limit")
        self.assertEqual(run.resource_stop["limit"], "wall")
        self.assertIsNone(owner.team)

    def test_historical_native_checkpoint_requires_new_resource_identity(self):
        selected = replace(config(), ecm_tiers=((11, 47, 0),))
        owner = NativeECMExecutor(settings=settings())
        self.addCleanup(owner.close)
        run = owner.run(1013, config=selected, budget=Budget(rss_bytes=None))
        old = copy.deepcopy(run.checkpoint)
        old.pop("resource_version")
        with self.assertRaises(ValueError):
            owner.run(1013, config=selected, budget=Budget(), checkpoint=old)


@unittest.skipUnless(
    (SOURCE / "build/pthread/manifest.json").exists(),
    "optional pthread build absent",
)
class CompiledTeamTests(unittest.TestCase):
    def setUp(self):
        self.environment = patch.dict(os.environ, ENV)
        self.environment.start()

    def tearDown(self):
        self.environment.stop()

    def test_unequal_private_ladders_all_schedules_and_team_sizes(self):
        # Three full-modulus widths; scalar lengths vary independently. This
        # verifies scheduling semantics, not a timing or worker-count study.
        exponents = (1, 3, 13, 257, 2**65 + 3, 2**127 + 7)
        for workers in (1, 2, 4):
            for schedule in ("dynamic", "static"):
                owner = team(workers, schedule)
                requests, expected, seen = [], {}, {}
                for index, scalar in enumerate(exponents):
                    n = (2 ** (200 + 66 * (index % 3))) - 189
                    point, a24 = (index + 2, 1), index + 2
                    requests.append(("ladder", scalar, point, n, a24))
                    expected[index] = core.scalar_multiply(
                        scalar, *point, n, a24
                    )
                try:

                    def publish(index, result):
                        self.assertNotIn(index, seen)
                        seen[index] = result
                        return False

                    owner.execute(requests, budget(), publish)
                    self.assertEqual(seen, expected)
                    self.assertEqual(owner.completed, len(requests))
                    self.assertTrue(
                        all(
                            p.lib.b4_step(p.context, 5, 128) == -1
                            for p in owner.packets
                        )
                    )
                finally:
                    owner.close()
                    owner.close()

    def test_simultaneous_success_stops_after_one_checked_publication(self):
        owner, seen = team(), []
        try:
            owner.execute(
                [("ladder", 1, (i + 2, 1), 1013, 2) for i in range(8)],
                budget(),
                lambda index, result: seen.append(index) or True,
            )
            self.assertEqual(seen, [0])
            self.assertEqual((owner.completed, owner.discarded), (1, 7))
        finally:
            owner.close()

    def test_cancelled_private_chunks_abort_and_can_be_reused(self):
        owner, calls = team(), 0

        def cancelled():
            nonlocal calls
            calls += 1
            # First setup has bounded transport/program charges too. Use a
            # counter that becomes armed only once native slices actually ran.
            return owner.slices >= 2

        try:
            with self.assertRaises(BudgetExhaustedError):
                owner.execute(
                    [("ladder", 2**4095 + 7, (2, 1), 2**1024 - 109, 2)],
                    budget(cancelled=cancelled),
                    lambda *_: False,
                )
            self.assertGreater(calls, 0)
            self.assertTrue(
                all(
                    p.lib.b4_step(p.context, 5, 128) == -1
                    for p in owner.packets
                )
            )
            seen = []
            owner.execute(
                [("ladder", 13, (2, 1), 1013, 2)],
                budget(),
                lambda index, result: seen.append(result) or False,
            )
            self.assertEqual(seen, [core.scalar_multiply(13, 2, 1, 1013, 2)])
        finally:
            owner.close()

    def test_duplicate_context_preflight_does_not_start_arithmetic(self):
        owner = team()
        try:
            owner._grow(1, 1013, budget())
            packet = owner.packets[0]
            packet.prepare_ladder(13, (2, 1), 1013, 2, budget())
            try:
                owner.contexts[0] = owner.contexts[1] = packet.context
                code = owner.lib.b4_pthread_step(
                    owner.owner,
                    owner.contexts,
                    2,
                    2,
                    0,
                    1,
                    5,
                    128,
                    owner.statuses,
                    owner.actual,
                )
                self.assertEqual(code, -1)
                self.assertEqual(owner.actual[0], 0)
                self.assertEqual(
                    packet.lib.b4_result(packet.context, packet.output), -1
                )
            finally:
                packet.end()
        finally:
            owner.close()

    def test_shared_verified_prac_records_match_scalar_frontier_recovery(self):
        from v2.ecm.chains import ChainPlan

        plan = ChainPlan(2000, "python-int", budget())
        n, a24 = 2**199 - 19, 2
        powers = [(p, row[0]) for p, row in list(plan.entries.items())[:16]]
        expected = [
            plan.execute(powers, (x, 1), n, a24, budget()) for x in (2, 3)
        ]
        entries = plan.prepare(powers, budget())
        owner, seen = team(), {}
        try:

            def publish(index, result):
                seen[index] = plan.recover(
                    entries, (index + 2, 1), *result, n, a24
                )
                return False

            owner.execute(
                [("records", entries, (x, 1), n, a24) for x in (2, 3)],
                budget(),
                publish,
            )
            self.assertEqual([seen[0], seen[1]], expected)
            self.assertIs(owner.packets[0].programs, owner.packets[1].programs)
            self.assertEqual(len(owner.programs.programs), len(entries))
        finally:
            owner.close()

    def test_partial_context_initialization_retains_ownership_and_unwinds(
        self,
    ):
        from v2.ecm.native_packet import NativePacket

        owner, calls = team(), 0
        allocate = NativePacket._allocate

        def fail_second(packet, *args):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise MemoryError("injected context initialization failure")
            return allocate(packet, *args)

        grant = budget()
        try:
            with patch.object(NativePacket, "_allocate", fail_second):
                with self.assertRaises(MemoryError):
                    owner.execute(
                        [("ladder", 13, (2, 1), 1013, 2)] * 2,
                        grant,
                        lambda *_: False,
                    )
            self.assertGreater(grant.used, 0)
            self.assertEqual(owner.slices, 0)
            self.assertEqual(len(owner.packets), 1)
            retained = owner.packets[0]
        finally:
            owner.close()
        self.assertIsNone(retained.context)

    def test_native_owner_uses_same_b12_assignments_and_exact_curve_outcomes(
        self,
    ):
        native = NativeECMExecutor(seed=19, settings=settings())
        control = ECMExecutor(seed=19, settings=settings(1))
        try:
            observed = native.run(
                143, config=config(), budget=budget(), fixed_work=True
            )
            expected = control.run(
                143, config=config("python"), budget=budget(), fixed_work=True
            )
            first = observed.checkpoint["payload"]
            second = expected.checkpoint["search"]["payload"]
            self.assertEqual(first["manifest"], second["manifest"])
            self.assertEqual(observed.reason, "curves_exhausted")
            self.assertEqual(set(first["completed"]), set(second["completed"]))
            for key, result in first["completed"].items():
                for field in ("id", "sigma", "factor"):
                    self.assertEqual(
                        result[field], second["completed"][key][field]
                    )
            self.assertEqual(observed.stats["charged_child_cpu"], 0)
        finally:
            native.close()
            control.close()

    def test_resume_keeps_charges_and_rejects_schedule_or_count_changes(self):
        owner = NativeECMExecutor(seed=17, settings=settings())
        cancelled = budget(cancelled=lambda: True)
        try:
            paused = owner.run(
                1013, config=config(), budget=cancelled, fixed_work=True
            )
            self.assertEqual(paused.reason, "cancelled")
            self.assertEqual(paused.checkpoint["payload"]["completed"], {})
            saved = copy.deepcopy(paused.checkpoint)
            resumed = owner.run(
                1013,
                config=config(),
                budget=budget(),
                checkpoint=saved,
                fixed_work=True,
            )
            self.assertEqual(resumed.reason, "curves_exhausted")
            self.assertGreaterEqual(
                resumed.checkpoint["payload"]["budget"]["work"],
                saved["payload"]["budget"]["work"],
            )
            for options in (
                dict(schedule="static"),
                dict(settings=settings(4)),
            ):
                other = NativeECMExecutor(
                    seed=17, **dict(dict(settings=settings()), **options)
                )
                try:
                    with self.assertRaises(ValueError):
                        other.run(
                            1013,
                            config=config(),
                            budget=budget(),
                            checkpoint=saved,
                            fixed_work=True,
                        )
                finally:
                    other.close()
        finally:
            owner.close()


if __name__ == "__main__":
    unittest.main()
