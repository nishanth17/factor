"""Whole-native production portfolio ownership and checked resume."""

import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from v2.common.arithmetic import BackendUnavailableError
from v2.ecm import native_packet
from v2.ecm.native_curve import NativeCurvePacket
from v2.ecm.native_packet import SOURCE
from v2.ecm.workers import ECMExecutor
from v2.execution.budget import Budget
from v2.execution.ecm_native import NativeECMExecutor
from v2.execution.ecm_openmp import OpenMPExecutor
from v2.execution.portfolio_config import PortfolioConfig
from v2.execution.resources import AllocationFailureError
from v2.execution.worker_contract import WorkerSettings
from v2.portfolio import factorize_bounded
from v2.tests.test_h3_resources import Provider
from v2.tests.test_native_correctness import deterministic_budget

N = 646469834424391669878723334605631644580135881618808383583597


def configuration(kernel="c-gmp-full"):
    return PortfolioConfig(
        trial_bound=2,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=((31, 233, 9),),
        ecm_chain_mode="auto",
        ecm_kernel=kernel,
        siqs=None,
        max_input_bits=1024,
    )


@unittest.skipUnless(
    all(
        (SOURCE / name).exists()
        for name in (
            "build/manifest.json",
            "build/pthread/manifest.json",
            "build/openmp/manifest.json",
        )
    ),
    "scalar/pthread/provided OpenMP builds required",
)
class NativeProductionTests(unittest.TestCase):
    def call(
        self, owner, *, n=N, checkpoint=None, budget=None, kernel="c-gmp-full"
    ):
        return factorize_bounded(
            n,
            seed=73,
            config=configuration(kernel),
            budget=budget or deterministic_budget(),
            checkpoint=checkpoint,
            ecm_executor=owner,
            qs_workers=None,
        )

    def test_benchmark_native_owners_honor_explicit_worker_settings(self):
        from v2.benchmarks.ecm.b4a.production_batch import ProductionSession

        for role, owner_type in (
            ("pthread", NativeECMExecutor),
            ("openmp", OpenMPExecutor),
        ):
            with self.subTest(role=role):
                request = dict(
                    role=role,
                    workers=6,
                    worker_settings=dict(
                        allow_wide_workers=True, parent_poll_seconds=0.02
                    ),
                )
                with patch.object(
                    owner_type, "__init__", return_value=None
                ) as init:
                    ProductionSession.make_owner(request)
                settings = init.call_args.kwargs["settings"]
                self.assertEqual(settings.workers, 6)
                self.assertTrue(settings.allow_wide_workers)
                self.assertEqual(settings.parent_poll_seconds, 0.02)

    def test_all_owners_publish_same_portfolio_result_and_retain_contexts(
        self,
    ):
        expected = None
        factories = (
            lambda: ECMExecutor(seed=73, settings=WorkerSettings(workers=1)),
            lambda: NativeECMExecutor(
                seed=73, settings=WorkerSettings(workers=4), capacity=32
            ),
            lambda: OpenMPExecutor(
                seed=73, settings=WorkerSettings(workers=4), capacity=32
            ),
        )
        for factory in factories:
            with factory() as owner:
                first = self.call(owner)
                self.assertEqual(first.result.reconstruct(), N)
                self.assertFalse(first.result.complete)
                if expected is None:
                    expected = first.result
                self.assertEqual(first.result, expected)
                team = getattr(owner, "team", None)
                if team is not None:
                    self.assertEqual(team.actual[0], 4)
                    if isinstance(owner, OpenMPExecutor):
                        self.assertEqual(
                            team.lib.b4_omp_participants(team.owner), 15
                        )
                second = self.call(owner)
                self.assertEqual(second.result, first.result)
                if team is not None:
                    self.assertIs(owner.team, team)
            owner.close()
            if team is not None:
                self.assertEqual(team.packets, [])
                self.assertEqual(team.blocks.used, 0)

    def test_portfolio_publishes_a_proper_native_factor(self):
        for kind in (NativeECMExecutor, OpenMPExecutor):
            with kind(seed=73, settings=WorkerSettings(workers=4)) as owner:
                run = self.call(owner, n=1009 * N)
                self.assertEqual(run.result.reconstruct(), 1009 * N)
                self.assertIn(1009, [int(p.value) for p in run.result.factors])
                self.assertTrue(
                    any(
                        event["stage"] == "ecm_workers"
                        and event["outcome"] == "factor"
                        for event in run.events
                    )
                )

    def test_private_contexts_grow_on_demand_and_survive_quota_changes(self):
        for kind in (NativeECMExecutor, OpenMPExecutor):
            with kind(
                seed=73, settings=WorkerSettings(workers=4), capacity=32
            ) as owner:
                for curves, expected in (
                    (3, 3),
                    (9, 9),
                    (3, 3),
                    (37, 32),
                    (9, 9),
                ):
                    config = replace(
                        configuration(), ecm_tiers=((31, 233, curves),)
                    )
                    run = factorize_bounded(
                        N,
                        seed=73,
                        config=config,
                        budget=deterministic_budget(),
                        ecm_executor=owner,
                        qs_workers=None,
                    )
                    self.assertEqual(run.result.reconstruct(), N)
                    self.assertFalse(run.result.complete)
                    packets = owner.team.packets
                    self.assertEqual(len(packets), expected)
                    self.assertEqual(
                        len(
                            {
                                int(
                                    owner.team.ffi.cast("uintptr_t", p.context)
                                )
                                for p in packets
                            }
                        ),
                        expected,
                    )
                    retained = list(packets)
                    # A changed numerical manifest deliberately replaces the
                    # team. Repeating that manifest must reuse its contexts.
                    repeated = factorize_bounded(
                        N,
                        seed=73,
                        config=config,
                        budget=deterministic_budget(),
                        ecm_executor=owner,
                        qs_workers=None,
                    )
                    self.assertEqual(repeated.result, run.result)
                    self.assertEqual(owner.team.packets, retained)
                team = owner.team
            self.assertEqual(team.packets, [])
            self.assertEqual(team.blocks.used, 0)

    def test_portfolio_interruption_reopens_same_backend_and_rejects_another(
        self,
    ):
        for kind in (NativeECMExecutor, OpenMPExecutor):
            owner = kind(
                seed=73, settings=WorkerSettings(workers=4), capacity=32
            )
            try:
                grant = Budget(
                    provider=Provider(),
                    cancelled=lambda: owner.team is not None
                    and owner.team.slices > 0,
                )
                paused = self.call(owner, budget=grant)
                self.assertEqual(paused.reason, "cancelled")
                checkpoint = json.loads(json.dumps(paused.checkpoint))
                self.assertGreater(checkpoint["payload"]["work_used"], 0)
            finally:
                owner.close()
            with kind(
                seed=73, settings=WorkerSettings(workers=4), capacity=32
            ) as resumed_owner:
                resumed = self.call(resumed_owner, checkpoint=checkpoint)
                self.assertEqual(resumed.result.reconstruct(), N)
                self.assertNotEqual(resumed.reason, "cancelled")
            other_kind = (
                OpenMPExecutor
                if kind is NativeECMExecutor
                else NativeECMExecutor
            )
            with other_kind(
                seed=73, settings=WorkerSettings(workers=4), capacity=32
            ) as other:
                with self.assertRaises(ValueError):
                    self.call(other, checkpoint=checkpoint)

    def test_extended_openmp_counts_preserve_results_and_participation(self):
        for count in (8, 10, 12, 14, 16):
            with self.subTest(workers=count):
                settings = WorkerSettings(
                    workers=count, allow_wide_workers=True
                )
                with OpenMPExecutor(seed=73, settings=settings) as owner:
                    run = self.call(owner)
                    self.assertEqual(run.result.reconstruct(), N)
                    self.assertEqual(owner.team.actual[0], count)
                    activity = owner.team.worker_activity()
                    self.assertEqual(
                        activity["entered_mask"], (1 << count) - 1
                    )
                    self.assertEqual(len(activity["workers"]), count)

    def test_extended_teams_execute_and_reuse_all_compute_ranks(self):
        from v2.benchmarks.ecm.b4a.production_batch import ProductionSession

        n = 440054034826859056921224347106145832973598799449799064289819
        for role in ("pthread", "openmp"):
            for count in (8, 10, 12, 14, 16):
                with self.subTest(role=role, workers=count):
                    session = ProductionSession()
                    request = dict(
                        role=role,
                        workers=count,
                        capacity=32,
                        # Static ownership guarantees a task for each rank in
                        # this deliberately short fixture. Dynamic scheduling
                        # may legally finish before a late rank takes a task.
                        schedule="static",
                        worker_settings=dict(allow_wide_workers=True),
                        n=n,
                        seed=73,
                        b1=127,
                        b2=4000,
                        curves=37,
                    )
                    try:
                        identities = None
                        for repeat in range(2):
                            result = session.call(request)
                            self.assertIsNone(result["resource_stop"])
                            self.assertEqual(result["completed_curves"], 37)
                            self.assertEqual(result["remaining"], [n])
                            activity = session.owner.team.worker_activity()
                            self.assertEqual(activity["actual"], count)
                            self.assertEqual(
                                activity["entered_mask"], (1 << count) - 1
                            )
                            ranks = activity["workers"]
                            current = [row["identity"] for row in ranks]
                            self.assertEqual(len(set(current)), count)
                            if identities is not None:
                                self.assertEqual(current, identities)
                            identities = current
                            self.assertTrue(
                                all(row["tasks"] > 0 for row in ranks)
                            )
                            self.assertEqual(
                                sum(row["curves"] for row in ranks),
                                37 * (repeat + 1),
                            )
                    finally:
                        session.close()

    def test_openmp_explicit_start_uses_the_selected_runtime(self):
        with OpenMPExecutor(
            seed=73, settings=WorkerSettings(workers=4)
        ) as owner:
            owner.start(budget=deterministic_budget())
            self.assertEqual(owner.team.variant, "openmp")
            self.assertEqual(owner.team.actual[0], 4)
            self.assertEqual(
                owner.team.lib.b4_omp_participants(owner.team.owner), 15
            )
            result = self.call(owner)
            self.assertEqual(result.result.reconstruct(), N)
            self.assertEqual(owner.team.variant, "openmp")

    def test_openmp_partial_initialization_is_reachable_by_close(self):
        with OpenMPExecutor(
            seed=73, settings=WorkerSettings(workers=4)
        ) as owner:
            original = NativeCurvePacket.available
            calls = [0]

            def fail(packet, *args):
                calls[0] += 1
                if calls[0] == 2:
                    raise MemoryError("production partial-init injection")
                return original(packet, *args)

            with patch.object(NativeCurvePacket, "available", fail):
                with self.assertRaises(AllocationFailureError) as raised:
                    self.call(owner)
                self.assertIsInstance(raised.exception.__cause__, MemoryError)
            team = owner.team
            self.assertGreater(len(team.packets), 0)
        self.assertEqual(team.packets, [])
        self.assertEqual(team.blocks.used, 0)

    def test_missing_openmp_is_lazy_for_cheap_calls_and_explicit_for_ecm(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(native_packet, "SOURCE", Path(directory)):
                with OpenMPExecutor(seed=73) as owner:
                    cheap = self.call(owner, n=9)
                    self.assertTrue(cheap.result.complete)
                    self.assertIsNone(owner.team)
                    with self.assertRaises(BackendUnavailableError):
                        self.call(owner)
