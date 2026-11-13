"""Actual pthread arithmetic, dispatch faults and synchronized cancellation."""

import json
import os
import signal
import threading
import unittest
from dataclasses import replace
from pathlib import Path

from v2.ecm import chains, core
from v2.ecm.native_packet import SOURCE, NativePacket
from v2.ecm.native_team import NativeTeam, TeamRefusedError
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.ecm_native import NativeECMExecutor
from v2.execution.portfolio_config import (
    _execution_context,
    close_execution_context,
)
from v2.execution.stage_jobs import stage_one_plan
from v2.execution.worker_contract import WorkerSettings
from v2.tests.test_h3_resources import Clock, Provider
from v2.tests.test_native_correctness import (
    FIXTURES,
    deterministic_budget,
    fixed_curve,
)
from v2.tests.test_native_ecm import config


@unittest.skipUnless(
    (SOURCE / "build/pthread/manifest.json").exists(),
    "pthread native build required",
)
class PthreadTests(unittest.TestCase):
    def owner(self, **options):
        owner = NativeTeam(**options)
        self.addCleanup(owner.close)
        return owner

    def test_fixed_fixtures_every_boundary_all_threads_and_schedules(self):
        rows = []
        for case in json.loads(FIXTURES.read_text())["cases"]:
            control = fixed_curve(case, "python")
            for workers in (1, 2, 4):
                for schedule in ("static", "dynamic"):
                    for chunk in (1, 3):
                        with self.subTest(
                            case=case["id"],
                            workers=workers,
                            schedule=schedule,
                            chunk=chunk,
                        ):
                            owner = NativeTeam(
                                workers=workers,
                                schedule=schedule,
                                chunk_size=chunk,
                            )
                            try:
                                actual = fixed_curve(case, "c-gmp", owner)
                                self.assertEqual(actual[:3], control[:3])
                                self.assertEqual(
                                    actual[0]["factor"], case["factor"]
                                )
                                self.assertGreater(actual[3], 0)
                                self.assertGreater(owner.slices, 0)
                                self.assertEqual(owner.lib.b4_pthread(), 1)
                                self.assertEqual(owner.lib.b4_openmp(), 0)
                                rows.append(
                                    dict(
                                        case=case["id"],
                                        workers=workers,
                                        schedule=schedule,
                                        chunk=chunk,
                                        boundaries_equal=True,
                                        native_slices=owner.slices,
                                    )
                                )
                            finally:
                                owner.close()
        output = os.environ.get("FACTOR_B4A_CORRECTNESS_OUTPUT")
        if output:
            with Path(output).open("a") as stream:
                stream.write(json.dumps(dict(pthread=rows)) + "\n")

    def test_uneven_chunks_empty_queue_and_reuse_preserve_identity(self):
        for workers in (1, 2, 4, 6, 8, 10):
            for schedule in ("static", "dynamic"):
                for chunk in (1, 3):
                    owner = self.owner(
                        workers=workers, schedule=schedule, chunk_size=chunk
                    )
                    for count in (0, 1, 2, 5, 7):
                        requests, expected, seen = [], {}, {}
                        for i in range(count):
                            n = 2 ** (200 + i % 3 * 66) - 189
                            scalar, point = 2 ** (i * 7) + 13, (i + 2, 1)
                            requests.append(("ladder", scalar, point, n, 2))
                            expected[i] = core.scalar_multiply(
                                scalar, *point, n, 2
                            )

                        def publish(index, result):
                            self.assertNotIn(index, seen)
                            seen[index] = result
                            return False

                        owner.execute(
                            requests, deterministic_budget(), publish
                        )
                        self.assertEqual(seen, expected)
                    owner.close()
                    owner.close()
                    with self.assertRaises(RuntimeError):
                        owner.execute(
                            [], deterministic_budget(), lambda *_: False
                        )

    def test_partial_creation_worker_and_preparation_failures_drain(self):
        owner = self.owner(workers=4)
        grant = deterministic_budget()
        owner._grow(1, 1013, grant)
        for create_at, job_at, error in (
            (0, -1, TeamRefusedError),
            (2, -1, TeamRefusedError),
            (-1, 0, RuntimeError),
        ):
            self.assertEqual(
                owner.lib.b4_pthread_inject(owner.owner, create_at, job_at, 0),
                0,
            )
            seen = []
            with self.assertRaises(error):
                owner.execute(
                    [("ladder", 13, (2, 1), 1013, 2)],
                    grant,
                    lambda i, r: seen.append((i, r)) or False,
                )
            self.assertEqual(seen, [])
            self.assertTrue(all(not p._lock.locked() for p in owner.packets))
            self.assertEqual(
                owner.lib.b4_pthread_inject(owner.owner, -1, -1, 0), 0
            )
            owner.execute(
                [("ladder", 13, (2, 1), 1013, 2)],
                grant,
                lambda i, r: seen.append((i, r)) or False,
            )
            self.assertEqual(
                seen, [(0, core.scalar_multiply(13, 2, 1, 1013, 2))]
            )
        with self.assertRaises(ValueError):
            owner.execute(
                [
                    ("ladder", 13, (2, 1), 1013, 2),
                    ("ladder", 2**4096, (2, 1), 1013, 2),
                ],
                grant,
                lambda *_: self.fail("invalid batch published"),
            )
        self.assertTrue(all(not p._lock.locked() for p in owner.packets))

    def test_cancel_inside_foreign_dispatch_joins_before_reuse(self):
        owner = self.owner(workers=4)
        owner._grow(1, 1013, deterministic_budget())
        self.assertEqual(
            owner.lib.b4_pthread_inject(owner.owner, -1, -1, 1), 0
        )
        errors, seen = [], []

        def dispatch():
            try:
                owner.execute(
                    [("ladder", 2**4095 + 3, (2, 1), 1013, 2)],
                    deterministic_budget(),
                    lambda i, r: seen.append((i, r)) or False,
                )
            except BaseException as error:
                errors.append(error)

        thread = threading.Thread(target=dispatch)
        thread.start()
        try:
            # A C condition handshake proves all workers entered the gate.
            # No assertion depends on native arithmetic taking a chosen time.
            owner.lib.b4_pthread_wait(owner.owner, 4)
            with self.assertRaises(RuntimeError):
                owner.close()
            previous = signal.signal(signal.SIGUSR1, lambda *_: owner.cancel())
            try:
                signal.raise_signal(signal.SIGUSR1)
            finally:
                signal.signal(signal.SIGUSR1, previous)
        finally:
            owner.lib.b4_pthread_release(owner.owner)
            thread.join(10)
        self.assertFalse(thread.is_alive())
        self.assertEqual(seen, [])
        self.assertEqual(len(errors), 1)
        self.assertIsInstance(errors[0], BudgetExhaustedError)
        self.assertEqual(str(errors[0]), "cancelled")
        self.assertEqual(
            owner.lib.b4_pthread_inject(owner.owner, -1, -1, 0), 0
        )
        owner.execute(
            [("ladder", 13, (2, 1), 1013, 2)],
            deterministic_budget(),
            lambda i, r: seen.append(r) or False,
        )
        self.assertEqual(seen, [core.scalar_multiply(13, 2, 1, 1013, 2)])

    def test_interrupted_team_checkpoint_repeated_resume_keeps_history(self):
        owner = NativeECMExecutor(seed=17, settings=WorkerSettings(workers=2))
        self.addCleanup(owner.close)
        selected = config(ecm_kernel="c-gmp", ecm_tiers=((11, 47, 3),))
        checkpoint = None
        for _ in range(2):
            start = 0 if owner.team is None else owner.team.slices
            wall, cpu = Clock(), Clock()

            def cancel_after_slice():
                if owner.team is not None and owner.team.slices > start:
                    wall.value, cpu.value = 2, 3
                    return True
                return False

            grant = Budget(
                wall_clock=wall,
                cpu_clock=cpu,
                provider=Provider(),
                cancelled=cancel_after_slice,
            )
            paused = owner.run(
                1013,
                config=selected,
                budget=grant,
                checkpoint=checkpoint,
                fixed_work=True,
            )
            self.assertEqual(paused.reason, "cancelled")
            saved = paused.checkpoint["payload"]["budget"]
            if checkpoint is not None:
                before = checkpoint["payload"]["budget"]
                self.assertGreaterEqual(saved["wall"], before["wall"] + 2)
                self.assertGreaterEqual(saved["cpu"], before["cpu"] + 3)
                self.assertGreater(saved["work"], before["work"])
            grant.restore_charges(**saved)
            once = grant.charges()
            grant.restore_charges(**saved)
            self.assertEqual(grant.charges(), once)
            checkpoint = json.loads(json.dumps(paused.checkpoint))
        resumed = owner.run(
            1013,
            config=selected,
            budget=deterministic_budget(),
            checkpoint=checkpoint,
            fixed_work=True,
        )
        control = NativeECMExecutor(
            seed=17, settings=WorkerSettings(workers=2)
        )
        self.addCleanup(control.close)
        expected = control.run(
            1013,
            config=selected,
            budget=deterministic_budget(),
            fixed_work=True,
        )
        self.assertEqual(resumed.reason, "curves_exhausted")
        self.assertEqual(
            resumed.checkpoint["payload"]["completed"],
            expected.checkpoint["payload"]["completed"],
        )
        self.assertEqual(resumed.stats["charged_child_cpu"], 0)
        self.assertEqual(resumed.checkpoint["payload"]["runs"], 3)

    def test_root_stops_after_real_team_work_and_before_publication(self):
        for cause in ("wall", "cpu", "rss", "cancelled"):
            owner = self.owner()
            wall, cpu, provider = Clock(), Clock(), Provider()
            flag = [False]
            grant = Budget(
                seconds=10,
                cpu_seconds=10,
                rss_bytes=500,
                wall_clock=wall,
                cpu_clock=cpu,
                provider=provider,
                cancelled=lambda: flag[0],
            )
            owner._grow(1, 1013, grant)
            original = owner.lib

            class Library:
                def __getattr__(self, name):
                    return getattr(original, name)

                def b4_pthread_step(self, *args):
                    result = original.b4_pthread_step(*args)
                    if cause == "wall":
                        wall.value = 10
                    elif cause == "cpu":
                        cpu.value = 10
                    elif cause == "rss":
                        provider.rows = [replace(provider.rows[0], rss=500)]
                    else:
                        flag[0] = True
                    return result

            owner.lib = Library()
            try:
                with self.assertRaisesRegex(
                    BudgetExhaustedError,
                    cause if cause == "cancelled" else cause + "_limit",
                ):
                    owner.execute(
                        [("ladder", 13, (2, 1), 1013, 2)],
                        grant,
                        lambda *_: self.fail("stopped work published"),
                    )
            finally:
                owner.lib = original
            self.assertEqual(owner.slices, 1)
            self.assertTrue(all(not p._lock.locked() for p in owner.packets))

    def test_factor_publication_discards_unfinished_peer(self):
        owner, seen = self.owner(workers=2), []
        owner.execute(
            [
                ("ladder", 1, (2, 1), 1013, 2),
                ("ladder", 2**4095 + 3, (3, 1), 1013, 2),
            ],
            deterministic_budget(),
            lambda i, r: seen.append(i) or True,
        )
        self.assertEqual(seen, [0])
        self.assertEqual(owner.discarded, 1)

    def test_prac_eligibility_and_ladder_boundaries(self):
        for n in (10**39 - 1, 10**39 + 1, 10**80 - 1, 10**80 + 1):
            for bound in (1999, 2000, 2001):
                for curves in (7, 8):
                    selected = config(
                        ecm_tiers=((bound, bound + 100, curves),),
                        ecm_chain_mode="auto",
                    )
                    grant = deterministic_budget()
                    context, programs = _execution_context({}, selected, grant)
                    try:
                        plan = stage_one_plan(
                            dict(kind="ecm", n=n, b1=bound),
                            programs or context,
                            grant,
                        )
                        eligible = (
                            chains.supports_modulus(n)
                            and bound == 2000
                            and curves >= 8
                        )
                        self.assertEqual(plan is not None, eligible)
                    finally:
                        close_execution_context(programs or context)

    def test_full_width_carry_and_exceptional_operands(self):
        packet = NativePacket()
        self.addCleanup(packet.close)
        for bits in (3, 64, 65, 128, 200, 266, 333, 1024):
            n = 2**bits - 1
            for point in ((n - 1, n - 2), (0, 1), (1, 0)):
                for scalar in (0, 1, 2, 13, 2**4095 + 1):
                    expected = core.scalar_multiply(scalar, *point, n, n - 1)
                    actual = packet.ladder(
                        scalar, point, n, n - 1, deterministic_budget()
                    )
                    # Infinity is projective; the accepted scalar-zero/zero-Z
                    # boundary is canonical (1,0) in both implementations.
                    self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
