"""Deterministic generation, startup, cancellation and final join checks."""

import threading
import unittest

from v2.ecm.native_packet import SOURCE
from v2.ecm.native_team import NativeTeam, TeamRefusedError
from v2.execution.budget import BudgetExhaustedError
from v2.tests.test_native_correctness import deterministic_budget


@unittest.skipUnless(
    (SOURCE / "build/pthread/manifest.json").exists(),
    "pthread native build required",
)
class PersistentTeamTests(unittest.TestCase):
    def counters(self, owner):
        return tuple(
            int(owner.lib.b4_pthread_counter(owner.owner, kind))
            for kind in range(4)
        )

    def test_generations_reuse_handles_and_close_joins_exactly_once(self):
        for workers in range(1, 11):
            for schedule in ("static", "dynamic"):
                owner = NativeTeam(workers=workers, schedule=schedule)
                self.addCleanup(owner.close)
                for count in (1, 5, 2):
                    seen = []
                    owner.execute(
                        [
                            ("ladder", 2**511 + 3, (i + 2, 1), 1013, 2)
                            for i in range(count)
                        ],
                        deterministic_budget(),
                        lambda i, result: seen.append(i) or False,
                    )
                    self.assertEqual(seen, list(range(count)))
                    created, joined, generations, wakes = self.counters(owner)
                    self.assertEqual((created, joined), (workers, 0))
                    self.assertEqual(generations, owner.slices)
                    self.assertEqual(wakes, workers * generations)
                owner.close()
                self.assertEqual(owner.lifecycle[:2], (workers, workers))
                saved = owner.lifecycle
                owner.close()
                self.assertEqual(owner.lifecycle, saved)

    def test_partial_startup_never_publishes_and_joins_created_subset(self):
        owner = NativeTeam(workers=10)
        self.addCleanup(owner.close)
        owner._grow(1, 1013, deterministic_budget())
        owner.lib.b4_pthread_inject(owner.owner, 9, -1, 0)
        with self.assertRaises(TeamRefusedError):
            owner.execute(
                [("ladder", 1, (2, 1), 1013, 2)],
                deterministic_budget(),
                lambda *_: self.fail("partial startup published"),
            )
        self.assertEqual(self.counters(owner), (9, 9, 0, 0))
        owner.lib.b4_pthread_inject(owner.owner, -1, -1, 0)
        owner.execute(
            [("ladder", 1, (2, 1), 1013, 2)],
            deterministic_budget(),
            lambda *_: False,
        )
        self.assertEqual(self.counters(owner), (19, 9, 1, 10))
        owner.close()
        self.assertEqual(owner.lifecycle, (19, 19, 1, 10))

    def test_cancelled_generation_drains_then_reuses_the_same_workers(self):
        owner = NativeTeam(workers=10)
        self.addCleanup(owner.close)
        owner._grow(1, 1013, deterministic_budget())
        owner.lib.b4_pthread_inject(owner.owner, -1, -1, 1)
        errors = []

        def dispatch():
            try:
                owner.execute(
                    [("ladder", 255, (2, 1), 1013, 2)],
                    deterministic_budget(),
                    lambda *_: self.fail("cancelled generation published"),
                )
            except BaseException as error:
                errors.append(error)

        caller = threading.Thread(target=dispatch)
        caller.start()
        try:
            owner.lib.b4_pthread_wait(owner.owner, 10)
            self.assertEqual(owner.lib.b4_pthread_shutdown(owner.owner), -1)
            owner.cancel()
        finally:
            owner.lib.b4_pthread_release(owner.owner)
            caller.join(10)
        self.assertFalse(caller.is_alive())
        self.assertEqual(len(errors), 1)
        self.assertIsInstance(errors[0], BudgetExhaustedError)
        self.assertEqual(self.counters(owner), (10, 0, 1, 10))
        owner.lib.b4_pthread_inject(owner.owner, -1, -1, 0)
        owner.execute(
            [("ladder", 1, (2, 1), 1013, 2)],
            deterministic_budget(),
            lambda *_: False,
        )
        self.assertEqual(self.counters(owner), (10, 0, 2, 20))
        owner.close()
        self.assertEqual(owner.lifecycle, (10, 10, 2, 20))
