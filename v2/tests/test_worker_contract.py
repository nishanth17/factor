"""Pure scheduler/failure lane, without arithmetic or process imports."""

import copy
import itertools
import json
import struct
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.execution.worker_contract import (
    Decoder,
    EngagementRule,
    LeaseRefusedError,
    Ledger,
    Manifest,
    Search,
    WorkerSettings,
    canonical,
    digest,
    frame,
)


def manifest(tiers=((10, 100, 4), (20, 200, 2)), n=35):
    return Manifest(
        n,
        7919,
        "original/cofactor",
        dict(ecm_tiers=tiers, backend="python-int"),
    )


def message(search, position, token, *, factor=None, reason="complete"):
    return dict(
        id=search.manifest.assignment(position)["id"],
        token=token,
        attempt=search.attempts[position],
        factor=factor,
        reason=reason,
    )


class IdentityTests(unittest.TestCase):
    def test_workers_and_completion_order_do_not_change_declared_set(self):
        expected = {(0, i) for i in range(4)} | {(1, i) for i in range(2)}
        ids = None
        for workers in (1, 2, 4):
            search = Search(
                manifest(), WorkerSettings(workers=workers), fixed_work=True
            )
            seen, token = set(), 0
            while len(search.completed) < 6:
                issued = []
                for slot in range(workers):
                    position = search.next(slot)
                    if position is not None:
                        token += 1
                        search.dispatch(position, token)
                        issued.append((position, token))
                self.assertTrue(issued)
                # Reverse completion to disagree with assignment order.
                for position, ticket in reversed(issued):
                    seen.add(search.manifest.jobs[position])
                    search.receive(
                        position, ticket, message(search, position, ticket)
                    )
            actual = [search.manifest.assignment(i) for i in range(6)]
            self.assertEqual(seen, expected)
            if ids is not None:
                self.assertEqual(actual, ids)
            ids = actual

    def test_identity_binds_lineage_bounds_backend_seed_and_tier(self):
        first = manifest()
        for changed in (
            Manifest(35, 7920, first.lineage, first.options),
            Manifest(35, 7919, "different", first.options),
            manifest(((10, 101, 4), (20, 200, 2))),
            manifest(n=77),
        ):
            self.assertNotEqual(first.identity, changed.identity)
        with self.assertRaises(ValueError):
            Manifest(
                35,
                0,
                "root",
                dict(ecm_tiers=((10, 100, 1),), backend="gmpy2-mpz"),
            )

    def test_sigma_collisions_do_not_skip_or_redraw_assignments(self):
        with patch("v2.execution.worker_contract.random.Random") as rng:
            rng.return_value.randrange.return_value = 6
            first, second = manifest().assignment(0), manifest().assignment(1)
        self.assertEqual(first["sigma"], second["sigma"])
        self.assertNotEqual(first["id"], second["id"])
        self.assertNotEqual(first["seed"], second["seed"])

    def test_manifest_is_an_immutable_copy(self):
        options = dict(ecm_tiers=[[10, 100, 2]], backend="python-int")
        plan = Manifest(35, 0, "root", options)
        options["ecm_tiers"][0][2] = 7
        snapshot = plan.options
        snapshot["ecm_tiers"][0][2] = 9
        self.assertEqual(plan.tiers, ((10, 100, 2),))
        self.assertEqual(plan.options["ecm_tiers"][0][2], 2)


class LifecycleTests(unittest.TestCase):
    def test_first_observed_factor_and_late_result_ownership(self):
        for order in itertools.permutations((0, 1)):
            search = Search(manifest(), WorkerSettings(workers=2))
            for position in (0, 1):
                search.dispatch(position, position + 1)
            for position in order:
                search.receive(
                    position,
                    position + 1,
                    message(
                        search, position, position + 1, factor=(5, 7)[position]
                    ),
                )
            self.assertEqual(search.factor, (5, 7)[order[0]])
            self.assertIsNone(search.next(0))
            self.assertEqual(set(search.completed), {0, 1})

    def test_cancelled_publication_keeps_result_without_publishing(self):
        search = Search(manifest(), WorkerSettings())
        search.dispatch(0, 1)
        with self.assertRaises(ValueError):
            search.snapshot()
        search.receive(0, 1, message(search, 0, 1, factor=5), publish=False)
        self.assertIsNone(search.factor)
        self.assertIn(0, search.completed)
        self.assertFalse(search.retired)

    def test_fault_before_ack_retries_same_curve_with_a_new_attempt(self):
        search = Search(manifest(), WorkerSettings())
        search.dispatch(0, 1)
        old = message(search, 0, 1)
        # Lost result: keep the attempt charge, discard only the live handle.
        del search.active[0]
        saved = json.loads(canonical(search.snapshot()))
        resumed = Search.restore(saved, search.manifest, search.settings)
        self.assertEqual(resumed.next(0), 0)
        resumed.dispatch(0, 2)
        with self.assertRaises(ValueError):
            resumed.receive(0, 2, old)
        resumed.receive(0, 2, message(resumed, 0, 2))
        after = Search.restore(
            resumed.snapshot(), resumed.manifest, resumed.settings
        )
        self.assertEqual(after.next(0), 1)
        self.assertEqual(after.attempts[0], 2)

    def test_duplicate_corrupt_stale_and_nonproper_results_do_not_ack(self):
        for factor in (0, 1, 35, 2, -5, True, 5.0):
            search = Search(manifest(), WorkerSettings())
            search.dispatch(0, 1)
            with self.assertRaises(ValueError):
                search.receive(0, 1, message(search, 0, 1, factor=factor))
            self.assertIn(0, search.active)
        search = Search(manifest(), WorkerSettings())
        search.dispatch(0, 1)
        valid = message(search, 0, 1)
        wrong = dict(valid, id="retired-cofactor")
        with self.assertRaises(ValueError):
            search.receive(0, 1, wrong)
        search.receive(0, 1, valid)
        with self.assertRaises(ValueError):
            search.receive(0, 1, valid)

    def test_checkpoint_corruption_count_change_and_attempt_ceiling(self):
        search = Search(manifest(), WorkerSettings(max_attempts=1))
        search.dispatch(0, 1)
        del search.active[0]
        saved = search.snapshot()
        with self.assertRaises(LeaseRefusedError):
            search.next(0)
        with self.assertRaises(ValueError):
            Search.restore(
                saved, search.manifest, replace(search.settings, workers=2)
            )
        corrupt = copy.deepcopy(saved)
        corrupt["payload"]["attempts"][0] = 0
        with self.assertRaises(ValueError):
            Search.restore(corrupt, search.manifest, search.settings)
        corrupt["sha256"] = digest(corrupt["payload"])
        corrupt["payload"]["factor"] = 5
        corrupt["sha256"] = digest(corrupt["payload"])
        with self.assertRaises(ValueError):
            Search.restore(corrupt, search.manifest, search.settings)


class ResourceTests(unittest.TestCase):
    def test_missing_report_is_not_free_and_refund_is_exactly_once(self):
        ledger = Ledger(100, 10)
        first = ledger.reserve("curve-a", 40, 4)
        second = ledger.reserve("curve-b", 40, 4)
        ledger.lose(first)
        ledger.settle(second, "curve-b", 13, 1.5)
        self.assertEqual(ledger.work, 53)
        self.assertEqual(ledger.cpu, 5.5)
        with self.assertRaises(KeyError):
            ledger.settle(second, "curve-b", 0, 0)
        ledger.charge(parent_cpu=4)
        with self.assertRaises(LeaseRefusedError):
            ledger.reserve("unfunded", 1, 1)
        self.assertEqual(ledger.work, 53)

    def test_bad_report_never_refunds_and_partial_grant_can_be_used(self):
        ledger = Ledger(7, 1)
        ticket = ledger.reserve("last", 7, 1)
        for identity, work, cpu in (
            ("wrong", 0, 0),
            ("last", 8, 0),
            ("last", 0, 2),
            ("last", -1, 0),
            ("last", 0, float("nan")),
        ):
            with self.assertRaises(ValueError):
                ledger.settle(ticket, identity, work, cpu)
            self.assertEqual((ledger.work, ledger.cpu), (7, 1))
        ledger.settle(ticket, "last", 7, 1)
        with self.assertRaises(LeaseRefusedError):
            Ledger(0, 0).reserve("zero", 1, 0)
        with self.assertRaises(MemoryError):
            WorkerSettings(memory_bytes=2**20)
        with self.assertRaises(MemoryError):
            Search(manifest(), WorkerSettings(checkpoint_bytes=4096))

    def test_frame_fragmentation_coalescing_corruption_and_oversize(self):
        raw = frame(dict(kind="fake", id=7), 1024)
        for boundary in range(len(raw)):
            decoder = Decoder(1024)
            self.assertIsNone(decoder.feed(raw[:boundary]))
            self.assertEqual(
                decoder.feed(raw[boundary:]), dict(kind="fake", id=7)
            )
        decoder = Decoder(1024)
        self.assertEqual(decoder.feed(raw + raw), dict(kind="fake", id=7))
        self.assertEqual(decoder.feed(b""), dict(kind="fake", id=7))
        for raw_bad in (
            struct.pack("!I", 1025),
            b"\x00\x00\x00\x00",
            struct.pack("!I", 1) + b"x",
        ):
            with self.assertRaises(ValueError):
                Decoder(1024).feed(raw_bad)
        with self.assertRaises(ValueError):
            frame(dict(payload="x" * 2000), 1024)

    def test_engagement_uses_only_declared_plan_and_pool_readiness(self):
        plan = manifest()
        rule = EngagementRule(plan.schedule, 5, 1, (0.5, 0.5))
        self.assertFalse(rule.engages(plan, ready=False))
        self.assertTrue(rule.engages(plan, ready=True))
        self.assertFalse(rule.engages(plan, ready=True, remaining=(0, 1)))
        self.assertFalse(EngagementRule().engages(plan, ready=True))
        for value in (None, -1, float("nan"), float("inf"), True):
            invalid = replace(rule, seconds_per_curve=(value, 1))
            self.assertFalse(invalid.engages(plan, ready=True))
        self.assertFalse(
            replace(rule, schedule="stale-source").engages(plan, ready=True)
        )
        # No factor labels, curve outcomes or worker/core type are arguments.
        self.assertEqual(
            rule.engages(plan, ready=False), rule.engages(plan, ready=False)
        )


if __name__ == "__main__":
    unittest.main()
