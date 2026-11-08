"""Pure scheduler/failure lane, without arithmetic or process imports."""

import copy
import itertools
import json
import struct
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.execution.worker_contract import (
    RESIDUE_VERSION,
    CostRow,
    Decoder,
    EngagementRule,
    LeaseRefusedError,
    Ledger,
    Manifest,
    Search,
    WorkerSettings,
    canonical,
    checked_residue,
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
        self.assertNotIn(0, search.completed)
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
        search.receive(0, 1, message(search, 0, 1, reason="error"))
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
    def test_missing_report_is_unknown_and_ack_is_exactly_once(self):
        ledger = Ledger()
        first = ledger.reserve("curve-a")
        second = ledger.reserve("curve-b")
        ledger.lose(first)
        ledger.settle(second, "curve-b", 13, 1.5)
        self.assertEqual((ledger.work, ledger.cpu), (13, 1.5))
        self.assertEqual(ledger.receipts[0]["status"], "lost")
        with self.assertRaises(KeyError):
            ledger.settle(second, "curve-b", 0, 0)
        ledger.charge(parent_cpu=4)
        ledger.reserve("next")
        self.assertEqual(ledger.work, 13)

    def test_bad_reports_do_not_change_observed_charges(self):
        ledger = Ledger()
        ticket = ledger.reserve("last")
        for identity, work, cpu in (
            ("wrong", 0, 0),
            ("last", -1, 0),
            ("last", 0, float("nan")),
        ):
            with self.assertRaises(ValueError):
                ledger.settle(ticket, identity, work, cpu)
            self.assertEqual((ledger.work, ledger.cpu), (0, 0))
        ledger.settle(ticket, "last", 8, 2)
        self.assertEqual((ledger.work, ledger.cpu), (8, 2))
        with self.assertRaises(TypeError):
            WorkerSettings(memory_bytes=2**20)
        with self.assertRaises(TypeError):
            WorkerSettings(checkpoint_bytes=4096)

    def test_frame_fragmentation_coalescing_corruption_and_oversize(self):
        raw = frame(dict(kind="fake", id=7))
        for boundary in range(len(raw)):
            decoder = Decoder()
            self.assertIsNone(decoder.feed(raw[:boundary]))
            self.assertEqual(
                decoder.feed(raw[boundary:]), dict(kind="fake", id=7)
            )
        decoder = Decoder()
        self.assertEqual(decoder.feed(raw + raw), dict(kind="fake", id=7))
        self.assertEqual(decoder.feed(b""), dict(kind="fake", id=7))
        for raw_bad in (
            b"\x00\x00\x00\x00",
            struct.pack("!I", 1) + b"x",
        ):
            with self.assertRaises(ValueError):
                Decoder().feed(raw_bad)
        self.assertEqual(
            Decoder().feed(frame(dict(payload="x" * 2000))),
            dict(payload="x" * 2000),
        )

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


def boundary(search, position, token):
    assignment = search.manifest.assignment(position)
    payload = dict(
        version=RESIDUE_VERSION,
        id=assignment["id"],
        metadata=search.manifest.metadata,
        n=assignment["n"],
        b1=assignment["b1"],
        b2=assignment["b2"],
        sigma=assignment["sigma"],
        value=[2, 3],
        a24=12,
        chain_chunks=0,
    )
    return dict(
        message(search, position, token, reason="boundary"),
        stage="stage_one",
        phase="stage_two_setup",
        sigma=assignment["sigma"],
        residue=dict(payload, sha256=digest(payload)),
    )


class StagedContractTests(unittest.TestCase):
    def test_queue_capacity_includes_inflight_and_preserves_tier_barriers(
        self,
    ):
        search = Search(
            manifest(),
            WorkerSettings(
                workers=2, stage_policy="staged", residue_capacity=2
            ),
            fixed_work=True,
        )
        search.dispatch(0, 1)
        search.dispatch(1, 2)
        self.assertIsNone(search.next())
        with self.assertRaises(LeaseRefusedError):
            search.dispatch(2, 3)
        search.receive(1, 2, boundary(search, 1, 2))
        self.assertEqual(search.next(), 1)
        search.receive(0, 1, boundary(search, 0, 1))
        search.dispatch(0, 3)
        done = dict(
            message(search, 0, 3, factor=5),
            stage="stage_two",
            attempt=search.continuation_attempts[0],
            boundary_sha256=search.boundaries[0]["sha256"],
        )
        search.receive(0, 3, done)
        self.assertEqual(search.next(), 2)
        self.assertNotIn(0, search.boundaries)
        self.assertNotIn(4, search.active)

    def test_checked_boundary_restart_retries_only_unfinished_continuation(
        self,
    ):
        search = Search(
            manifest(),
            WorkerSettings(stage_policy="staged", residue_capacity=1),
        )
        search.dispatch(0, 1)
        search.receive(0, 1, boundary(search, 0, 1))
        search.dispatch(0, 2)
        del search.active[
            0
        ]  # Crash: lease remains charged in the outer ledger.
        resumed = Search.restore(
            json.loads(canonical(search.snapshot())),
            search.manifest,
            search.settings,
        )
        self.assertEqual(resumed.next(), 0)
        self.assertEqual(resumed.stage(0), "stage_two")
        self.assertEqual(resumed.attempts[0], 1)
        resumed.dispatch(0, 3)
        self.assertEqual(resumed.continuation_attempts[0], 2)
        bad = dict(
            message(resumed, 0, 3),
            stage="stage_two",
            attempt=2,
            boundary_sha256="0" * 64,
        )
        with self.assertRaises(ValueError):
            resumed.receive(0, 3, bad)
        self.assertIn(0, resumed.active)

    def test_bad_residue_never_acknowledges_or_consumes_queue_capacity(self):
        for changes in (
            dict(value=[2, 5]),
            dict(value=[True, 3]),
            dict(a24=35),
            dict(b1=11),
            dict(sigma=17),
            dict(baby=[]),
            dict(value=[2]),
        ):
            search = Search(manifest(), WorkerSettings(stage_policy="staged"))
            search.dispatch(0, 1)
            msg = boundary(search, 0, 1)
            raw = dict(msg["residue"], **changes)
            raw["sha256"] = digest(
                {k: v for k, v in raw.items() if k != "sha256"}
            )
            with self.assertRaises(ValueError):
                search.receive(0, 1, dict(msg, residue=raw))
            self.assertIn(0, search.active)
            self.assertFalse(search.boundaries)
        search = Search(manifest(), WorkerSettings())
        search.dispatch(0, 1)
        with self.assertRaises(ValueError):
            search.receive(0, 1, boundary(search, 0, 1))
        with self.assertRaises(ValueError):
            checked_residue({}, manifest().assignment(0), manifest().metadata)

    def test_restart_rejects_unpaid_boundaries_and_reset_charges(self):
        search = Search(
            manifest(),
            WorkerSettings(stage_policy="staged", residue_capacity=1),
        )
        ledger = Ledger()
        token = ledger.reserve(search.manifest.assignment(0)["id"])
        search.dispatch(0, token)
        msg = boundary(search, 0, token)
        search.receive(0, token, msg)
        ledger.settle(token, msg["id"], 7, 0.2)
        ledger.receipts[0].update(
            stage="stage_one",
            attempt=1,
            reason="boundary",
            boundary_sha256=msg["residue"]["sha256"],
        )
        search.check_leases(ledger.receipts, dict(work=7, cpu=0.2))
        with self.assertRaises(ValueError):
            search.check_leases(ledger.receipts, dict(work=0, cpu=0.2))
        with self.assertRaises(ValueError):
            search.check_leases(ledger.receipts, dict(work=7, cpu=0))
        with self.assertRaises(ValueError):
            search.check_leases([], dict(work=7, cpu=0.2))
        with self.assertRaises(ValueError):
            search.check_leases(ledger.receipts * 2, dict(work=14, cpu=0.4))

    def test_old_checkpoint_version_and_future_tier_issue_are_rejected(self):
        search = Search(manifest(), WorkerSettings())
        saved = search.snapshot()
        saved["payload"]["version"] = "b12-f4-v1"
        saved["sha256"] = digest(saved["payload"])
        with self.assertRaises(ValueError):
            Search.restore(saved, search.manifest, search.settings)
        saved = search.snapshot()
        saved["payload"]["attempts"][4] = 1
        saved["sha256"] = digest(saved["payload"])
        with self.assertRaises(ValueError):
            Search.restore(saved, search.manifest, search.settings)


class OptimizationContractTests(unittest.TestCase):
    def test_math_identity_stays_v1_and_metadata_is_modulus_independent(self):
        first, second = manifest(n=35), manifest(n=1009)
        expected = digest(
            ["b12-f4-v1", first.n, first.seed, first.lineage, first.options]
        )
        self.assertEqual(first.identity, expected)
        self.assertEqual(first.metadata, second.metadata)
        self.assertNotEqual(first.schedule, second.schedule)
        assignment = first.assignment(0)
        assignment["seed"] = 0
        self.assertNotEqual(first.assignment(0)["seed"], 0)

    def test_explicit_wider_workers_and_simultaneous_cache_ownership(self):
        with self.assertRaises(ValueError):
            WorkerSettings(workers=6)
        self.assertEqual(
            WorkerSettings(workers=10, allow_wide_workers=True).workers, 10
        )
        settings = WorkerSettings(
            workers=10,
            allow_wide_workers=True,
            metadata_slots=2,
            decoded_cache_bytes=2**20,
        )
        self.assertGreaterEqual(settings.owned_bytes, 0)
        for changes in (
            dict(packet_curves=5),
            dict(poll_every=65),
            dict(metadata_slots=3),
            dict(parent_poll_seconds=1),
            dict(stage_policy="unknown"),
        ):
            with self.assertRaises(ValueError):
                adapt_resource_configuration(WorkerSettings, **changes)

    def test_cost_table_distinguishes_cold_pool_and_metadata_warmth(self):
        plan = manifest()
        row = CostRow(
            plan.schedule,
            10,
            1,
            (1, 1),
            pool_seconds=5,
            seconds_per_stage_two=(0.2, 0.2),
        )
        rule = EngagementRule(cost_table=(row,))
        self.assertFalse(rule.engages(plan, ready=False))
        self.assertTrue(rule.engages(plan, ready=False, pool_started=True))
        self.assertTrue(rule.engages(plan, ready=True))
        self.assertFalse(
            rule.engages(
                plan, ready=True, remaining=(0, 0), remaining_stage_two=(2, 0)
            )
        )
        no_stage_cost = EngagementRule(plan.schedule, 1, 1, (1, 1))
        self.assertFalse(
            no_stage_cost.engages(
                plan, ready=True, remaining=(0, 0), remaining_stage_two=(1, 0)
            )
        )
        invalid = EngagementRule(
            cost_table=(
                {
                    "schedule": plan.schedule,
                    "seconds_per_curve": [float("nan")],
                },
            )
        )
        self.assertFalse(invalid.engages(plan, ready=True))
        self.assertFalse(
            EngagementRule(cost_table=(row, row)).engages(plan, ready=True)
        )
        self.assertFalse(
            EngagementRule(cost_table=(row,) * 33).engages(plan, ready=True)
        )

    def test_decoder_recognizes_buffered_and_invalid_following_frames(self):
        raw = frame(dict(a=1))
        decoder = Decoder()
        self.assertEqual(decoder.feed(raw + raw), dict(a=1))
        self.assertTrue(decoder.complete)
        self.assertEqual(decoder.feed(b""), dict(a=1))
        self.assertFalse(decoder.complete)
        self.assertEqual(
            decoder.feed(raw + struct.pack("!I", 2048)), dict(a=1)
        )
        self.assertFalse(decoder.complete)
        self.assertIsNone(decoder.feed(b""))


class DurablePositionTests(unittest.TestCase):
    def test_double_digit_positions_survive_actual_json_round_trip(self):
        plan = manifest(((10, 100, 16),))
        settings = WorkerSettings()
        search = Search(plan, settings, fixed_work=True)
        for position in range(16):
            token = position + 1
            search.dispatch(position, token)
            search.receive(position, token, message(search, position, token))
        encoded = json.dumps(search.snapshot())
        restored = Search.restore(
            json.loads(encoded), plan, settings, fixed_work=True
        )
        self.assertEqual(set(restored.completed), set(range(16)))
        self.assertEqual(restored.attempts, search.attempts)

    def test_noncanonical_position_is_rejected_after_resealing(self):
        plan = manifest(((10, 100, 16),))
        settings = WorkerSettings()
        search = Search(plan, settings, fixed_work=True)
        search.dispatch(0, 1)
        search.receive(0, 1, message(search, 0, 1))
        saved = search.snapshot()
        saved["payload"]["completed"]["00"] = saved["payload"][
            "completed"
        ].pop("0")
        saved["sha256"] = digest(saved["payload"])
        with self.assertRaisesRegex(ValueError, "noncanonical"):
            Search.restore(saved, plan, settings, fixed_work=True)

    def test_sixteen_slot_traversal_round_trip_keeps_remaining_identity(self):
        settings = WorkerSettings(workers=16, allow_wide_workers=True)
        plan = manifest(((10, 100, 32),))
        search = Search(plan, settings, fixed_work=True)
        issued = []
        for slot in range(16):
            position = search.next(slot)
            token = slot + 1
            search.dispatch(position, token)
            issued.append((position, token))
        for position, token in reversed(issued):
            search.receive(position, token, message(search, position, token))
        restored = Search.restore(
            json.loads(json.dumps(search.snapshot())),
            plan,
            settings,
            fixed_work=True,
        )
        self.assertEqual(restored.next(0), 16)
        self.assertEqual(set(restored.completed), set(range(16)))
        self.assertEqual(restored.manifest.assignment(16), plan.assignment(16))
        for count in (10, 12, 14, 16):
            with self.assertRaises(ValueError):
                WorkerSettings(workers=count)
        with self.assertRaises(ValueError):
            WorkerSettings(workers=17, allow_wide_workers=True)
        self.assertEqual(
            WorkerSettings(workers=16, allow_wide_workers=True).workers, 16
        )


if __name__ == "__main__":
    unittest.main()
