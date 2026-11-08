"""Independent metadata traversal and exact checked-worker contracts."""

import copy
import json
import multiprocessing
import unittest
from dataclasses import replace
from math import isqrt
from types import SimpleNamespace

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.execution.budget import BudgetExhaustedError
from v2.execution.ecm_stage_two import StageTwoSchedule, advance_prepared_job
from v2.execution.ecm_worker import CurveEngine
from v2.execution.ecm_workers import ECMExecutor, arithmetic_options
from v2.execution.stage_jobs import advance_job, new_job
from v2.execution.worker_contract import (
    Ledger,
    Manifest,
    Search,
    WorkerSettings,
    digest,
)
from v2.tests.resource_stops import ActionBudget, ProspectiveActionBudget
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_ecm_workers import allowance, configuration, task


def trial_primes(lo, hi):
    return [
        p
        for p in range(max(2, lo), hi)
        if all(p % d for d in range(2, isqrt(p) + 1))
    ]


def record(result):
    return tuple(
        result[k] for k in ("id", "sigma", "factor", "phase", "chain_chunks")
    )


class PreparedIndexTests(unittest.TestCase):
    def test_independent_giant_walk_matches_every_prime_and_batch_boundary(
        self,
    ):
        for b1, b2, size in ((7, 257, 17), (20, 400, 32), (51, 2000, 97)):
            origin = b1 if b1 % 2 else b1 - 1
            distance = min(isqrt(b2), (origin - 1) // 2)
            schedule = StageTwoSchedule(size, 2**20)
            for repetition in range(2):
                seen = []
                center = origin
                for lo in range(b1 + 1, b2 + 1, 2 * size):
                    hi = min(b2 + 1, lo + 2 * size)
                    values = trial_primes(lo, hi)
                    cursor = dict(left=lo, next=hi, values=values, index=0)
                    while cursor["index"] < len(values):
                        start = cursor["index"]
                        while values[start] > center + 2 * distance:
                            center += 2 * distance
                        quota = (3, 11, 128)[start % 3]
                        expected = []
                        for p in values[start:]:
                            if (
                                p > center + 2 * distance
                                or len(expected) == quota
                            ):
                                break
                            # Enumerate the even baby multiples independently.
                            offsets = [
                                i
                                for i in range(1, distance + 1)
                                if center + 2 * i == p
                            ]
                            self.assertEqual(len(offsets), 1)
                            expected.append(offsets[0])
                        end, actual = schedule.take(
                            cursor,
                            b1,
                            b2,
                            distance,
                            center,
                            quota,
                            allowance(),
                        )
                        self.assertEqual(list(actual[start:end]), expected)
                        self.assertEqual(end, start + len(expected))
                        seen.extend(values[start:end])
                        cursor["index"] = end
                self.assertEqual(seen, trial_primes(b1 + 1, b2 + 1))
            self.assertGreater(schedule.hits, 0)
            self.assertGreaterEqual(schedule.used_bytes, 0)

    def test_group_reads_reuse_one_immutable_block_without_slice_allocation(
        self,
    ):
        schedule = StageTwoSchedule(16, 65536)
        cursor = dict(left=8, next=40, values=trial_primes(8, 40), index=0)
        short_end, short = schedule.take(cursor, 7, 45, 3, 7, 1, allowance())
        full_end, full = schedule.take(cursor, 7, 45, 3, 7, 128, allowance())
        self.assertIs(short, full)
        self.assertEqual(short_end, 1)
        self.assertEqual(full_end, 2)
        self.assertEqual(short[:full_end], (2, 3))

    def test_interrupted_build_never_publishes_partial_metadata(self):
        schedule = StageTwoSchedule(16, 65536)
        cursor = dict(left=8, next=40, values=trial_primes(8, 40), index=0)
        polls = []
        budget = Budget(
            work_limit=10000,
            seconds=None,
            cpu_seconds=None,
            cancelled=lambda: polls.append(1) or len(polls) >= 3,
        )
        before = schedule.used_bytes
        with self.assertRaisesRegex(BudgetExhaustedError, "cancelled"):
            schedule.take(cursor, 7, 45, 3, 7, 128, budget)
        self.assertFalse(schedule.blocks)
        self.assertEqual(schedule.build_bytes, 0)
        self.assertEqual(schedule.used_bytes, before)
        self.assertEqual(cursor["index"], 0)
        self.assertGreater(budget.used, 0)
        self.assertIsNotNone(
            schedule.take(cursor, 7, 45, 3, 7, 128, allowance())
        )

    def test_funded_scratch_refuses_retention_without_partial_allocation(self):
        schedule = StageTwoSchedule(16, 8192)
        cursor = dict(left=8, next=40, values=trial_primes(8, 40), index=0)
        self.assertIsNone(schedule.take(cursor, 7, 45, 3, 7, 128, allowance()))
        self.assertEqual(schedule.used_bytes, 8192)
        self.assertEqual(schedule.refusals, 1)
        self.assertFalse(schedule.blocks)
        with self.assertRaises(MemoryError):
            StageTwoSchedule(16, 8191)

    def test_exact_term_polynomial_and_nonunit_recovery(self):
        config = configuration(ecm_tiers=((7, 45, 1),), gcd_batch=2)
        initial = new_job("ecm", 35, 7, 7, 45)
        initial.update(
            phase="stage_two",
            distance=3,
            center=7,
            baby=[None, [2, 1], [3, 5], [4, 7]],
            giant=[1, 0],
            previous=[1, 1],
            value=[1, 1],
            a24=1,
            cursor=dict(left=8, next=14, hi=46, values=[11, 13], index=0),
        )
        reference, prepared = copy.deepcopy(initial), copy.deepcopy(initial)
        old_context = SimpleNamespace()
        new_context = SimpleNamespace(
            stage_two_schedule=StageTwoSchedule(16, 65536)
        )
        advance_job(reference, allowance(), old_context, config)
        advance_prepared_job(prepared, allowance(), new_context, config)
        self.assertEqual(prepared["terms"], [5, 7])
        self.assertEqual(prepared["product"], (5 * 7) % 35)
        self.assertEqual(prepared, reference)
        for _ in range(5):
            if reference["done"]:
                break
            advance_job(reference, allowance(), old_context, config)
            advance_prepared_job(prepared, allowance(), new_context, config)
            self.assertEqual(prepared, reference)
        self.assertTrue(prepared["done"])
        self.assertEqual(prepared["factor"], 5)

    def test_cursor_semantics_and_changed_provider_are_rejected(self):
        cache = StageTwoSchedule(16, 65536)
        cursor = dict(left=8, next=40, values=trial_primes(8, 40), index=0)
        cache.take(cursor, 7, 45, 3, 7, 128, allowance())
        with self.assertRaises(ValueError):
            cache.take(cursor, 7, 45, 3, 13, 128, allowance())
        cursor["values"] = cursor["values"][:-1]
        with self.assertRaises(ValueError):
            cache.take(cursor, 7, 45, 3, 7, 128, allowance())


class CheckedBudgetTests(unittest.TestCase):
    def test_prospective_matches_budget_oracle_with_large_exact_counters(self):
        for used in (0, 2**63 - 1, 10**30):
            plain = Budget(
                work_limit=10**40, used=used, seconds=None, cpu_seconds=None
            )
            fast = ProspectiveActionBudget(
                poll_every=1,
                work_limit=10**40,
                used=used,
                seconds=None,
                cpu_seconds=None,
            )
            for amount in (0, 1, 2, 128, 2**63 - 1, 0):
                plain.consume(amount)
                fast.consume(amount)
                self.assertEqual(fast.used, plain.used)
        for cls in (ActionBudget, ProspectiveActionBudget):
            budget = cls(
                poll_every=64, work_limit=3, seconds=None, cpu_seconds=None
            )
            budget.consume(3)
            with self.assertRaisesRegex(BudgetExhaustedError, "work_limit"):
                budget.consume(1)
            self.assertEqual(budget.used, 3)
            for bad in (True, -1, None, 1.5, 2**63):
                with self.assertRaises(ValueError):
                    budget.consume(bad)
                self.assertEqual(budget.used, 3)

    def test_deferred_boundaries_keep_cancellation_and_expiry_precedence(self):
        for cls in (ActionBudget, ProspectiveActionBudget):
            stopped = [False]
            budget = cls(
                poll_every=4,
                work_limit=100,
                seconds=None,
                cpu_seconds=None,
                cancelled=lambda: stopped[0],
            )
            budget.consume(1)
            stopped[0] = True
            with self.assertRaisesRegex(BudgetExhaustedError, "cancelled"):
                for _ in range(4):
                    budget.consume(1)
            self.assertLessEqual(budget.used, 4)
            budget.used = 100
            with self.assertRaisesRegex(BudgetExhaustedError, "cancelled"):
                budget.consume(1)
            for limits, reason in (
                (dict(seconds=0), "wall_limit"),
                (dict(cpu_seconds=0), "cpu_limit"),
            ):
                expired = cls(poll_every=64, work_limit=0, **limits)
                with self.assertRaisesRegex(BudgetExhaustedError, reason):
                    expired.consume(1)
                self.assertEqual(expired.used, 0)


class HeadroomWorkerTests(unittest.TestCase):
    def test_candidates_preserve_curves_and_reuse_one_context(
        self,
    ):
        config = configuration(ecm_tiers=((20, 400, 4),))
        expected = None
        for options in (
            {},
            dict(budget_mode="prospective"),
            dict(bind_context=True),
            dict(stage_two_cache_bytes=65536),
            dict(stage_two_cache_bytes=8192),
            dict(
                stage_two_cache_bytes=65536,
                bind_context=True,
                budget_mode="prospective",
            ),
        ):
            settings = adapt_resource_configuration(
                WorkerSettings,
                worker_imports="ecm",
                decoded_cache_bytes=65536,
                **options,
            )
            engine = CurveEngine(settings)
            plan = Manifest(
                1000003, 7919, "root", arithmetic_options(config, settings)
            )
            rows = []
            for repeat in range(2):
                actual = [
                    record(
                        engine.execute(task(plan, i), cancelled=lambda: False)
                    )
                    for i in range(4)
                ]
                self.assertEqual(actual, rows or actual)
                rows = actual
            self.assertEqual(engine.preparations, 1)
            if options.get("stage_two_cache_bytes") == 65536:
                self.assertGreater(engine.stage_two_hits, 0)
                built = engine.stage_two_build_bytes
                engine.execute(task(plan, 0), cancelled=lambda: False)
                self.assertEqual(engine.stage_two_build_bytes, built)
            if options.get("stage_two_cache_bytes") == 8192:
                self.assertGreater(engine.stage_two_refusals, 0)
            self.assertEqual(rows, expected or rows)
            expected = rows

    def test_resource_bounds_do_not_change_search_identity(self):
        config = configuration()
        base = WorkerSettings()
        changed = replace(
            base,
            stage_two_cache_bytes=65536,
            bind_context=True,
            budget_mode="prospective",
        )
        self.assertEqual(
            changed.owned_bytes - base.owned_bytes, 2 * (65536 + 4096)
        )
        plans = [
            Manifest(1000003, 7, "root", arithmetic_options(config, s))
            for s in (base, changed)
        ]
        self.assertEqual(plans[0].identity, plans[1].identity)
        for values in (
            dict(bind_context=1),
            dict(budget_mode="bad"),
            dict(packet_policy="bad"),
        ):
            with self.assertRaises(ValueError):
                adapt_resource_configuration(WorkerSettings, **values)
        with self.assertRaises(MemoryError):
            arithmetic_options(
                configuration(segment_size=4096),
                replace(base, stage_two_cache_bytes=65536),
            )

    def test_new_wire_counters_remain_bounded_and_untrusted(self):
        config = configuration()
        settings = WorkerSettings(stage_two_cache_bytes=65536)
        plan = Manifest(
            1000003, 7, "root", arithmetic_options(config, settings)
        )
        valid = CurveEngine(settings).execute(
            task(plan, 0), cancelled=lambda: False
        )
        for field in (
            "stage_two_hits",
            "stage_two_misses",
            "stage_two_refusals",
            "stage_two_build_bytes",
        ):
            for value in (True, -1, 2**63):
                search = Search(plan, settings, fixed_work=True)
                search.dispatch(0, 1)
                damaged = dict(valid, **{field: value})
                with self.assertRaises(ValueError):
                    search.receive(0, 1, damaged)
        search = Search(plan, settings, fixed_work=True)
        search.dispatch(0, 1)
        with self.assertRaises(ValueError):
            search.receive(0, 1, dict(valid, stage_two_points=[1, 2]))

    def test_old_resource_version_is_rejected_with_a_valid_digest(self):
        settings = WorkerSettings()
        plan = Manifest(
            1000003, 7, "root", arithmetic_options(configuration(), settings)
        )
        checkpoint = Search(plan, settings, fixed_work=True).snapshot()
        checkpoint["payload"]["version"] = "b12-f4-v5"
        checkpoint["sha256"] = digest(checkpoint["payload"])
        with self.assertRaisesRegex(ValueError, "unsupported"):
            Search.restore(checkpoint, plan, settings, fixed_work=True)

    def test_tail_packets_keep_singletons_for_remaining_quota(self):
        for curves, policy, expected in (
            (8, "tail_one", 2),
            (4, "tail_one", 1),
            (4, "fixed", 2),
        ):
            settings = WorkerSettings(
                workers=2, packet_curves=2, packet_policy=policy
            )
            config = configuration(ecm_tiers=((20, 200, curves),))
            plan = Manifest(
                1000003, 7, "root", arithmetic_options(config, settings)
            )
            search = Search(plan, settings, fixed_work=True)
            with ECMExecutor(settings=settings) as owner:
                _, issued = owner._packet(
                    plan,
                    search,
                    Ledger(10**8, 100),
                    10**20,
                    dict(index=0, metadata_cache=[plan.metadata]),
                    peers=2,
                )
            self.assertEqual(len(issued), expected)

    def test_spawn_checked_restart_and_no_children_with_combined_candidates(
        self,
    ):
        before = {p.pid for p in multiprocessing.active_children()}
        config = configuration(ecm_tiers=((20, 400, 4),))
        settings = WorkerSettings(
            workers=2,
            worker_imports="ecm",
            decoded_cache_bytes=65536,
            stage_two_cache_bytes=65536,
            bind_context=True,
            budget_mode="prospective",
            packet_curves=2,
            packet_policy="tail_one",
        )
        with ECMExecutor(
            seed=7919, settings=settings, force_workers=True
        ) as owner:
            first = owner.run(
                1000003, config=config, budget=allowance(), fixed_work=True
            )
        self.assertEqual(first.reason, "curves_exhausted")
        self.assertGreater(first.stats["stage_two_hits"], 0)
        checkpoint = json.loads(json.dumps(first.checkpoint))
        with ECMExecutor(
            seed=7919, settings=settings, force_workers=True
        ) as owner:
            resumed = owner.run(
                1000003,
                config=config,
                budget=allowance(),
                checkpoint=checkpoint,
                fixed_work=True,
            )
        self.assertEqual(
            checkpoint["search"]["payload"]["completed"],
            resumed.checkpoint["search"]["payload"]["completed"],
        )
        with ECMExecutor(
            seed=7919,
            settings=replace(settings, workers=1),
            force_workers=True,
        ) as owner:
            with self.assertRaises(ValueError):
                owner.run(
                    1000003,
                    config=config,
                    budget=allowance(),
                    checkpoint=checkpoint,
                    fixed_work=True,
                )
        self.assertEqual(
            {p.pid for p in multiprocessing.active_children()}, before
        )
