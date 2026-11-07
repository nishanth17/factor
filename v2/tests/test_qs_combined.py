"""Independent current-worker flags retain exact paid traversal and resume."""

import copy
import hashlib
import json
import multiprocessing
import os
import unittest
from dataclasses import replace
from itertools import product
from pathlib import Path
from unittest.mock import patch

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.work_budget import PollingBudget
from v2.qs.current_parallel import CurrentParallelSIQSJob, QSWorkerConfig
from v2.qs.packet_lease import PacketLeaseBudget
from v2.qs.parallel import CollectionPool
from v2.qs.worker_packet import canonical, collect_packet
from v2.qs.worker_state import (
    WorkerState,
    compact_task,
    context_key,
    expand_task,
)
from v2.qs.worker_wire import PackedRows, decode, legacy_length
from v2.tests.test_current_qs_workers import (
    N,
    allowance,
    configuration,
    prepared,
    receive,
    settings,
)


def crash_packet(task, state=None):
    """Terminate only the worker deliberately assigned this test lease."""
    os._exit(7)


class PacketLeaseTests(unittest.TestCase):
    def test_mixed_reservations_match_incumbent_poll_prefixes(self):
        traces = []
        for cls in (PollingBudget, PacketLeaseBudget):
            raw = Budget(work_limit=100, seconds=None, cpu_seconds=None)
            trace = []
            view = cls(raw, poll=lambda b: trace.append(b.used))
            view.consume(0)
            view.consume(7)
            view.consume_span(50)
            view.consume_actions(4, 3)
            view.consume_span(10)
            view.consume(0)
            before = view.used
            with self.assertRaises(BudgetExhaustedError):
                view.consume(100)
            self.assertEqual(view.used, before)
            traces.append((view.used, trace, view.reason))
        self.assertEqual(traces[0], traces[1])

    def test_external_stop_allows_at_most_63_more_atomic_actions(self):
        stopped = [False]
        raw = Budget(
            work_limit=1000,
            seconds=None,
            cpu_seconds=None,
            cancelled=lambda: stopped[0],
        )
        view = PacketLeaseBudget(raw)
        view.consume()
        stopped[0] = True
        for _ in range(63):
            view.consume()
        with self.assertRaises(BudgetExhaustedError):
            view.consume()
        self.assertEqual(
            (view.used, raw.used, view.reason), (64, 64, "cancelled")
        )

    def test_forced_poll_and_group_boundaries_do_not_skip_stop(self):
        for operation in (
            lambda b: b.consume(0),
            lambda b: b.consume_actions(64, 1),
            lambda b: b.consume_span(64),
        ):
            raw = Budget(work_limit=1000, seconds=None, cpu_seconds=None)
            view = PacketLeaseBudget(raw)
            view.consume()
            raw.cancelled = lambda: True
            with self.assertRaises(BudgetExhaustedError):
                operation(view)
            self.assertEqual(view.used, 1)

    def test_private_wall_cpu_and_work_refusals_remain_exact(self):
        for kwargs, reason in (
            (dict(seconds=0), "wall_limit"),
            (dict(cpu_seconds=0), "cpu_limit"),
            (dict(work_limit=0), "work_limit"),
        ):
            defaults = dict(seconds=None, cpu_seconds=None, work_limit=100)
            defaults.update(kwargs)
            budget = PacketLeaseBudget(Budget(**defaults))
            with self.assertRaises(BudgetExhaustedError):
                budget.consume()
            self.assertEqual((budget.used, budget.reason), (0, reason))

    def test_leased_packet_matches_incumbent_atoms_coverage_and_work(self):
        for extended in (False, True):
            reference = prepared(extended)
            reference.run(fixed_work=True)
            leased = prepared(
                extended,
                worker_settings=replace(settings(), lease_budget=True),
            )
            leased.run(fixed_work=True)
            self.assertEqual(
                dict(reference.engine.collector._atoms),
                dict(leased.engine.collector._atoms),
            )
            for key in ("scanned", "worker_work", "useful_rows"):
                self.assertEqual(reference.stats[key], leased.stats[key])

    def test_pending_resume_pins_lease_policy_and_reconciles_once(self):
        job = prepared(worker_settings=replace(settings(), lease_budget=True))
        task = job._task(0, 10**8, 10, 10)
        report = collect_packet(task)
        receive(job, task, report)
        checkpoint = job.checkpoint()
        restored = CurrentParallelSIQSJob.from_checkpoint(
            checkpoint, budget=allowance()
        )
        self.assertEqual(restored.stats["lease_reconciliations"], 1)
        with self.assertRaises(ValueError):
            restored._receive(report, task, restored.receipts[0])
        with self.assertRaises(ValueError):
            CurrentParallelSIQSJob.from_checkpoint(
                checkpoint, budget=allowance(), workers=settings()
            )
        restored.run(fixed_work=True)
        reference = prepared()
        reference.run(fixed_work=True)
        self.assertEqual(
            dict(restored.engine.collector._atoms),
            dict(reference.engine.collector._atoms),
        )

    def test_stale_guard_rejected_and_lost_report_keeps_full_grants(self):
        job = prepared(worker_settings=replace(settings(), lease_budget=True))
        task = job._task(0, 10**8, 10, 10)
        report = collect_packet(task)
        altered = copy.deepcopy(report)
        altered["lease_policy"] = "peer"
        with self.assertRaises(ValueError):
            job._validate_report(altered, task)
        job.budget.consume(task["work"])
        before, cpu = job.budget.used, job.budget.prior_cpu
        receipt = dict(status="issued")
        job._lost(task, receipt)
        self.assertEqual(job.budget.used, before)
        self.assertEqual(receipt["charged_work"], task["work"])
        self.assertEqual(job.budget.prior_cpu, cpu + task["cpu"])
        job._lost(task, receipt)
        self.assertEqual(job.budget.prior_cpu, cpu + task["cpu"])

    def test_publication_stop_retains_a_completed_paid_prefix(self):
        job = prepared(worker_settings=replace(settings(), lease_budget=True))
        task = job._task(0, 10**8, 10, 10)
        reconcile = PacketLeaseBudget.reconcile

        def stop_at_publication(view):
            result = reconcile(view)
            view.budget.cancelled = lambda: True
            return result

        with patch.object(PacketLeaseBudget, "reconcile", stop_at_publication):
            report = collect_packet(task)
        self.assertEqual(report["publication_stop"], "cancelled")
        self.assertEqual(report["reason"], "complete")
        checked, used = job._validate_report(report, task)
        self.assertEqual(used, report["work"])
        self.assertTrue(checked["rows"])

    def test_spawned_lease_path_matches_local_and_retires_owned_children(self):
        owned_before = {p.pid for p in multiprocessing.active_children()}
        reference = prepared(
            worker_settings=replace(settings(), lease_budget=True)
        )
        reference.run(fixed_work=True)
        spawned = prepared(
            worker_settings=replace(settings("process", 2), lease_budget=True)
        )
        with CollectionPool("process", 2) as pool:
            result = spawned.run(pool=pool, fixed_work=True)
        self.assertTrue(pool.closed)
        self.assertTrue(result.stats["schedule_complete"])
        self.assertEqual(
            dict(reference.engine.collector._atoms),
            dict(spawned.engine.collector._atoms),
        )
        self.assertEqual(
            {p.pid for p in multiprocessing.active_children()}, owned_before
        )

    def test_crash_before_acknowledgement_keeps_prepaid_grants(self):
        job = prepared(
            worker_settings=replace(settings("process", 2), lease_budget=True)
        )
        with CollectionPool("process", 2) as pool:
            with patch("v2.qs.current_parallel.execute_packet", crash_packet):
                result = job.run(pool=pool, fixed_work=True)
        self.assertEqual(result.reason, "worker_failure")
        self.assertTrue(pool.closed)
        self.assertTrue(job.stats["lost_work"])
        self.assertGreaterEqual(job.budget.used, job.stats["lost_work"])
        self.assertTrue(all(r["status"] == "lost" for r in job.receipts))
        self.assertTrue(all(not p.is_alive() for p in pool._owned_processes))

    def test_retiring_worker_after_acknowledgement_preserves_resume(self):
        job = prepared(
            worker_settings=replace(settings("process", 2), lease_budget=True)
        )
        with CollectionPool("process", 2) as pool:
            job.run(pool=pool, max_assignments=1, fixed_work=True)
            checkpoint = job.checkpoint()
            pool.abort()
        restored = CurrentParallelSIQSJob.from_checkpoint(
            checkpoint, budget=allowance()
        )
        restored.run(fixed_work=True)
        reference = prepared()
        reference.run(fixed_work=True)
        self.assertEqual(
            dict(reference.engine.collector._atoms),
            dict(restored.engine.collector._atoms),
        )

    def test_all_off_settings_keep_the_frozen_legacy_record(self):
        path = Path(__file__).resolve().parents[1] / (
            "benchmarks/inputs/a7a/incumbent_worker_settings.json"
        )
        legacy = json.loads(path.read_text())
        current = QSWorkerConfig(**legacy)
        self.assertEqual(current.record(), legacy)
        self.assertNotIn("lease_budget", current.record())
        self.assertTrue(
            replace(current, lease_budget=True).record()["lease_budget"]
        )


class PersistentStateTests(unittest.TestCase):
    def test_retained_pool_rejects_disabling_or_repricing_live_state(self):
        job = prepared(worker_settings=replace(settings(), reuse_state=True))
        with CollectionPool("serial", 1) as pool:
            job.run(pool=pool, fixed_work=True)
            for worker_settings in (
                settings(),
                replace(
                    settings(), reuse_state=True, state_cache_bytes=32 * 2**20
                ),
            ):
                other = prepared(worker_settings=worker_settings)
                with self.assertRaises(ValueError):
                    other.run(pool=pool, fixed_work=True)

    def payload(self, job):
        base = job.source.base
        payload = dict(
            job=job.identity,
            execution=job.execution_identity,
            config=job.config,
            settings=job.workers,
            base=(base.n, base.multiplier, base.bound, base.entries),
            work=10**6,
        )
        payload["key"] = context_key(payload)
        return payload

    def test_compact_deltas_require_the_installed_identity(self):
        job = prepared(worker_settings=replace(settings(), reuse_state=True))
        payload = self.payload(job)
        runtime = WorkerState(payload, allowance())
        task = compact_task(job._task(0, 10**8, 10, 10), runtime.key)
        self.assertTrue({"base", "config", "settings"}.isdisjoint(task))
        expanded = expand_task(task, {"_qs_context": runtime})
        self.assertIs(expanded["_runtime"], runtime)
        self.assertIs(expanded["config"], job.config)
        with self.assertRaises(ValueError):
            expand_task(task, {})
        with self.assertRaises(ValueError):
            expand_task(
                dict(task, context_key="stale"), {"_qs_context": runtime}
            )
        with self.assertRaises(ValueError):
            WorkerState(dict(payload, key="stale"), allowance())

    def test_cache_is_an_explicit_refusable_owned_allowance(self):
        job = prepared(
            worker_settings=replace(
                settings(), reuse_state=True, state_cache_bytes=4096
            )
        )
        with self.assertRaises(MemoryError):
            WorkerState(self.payload(job), allowance())
        result = job.run(fixed_work=True)
        self.assertEqual(result.reason, "memory_limit")
        self.assertEqual(job.stats["state_install_work"], 0)

    def test_state_alone_and_with_leases_preserves_exact_serial_atoms(self):
        for extended in (False, True):
            reference = prepared(extended)
            reference.run(fixed_work=True)
            for leased in (False, True):
                for mode, count in (("serial", 1), ("process", 2)):
                    job = prepared(
                        extended,
                        worker_settings=replace(
                            settings(mode, count),
                            reuse_state=True,
                            lease_budget=leased,
                        ),
                    )
                    with CollectionPool(mode, count) as pool:
                        result = job.run(pool=pool, fixed_work=True)
                    self.assertTrue(result.stats["schedule_complete"])
                    self.assertEqual(
                        dict(reference.engine.collector._atoms),
                        dict(job.engine.collector._atoms),
                    )
                    self.assertEqual(
                        reference.stats["scanned"], job.stats["scanned"]
                    )
                    self.assertEqual(
                        job.stats["state_context_hits"],
                        job.stats["received_packets"],
                    )
                    if mode == "serial":
                        self.assertGreater(job.stats["state_family_reuses"], 0)
                        self.assertGreater(
                            job.stats["state_collector_reuses"], 0
                        )

    def test_quiescent_resume_rebuilds_state_without_replaying_paid_rows(self):
        job = prepared(worker_settings=replace(settings(), reuse_state=True))
        with CollectionPool("serial", 1) as pool:
            job.run(pool=pool, max_assignments=1, fixed_work=True)
            checkpoint = job.checkpoint()
        restored = CurrentParallelSIQSJob.from_checkpoint(
            checkpoint, budget=allowance()
        )
        before = restored.stats["state_install_work"]
        restored.run(fixed_work=True)
        self.assertGreater(restored.stats["state_install_work"], before)
        reference = prepared()
        reference.run(fixed_work=True)
        self.assertEqual(
            dict(restored.engine.collector._atoms),
            dict(reference.engine.collector._atoms),
        )

    def test_install_failure_cannot_refund_prepaid_context_work(self):
        job = prepared(worker_settings=replace(settings(), reuse_state=True))
        before = job.budget.used
        with CollectionPool("serial", 1) as pool:
            with patch(
                "v2.qs.worker_state.install_context",
                side_effect=RuntimeError("lost acknowledgement"),
            ):
                with self.assertRaises(RuntimeError):
                    job.run(pool=pool, fixed_work=True)
            self.assertTrue(pool.closed)
        self.assertGreater(job.stats["state_install_work"], 0)
        self.assertGreaterEqual(
            job.budget.used, before + job.stats["state_install_work"]
        )
        self.assertGreaterEqual(
            job.budget.prior_cpu, job.workers.assignment_cpu
        )

    def test_install_crash_retires_only_owned_workers_without_free_setup(self):
        owned_before = {p.pid for p in multiprocessing.active_children()}
        job = prepared(
            worker_settings=replace(settings("process", 2), reuse_state=True)
        )
        before = job.budget.used
        with CollectionPool("process", 2) as pool:
            with patch("v2.qs.worker_state.install_context", crash_packet):
                with self.assertRaises(Exception):
                    job.run(pool=pool, fixed_work=True)
        self.assertTrue(pool.closed)
        self.assertTrue(all(not p.is_alive() for p in pool._owned_processes))
        self.assertGreaterEqual(
            job.budget.used, before + job.stats["state_install_work"]
        )
        self.assertEqual(
            {p.pid for p in multiprocessing.active_children()}, owned_before
        )

    def test_pending_report_rejects_a_different_context_before_refund(self):
        job = prepared(worker_settings=replace(settings(), reuse_state=True))
        task = job._task(0, 10**8, 10, 10)
        report = collect_packet(task)
        report["context_key"] = "different context"
        before = job.budget.used
        with self.assertRaises(ValueError):
            job._validate_report(report, task)
        self.assertEqual(job.budget.used, before)


class PackedTransportTests(unittest.TestCase):
    def test_binary_roundtrip_and_exact_legacy_accounting(self):
        rows = []
        for gray in (0, 127, 128, 65535):
            for position in (-65536, -1, 0, 1, 65536):
                row = [
                    gray,
                    position,
                    -1 if position < 0 else 1,
                    [[2, 128], [999983, 1]],
                    1,
                    [1000000007, 1000000007],
                    2**96 + 1,
                    True,
                ]
                self.assertEqual(legacy_length(row), len(canonical(row)))
                rows.append(row)
        writer = PackedRows(65536)
        for row in rows:
            writer.append(row)
        raw = writer.finish()
        self.assertLess(len(raw), len(canonical(rows)))
        self.assertEqual(
            decode(
                raw,
                max_rows=64,
                max_exponents=8,
                max_bits=256,
                legacy_limit=65536,
            ),
            rows,
        )
        writer.reset()
        writer.append([0, 0, 1, [], 1, [], None, False])
        self.assertEqual(raw[:4], b"QSP1")
        self.assertEqual(int.from_bytes(raw[4:8], "little"), len(rows))

    def test_decoder_refuses_truncation_noncanonical_values_and_expansion(
        self,
    ):
        writer = PackedRows(4096)
        writer.append([0, 0, 1, [], 1, [], None, False])
        raw = writer.finish()
        corruptions = (
            raw[:-1],
            raw + b"x",
            b"bad!" + raw[4:],
            b"QSP1" + (1000).to_bytes(4, "little") + raw[8:],
            raw[:8] + b"\x80\x00" + raw[9:],
        )
        for payload in corruptions:
            with self.assertRaises(ValueError):
                decode(
                    payload,
                    max_rows=2,
                    max_exponents=3,
                    max_bits=128,
                    legacy_limit=4096,
                )
        with self.assertRaises(ValueError):
            decode(
                raw, max_rows=2, max_exponents=3, max_bits=128, legacy_limit=20
            )

    def test_packed_packet_preserves_legacy_cursor_signature_and_work(self):
        for limit in (4096, 65536):
            jobs = [
                prepared(
                    worker_settings=replace(
                        settings(), packet_bytes=limit, packed_transport=packed
                    )
                )
                for packed in (False, True)
            ]
            reports = [collect_packet(j._task(0, 10**8, 10, 10)) for j in jobs]
            checked = [
                j._validate_report(r, j._task(0, 10**8, 10, 10))[0]
                for j, r in zip(jobs, reports)
            ]
            self.assertEqual(checked[0]["rows"], checked[1]["rows"])
            for field in ("cursor", "scanned", "work", "reason"):
                self.assertEqual(reports[0][field], reports[1][field])
            self.assertLessEqual(
                len(reports[1]["packet"]), len(reports[0]["packet"])
            )

    def test_packed_pending_checkpoint_and_l_r_combinations(self):
        reference = prepared(True)
        reference.run(fixed_work=True)
        for leased, reused in (
            (False, False),
            (True, False),
            (False, True),
            (True, True),
        ):
            job = prepared(
                True,
                worker_settings=replace(
                    settings(),
                    packed_transport=True,
                    lease_budget=leased,
                    reuse_state=reused,
                ),
            )
            task = job._task(0, 10**8, 10, 10)
            receive(job, task, collect_packet(task))
            restored = CurrentParallelSIQSJob.from_checkpoint(
                job.checkpoint(), budget=allowance()
            )
            restored.run(fixed_work=True)
            self.assertEqual(
                dict(restored.engine.collector._atoms),
                dict(reference.engine.collector._atoms),
            )

    def test_self_consistent_wire_corruption_still_fails_exact_math(self):
        job = prepared(
            worker_settings=replace(settings(), packed_transport=True)
        )
        task = job._task(0, 10**8, 10, 10)
        report = collect_packet(task)
        checked, _ = job._validate_report(report, task)
        checked["rows"][0][2] *= -1
        writer = PackedRows(job.workers.packet_bytes)
        for row in checked["rows"]:
            writer.append(row)
        report["packet"] = writer.finish()
        report["sha256"] = hashlib.sha256(report["packet"]).hexdigest()
        receive(job, task, report)
        with self.assertRaises(ValueError):
            job._merge()
        self.assertFalse(job.engine.collector._atoms)

    def test_spawned_packed_combination_cleans_its_transient_buffer(self):
        reference = prepared(True)
        reference.run(fixed_work=True)
        job = prepared(
            True,
            worker_settings=replace(
                settings("process", 2),
                packed_transport=True,
                lease_budget=True,
                reuse_state=True,
            ),
        )
        with CollectionPool("process", 2) as pool:
            result = job.run(pool=pool, fixed_work=True)
        self.assertTrue(result.stats["schedule_complete"])
        self.assertEqual(
            dict(job.engine.collector._atoms),
            dict(reference.engine.collector._atoms),
        )


class VerificationTests(unittest.TestCase):
    def test_each_verification_choice_alone_and_with_l_r_t_matches_oracles(
        self,
    ):
        for extended in (False, True):
            reference = prepared(extended)
            reference.run(fixed_work=True)
            for mode in ("worker", "batch"):
                for combined in (False, True):
                    for backend, count in (("serial", 1), ("process", 2)):
                        job = prepared(
                            extended,
                            worker_settings=replace(
                                settings(backend, count),
                                verification=mode,
                                lease_budget=combined,
                                reuse_state=combined,
                                packed_transport=combined,
                            ),
                        )
                        with CollectionPool(backend, count) as pool:
                            result = job.run(pool=pool, fixed_work=True)
                        self.assertTrue(result.stats["schedule_complete"])
                        self.assertEqual(
                            dict(job.engine.collector._atoms),
                            dict(reference.engine.collector._atoms),
                        )
                        self.assertEqual(
                            job.stats["scanned"], reference.stats["scanned"]
                        )

    def test_worker_claim_cannot_bypass_independent_checked_insertion(self):
        job = prepared(
            True, worker_settings=replace(settings(), verification="worker")
        )
        task = job._task(0, 10**8, 10, 10)
        report = collect_packet(task)
        rows = json.loads(report["packet"])
        rows[0][2] *= -1
        report["packet"] = canonical(rows).encode("ascii")
        report["sha256"] = hashlib.sha256(report["packet"]).hexdigest()
        receive(job, task, report)
        with self.assertRaises(ValueError):
            job._merge()
        self.assertFalse(job.engine.collector._atoms)

    def test_batch_checks_all_members_before_their_graph_insertion(self):
        job = prepared(
            True, worker_settings=replace(settings(), verification="batch")
        )
        task = job._task(0, 10**8, 10, 10)
        report = collect_packet(task)
        rows = json.loads(report["packet"])
        self.assertGreater(len(rows), 1)
        rows[1][2] *= -1
        report["packet"] = canonical(rows).encode("ascii")
        report["sha256"] = hashlib.sha256(report["packet"]).hexdigest()
        receive(job, task, report)
        with self.assertRaises(ValueError):
            job._merge()
        self.assertFalse(job.engine.collector._atoms)

    def test_cancelled_batch_retains_a_paid_ordered_prefix_for_resume(self):
        job = prepared(
            True, worker_settings=replace(settings(), verification="batch")
        )
        task = job._task(0, 10**8, 10, 10)
        receive(job, task, collect_packet(task))
        admit = job.engine.collector._admit

        def cancel_after_first(atom, stats):
            refusal = admit(atom, stats)
            job.budget.cancelled = lambda: True
            return refusal

        with patch.object(job.engine.collector, "_admit", cancel_after_first):
            with self.assertRaises(BudgetExhaustedError):
                job._merge()
        self.assertEqual(job.pending[0]["admitted"], 1)
        checkpoint = job.checkpoint()
        restored = CurrentParallelSIQSJob.from_checkpoint(
            checkpoint, budget=allowance()
        )
        with self.assertRaises(ValueError):
            CurrentParallelSIQSJob.from_checkpoint(
                checkpoint,
                budget=allowance(),
                workers=replace(job.workers, verification="worker"),
            )
        restored.run(fixed_work=True)
        reference = prepared(True)
        reference.run(fixed_work=True)
        self.assertEqual(
            dict(restored.engine.collector._atoms),
            dict(reference.engine.collector._atoms),
        )


class SchedulingTests(unittest.TestCase):
    def test_all_fixed_order_flag_combinations_match_the_exact_oracle(self):
        reference = prepared(True)
        reference.run(fixed_work=True)
        expected = dict(reference.engine.collector._atoms)
        for leased, reused, packed, ahead, whole in product(
            (False, True), repeat=5
        ):
            for verification in ("parent", "worker", "batch"):
                worker_settings = replace(
                    settings(),
                    assignment_work=2 * 10**8,
                    lease_budget=leased,
                    reuse_state=reused,
                    packed_transport=packed,
                    lookahead=2 if ahead else 0,
                    whole_family=whole,
                    verification=verification,
                )
                job = prepared(True, worker_settings=worker_settings)
                result = job.run(fixed_work=True)
                self.assertTrue(result.stats["schedule_complete"])
                self.assertEqual(dict(job.engine.collector._atoms), expected)
                self.assertEqual(
                    job.stats["scanned"], reference.stats["scanned"]
                )

    def test_four_and_six_worker_combinations_cover_more_packets_than_slots(
        self,
    ):
        config = replace(configuration(True), family_count=6)
        reference = CurrentParallelSIQSJob(
            N,
            seed=29,
            config=config,
            budget=allowance(),
            workers=replace(settings(), polynomials_per_task=1),
        )
        reference.run(fixed_work=True)
        for count in (4, 6):
            worker_settings = replace(
                settings("process", count),
                polynomials_per_task=1,
                lookahead=count,
                lease_budget=True,
                reuse_state=True,
                packed_transport=True,
                verification="worker",
            )
            job = CurrentParallelSIQSJob(
                N,
                seed=29,
                config=config,
                budget=allowance(),
                workers=worker_settings,
            )
            with CollectionPool("process", count) as pool:
                result = job.run(pool=pool, fixed_work=True)
            self.assertTrue(result.stats["schedule_complete"])
            self.assertEqual(
                dict(job.engine.collector._atoms),
                dict(reference.engine.collector._atoms),
            )

    def test_small_whole_family_lease_keeps_its_paid_partial_checkpoint(self):
        job = prepared(worker_settings=replace(settings(), whole_family=True))
        first = job.run(fixed_work=True)
        self.assertEqual(first.reason, "work_limit")
        self.assertGreater(first.stats["scanned"], 0)
        restored = CurrentParallelSIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )
        for _ in range(4):
            result = restored.run(fixed_work=True)
            if result.stats["schedule_complete"]:
                break
        self.assertTrue(result.stats["schedule_complete"])
        reference = prepared()
        reference.run(fixed_work=True)
        self.assertEqual(
            dict(restored.engine.collector._atoms),
            dict(reference.engine.collector._atoms),
        )

    def test_lookahead_and_family_granularity_keep_the_polynomial_stream(self):
        reference = prepared(True)
        reference.run(fixed_work=True)
        for whole in (False, True):
            for combined in (False, True):
                job = prepared(
                    True,
                    worker_settings=replace(
                        settings("process", 2),
                        lookahead=2,
                        whole_family=whole,
                        lease_budget=combined,
                        reuse_state=combined,
                        packed_transport=combined,
                        verification="worker" if combined else "parent",
                    ),
                )
                with CollectionPool("process", 2) as pool:
                    result = job.run(pool=pool, fixed_work=True)
                self.assertTrue(result.stats["schedule_complete"])
                self.assertEqual(
                    job.stats["scanned"], reference.stats["scanned"]
                )
                self.assertEqual(
                    dict(job.engine.collector._atoms),
                    dict(reference.engine.collector._atoms),
                )

    def test_paid_pending_packets_beyond_one_wave_survive_resume(self):
        job = prepared(
            worker_settings=replace(
                settings("process", 2),
                lookahead=2,
                packed_transport=True,
            )
        )
        task = job._task(3, 10**8, 10, 10)
        receive(job, task, collect_packet(task))
        self.assertIn(3, job.pending)
        restored = CurrentParallelSIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )
        self.assertIn(3, restored.pending)
        restored.run(fixed_work=True)
        reference = prepared()
        reference.run(fixed_work=True)
        self.assertEqual(
            dict(restored.engine.collector._atoms),
            dict(reference.engine.collector._atoms),
        )

    def test_fixed_assignment_quota_does_not_turn_into_lookahead_work(self):
        job = prepared(
            worker_settings=replace(settings("process", 2), lookahead=4)
        )
        with CollectionPool("process", 2) as pool:
            result = job.run(pool=pool, max_assignments=1, fixed_work=True)
        self.assertEqual(result.stats["assignments"], 1)
        self.assertEqual(len(job.receipts), 1)
        self.assertEqual(job.stats["scanned"], 2 * 129)

    def test_whole_family_keeps_the_configured_prefix_not_extra_polynomials(
        self,
    ):
        job = prepared(
            worker_settings=replace(
                settings(), whole_family=True, assignment_work=2 * 10**8
            )
        )
        self.assertEqual(job.packet_width, job.config.gray_limit)
        self.assertEqual(job.assignment_count, len(job.source.assignments))
        result = job.run(fixed_work=True)
        self.assertEqual(result.stats["scanned"], 8 * 129)


class ArrivalPolicyTests(unittest.TestCase):
    def test_ordered_policy_cannot_be_bypassed_by_an_explicit_index(self):
        job = prepared(worker_settings=replace(settings(), lookahead=3))
        self.publication(job, 1)
        with self.assertRaises(ValueError):
            job._merge(1)

    def publication(self, job, index):
        """Replay a complete publication to isolate admission from solving."""
        task = job._task(index, 2 * 10**8, 10, 10)
        task["fixed_work"] = True
        report = collect_packet(task)
        self.assertEqual(report["reason"], "complete")
        report["divisor"] = None
        receive(job, task, report)

    def job(self):
        job = prepared(
            True,
            worker_settings=replace(
                settings("serial", 1),
                arrival_order=True,
                lookahead=3,
                packed_transport=True,
                lease_budget=True,
                verification="batch",
                assignment_work=2 * 10**8,
            ),
        )
        job.fixed_work = False
        return job

    def test_arrival_is_rejected_as_a_fixed_work_comparison(self):
        job = self.job()
        job.fixed_work = None
        with self.assertRaises(ValueError):
            job.run(fixed_work=True)

    def test_completed_ahead_witness_restores_until_frontier_catches_up(self):
        job = self.job()
        self.publication(job, 1)
        self.assertEqual(job._ready_assignment(), 1)
        with patch.object(job, "_solve"):
            job._merge(1)
        self.assertEqual(job.next_assignment, 0)
        self.assertIn(1, job.ahead)
        restored = CurrentParallelSIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )
        self.assertIn(1, restored.ahead)
        self.publication(restored, 0)
        with patch.object(restored, "_solve"):
            restored._merge(0)
        self.assertEqual(restored.next_assignment, 2)
        self.assertFalse(restored.ahead)
        self.assertEqual(restored.stats["assignments"], 2)

    def test_partial_admission_keeps_arrival_priority_across_resume(self):
        job = self.job()
        self.publication(job, 1)
        admit = job.engine.collector._admit

        def cancel(atom, stats):
            result = admit(atom, stats)
            job.budget.cancelled = lambda: True
            return result

        with patch.object(job.engine.collector, "_admit", cancel):
            with self.assertRaises(BudgetExhaustedError):
                job._merge(1)
        self.assertEqual(job.pending[1]["admitted"], 1)
        restored = CurrentParallelSIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )
        self.publication(restored, 0)
        self.assertEqual(restored._ready_assignment(), 1)

    def test_corrupt_arrival_history_cannot_hide_an_unpaid_frontier(self):
        job = self.job()
        self.publication(job, 1)
        with patch.object(job, "_solve"):
            job._merge(1)
        checkpoint = job.checkpoint()
        payload = json.loads(checkpoint["blob"])
        payload["receipts"][0]["arrival_seq"] = 10
        blob = canonical(payload)
        altered = dict(
            checkpoint,
            blob=blob,
            sha256=hashlib.sha256(blob.encode("ascii")).hexdigest(),
        )
        with self.assertRaises(ValueError):
            CurrentParallelSIQSJob.from_checkpoint(altered, budget=allowance())

    def test_spawned_trajectory_variant_returns_a_proper_factor(self):
        worker_settings = replace(
            settings("process", 2),
            arrival_order=True,
            lookahead=2,
            lease_budget=True,
            reuse_state=True,
            packed_transport=True,
            verification="worker",
        )
        job = CurrentParallelSIQSJob(
            N,
            seed=29,
            config=configuration(True),
            budget=allowance(),
            workers=worker_settings,
        )
        with CollectionPool("process", 2) as pool:
            result = job.run(pool=pool)
        self.assertIsNotNone(result.divisor)
        self.assertEqual(result.divisor * result.cofactor, N)
        self.assertTrue(all(not p.is_alive() for p in pool._owned_processes))


if __name__ == "__main__":
    unittest.main()
