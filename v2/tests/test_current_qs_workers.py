"""Small exact oracles for current-collector packets and failure ownership."""

import copy
import hashlib
import json
import multiprocessing
import os
import time
import unittest
from collections import defaultdict
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import Budget
from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration as resource_configuration,
)
from v2.benchmarks.infrastructure.resource_control import (
    historical_replace as replace,
)
from v2.qs.current_checkpoint import decode_config, encode_config
from v2.qs.current_parallel import CurrentParallelSIQSJob, QSWorkerConfig
from v2.qs.extended_policy import BoundedYieldSieveConfig
from v2.qs.families import _checksum
from v2.qs.parallel import CollectionPool
from v2.qs.relations import verify_combined
from v2.qs.sieve_collector import GrowingDoubleLargeSieveConfig
from v2.qs.siqs import SIQSConfig, make_family
from v2.qs.worker_packet import (
    canonical,
    collect_packet,
    decode_atom,
    polynomial_at,
)

N = 4001 * 5003


def allowance(work=2 * 10**9):
    return Budget(work_limit=work, seconds=15, cpu_seconds=60)


def configuration(extended=False, **kwargs):
    collector = (
        BoundedYieldSieveConfig if extended else GrowingDoubleLargeSieveConfig
    )(
        block_width=32,
        score_policy="powers",
        division="bucket",
        residual_bound=1000,
        large_prime_bound=1000,
        large_product_bound=1000000,
        sparse_queue_bytes=2**20,
        **(
            dict(
                batch_policy="off",
                split_policy="rho-first-v3",
                recompute_prime_limit=0,
            )
            if extended
            else {}
        ),
    )
    return resource_configuration(
        SIQSConfig,
        base_bound=100,
        factor_count=3,
        pool_size=8,
        family_count=2,
        half_width=64,
        max_half_width=64,
        assignment_policy="windowed_pair",
        polynomials_per_family=4,
        sign_policy="swap_high",
        polynomial_policy="parity2",
        root_update_policy="batch",
        batch_width=32,
        memory_bytes=64 * 2**20,
        checkpoint_bytes=4096,
        collector=collector,
        **kwargs,
    )


def settings(mode="serial", workers=1, **kwargs):
    return resource_configuration(
        QSWorkerConfig,
        mode=mode,
        workers=workers,
        packet_bytes=65536,
        max_batch_atoms=512,
        parent_work_reserve=1000,
        assignment_work=10**8,
        polynomials_per_task=2,
        **kwargs,
    )


def prepared(extended=False, worker_settings=None):
    job = CurrentParallelSIQSJob(
        N,
        seed=29,
        config=configuration(extended),
        workers=worker_settings or settings(),
        budget=allowance(),
    )
    job.fixed_work = True
    job._setup()
    return job


def packet(job, index=0):
    task = job._task(index, 10**8, 10, 10)
    return task, collect_packet(task)


def receive(job, task, report):
    job.budget.consume(task["work"])
    job.attempts[task["id"]] = task["attempt"]
    receipt = dict(
        id=task["id"],
        attempt=task["attempt"],
        token=task["token"],
        start=task["cursor"],
        work=task["work"],
        cpu_grant=task["cpu"],
        status="issued",
        charged_work=task["work"],
        split_grant=task["split_grant"],
    )
    job.receipts.append(receipt)
    job._receive(report, task, receipt)
    return receipt


def corrupt(checkpoint, change):
    result = copy.deepcopy(checkpoint)
    payload = json.loads(result["blob"])
    change(payload)
    result["blob"] = canonical(payload)
    result["sha256"] = _checksum(payload)
    return result


def crash_packet(task):
    os._exit(23)


def missing_packet(task):
    time.sleep(5)


class PacketTests(unittest.TestCase):
    def test_crt_oracle_matches_every_gray_and_multiplier_case(self):
        for multiplier in (1, 2, 3, 5):
            job = CurrentParallelSIQSJob(
                N,
                seed=29,
                config=configuration(multiplier=multiplier),
                workers=settings(),
                budget=allowance(),
            )
            job._setup()
            if job.source.base is None:
                continue
            base, primes = job.source.base, job.source.assignments[0]
            family = resource_configuration(
                make_family,
                base,
                primes,
                config=job.config,
                budget=job.budget,
                memory_bytes=64 * 2**20,
            )
            for ordinal in range(job.config.gray_limit):
                step = family.next()
                expected = polynomial_at(
                    base, primes, ordinal, job.config, job.budget
                )
                self.assertEqual(step.polynomial, expected)
                for entry, roots in zip(base.entries, step.roots):
                    actual = {
                        x
                        for x in range(entry.prime)
                        if expected.value(x) % entry.prime == 0
                    }
                    oracle = (
                        set(range(entry.prime))
                        if roots.all_positions
                        else set(roots.roots)
                    )
                    self.assertEqual(actual, oracle)

    def test_packet_atoms_pass_independent_trial_division(self):
        for extended in (False, True):
            job = prepared(extended)
            task, report = packet(job)
            self.assertEqual(report["reason"], "complete")
            self.assertEqual(report["scanned"], 2 * 129)
            rows = json.loads(report["packet"])
            self.assertTrue(rows)
            for row in rows:
                poly = polynomial_at(
                    job.source.base,
                    task["primes"],
                    row[0],
                    job.config,
                    job.budget,
                )
                atom = decode_atom(
                    row,
                    poly,
                    job.config.collector,
                    job.source.base,
                    job.budget,
                )
                value = abs(atom.u**2 - job.source.base.n_prime)
                factors = []
                for p in job.source.base.primes:
                    exponent = 0
                    while value % p == 0:
                        value //= p
                        exponent += 1
                    if exponent:
                        factors.append((p, exponent))
                self.assertEqual(tuple(atom.exponents), tuple(factors))
                self.assertEqual(
                    value, atom.residual_product * atom.square_correction**2
                )

    def test_serial_stream_and_spawned_worker_counts_keep_exact_coverage(self):
        expected = None
        owned_before = {p.pid for p in multiprocessing.active_children()}
        for mode, workers in (
            ("serial", 1),
            ("process", 1),
            ("process", 2),
            ("process", 4),
        ):
            job = prepared(worker_settings=settings(mode, workers))
            with CollectionPool(mode, workers) as pool:
                result = job.run(pool=pool, fixed_work=True)
            self.assertEqual(result.divisor * result.cofactor, N)
            self.assertTrue(result.stats["schedule_complete"])
            self.assertEqual(result.stats["scanned"], 8 * 129)
            signature = dict(job.engine.collector._atoms)
            if expected is None:
                expected = signature
            else:
                self.assertEqual(signature, expected)
        self.assertEqual(
            {p.pid for p in multiprocessing.active_children()}, owned_before
        )

    def test_wider_counts_cover_more_packets_than_slots_and_resume(self):
        expected = None
        owned_before = {p.pid for p in multiprocessing.active_children()}
        for count in (1, 3, 6, 8, 10):
            mode = "serial" if count == 1 else "process"
            job = CurrentParallelSIQSJob(
                N,
                seed=29,
                config=replace(configuration(), family_count=6),
                workers=replace(settings(mode, count), polynomials_per_task=1),
                budget=allowance(),
            )
            job._setup()
            job = CurrentParallelSIQSJob.from_checkpoint(
                job.checkpoint(), budget=allowance()
            )
            self.assertEqual(job.workers.workers, count)
            self.assertGreater(job.assignment_count, 10)
            with CollectionPool(mode, count) as pool:
                result = job.run(pool=pool, fixed_work=True)

            self.assertTrue(result.stats["schedule_complete"])
            self.assertEqual(
                result.stats["scanned"], job.assignment_count * 129
            )
            signature = dict(job.engine.collector._atoms)
            if expected is None:
                expected = signature
            else:
                self.assertEqual(signature, expected)
        self.assertEqual(
            {p.pid for p in multiprocessing.active_children()}, owned_before
        )

    def test_cross_packet_partials_and_every_combined_row_verify(self):
        job = prepared()
        per_packet = []
        for index in range(job.assignment_count):
            task, report = packet(job, index)
            receive(job, task, report)
            per_packet.append(
                set(
                    a.relation_id for a in job.engine.collector._atoms.values()
                )
            )
            job._merge()
        collector = job.engine.collector
        self.assertTrue(collector._combined)
        self.assertTrue(
            any(
                len(
                    {
                        job.contexts[collector._atoms[i].polynomial.identity][
                            0
                        ]
                        for i in row.atom_ids
                    }
                )
                > 1
                for row in collector._combined
            )
        )
        for row in collector._combined:
            self.assertTrue(
                verify_combined(
                    row,
                    job.source.base,
                    collector._atoms,
                    budget=allowance(),
                )
            )
        # Re-admission is idempotent, including graph edges.
        before = (len(collector._atoms), len(collector._rows))
        stats = defaultdict(int)
        for atom in list(collector._atoms.values()):
            collector._admit(atom, stats)
        self.assertEqual(before, (len(collector._atoms), len(collector._rows)))

    def test_packet_limit_keeps_prefix_then_continues_without_loss(self):
        job = prepared(worker_settings=replace(settings(), max_batch_atoms=2))
        result = job.run(fixed_work=True)
        self.assertTrue(result.stats["schedule_complete"])
        self.assertGreater(len(job.receipts), job.assignment_count)
        reference = prepared()
        reference.run(fixed_work=True)
        self.assertEqual(
            dict(job.engine.collector._atoms),
            dict(reference.engine.collector._atoms),
        )
        self.assertEqual(job.stats["scanned"], reference.stats["scanned"])

    def test_current_wide_windows_are_not_restricted_to_historical_caps(self):
        config = replace(
            configuration(), half_width=65536, max_half_width=65536
        )
        job = CurrentParallelSIQSJob(
            N, config=config, workers=settings(), budget=allowance()
        )
        self.assertEqual(job.config.half_width, 65536)

    def test_bounded_c9a_policy_roundtrip_preserves_every_field(self):
        config = configuration(True)
        self.assertEqual(decode_config(encode_config(config)), config)
        job = prepared(True)
        result = job.run(fixed_work=True)
        self.assertIn(result.reason, ("families_exhausted", "factor_found"))
        restored = CurrentParallelSIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )
        self.assertEqual(restored.config, config)
        self.assertEqual(
            len(restored.engine.collector._rows),
            len(job.engine.collector._rows),
        )


class OwnershipTests(unittest.TestCase):
    def test_bad_packets_cannot_refund_or_enter_parent_store(self):
        for field, value in (
            ("job", "stale-ecm"),
            ("kind", "ecm"),
            ("id", 8),
            ("work", -1),
            ("scanned", -1),
            ("cursor", [7, 0]),
            ("packet", b"[]"),
        ):
            job = prepared()
            task, report = packet(job)
            report[field] = value
            used = job.budget.used
            with self.assertRaises(ValueError):
                job._validate_report(report, task)
            self.assertEqual(job.budget.used, used)
            self.assertFalse(job.engine.collector._atoms)

    def test_bad_atom_context_and_exact_payload_are_rejected(self):
        job = prepared()
        task, report = packet(job)
        rows = json.loads(report["packet"])
        rows[0][2] *= -1
        report["packet"] = canonical(rows).encode()
        report["sha256"] = hashlib.sha256(report["packet"]).hexdigest()
        receive(job, task, report)
        with self.assertRaises(ValueError):
            job._merge()
        self.assertFalse(job.engine.collector._atoms)

    def test_duplicate_report_is_not_refunded_twice(self):
        job = prepared()
        task, report = packet(job)
        receipt = receive(job, task, report)
        used = job.budget.used
        with self.assertRaises(ValueError):
            job._receive(report, task, receipt)
        self.assertEqual(job.budget.used, used)

    def test_root_stop_refuses_dispatch_without_issuing_packets(self):
        job = prepared()
        job.budget.cancelled = lambda: True
        result = job.run(fixed_work=True)
        self.assertEqual(result.reason, "cancelled")
        self.assertFalse(job.receipts)

    def test_quota_refusal_and_cancellation_keep_reported_prefix(self):
        job = prepared()
        task, _ = packet(job)
        task["work"] = 140000
        report = collect_packet(task)
        self.assertEqual(report["reason"], "complete")
        self.assertGreater(report["work"], task["work"])
        cancelled = collect_packet(task, cancelled=lambda: True)
        self.assertEqual(cancelled["reason"], "cancelled")
        self.assertEqual(cancelled["scanned"], 0)

    def test_pool_reuse_and_close_guards(self):
        with CollectionPool("process", 2) as pool:
            for _ in range(2):
                job = prepared(worker_settings=settings("process", 2))
                job.run(pool=pool, fixed_work=True)
                self.assertFalse(pool.lock.locked())
                self.assertFalse(pool.closed)
        self.assertTrue(pool.closed)
        with self.assertRaises(ValueError):
            pool.begin(allowance())

    def test_worker_crash_retains_lease_and_releases_owned_children(self):
        job = prepared(worker_settings=settings("process", 2))
        before = {p.pid for p in multiprocessing.active_children()}
        with CollectionPool("process", 2) as pool:
            with patch("v2.qs.current_parallel.execute_packet", crash_packet):
                result = job.run(pool=pool, fixed_work=True)
            self.assertEqual(result.reason, "worker_failure")
            self.assertTrue(pool.closed)
        self.assertTrue(any(r["status"] == "lost" for r in job.receipts))
        self.assertGreater(job.stats["lost_packets"], 0)
        self.assertEqual(
            {p.pid for p in multiprocessing.active_children()}, before
        )

    def test_missing_report_has_a_bounded_owned_pool_retirement(self):
        job = prepared(
            worker_settings=replace(settings("process", 2), drain_seconds=0.02)
        )
        job.budget.seconds = 0.05
        with CollectionPool("process", 2) as pool:
            with patch(
                "v2.qs.current_parallel.execute_packet", missing_packet
            ):
                result = job.run(pool=pool, fixed_work=True)
            self.assertEqual(result.reason, "worker_failure")
            self.assertTrue(pool.closed)
        self.assertTrue(all(r["status"] == "lost" for r in job.receipts))

    def test_abort_keeps_created_children_after_executor_table_is_lost(self):
        before = {p.pid for p in multiprocessing.active_children()}
        with CollectionPool("process", 2) as pool:
            for _ in range(2):
                pool.executor.submit(missing_packet, {})
            processes = tuple(pool._owned_processes)
            self.assertEqual(len(processes), 2)
            pool.executor._processes = None
            pool.abort()
            self.assertTrue(all(not p.is_alive() for p in processes))
        self.assertEqual(
            {p.pid for p in multiprocessing.active_children()}, before
        )

    def test_portfolio_recursive_reconstruction_and_worker_marker_resume(self):
        from v2.portfolio import PortfolioConfig, factorize_bounded

        config = resource_configuration(
            PortfolioConfig,
            trial_bound=100,
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=(),
            siqs=configuration(),
            memory_bytes=128 * 2**20,
        )
        workers = settings("process", 2)
        with CollectionPool("process", 2) as pool:
            run = factorize_bounded(
                97 * N,
                seed=29,
                config=config,
                qs_workers=workers,
                qs_pool=pool,
                budget=allowance(),
            )
        self.assertEqual(run.result.reconstruct(), 97 * N)
        self.assertTrue(run.result.complete)
        self.assertTrue(run.result.proven)
        self.assertTrue(any(e["stage"] == "siqs" for e in run.events))
        paused = factorize_bounded(
            N,
            seed=29,
            config=config,
            qs_workers=workers,
            budget=Budget(work_limit=50, seconds=5, cpu_seconds=20),
        )
        self.assertFalse(paused.result.complete)
        with self.assertRaisesRegex(ValueError, "worker settings changed"):
            factorize_bounded(
                N,
                config=config,
                checkpoint=paused.checkpoint,
                budget=allowance(),
            )
        resumed = factorize_bounded(
            N,
            config=config,
            qs_workers=workers,
            checkpoint=paused.checkpoint,
            budget=allowance(),
        )
        self.assertEqual(resumed.result.reconstruct(), N)
        self.assertTrue(resumed.result.complete)

    def test_pending_checkpoint_cursor_and_corruption(self):
        job = prepared()
        task, report = packet(job)
        receive(job, task, report)
        checkpoint = job.checkpoint()
        restored = CurrentParallelSIQSJob.from_checkpoint(
            checkpoint, budget=allowance()
        )
        self.assertEqual(restored.pending[0]["rows"], job.pending[0]["rows"])
        used = restored.budget.used
        restored.run(fixed_work=True)
        self.assertGreater(restored.budget.used, used)
        for mutation in (
            lambda p: p.update(next_assignment=99),
            lambda p: p["receipts"][0].update(charged_work=10**15),
            lambda p: p["pending"][0][1].update(admitted=99999),
            lambda p: p.update(contexts={"bad": [0, 0]}),
        ):
            with self.assertRaises(ValueError):
                CurrentParallelSIQSJob.from_checkpoint(
                    corrupt(checkpoint, mutation), budget=allowance()
                )
        with self.assertRaisesRegex(ValueError, "settings changed"):
            CurrentParallelSIQSJob.from_checkpoint(
                checkpoint, budget=allowance(), workers=settings("process", 2)
            )

    def test_paused_merge_restores_exact_store_and_charges_rebuild(self):
        job = prepared()
        task, report = packet(job)
        receive(job, task, report)
        original = job.engine.collector._admit
        calls = []

        def pause(atom, stats):
            calls.append(atom)
            if len(calls) == 3:
                raise MemoryError("pause publication")
            return original(atom, stats)

        with patch.object(job.engine.collector, "_admit", pause):
            with self.assertRaises(MemoryError):
                job._merge()
        self.assertEqual(job.pending[0]["admitted"], 2)
        restored = CurrentParallelSIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )
        self.assertEqual(restored.pending[0]["admitted"], 2)
        restored.run(fixed_work=True)
        reference = prepared()
        reference.run(fixed_work=True)
        self.assertEqual(
            dict(restored.engine.collector._atoms),
            dict(reference.engine.collector._atoms),
        )

    def test_unknown_native_or_nested_backends_are_refused(self):
        with self.assertRaises(ValueError):
            settings("thread", 2)
        for count in (0, -1, 11, True):
            with self.assertRaises(ValueError):
                QSWorkerConfig(workers=count)
        job = CurrentParallelSIQSJob(
            N, config=configuration(), workers=settings(), budget=allowance()
        )
        self.assertEqual(
            job.run(max_assignments=1, fixed_work=True).reason, "paused"
        )

    def test_removed_worker_grants_are_rejected_by_name(self):
        from v2.execution.resources import RemovedResourceLimitError

        for name in (
            "memory_bytes",
            "worker_memory_bytes",
            "checkpoint_bytes",
            "context_limit",
            "receipt_limit",
        ):
            with self.subTest(name=name):
                with self.assertRaises(RemovedResourceLimitError) as caught:
                    QSWorkerConfig(**{name: 32 * 2**30})
                self.assertEqual(caught.exception.argument, name)


class SharedCPUReadingTests(unittest.TestCase):
    def test_snapshot_observes_live_updates_and_new_run_baselines(self):
        from v2.qs.worker_packet import _process_cpu_reader

        context = multiprocessing.get_context("spawn")
        state = dict(
            cpu=context.Array("d", (3.0, 5.0, 2.0)),
            cpu_base=context.Array("d", (1.0, 1.5, 3.0)),
        )
        used = _process_cpu_reader(state)
        self.assertEqual(used(), 5.5)
        state["cpu"][0], state["cpu"][2] = 4.0, 6.0
        self.assertEqual(used(), 9.5)

        for index, value in enumerate(state["cpu"]):
            state["cpu_base"][index] = value
        resumed = _process_cpu_reader(state)
        self.assertEqual(resumed(), 0)
        state["cpu"][1] += 2.0
        self.assertEqual(resumed(), 2.0)

    def test_adapter_still_checks_peer_cpu_parent_cpu_and_external_stop(self):
        from v2.qs import parallel
        from v2.qs.worker_packet import execute_packet

        context = multiprocessing.get_context("spawn")

        state = dict(
            cpu=context.Array("d", (0.0, 0.0)),
            cpu_base=context.Array("d", (0.0, 0.0)),
            rss=context.Array("q", 2),
            stop=context.Event(),
            parent_cpu=context.Value("d", 0),
            cpu_limit=context.Value("d", 3),
        )

        def probe(task, *, cancelled, publish_cpu):
            self.assertFalse(cancelled())
            state["cpu"][1] = 4.0
            self.assertTrue(cancelled())
            state["cpu"][1] = 0.0
            state["parent_cpu"].value = 4.0
            self.assertTrue(cancelled())
            state["parent_cpu"].value = 0.0
            state["stop"].set()
            self.assertTrue(cancelled())
            return "checked"

        with (
            patch.object(parallel, "_PROCESS_STATE", (state, 0)),
            patch("v2.qs.worker_packet.collect_packet", side_effect=probe),
        ):
            self.assertEqual(execute_packet({}), "checked")


class PrefixGeometryTests(unittest.TestCase):
    def test_packet_sizes_preserve_requested_positions_and_exact_atoms(self):
        from v2.benchmarks.qs.a7a.cell import prefix

        expected = None
        for width in (1, 2, 4, 16, 64):
            budget = allowance()
            budget.a7a_seed = 29
            # The tiny ownership fixtures fund two-polynomial packets.
            # Geometry equality needs one finite grant that also funds a
            # whole four-polynomial fixture family, without a lease pause.
            workers = replace(
                settings(),
                polynomials_per_task=width,
                assignment_work=10**9,
            )
            row = prefix(
                N,
                configuration(),
                workers,
                budget,
                999,
                None,
                polynomials=8,
            )
            self.assertEqual(row["stats"]["scanned"], 8 * 129)
            self.assertGreater(row["stats"]["admitted_atoms"], 0)
            if expected is None:
                expected = row["signature"]
            self.assertEqual(row["signature"], expected)

    def test_unaligned_or_unbounded_quota_refuses_before_consuming_work(self):
        from v2.benchmarks.qs.a7a.cell import prefix

        for quota in (0, True, 8193, 9, 3):
            budget = allowance()
            budget.a7a_seed = 29
            with self.assertRaises(ValueError):
                prefix(
                    N,
                    configuration(),
                    settings(),
                    budget,
                    1,
                    None,
                    polynomials=quota,
                )
            self.assertEqual(budget.used, 0)


class SplitLeaseTests(unittest.TestCase):
    def test_retained_helpers_share_the_root_rss_limit(self):
        from v2.execution.budget import Budget as RootBudget
        from v2.execution.resources import ProcessSample
        from v2.tests.test_h3_resources import Clock, Provider

        clock, provider = Clock(), Provider()
        budget = RootBudget(
            provider=provider, wall_clock=clock, cpu_clock=clock, rss_bytes=500
        )
        job = CurrentParallelSIQSJob(
            N, config=configuration(), workers=settings(), budget=budget
        )
        provider.rows.append(ProcessSample(999999, 1, 450, 0))
        clock.value = 1
        stopped = job.run(fixed_work=True)
        self.assertEqual(stopped.reason, "rss_limit")
        self.assertEqual(stopped.resource_stop["observed"], 550)
        self.assertFalse(job.receipts)

    def test_unused_helpers_need_no_hypothetical_reservation(self):
        from v2.execution.budget import Budget as RootBudget
        from v2.tests.test_h3_resources import Clock, Provider

        clock, provider = Clock(), Provider()
        budget = RootBudget(
            provider=provider, wall_clock=clock, cpu_clock=clock, rss_bytes=500
        )
        job = CurrentParallelSIQSJob(
            N, config=configuration(), workers=settings(), budget=budget
        )
        result = job.run(max_assignments=1, fixed_work=True)
        self.assertEqual(result.reason, "paused")
        self.assertEqual(result.stats["observed_rss_bytes"], 100)

    def test_cumulative_split_quota_is_not_reset_per_packet_or_restore(self):
        config = configuration()
        config = replace(
            config, collector=replace(config.collector, split_call_limit=2)
        )
        job = CurrentParallelSIQSJob(
            N, config=config, workers=settings(), budget=allowance()
        )
        job.run(fixed_work=True)
        self.assertLessEqual(job.split_used, 2)
        self.assertEqual(job.split_used, job.engine.collector._split_calls)
        restored = CurrentParallelSIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )
        self.assertEqual(restored.split_used, job.split_used)

    def test_lost_split_reservation_cannot_be_refunded(self):
        job = prepared()
        task, _ = packet(job)
        task["split_grant"] = 7
        receipt = dict(status="issued")
        job._lost(task, receipt)
        self.assertEqual(job.split_used, 7)
        self.assertEqual(receipt["split_calls"], 7)
