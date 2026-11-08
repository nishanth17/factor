"""Deterministic resource ownership; no elapsed-time assertions."""

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.resources import (
    ProcessSample,
    RemovedResourceLimitError,
    atomic_save,
)


class Clock:
    def __init__(self):
        self.value = 0.0

    def __call__(self):
        return self.value


class Provider:
    def __init__(self):
        self.rows = [ProcessSample(os.getpid(), 1, 100, 0)]

    def physical_memory(self):
        return 1000

    def sample(self, roots):
        return list(self.rows)


class ResourceTests(unittest.TestCase):
    def setUp(self):
        self.wall, self.cpu, self.provider = Clock(), Clock(), Provider()

    def budget(self, **kwargs):
        return Budget(
            wall_clock=self.wall,
            cpu_clock=self.cpu,
            provider=self.provider,
            **kwargs,
        )

    def test_defaults_and_work_counter(self):
        budget = self.budget()
        self.assertIsNone(budget.seconds)
        self.assertIsNone(budget.cpu_seconds)
        self.assertEqual(budget.rss_bytes, 800)
        budget.consume(10**100)
        self.assertEqual(budget.used, 10**100)

    def test_qs_packet_capacity_continues_identical_assignment(self):
        from dataclasses import replace

        from v2.qs.parallel import (
            CollectionPool,
            ParallelConfig,
            ParallelSIQSJob,
            _collect,
        )

        config = ParallelConfig(family_count=2, pool_size=8)
        job = ParallelSIQSJob(4001 * 5003, config=config, budget=self.budget())
        job._setup()
        results = []
        for capacity in (1, config.max_batch_atoms):
            selected = replace(config, max_batch_atoms=capacity)
            task = (
                job.base,
                job._assignment_primes(0),
                0,
                selected,
                None,
                True,
            )
            with CollectionPool() as pool:
                result = _collect(task, pool.state)
            packets = (*result["packets"], result["atoms"])
            self.assertTrue(all(len(packet) <= capacity for packet in packets))
            self.assertEqual(result["reason"], "complete")
            results.append(
                (
                    result,
                    tuple(
                        atom.relation_id
                        for packet in packets
                        for atom in packet
                    ),
                )
            )
        self.assertTrue(results[0][0]["packets"])
        self.assertEqual(results[0][1], results[1][1])
        for field in ("work", "scanned", "divisor"):
            self.assertEqual(results[0][0][field], results[1][0][field])

    def test_removed_work_argument(self):
        with self.assertRaisesRegex(RemovedResourceLimitError, "work_limit"):
            self.budget(work_limit=10**20)

    def test_precedence_and_payload(self):
        budget = self.budget(seconds=1, cpu_seconds=1, rss_bytes=100)
        self.wall.value = self.cpu.value = 2
        budget.cancelled = lambda: True
        with self.assertRaisesRegex(BudgetExhaustedError, "cancelled"):
            budget.consume()
        budget.cancelled = None
        for name, limit in (("wall", 1), ("cpu", 1), ("rss", 100)):
            with self.assertRaisesRegex(BudgetExhaustedError, name + "_limit"):
                budget.publication()
            self.assertEqual(budget.stop["limit"], name)
            self.assertEqual(budget.stop["limit_value"], limit)
            self.assertEqual(budget.stop["scope"], "factorization")
            self.assertGreaterEqual(budget.stop["observed"], limit)
            setattr(
                budget,
                {"wall": "seconds", "cpu": "cpu_seconds", "rss": "rss_bytes"}[
                    name
                ],
                None,
            )

    def test_sampled_current_tree_rss(self):
        budget = self.budget(rss_bytes=500)
        self.provider.rows.append(ProcessSample(123, 2, 450, 1))
        with self.assertRaisesRegex(BudgetExhaustedError, "rss_limit"):
            budget.publication()
        self.assertEqual(budget.stop["observed"], 550)
        self.provider.rows.pop()
        budget.publication()
        self.assertEqual(budget.observed_rss, 100)

    def test_live_ack_exit_no_double_charge(self):
        budget = self.budget()
        self.provider.rows.append(ProcessSample(123, 2, 200, 3))
        budget.publication()
        self.assertEqual(budget.cpu_used, 3)
        budget.observe_child((123, 2), 3)
        budget.observe_child((123, 2), 5)
        budget.observe_child((123, 2), 5)
        self.assertEqual(budget.cpu_used, 5)
        self.provider.rows.pop()
        budget.publication()
        self.assertEqual(budget.cpu_used, 5)
        budget.observe_child((123, 3), 2)
        self.assertEqual(budget.cpu_used, 7)

    def test_close_start_excludes_idle(self):
        budget = self.budget()
        self.wall.value, self.cpu.value = 4, 2
        budget.close()
        self.wall.value, self.cpu.value = 104, 52
        self.assertEqual((budget.wall_used, budget.cpu_used), (4, 2))
        budget.start()
        self.wall.value, self.cpu.value = 107, 53
        budget.close().close()
        self.assertEqual((budget.wall_used, budget.cpu_used), (7, 3))

    def test_saved_history_and_retained_owner(self):
        budget = self.budget(cpu_seconds=60)
        budget.restore_charges(work=10, wall=20, cpu=40)
        self.cpu.value = 15
        self.assertEqual(budget.cpu_used, 55)
        budget.restore_charges(work=10, wall=20, cpu=40)
        self.assertEqual(budget.cpu_used, 55)
        fresh = self.budget(cpu_seconds=60)
        fresh.restore_charges(work=10, wall=20, cpu=40)
        self.assertEqual(fresh.cpu_used, 40)

    def test_saved_prefix_adds_fresh_startup_and_keeps_unsaved_costs(self):
        budget = self.budget(cpu_seconds=100)
        budget.consume(5)
        self.wall.value, self.cpu.value = 10, 7
        budget.close()
        self.wall.value, self.cpu.value = 110, 57
        budget.start()
        budget.restore_charges(work=100, wall=20, cpu=40)
        self.assertEqual(
            (budget.used, budget.wall_used, budget.cpu_used), (105, 30, 47)
        )
        saved = budget.charges()
        self.wall.value, self.cpu.value = 112, 60
        budget.close()
        budget.start()
        budget.restore_charges(
            work=saved["work"], wall=saved["wall"], cpu=saved["cpu"]
        )
        self.assertEqual((budget.wall_used, budget.cpu_used), (32, 50))
        budget.restore_charges(work=100, wall=20, cpu=40)
        self.assertEqual(
            (budget.used, budget.wall_used, budget.cpu_used), (105, 32, 50)
        )

    def test_fresh_stale_resume_inherits_only_recorded_history(self):
        for _ in range(3):
            owner = self.budget()
            owner.restore_charges(work=10, wall=20, cpu=40)
            self.assertEqual(
                (owner.used, owner.wall_used, owner.cpu_used), (10, 20, 40)
            )
            owner.close()

    def test_recorded_encoding_overhead_is_not_added_twice(self):
        from v2.execution.resources import record_checkpoint

        owner = self.budget()
        owner.consume(5)
        self.wall.value, self.cpu.value = 12, 10
        record_checkpoint(
            owner,
            {
                "payload": {"work_used": 5, "wall_used": 10, "cpu_used": 7},
                "serialization_overhead": [2, 3],
            },
        )
        owner.close()
        self.wall.value, self.cpu.value = 112, 60
        owner.start()
        owner.restore_charges(work=5, wall=12, cpu=10)
        self.assertEqual(
            (owner.used, owner.wall_used, owner.cpu_used), (5, 12, 10)
        )

    def test_family_resource_prefix_is_remembered(self):
        from v2.execution.resources import record_checkpoint

        owner = self.budget()
        owner.consume(5)
        self.wall.value, self.cpu.value = 12, 10
        record_checkpoint(
            owner,
            {
                "payload": {
                    "resources": {
                        "work_used": 5,
                        "wall_used": 12,
                        "cpu_used": 10,
                    }
                }
            },
        )
        owner.close().start()
        owner.restore_charges(work=5, wall=12, cpu=10)
        self.assertEqual(
            (owner.used, owner.wall_used, owner.cpu_used), (5, 12, 10)
        )

    def test_invalid_resource_observations_are_provider_errors(self):
        from v2.execution.resources import ResourceProviderError

        budget = self.budget()
        self.wall.value = float("nan")
        with self.assertRaises(ResourceProviderError):
            budget.consume()
        self.wall.value = -1
        with self.assertRaises(ResourceProviderError):
            budget.consume()
        self.wall.value = 0
        self.provider.rows[0] = ProcessSample(os.getpid(), 1, -1, 0)
        with self.assertRaises(ResourceProviderError):
            budget.publication()

    def test_atomic_optional_oversized_and_interrupted_save(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.json"
            old = {"published": "old"}
            atomic_save(path, old)
            large = {"published": "x" * (2**20 + 1)}
            self.assertIsNone(atomic_save(path, large))
            self.assertEqual(json.loads(path.read_text()), large)
            atomic_save(path, old)
            event = atomic_save(path, large, checkpoint_bytes=2**20)
            self.assertEqual(event["reason"], "checkpoint_size_skip")
            self.assertGreater(event["requested"], event["limit_value"])
            self.assertEqual(json.loads(path.read_text()), old)
            with patch("os.replace", side_effect=OSError("injected")):
                with self.assertRaisesRegex(OSError, "injected"):
                    atomic_save(path, large)
            self.assertEqual(json.loads(path.read_text()), old)
            self.assertEqual(list(Path(directory).iterdir()), [path])


class StorageMigrationTests(unittest.TestCase):
    def test_internal_limit_arguments_are_named_errors(self):
        from v2.portfolio import PortfolioConfig
        from v2.qs import SieveConfig, SIQSConfig
        from v2.qs.extended_policy import BoundedYieldSieveConfig

        for factory, argument in (
            (SieveConfig, "max_atoms"),
            (SieveConfig, "max_relations"),
            (SieveConfig, "memory_bytes"),
            (SIQSConfig, "checkpoint_bytes"),
            (PortfolioConfig, "memory_bytes"),
            (BoundedYieldSieveConfig, "spill_bytes"),
        ):
            with self.subTest(factory=factory, argument=argument):
                with self.assertRaisesRegex(
                    RemovedResourceLimitError, argument
                ):
                    factory(**{argument: None})

    def test_removed_positional_limits_are_named_before_validation(self):
        from v2.execution.worker_contract import (
            WorkerSettings,
            checked_residue,
        )
        from v2.qs.sieve_collector import SieveConfig
        from v2.qs.sss import SSSConfig

        for factory, count, argument in (
            (WorkerSettings, 2, "assignment_work"),
            (SieveConfig, 11, "max_relations"),
            (SSSConfig, 7, "max_candidates"),
            (checked_residue, 4, "limit"),
        ):
            with self.subTest(factory=factory.__name__):
                with self.assertRaises(RemovedResourceLimitError) as failed:
                    factory(*([None] * count))
                self.assertEqual(failed.exception.argument, argument)

    def test_graph_65536_is_not_an_admission_ceiling(self):
        from v2.qs.large_primes import LargePrimeForest, graph_reserve

        class BoundaryDict(dict):
            def __len__(self):
                return 65536

        graph = LargePrimeForest()
        graph.edges = BoundaryDict()
        plan = graph.plan_link(101, 103, "boundary", Budget(rss_bytes=None))
        graph.commit_link(plan)
        self.assertIn("boundary", graph.edges)
        self.assertGreater(graph_reserve(65537), graph_reserve(65536))

    def test_matrix_estimate_is_diagnostic(self):
        from v2.qs.linear_algebra import filter_matrix

        with patch(
            "v2.qs.linear_algebra.matrix_workspace", return_value=10**30
        ):
            matrix = filter_matrix((3, 3), budget=Budget(rss_bytes=None))
        self.assertEqual(matrix.original_rows, (3, 3))
        self.assertEqual(matrix.workspace_bytes, 10**30)

    def test_old_checkpoint_refused_and_saved_resume(self):
        from v2.execution.resources import UnsupportedResourceVersionError
        from v2.qs.siqs import SIQSConfig, SIQSJob

        job = SIQSJob(
            1009 * 1013, config=SIQSConfig(), budget=Budget(rss_bytes=None)
        )
        job.run(max_blocks=1)
        snapshot = job.checkpoint()
        old = dict(snapshot)
        del old["resource_version"]
        with self.assertRaises(UnsupportedResourceVersionError):
            SIQSJob.from_checkpoint(old, budget=Budget(rss_bytes=None))
        restored = SIQSJob.from_checkpoint(
            snapshot, budget=Budget(rss_bytes=None)
        )
        self.assertEqual(
            restored.engine.next_position, job.engine.next_position
        )
        self.assertGreaterEqual(restored.budget.used, job.budget.used)

    def test_stage_share_hands_off_with_unfunded_reserve(self):
        import contextlib
        import io

        from v2.execution.allocation import ECMAllocation
        from v2.portfolio import PortfolioConfig, factorize_bounded
        from v2.qs.siqs import SIQSConfig

        output = io.StringIO()
        config = PortfolioConfig(
            allocation=ECMAllocation(
                "pretest", pretest_work=0, fallback_work=10**100
            ),
            siqs=SIQSConfig(mode="qs", base_bound=100),
            ecm_chain_mode="off",
        )
        with (
            contextlib.redirect_stdout(output),
            contextlib.redirect_stderr(output),
        ):
            run = factorize_bounded(
                1009 * 1013,
                config=config,
                budget=Budget(rss_bytes=None),
                ecm_executor=None,
            )
        self.assertEqual(output.getvalue(), "")
        self.assertEqual(run.result.reconstruct(), 1009 * 1013)
        self.assertNotIn(
            run.reason, ("pretest_exhausted", "insufficient_fallback_work")
        )


class ChildReaperTests(unittest.TestCase):
    def test_live_ack_exit_and_crash_use_one_total(self):
        from types import SimpleNamespace

        provider = Provider()
        provider.rows.append(ProcessSample(123, 10, 50, 2))
        budget = Budget(
            provider=provider, wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        process = SimpleNamespace(
            pid=123, _popen=SimpleNamespace(pid=123, returncode=None)
        )
        budget.track_process(process)
        budget.acknowledge(123, 3)
        usage = SimpleNamespace(ru_utime=4, ru_stime=1)
        with patch("os.wait4", return_value=(123, 9, usage)) as wait:
            self.assertEqual(process._popen.poll(), -9)
            self.assertEqual(process._popen.poll(), -9)
            wait.assert_called_once_with(123, os.WNOHANG)
        # The initial retained CPU was a baseline, so only 5 - 2 is new.
        self.assertEqual(budget.cpu_used, 3)
        budget.acknowledge(123, 5)
        self.assertEqual(budget.cpu_used, 3)

    def test_retained_process_rebinds_and_excludes_idle_cpu(self):
        from types import SimpleNamespace

        provider = Provider()
        provider.rows.append(ProcessSample(123, 10, 50, 2))
        first = Budget(
            provider=provider, wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        process = SimpleNamespace(
            pid=123, _popen=SimpleNamespace(pid=123, returncode=None)
        )
        first.track_process(process)
        first.close()
        provider.rows[-1] = ProcessSample(123, 10, 50, 7)
        second = Budget(
            provider=provider, wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        second.track_process(process)
        usage = SimpleNamespace(ru_utime=8, ru_stime=1)
        with patch("os.wait4", return_value=(123, 0, usage)):
            self.assertEqual(process._popen.poll(0), 0)
        self.assertEqual(first.cpu_used, 0)
        self.assertEqual(second.cpu_used, 2)


class IncidentCoverageTests(unittest.TestCase):
    def test_75_80_metadata_is_not_a_memory_grant(self):
        from v2.qs.siqs import SIQSConfig

        for estimate in (17541809920, 17721011200):
            with self.subTest(estimate=estimate):
                with patch(
                    "v2.qs.siqs.root_update_workspace", return_value=estimate
                ):
                    config = SIQSConfig(
                        assignment_policy="flyer",
                        factor_count=8,
                        pool_size=64,
                        family_count=100000,
                        polynomials_per_family=64,
                        root_update_policy="batch",
                    )
                with patch(
                    "v2.qs.siqs.root_update_workspace", return_value=estimate
                ):
                    self.assertGreater(config.metadata_reserve, 12884901888)
        provider = Provider()
        budget = Budget(
            provider=provider, wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        budget.publication()
        provider.rows[0] = ProcessSample(os.getpid(), 1, 801, 0)
        with self.assertRaises(BudgetExhaustedError) as stopped:
            budget.publication()
        self.assertEqual(
            stopped.exception.payload,
            dict(
                reason="rss_limit",
                limit="rss",
                observed=801,
                limit_value=800,
                scope="factorization",
            ),
        )

    def test_unattributed_incident_candidate_sites_have_no_byte_veto(self):
        import ast

        import v2.qs

        directory = Path(v2.qs.__file__).parent
        names = (
            "sieve_collector",
            "large_prime_collector",
            "extended_collector",
            "opaque_conversion",
            "yield_memory",
            "checkpoint",
            "large_prime_checkpoint",
            "extended_checkpoint",
            "linear_algebra",
            "factor_base",
            "families",
            "assignment_stream",
            "root_updates",
            "multiplier",
            "norm_multiplier",
            "relations",
            "recovery_batch",
            "smooth_batch",
            "sss",
            "sss_checkpoint",
            "parallel",
        )
        for name in names:
            with self.subTest(candidate=name):
                tree = ast.parse((directory / (name + ".py")).read_text())
                for node in ast.walk(tree):
                    if isinstance(node, ast.If):
                        predicate = ast.unparse(node.test)
                        self.assertFalse(
                            any(
                                isinstance(item, ast.Attribute)
                                and item.attr == "memory_limit"
                                for item in ast.walk(node.test)
                            )
                        )
                        self.assertNotIn("config.memory_bytes", predicate)

    def test_actual_cpu_does_not_preapprove_packet_grants(self):
        from v2.ecm.workers import ECMExecutor, arithmetic_options
        from v2.execution.worker_contract import (
            Ledger,
            Manifest,
            Search,
            WorkerSettings,
        )
        from v2.portfolio import PortfolioConfig

        settings = WorkerSettings(workers=4, packet_curves=3)
        config = PortfolioConfig(
            ecm_tiers=((30, 100, 12),), ecm_chain_mode="off"
        )
        manifest = Manifest(
            1009 * 1013, 7, "incident", arithmetic_options(config, settings)
        )
        search = Search(manifest, settings, fixed_work=True)
        ledger = Ledger()
        executor = ECMExecutor(settings=settings)
        from v2.ecm.worker import CurveEngine

        for index in range(4):
            packet, issued = executor._packet(
                manifest,
                search,
                ledger,
                None,
                dict(index=index, metadata_cache=[]),
                peers=4 - index,
            )
            self.assertEqual(len(issued), 3)
            self.assertNotIn("cpu", packet["jobs"][0])
            self.assertNotIn("work", packet["jobs"][0])
            engine = CurveEngine(settings)
            for job, (position, token) in zip(packet["jobs"], issued):
                task = dict(
                    job,
                    schedule=manifest.metadata,
                    options=packet.get("options"),
                    packet=packet["packet"],
                )
                result = engine.execute(task, cancelled=lambda: False)
                self.assertEqual(result["reason"], "complete")
                search.receive(position, token, result)
                ledger.settle(token, result["id"], result["work"], 0)
        self.assertEqual(len(search.completed), 12)
        self.assertFalse(ledger.active)
        cpu = Clock()
        budget = Budget(
            cpu_seconds=440,
            cpu_clock=cpu,
            wall_clock=lambda: 0,
            provider=Provider(),
        )
        cpu.value = 100
        budget.consume(1)
        cpu.value = 440
        with self.assertRaisesRegex(BudgetExhaustedError, "cpu_limit"):
            budget.consume(1)

    def test_cancellation_before_publication_keeps_store_valid(self):
        from collections import Counter

        from v2.qs.factor_base import build_factor_base
        from v2.qs.polynomial import qs_polynomial
        from v2.qs.reference_collector import collect_block
        from v2.qs.sieve_collector import SieveCollector
        from v2.qs.storage import require_published

        owner = Budget(
            provider=Provider(), wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        base = build_factor_base(
            1009 * 1013, bound=100, budget=owner
        ).factor_base
        polynomial = qs_polynomial(base)
        atom = collect_block(polynomial, base, 0, 1, budget=owner).relations[0]
        collector = SieveCollector(polynomial, base, budget=owner)
        original = owner.publication
        owner.publication = lambda: (_ for _ in ()).throw(
            BudgetExhaustedError("cancelled")
        )
        with self.assertRaises(BudgetExhaustedError):
            collector._admit(atom, Counter())
        require_published(collector)
        self.assertFalse(collector._atoms)
        owner.publication = original
        collector._admit(atom, Counter())
        require_published(collector)
        self.assertEqual(tuple(collector._atoms), (atom.relation_id,))

    def test_allocation_interruption_poison_is_not_an_rss_stop(self):
        from collections import Counter

        from v2.qs.factor_base import build_factor_base
        from v2.qs.polynomial import qs_polynomial
        from v2.qs.reference_collector import collect_block
        from v2.qs.sieve_collector import SieveCollector, SieveConfig
        from v2.qs.storage import RelationPublicationError, require_published

        budget = Budget(
            provider=Provider(), wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        base = build_factor_base(
            1009 * 1013, bound=100, budget=budget
        ).factor_base
        polynomial = qs_polynomial(base)
        atom = collect_block(polynomial, base, 0, 1, budget=budget).relations[
            0
        ]
        owner = SieveCollector(
            polynomial, base, config=SieveConfig(), budget=budget
        )

        class BrokenList(list):
            def append(self, value):
                raise MemoryError("injected publication")

        owner._full = BrokenList()
        with self.assertRaises(RelationPublicationError) as failed:
            owner._admit(atom, Counter())
        self.assertEqual(failed.exception.reason, "allocation_failure")
        self.assertIsInstance(failed.exception.__cause__, MemoryError)
        self.assertIsNone(budget.stop)
        with self.assertRaises(RelationPublicationError):
            require_published(owner)


class ResumeOrdinalTests(unittest.TestCase):
    def test_cancellation_resumes_do_not_consume_failure_retry_effort(self):
        from v2.execution.worker_contract import (
            Manifest,
            Search,
            WorkerSettings,
        )

        settings = WorkerSettings()
        manifest = Manifest(
            35,
            7,
            "repeated-resume",
            dict(ecm_tiers=((10, 30, 1),), backend="python-int"),
        )
        search = Search(manifest, settings)
        for ordinal in range(1, 7):
            self.assertEqual(search.next(), 0)
            search.dispatch(0, ordinal)
            search.receive(
                0,
                ordinal,
                dict(
                    id=manifest.assignment(0)["id"],
                    token=ordinal,
                    attempt=ordinal,
                    reason="cancelled",
                    factor=None,
                    stage="curve",
                    work=1,
                ),
            )
            search = Search.restore(search.snapshot(), manifest, settings)
        self.assertEqual(search.attempts, [6])
        self.assertEqual(search.failed_attempts, [0])

    def test_numerical_failure_retry_effort_still_has_its_original_value(self):
        from v2.execution.worker_contract import (
            LeaseRefusedError,
            Manifest,
            Search,
            WorkerSettings,
        )

        settings = WorkerSettings()
        manifest = Manifest(
            35,
            7,
            "failure-retries",
            dict(ecm_tiers=((10, 30, 1),), backend="python-int"),
        )
        search = Search(manifest, settings)
        for ordinal in range(1, settings.max_attempts + 1):
            search.dispatch(0, ordinal)
            search.receive(
                0,
                ordinal,
                dict(
                    id=manifest.assignment(0)["id"],
                    token=ordinal,
                    attempt=ordinal,
                    reason="error",
                    factor=None,
                    stage="curve",
                    work=1,
                ),
            )
        with self.assertRaisesRegex(LeaseRefusedError, "attempt_limit"):
            search.next()


class DrainPublicationTests(unittest.TestCase):
    def test_drained_result_settles_without_publishing_a_completed_cursor(
        self,
    ):
        from v2.execution.worker_contract import (
            Manifest,
            Search,
            WorkerSettings,
        )

        manifest = Manifest(
            35,
            7,
            "drain",
            dict(ecm_tiers=((10, 30, 1),), backend="python-int"),
        )
        search = Search(manifest, WorkerSettings())
        search.dispatch(0, 1)
        search.receive(
            0,
            1,
            dict(
                id=manifest.assignment(0)["id"],
                token=1,
                attempt=1,
                reason="complete",
                factor=5,
                stage="curve",
                work=1,
            ),
            publish=False,
        )
        self.assertEqual(search.active, {})
        self.assertEqual(search.completed, {})
        self.assertIsNone(search.factor)
        self.assertEqual(search.next(), 0)


class ProcessIdentityRaceTests(unittest.TestCase):
    def test_report_before_first_live_sample_reconciles_with_birth_and_exit(
        self,
    ):
        from types import SimpleNamespace

        provider = Provider()
        budget = Budget(
            provider=provider, wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        process = SimpleNamespace(
            pid=123, _popen=SimpleNamespace(pid=123, returncode=None)
        )
        budget.track_process(process)
        budget.acknowledge(123, 2)
        provider.rows.append(ProcessSample(123, 99, 50, 3))
        budget.publication()
        self.assertEqual(budget.cpu_used, 3)
        usage = SimpleNamespace(ru_utime=4, ru_stime=1)
        with patch("os.wait4", return_value=(123, 0, usage)):
            process._popen.poll()
        self.assertEqual(budget.cpu_used, 5)
        self.assertNotIn(123, budget.roots)

    def test_reused_pid_with_no_live_sample_does_not_reuse_old_cpu_counter(
        self,
    ):
        from types import SimpleNamespace

        provider = Provider()
        budget = Budget(
            provider=provider, wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        for cpu in (3, 1):
            process = SimpleNamespace(
                pid=123, _popen=SimpleNamespace(pid=123, returncode=None)
            )
            budget.track_process(process)
            usage = SimpleNamespace(ru_utime=cpu, ru_stime=0)
            with patch("os.wait4", return_value=(123, 0, usage)):
                process._popen.poll()
        self.assertEqual(budget.cpu_used, 4)


class ExecutorObservationTests(unittest.TestCase):
    def test_process_is_registered_before_manager_visibility(self):
        from types import SimpleNamespace
        from unittest.mock import Mock

        from v2.qs.parallel import ObservedProcessPoolExecutor

        process = object()
        owner = SimpleNamespace(track_process=Mock())
        executor = object.__new__(ObservedProcessPoolExecutor)
        executor._processes = {}
        executor.resource_owner = owner

        def spawn(instance):
            instance._processes[123] = process

        with patch(
            "concurrent.futures.ProcessPoolExecutor._spawn_process", new=spawn
        ):
            executor._spawn_process()
        owner.track_process.assert_called_once_with(process)


class SharedEntryTests(unittest.TestCase):
    def test_nested_entry_keeps_root_and_stage_view(self):
        from v2.execution.allocation import (
            ECMAllocation,
            HandoffRequiredError,
            PretestBudget,
        )
        from v2.execution.resources import accounted_entry, active_owner
        from v2.qs.factor_base import build_factor_base

        owner = Budget(
            provider=Provider(), wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        policy = ECMAllocation("pretest", pretest_work=1)
        view = PretestBudget(owner, policy, fallback=True)

        @accounted_entry
        def inner(*, budget):
            self.assertIs(budget, view)
            self.assertIs(active_owner(), owner)
            with self.assertRaises(HandoffRequiredError):
                budget.consume(2)
            build_factor_base(1009 * 1013, bound=40)

        @accounted_entry
        def outer(*, budget):
            self.assertIs(budget, owner)
            inner(budget=view)

        outer(budget=owner)
        self.assertGreater(owner.used, 1)
        self.assertIsNone(owner.reason)
        self.assertFalse(owner._active)

    def test_nested_entry_does_not_create_an_independent_limit(self):
        from v2.execution.resources import accounted_entry

        owner = Budget(
            provider=Provider(), wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        other = Budget(
            seconds=0,
            provider=Provider(),
            wall_clock=lambda: 0,
            cpu_clock=lambda: 0,
        )

        @accounted_entry
        def inner(*, budget):
            self.assertIs(budget, owner)
            budget.consume(7)

        @accounted_entry
        def outer(*, budget):
            inner(budget=other)

        outer(budget=owner)
        self.assertEqual(owner.used, 7)
        self.assertEqual(other.used, 0)

    def test_qs_context_close_owns_only_its_active_interval(self):
        from types import SimpleNamespace

        from v2.qs.parallel import CollectionPool

        wall, cpu = Clock(), Clock()
        owner = Budget(provider=Provider(), wall_clock=wall, cpu_clock=cpu)
        owner.close()
        wall.value, cpu.value = 100, 50

        def shutdown():
            wall.value += 2
            cpu.value += 1

        pool = object.__new__(CollectionPool)
        pool.mode, pool.executor = "serial", None
        pool.state = {"stop": SimpleNamespace(set=shutdown)}
        pool.close(budget=owner)
        self.assertEqual((owner.wall_used, owner.cpu_used), (2, 1))
        self.assertTrue(pool.closed)


class CandidateBoundaryTests(unittest.TestCase):
    def test_candidate_stores_grow_and_keep_state_on_root_rss_stop(self):
        from v2.qs.checkpoint import _restore_store, _store
        from v2.qs.extended_policy import (
            BoundedYieldSieveConfig,
            ExtendedYieldSieveConfig,
        )
        from v2.qs.factor_base import build_factor_base
        from v2.qs.polynomial import qs_polynomial
        from v2.qs.sieve_collector import (
            CombinedYieldSieveConfig,
            DoubleLargeSieveConfig,
            GrowingDoubleLargeSieveConfig,
            SieveCollector,
            SieveConfig,
            YieldSieveConfig,
        )

        for factory in (
            SieveConfig,
            DoubleLargeSieveConfig,
            GrowingDoubleLargeSieveConfig,
            YieldSieveConfig,
            CombinedYieldSieveConfig,
            ExtendedYieldSieveConfig,
            BoundedYieldSieveConfig,
        ):
            with self.subTest(candidate=factory.__name__):
                clock, provider = Clock(), Provider()
                budget = Budget(
                    provider=provider, wall_clock=clock, cpu_clock=lambda: 0
                )
                base = build_factor_base(
                    1009 * 1013, bound=100, budget=budget
                ).factor_base
                polynomial = qs_polynomial(base)
                config = (
                    factory()
                    if factory is SieveConfig
                    else factory(
                        large_prime_bound=1000, large_product_bound=1000000
                    )
                )
                collector = SieveCollector(
                    polynomial, base, config=config, budget=budget
                )
                collector._workspace = 10**30
                before = _store(collector)
                provider.rows[0] = ProcessSample(os.getpid(), 1, 900, 0)
                clock.value = 1
                stopped = collector.collect(-2, 3)
                self.assertEqual(stopped.reason, "rss_limit")
                self.assertEqual(stopped.resource_stop["observed"], 900)
                self.assertEqual(_store(collector), before)
                provider.rows[0] = ProcessSample(os.getpid(), 1, 100, 0)
                clock.value = 2
                result = collector.collect(stopped.next_position, 3)
                self.assertIn(result.reason, ("complete", "factor_found"))
                published = json.loads(json.dumps(_store(collector)))
                restored = SieveCollector(
                    polynomial, base, config=config, budget=budget
                )
                restored._workspace = 10**30
                _restore_store(published, restored, budget)
                self.assertEqual(
                    json.loads(json.dumps(_store(restored))), published
                )

    def test_encoding_allocation_and_unsupported_values_keep_old_save(self):
        from v2.execution.resources import (
            AllocationFailureError,
            UnsupportedEncodingError,
        )

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.json"
            atomic_save(path, {"published": "old"})
            with patch(
                "v2.execution.resources.json.dumps",
                side_effect=MemoryError("encoder allocation"),
            ):
                with self.assertRaises(AllocationFailureError) as failed:
                    atomic_save(path, {"published": "new"})
            self.assertIsInstance(failed.exception.__cause__, MemoryError)
            self.assertEqual(failed.exception.owner, "checkpoint_save")
            with self.assertRaises(UnsupportedEncodingError):
                atomic_save(path, {"unsupported": object()})
            self.assertEqual(
                json.loads(path.read_text()), {"published": "old"}
            )

    def test_allocation_failure_keeps_owner_through_process_serialization(
        self,
    ):
        import pickle

        from v2.execution.resources import AllocationFailureError

        original = AllocationFailureError("matrix elimination", "QSJob")
        decoded = pickle.loads(pickle.dumps(original))
        self.assertEqual(decoded.payload, original.payload)
        self.assertEqual(str(decoded), str(original))


class DeadlineContinuationTests(unittest.TestCase):
    def test_stage_wall_share_hands_remaining_root_to_qs(self):
        from v2 import portfolio
        from v2.execution.allocation import ECMAllocation
        from v2.qs.siqs import SIQSConfig

        clock = Clock()
        budget = Budget(
            seconds=100,
            provider=Provider(),
            wall_clock=clock,
            cpu_clock=lambda: 0,
        )
        config = portfolio.PortfolioConfig(
            trial_bound=2,
            rho_attempts=0,
            pm1_attempts=0,
            fermat_steps=0,
            ecm_tiers=((2, 2, 1),),
            allocation=ECMAllocation(
                "pretest", pretest_work=10**20, pretest_seconds=35
            ),
            siqs=SIQSConfig(
                mode="qs", base_bound=100, half_width=128, max_half_width=128
            ),
        )
        advance = portfolio.advance_job

        def expire_share(*args, **kwargs):
            clock.value = 36
            return advance(*args, **kwargs)

        with patch("v2.portfolio.advance_job", side_effect=expire_share):
            run = portfolio.factorize_bounded(
                1009 * 1013, config=config, budget=budget, ecm_executor=None
            )
        self.assertTrue(run.result.complete)
        self.assertEqual(run.result.reconstruct(), 1009 * 1013)
        self.assertIsNone(run.resource_stop)
        self.assertTrue(
            any(event.get("outcome") == "pretest_wall" for event in run.events)
        )
        self.assertEqual(budget.wall_used, 36)

    def test_deadline_checkpoint_resumes_fresh_and_retained_history(self):
        from v2.qs.sieve_collector import SieveCollector, SieveConfig
        from v2.qs.siqs import SIQSConfig, SIQSJob

        clock = Clock()
        budget = Budget(
            seconds=3,
            provider=Provider(),
            wall_clock=clock,
            cpu_clock=lambda: 0,
        )
        config = SIQSConfig(
            mode="qs",
            base_bound=100,
            half_width=16,
            max_half_width=16,
            batch_width=17,
            collector=SieveConfig(block_width=17),
        )
        job = SIQSJob(1009 * 1013, config=config, budget=budget)
        collect = SieveCollector._collect_progress

        def expire_after_publication(instance, *args):
            result = collect(instance, *args)
            clock.value = 4
            return result

        with patch.object(
            SieveCollector, "_collect_progress", new=expire_after_publication
        ):
            stopped = job.run()
        self.assertEqual(stopped.reason, "wall_limit")
        self.assertEqual(stopped.resource_stop["observed"], 4)
        self.assertEqual(job.engine.next_position, 1)
        saved = json.loads(json.dumps(job.checkpoint()))
        clock.value = 104
        budget.seconds = 10
        fresh = Budget(
            seconds=10,
            provider=Provider(),
            wall_clock=clock,
            cpu_clock=lambda: 0,
        )
        for owner in (fresh, budget):
            restored = SIQSJob.from_checkpoint(saved, budget=owner)
            result = restored.run()
            self.assertEqual(result.reason, "factor_found")
            self.assertEqual(result.divisor * result.cofactor, job.n)
            self.assertEqual(owner.wall_used, 4)


class FinalCpuTests(unittest.TestCase):
    def test_final_self_cpu_does_not_add_descendant_usage_twice(self):
        from types import SimpleNamespace

        provider = Provider()
        owner = Budget(
            provider=provider, wall_clock=lambda: 0, cpu_clock=lambda: 0
        )
        provider.rows.extend(
            (ProcessSample(123, 10, 100, 3), ProcessSample(456, 11, 100, 7))
        )
        process = SimpleNamespace(
            pid=123, _popen=SimpleNamespace(pid=123, returncode=None)
        )
        owner.track_process(process)
        provider.peek_exit = lambda pid, flag: ProcessSample(pid, 10, 0, 5)
        rolled = SimpleNamespace(ru_utime=10, ru_stime=2)
        with patch("os.wait4", return_value=(123, 0, rolled)):
            self.assertEqual(process._popen.poll(), 0)
        self.assertEqual(owner.cpu_used, 12)
