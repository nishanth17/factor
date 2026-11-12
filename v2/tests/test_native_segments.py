"""Coarse native exactness, observed limits, fallback and owned lifetimes."""

import json
import os
import subprocess
import sys
import unittest
from contextlib import contextmanager
from dataclasses import replace
from functools import wraps
from pathlib import Path

from v2.ecm.native_curve_team import NativeCurveTeam, OpenMPCurveTeam
from v2.ecm.native_packet import SOURCE
from v2.ecm.native_resources import native_bulk_allowed
from v2.ecm.native_segment import NativeSegment
from v2.execution import stage_jobs
from v2.execution.allocation import ECMAllocation, PretestBudget
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.ecm_native import _BudgetView
from v2.execution.portfolio_config import (
    PortfolioConfig,
    _execution_context,
    close_execution_context,
)
from v2.tests.test_h3_resources import Provider
from v2.tests.test_native_correctness import (
    FIXTURES,
    deterministic_budget,
    selected,
)


def isolated(function):
    """Avoid inheriting another test's still-live spawn resource tracker."""

    @wraps(function)
    def run(self):
        # Discovery imports this file as test_native_segments. The isolated
        # interpreter starts at the repository root and needs the package path.
        case_id = (
            f"v2.tests.test_native_segments.{type(self).__qualname__}."
            f"{self._testMethodName}"
        )
        if os.environ.get("FACTOR_NATIVE_SEGMENT_CHILD") != case_id:
            env = dict(
                os.environ,
                FACTOR_NATIVE_SEGMENT_CHILD=case_id,
                OMP_WAIT_POLICY="PASSIVE",
                KMP_BLOCKTIME="0",
                OMP_DYNAMIC="FALSE",
            )
            subprocess.run(
                [sys.executable, "-B", "-m", "unittest", "-v", case_id],
                env=env,
                cwd=Path(__file__).resolve().parents[2],
                check=True,
                timeout=30,
            )
            return
        function(self)

    return run


@contextmanager
def packet_scope(*, callback=None):
    corpus = FIXTURES.parents[1] / "corpora/b4a_production_3p_v1.json"
    case = next(
        row
        for row in json.loads(corpus.read_text())["fixtures"]
        if row["population"] == "failed" and row["digits"] == 60
    )
    budget = Budget(
        seconds=10, cpu_seconds=20, rss_bytes=8 * 2**30, cancelled=callback
    )
    config = replace(
        selected("c-gmp-full", 2000, 81000),
        ecm_chain_mode="auto",
        ecm_tiers=((2000, 81000, 9),),
    )
    context, programs = _execution_context({}, config, budget)
    source = programs or context
    packet = source.native_packet
    job = stage_jobs.new_job("ecm", case["n"], 73, 2000, 81000)
    try:
        yield packet, job, budget, source, config
    finally:
        close_execution_context(source)
        budget.close()


@unittest.skipUnless(
    sys.platform == "darwin" and (SOURCE / "build/manifest.json").exists(),
    "native macOS build required",
)
class NativeSegmentTests(unittest.TestCase):
    @isolated
    def test_interruption_before_or_after_foreign_entry_settles_exact_work(
        self,
    ):
        with packet_scope() as (packet, job, budget, source, config):
            stage_jobs.advance_job(job, budget, source, config)
            original = packet.lib

            class InterruptedLibrary:
                def __init__(self, execute):
                    self.execute = execute

                def __getattr__(self, name):
                    if name != "b4_segment":
                        return getattr(original, name)

                    def stop(*args):
                        if self.execute:
                            original.b4_segment(*args)
                        raise BudgetExhaustedError("cancelled")

                    return stop

            try:
                charged, used = packet.charged, budget.used
                packet.lib = InterruptedLibrary(False)
                with self.assertRaises(BudgetExhaustedError):
                    stage_jobs.advance_job(job, budget, source, config)
                self.assertEqual(packet.charged, charged)
                self.assertEqual(budget.used, used)
                self.assertEqual(packet._segment.returns.get("not_entered"), 1)

                packet.lib = InterruptedLibrary(True)
                with self.assertRaises(BudgetExhaustedError):
                    stage_jobs.advance_job(job, budget, source, config)
                self.assertEqual(packet.charged, int(packet._read().work))
                self.assertEqual(
                    packet.blocks.charged_metadata,
                    int(packet.lib.b4_index_stat(packet.blocks.index, 2)),
                )
                self.assertGreater(budget.used, used)
            finally:
                packet.lib = original

    @isolated
    def test_plan_hint_stays_lazy_and_rejects_drift(self):
        with packet_scope() as (packet, job, budget, source, config):
            packet.bind_job(job, budget, source, config)
            index = packet.blocks.ensure_index(budget)
            self.assertEqual(packet.chain_hint(budget), 2)
            packet.attach_index(index, budget)
            self.assertIsNone(packet.blocks.metadata)
            metadata = packet.blocks.native_metadata()
            metadata._catalogue(True, budget)
            self.assertFalse(metadata.plans)
            self.assertEqual(packet.chain_hint(budget), 2)
            packet.attach_index(index, budget)
            self.assertIsNone(packet.plan)
            packet.plan = packet._plan(job, source, budget)
            self.assertEqual(packet.chain_hint(budget), 1)
            self.assertEqual(
                packet.lib.b4_curve_attach_index_hint(
                    packet.context, index, 1
                ),
                0,
            )
            self.assertEqual(
                packet.lib.b4_curve_attach_index_hint(
                    packet.context, index, 0
                ),
                -1,
            )
            self.assertEqual(
                packet.lib.b4_curve_attach_index_hint(
                    packet.context, index, 2
                ),
                -1,
            )

    @isolated
    def test_native_owned_blocks_are_protected_until_context_drain(self):
        with packet_scope() as (packet, job, budget, source, config):
            for _ in range(10000):
                stage_jobs.advance_job(job, budget, source, config)
                if (
                    packet.blocks.index is not None
                    and packet.lib.b4_index_stat(packet.blocks.index, 3)
                ):
                    break
            else:
                self.fail("no native-owned metadata was generated")
            index = packet.blocks.index
            self.assertEqual(packet.lib.b4_index_destroy(index), -1)
            self.assertGreater(packet.blocks.charged_metadata, 0)
            self.assertGreater(packet._segment.observations, 0)
            packet.sync(allow_result=job["done"])
            self.assertGreater(len(job["cursor"]["values"]), 0)
            packet.discard(sync=False)
            packet.blocks.close()
            self.assertEqual(packet.blocks.used, 0)
            self.assertIsNone(packet.blocks.index)

    @isolated
    def test_actual_native_rss_observation_refuses_before_arithmetic(self):
        with packet_scope() as (packet, job, budget, source, config):
            packet.bind_job(job, budget, source, config)
            segment = NativeSegment(packet)
            try:
                self.assertTrue(segment.prepare(budget))
                budget.rss_bytes = segment.policy.rss = 1
                status = packet.lib.b4_segment(
                    packet.context, 128, segment.policy, segment.info
                )
                self.assertEqual(status, 0)
                self.assertEqual(segment.info.reason, 5)
                self.assertGreater(segment.info.observations, 0)
                self.assertGreater(segment.info.rss, 1)
                self.assertEqual(segment.info.work[0], 0)
                segment.settle([packet], budget)
                # Simulate RSS falling before Python publication. The original
                # observed refusal must survive a later zero-RSS sample.
                provider = Provider()
                provider.rows = [
                    replace(
                        provider.rows[0],
                        birth=budget.roots[os.getpid()],
                        rss=0,
                    )
                ]
                budget.provider = provider
                with self.assertRaises(BudgetExhaustedError) as stopped:
                    segment.check(budget, packet)
                self.assertEqual(str(stopped.exception), "rss_limit")
                self.assertEqual(
                    budget.stop["observed"], int(segment.info.rss)
                )
            finally:
                budget.rss_bytes = 8 * 2**30
                segment.close()

    @isolated
    def test_real_child_uses_full_tree_and_original_slice_then_reuses(self):
        with packet_scope() as (packet, job, budget, source, config):
            child = subprocess.Popen(
                [
                    sys.executable,
                    "-B",
                    "-c",
                    "import sys;print('ready',flush=True);sys.stdin.read()",
                ],
                stdout=subprocess.PIPE,
                stdin=subprocess.PIPE,
                text=True,
            )
            try:
                self.assertEqual(child.stdout.readline().strip(), "ready")
                stage_jobs.advance_job(job, budget, source, config)
                self.assertGreater(packet.return_reasons.get("children", 0), 0)
                self.assertGreater(packet.calls, 1)
            finally:
                child.stdin.close()
                try:
                    child.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    child.terminate()
                    try:
                        child.wait(timeout=1)
                    except subprocess.TimeoutExpired:
                        child.kill()
                        child.wait(timeout=1)
                child.stdout.close()
            for _ in range(10000):
                if job["done"]:
                    break
                stage_jobs.advance_job(job, budget, source, config)
            self.assertTrue(job["done"])
            self.assertGreater(packet._segment.waves, 0)

    @isolated
    def test_custom_clocks_providers_and_callbacks_keep_exact_fallback(self):
        budget = deterministic_budget()
        try:
            self.assertFalse(native_bulk_allowed(budget))
        finally:
            budget.close()
        budget = Budget(seconds=10, cpu_seconds=20, rss_bytes=8 * 2**30)
        try:
            reservation = PretestBudget(
                budget,
                ECMAllocation("pretest", pretest_work=0, pretest_seconds=1),
                fallback=False,
            )
            self.assertTrue(native_bulk_allowed(_BudgetView(budget)))
            self.assertFalse(native_bulk_allowed(reservation))
            self.assertFalse(native_bulk_allowed(_BudgetView(reservation)))
        finally:
            budget.close()
        count = [0]

        def callback():
            count[0] += 1
            return False

        with packet_scope(callback=callback) as (
            packet,
            job,
            budget,
            source,
            config,
        ):
            for _ in range(10000):
                stage_jobs.advance_job(job, budget, source, config)
                if job["done"]:
                    break
            self.assertTrue(job["done"])
            self.assertGreater(count[0], 0)
            self.assertIsNone(packet.blocks.index)
            self.assertIsNone(packet.blocks.metadata)

    @isolated
    def test_team_fixed_work_canonical_jobs_and_counters_match_legacy(self):
        corpus = FIXTURES.parents[1] / "corpora/b4a_production_3p_v1.json"
        n = next(
            row["n"]
            for row in json.loads(corpus.read_text())["fixtures"]
            if row["population"] == "failed" and row["digits"] == 60
        )

        def run(owner_type, coarse, b1, b2):
            budget = Budget(seconds=10, cpu_seconds=30, rss_bytes=8 * 2**30)
            config = PortfolioConfig(
                ecm_kernel="c-gmp-full",
                ecm_chain_mode="auto",
                ecm_tiers=((b1, b2, 9),),
                siqs=None,
            )
            context, programs = _execution_context({}, config, budget)
            source = programs or context
            team = owner_type(workers=2, capacity=9)
            team.coarse_segments = coarse
            active = {
                i: stage_jobs.new_job("ecm", n, 73 + i, b1, b2)
                for i in range(9)
            }
            rows = {}
            try:
                for _ in range(10000):
                    ready = team.advance_jobs(
                        list(active.items()), budget, source, config
                    )
                    for position, job in ready:
                        packet = team.bindings[position]
                        packet.sync(allow_result=True)
                        rows[position] = (
                            job,
                            dict(packet.statistics),
                            int(packet._read().work),
                        )
                        active.pop(position)
                    if not active:
                        break
                self.assertFalse(active)
                if coarse:
                    self.assertGreater(team._segment.observations, 0)
                return rows, team.foreign_returns
            finally:
                team.close()
                close_execution_context(source)
                budget.close()

        for owner_type in (NativeCurveTeam, OpenMPCurveTeam):
            for b1, b2 in ((127, 4000), (2000, 81000), (11000, 1100000)):
                with self.subTest(backend=owner_type.variant, b1=b1):
                    expected, old_calls = run(owner_type, False, b1, b2)
                    actual, new_calls = run(owner_type, True, b1, b2)
                    self.assertEqual(actual, expected)
                    self.assertLess(new_calls, old_calls)

    @isolated
    def test_scalar_primary_anchors_match_canonical_legacy_jobs(self):
        corpus = FIXTURES.parents[1] / "corpora/b4a_production_3p_v1.json"
        n = next(
            row["n"]
            for row in json.loads(corpus.read_text())["fixtures"]
            if row["population"] == "failed" and row["digits"] == 80
        )

        def run(coarse, b1, b2):
            budget = Budget(seconds=10, cpu_seconds=20, rss_bytes=8 * 2**30)
            config = PortfolioConfig(
                ecm_kernel="c-gmp-full",
                ecm_chain_mode="auto",
                ecm_tiers=((b1, b2, 9),),
                siqs=None,
            )
            context, programs = _execution_context({}, config, budget)
            source = programs or context
            packet = source.native_packet
            packet.coarse_segments = coarse
            rows = []
            try:
                for seed in range(73, 76):
                    job = stage_jobs.new_job("ecm", n, seed, b1, b2)
                    for _ in range(100000):
                        stage_jobs.advance_job(job, budget, source, config)
                        if job["done"]:
                            break
                    self.assertTrue(job["done"])
                    packet.sync(allow_result=True)
                    rows.append(
                        (
                            job,
                            dict(packet.statistics),
                            int(packet._read().work),
                        )
                    )
                if coarse:
                    self.assertGreater(packet._segment.observations, 0)
                return rows, packet.calls
            finally:
                close_execution_context(source)
                budget.close()

        for b1, b2 in ((2000, 81000), (11000, 1100000)):
            with self.subTest(b1=b1):
                expected, old_calls = run(False, b1, b2)
                actual, new_calls = run(True, b1, b2)
                self.assertEqual(actual, expected)
                if b1 == 2000:
                    self.assertTrue(
                        any(row[0].get("chain_chunks", 0) for row in actual)
                    )
                self.assertLess(new_calls, old_calls)

    @isolated
    def test_production_quotas_retain_plan_references_and_metadata(self):
        from v2.benchmarks.ecm.b4a.production_batch import ProductionSession

        corpus = FIXTURES.parents[1] / "corpora/b4a_production_3p_v1.json"
        n = next(
            row["n"]
            for row in json.loads(corpus.read_text())["fixtures"]
            if row["population"] == "failed" and row["digits"] == 80
        )
        for role in ("scalar", "pthread", "openmp"):
            for b1, b2, curves in ((2000, 81000, 9), (11000, 1100000, 37)):
                with self.subTest(role=role, b1=b1):
                    request = dict(
                        role=role,
                        n=n,
                        b1=b1,
                        b2=b2,
                        curves=curves,
                        seed=73,
                        workers=1
                        if role == "scalar"
                        else 4
                        if b1 == 2000
                        else 8,
                        capacity=32,
                        schedule="dynamic",
                        chunk_size=1,
                        worker_settings={"allow_wide_workers": True},
                        seconds=20,
                        cpu_seconds=120,
                    )
                    session = ProductionSession()
                    try:
                        first = session.call(request)
                        second = session.call(request)
                        for row in (first, second):
                            self.assertIsNone(row["resource_stop"])
                            self.assertEqual(row["completed_curves"], curves)
                            self.assertEqual(row["reason"], "exhausted")
                        self.assertEqual(
                            first["manifests"], second["manifests"]
                        )
                        self.assertEqual(
                            second["state"]["return_reasons"].get(
                                "metadata", 0
                            )
                            - second["before_state"]["return_reasons"].get(
                                "metadata", 0
                            ),
                            0,
                        )
                    finally:
                        session.close()
