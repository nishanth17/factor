"""Independent oracles for bounded immutable ECM metadata reuse."""

import itertools
import json
import multiprocessing
import os
import resource
import subprocess
import sys
import threading
import time
import unittest
from dataclasses import replace
from math import isqrt
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.ecm_worker import CurveEngine
from v2.execution.ecm_worker_state import DecodedPrograms
from v2.execution.ecm_workers import ECMExecutor, arithmetic_options
from v2.execution.worker_contract import Manifest, WorkerSettings
from v2.tests.test_ecm_workers import allowance, configuration, task


class PrimeSource:
    segment_size = 8
    base_primes = ()

    def __init__(self):
        self.calls = 0

    def prime_segment(self, lo, hi):
        self.calls += 1
        return [
            n
            for n in range(max(2, lo), hi)
            if all(n % d for d in range(2, isqrt(n) + 1))
        ]


class ScheduleReuseTests(unittest.TestCase):
    def test_streamed_prime_power_oracle_and_immutable_hits(self):
        source = PrimeSource()
        cache = DecodedPrograms(source, 2**17)
        first = cache.program_segment(2, 11, allowance(), bound=10)
        self.assertEqual(first, (2, 3, 5, 7))
        self.assertEqual(cache.power_values(2, 11, 10, 0, 4), (8, 9, 5, 7))
        self.assertIs(
            cache.program_segment(2, 11, allowance(), bound=10), first
        )
        self.assertEqual(source.calls, 1)
        with self.assertRaises(TypeError):
            first[0] = 13
        with self.assertRaises(ValueError):
            cache.power_values(2, 11, 10, 0, 5)

    def test_a_schedule_longer_than_64_blocks_is_reused_in_full(self):
        source = PrimeSource()
        cache = DecodedPrograms(source, 2**17)
        expected = []
        for left in range(2, 642, 8):
            expected.append(cache.program_segment(left, left + 8, allowance()))
        self.assertEqual(len(cache.blocks), 80)
        for left, values in zip(range(2, 642, 8), expected):
            self.assertIs(
                cache.program_segment(left, left + 8, allowance()), values
            )
        self.assertEqual(source.calls, 80)
        self.assertEqual(cache.hits, 80)
        self.assertLessEqual(cache.used_bytes, cache.memory_bytes)

    def test_insufficient_capacity_keeps_a_reusable_prefix_not_scan_thrash(
        self,
    ):
        source = PrimeSource()
        cache = DecodedPrograms(source, 1024)
        first = cache.program_segment(2, 11, allowance(), bound=10)
        for _ in range(3):
            tail = cache.program_segment(11, 19, allowance())
            self.assertEqual(tail, [11, 13, 17])
            self.assertIs(
                cache.program_segment(2, 11, allowance(), bound=10), first
            )
        self.assertEqual(source.calls, 4)
        self.assertEqual(cache.refusals, 3)
        self.assertEqual(cache.hits, 3)
        self.assertLessEqual(cache.used_bytes, 1024)

    def test_failed_build_does_not_publish_a_partial_entry(self):
        source = PrimeSource()
        cache = DecodedPrograms(source, 2**17)
        with self.assertRaises(BudgetExhaustedError):
            cache.program_segment(2, 11, allowance(work_limit=8), bound=10)
        self.assertEqual(source.calls, 1)
        self.assertFalse(cache.blocks)
        self.assertEqual(cache.used_bytes, 0)
        self.assertEqual(
            cache.program_segment(2, 11, allowance(), bound=10), (2, 3, 5, 7)
        )

    def test_deep_streamed_curves_match_uncached_math_and_reuse_across_moduli(
        self,
    ):
        config = configuration(ecm_tiers=((3000, 9000, 3),))
        settings = WorkerSettings(decoded_cache_bytes=2**20)
        cached, ordinary = CurveEngine(settings), CurveEngine()
        outputs = []
        for n in (2**61 - 1, 2**89 - 1):
            plan = Manifest(
                n, 7919, "root", arithmetic_options(config, settings)
            )
            self.assertEqual(plan.options["ecm_chain_mode"], "off")
            for position in range(3):
                message = task(plan, position)
                left = cached.execute(message, cancelled=lambda: False)
                right = ordinary.execute(message, cancelled=lambda: False)
                fields = (
                    "id",
                    "factor",
                    "reason",
                    "phase",
                    "sigma",
                    "chain_chunks",
                )
                self.assertEqual(
                    [left[k] for k in fields], [right[k] for k in fields]
                )
                outputs.append(left)
        self.assertTrue(any(x["phase"] == "stage_two" for x in outputs))
        self.assertEqual(cached.preparations, 1)
        self.assertGreater(cached.decode_hits, 0)
        self.assertLessEqual(
            cached.programs.used_bytes, settings.decoded_cache_bytes
        )


class ManifestReuseTests(unittest.TestCase):
    def test_reuse_keeps_identity_and_does_not_reuse_search_results(self):
        settings = WorkerSettings(manifest_cache_bytes=2**20)
        config = configuration()
        with ECMExecutor(settings=settings) as owner:
            first, hit = owner._manifest(35, config, "root")
            self.assertFalse(hit)
            assignment = first.assignment(0)
            assignment["seed"] = 0
            again, hit = owner._manifest(35, config, "root")
            self.assertTrue(hit)
            self.assertIs(first, again)
            self.assertNotEqual(again.assignment(0)["seed"], 0)
            for _ in range(2):
                run = owner.run(
                    35, config=config, budget=allowance(), fixed_work=True
                )
                self.assertEqual(run.stats["completed_before"], 0)
                self.assertEqual(run.stats["manifest_cache_hits"], 1)
        self.assertFalse(owner.manifests)
        self.assertEqual(owner.manifest_bytes, 0)

    def test_cache_keys_cover_modulus_lineage_and_all_math_options(self):
        owner = ECMExecutor(
            settings=WorkerSettings(manifest_cache_bytes=2**20)
        )
        config = configuration()
        reference, _ = owner._manifest(35, config, "a")
        for n, cfg, lineage in (
            (77, config, "a"),
            (35, config, "b"),
            (35, replace(config, gcd_batch=32), "a"),
        ):
            changed, hit = owner._manifest(n, cfg, lineage)
            self.assertFalse(hit)
            self.assertNotEqual(reference.identity, changed.identity)
        self.assertIs(owner._manifest(35, config, "a")[0], reference)
        owner.close()

    def test_refusal_eviction_and_ownership_stay_finite(self):
        config = configuration()
        refused = ECMExecutor(
            settings=WorkerSettings(manifest_cache_bytes=4096)
        )
        first, _ = refused._manifest(35, config, "root")
        second, hit = refused._manifest(35, config, "root")
        self.assertFalse(hit)
        self.assertIsNot(first, second)
        self.assertEqual(first.identity, second.identity)
        self.assertFalse(refused.manifests)
        refused.close()
        settings = WorkerSettings(manifest_cache_bytes=200000)
        with ECMExecutor(settings=settings) as owner:
            for n in range(35, 115):
                owner._manifest(n, config, "root")
                self.assertLessEqual(
                    owner.manifest_bytes, settings.manifest_cache_bytes
                )
                self.assertLessEqual(len(owner.manifests), 16)
            self.assertEqual(len(owner.manifests), 2)
        self.assertEqual(
            settings.owned_bytes - WorkerSettings().owned_bytes, 200000
        )

    def test_reused_manifest_does_not_repeat_seed_generation(self):
        config = configuration()
        with ECMExecutor(
            settings=WorkerSettings(manifest_cache_bytes=2**20)
        ) as owner:
            plan, _ = owner._manifest(35, config, "root")
            plan.assignment(0)
            with patch(
                "v2.execution.worker_contract.random.Random",
                side_effect=AssertionError,
            ):
                self.assertEqual(
                    owner._manifest(35, config, "root")[0].assignment(0),
                    plan.assignment(0),
                )


class WorkerKernelAndImportTests(unittest.TestCase):
    def test_reduced_kernels_match_polynomial_oracles_even_at_nonunits(self):
        from v2.ecm import chains, core

        for n in (9, 15, 35, 101, 10**59 + 7, 10**99 + 267):
            values = (0, 1, n - 1)
            for coordinates in itertools.product(values, repeat=6):
                self.assertEqual(
                    chains.point_add(*coordinates, n),
                    core.point_add(*coordinates, n),
                )
            for x, z, a24 in itertools.product(values, repeat=3):
                self.assertEqual(
                    chains.point_double(x, z, n, a24),
                    core.point_double(x, z, n, a24),
                )

    def test_slim_and_reduced_curves_preserve_exact_work_and_results(self):
        config = configuration(ecm_tiers=((20, 2000, 8),))
        legacy = CurveEngine()
        reduced = CurveEngine(
            WorkerSettings(worker_imports="ecm", stage_two_kernel="reduced")
        )
        stages = set()
        for n in (35, 1009 * 1013, 2**61 - 1):
            plan = Manifest(
                n, 7919, "root", arithmetic_options(config, WorkerSettings())
            )
            for position in range(8):
                message = task(plan, position)
                left = legacy.execute(message, cancelled=lambda: False)
                right = reduced.execute(message, cancelled=lambda: False)
                for key in (
                    "id",
                    "factor",
                    "reason",
                    "phase",
                    "sigma",
                    "chain_chunks",
                    "work",
                ):
                    self.assertEqual(left[key], right[key], (n, position, key))
                stages.add(left["phase"])
        self.assertIn("setup", stages)
        self.assertIn("stage_two", stages)

    def test_ecm_config_import_is_independent_of_relation_modules(self):
        if sys.implementation.name != "pypy":
            self.skipTest("PyPy implementing Python3.11 is required")
        code = """import json,sys
from dataclasses import asdict
from v2.execution.portfolio_config import PortfolioConfig
c=PortfolioConfig(ecm_tiers=((2000,200000,128),),siqs=None)
assert not any(x.startswith('v2.qs') for x in sys.modules)
assert PortfolioConfig.__module__=='v2.portfolio'
from v2.portfolio import PortfolioConfig as Public
assert Public is PortfolioConfig
assert asdict(Public(**asdict(c)))==asdict(c)
print(json.dumps({'passed':True}))
"""
        result = subprocess.run(
            [sys.executable, "-B", "-c", code],
            capture_output=True,
            text=True,
            timeout=8,
            env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
            check=True,
        )
        self.assertTrue(json.loads(result.stdout)["passed"])

    def test_new_options_are_bounded_and_do_not_move_curve_streams(self):
        for changes in (
            {"worker_imports": "bad"},
            {"stage_two_kernel": "bad"},
            {"manifest_cache_bytes": 16 * 2**20 + 1},
        ):
            with self.assertRaises(ValueError):
                adapt_resource_configuration(WorkerSettings, **changes)
        config = configuration()
        base = Manifest(
            35, 7919, "root", arithmetic_options(config, WorkerSettings())
        )
        settings = WorkerSettings(
            worker_imports="ecm",
            stage_two_kernel="reduced",
            decoded_cache_bytes=2**20,
            manifest_cache_bytes=2**20,
        )
        changed = Manifest(
            35, 7919, "root", arithmetic_options(config, settings)
        )
        self.assertEqual(base.identity, changed.identity)
        self.assertEqual(
            [base.assignment(i) for i in range(4)],
            [changed.assignment(i) for i in range(4)],
        )


class CommonCurveReplayTests(unittest.TestCase):
    def test_parameter_replay_binds_new_options_but_keeps_the_original_seed(
        self,
    ):
        from v2.benchmarks.ecm.b12_f4.kernel_replay import replay_assignment

        config = configuration()
        plan = Manifest(
            35, 7919, "root", arithmetic_options(config, WorkerSettings())
        )
        curve = plan.assignment(0)
        first = replay_assignment(curve, "gcd-128")
        second = replay_assignment(curve, "gcd-256")
        self.assertNotEqual(first["id"], second["id"])
        self.assertNotEqual(first["id"], curve["id"])
        self.assertEqual(first["seed"], second["seed"])
        self.assertEqual(first["sigma"], curve["sigma"])
        self.assertEqual(first["n"], curve["n"])
        self.assertEqual(plan.assignment(0), curve)


class ReuseProcessTests(unittest.TestCase):
    def test_slim_cached_reduced_pool_matches_serial_and_closes(self):
        from v2.benchmarks.ecm.b12_f4.profile_runner import signature

        before = {p.pid for p in multiprocessing.active_children()}
        config = configuration(ecm_tiers=((3000, 9000, 8),))
        settings = WorkerSettings(
            workers=2,
            worker_imports="ecm",
            stage_two_kernel="reduced",
            decoded_cache_bytes=2**20,
            packed_cache_bytes=2**20,
            manifest_cache_bytes=2**20,
            packet_curves=2,
            poll_every=16,
        )
        with ECMExecutor(
            seed=7919, settings=settings, force_workers=True
        ) as pool:
            with ECMExecutor(seed=7919, force_workers=False) as serial:
                for n in (2**61 - 1, 2**61 - 1, 2**89 - 1, 2**61 - 1):
                    actual = pool.run(
                        n, config=config, budget=allowance(), fixed_work=True
                    )
                    expected = serial.run(
                        n, config=config, budget=allowance(), fixed_work=True
                    )
                    self.assertEqual(
                        signature(actual.checkpoint, n, 8),
                        signature(expected.checkpoint, n, 8),
                    )
                self.assertEqual(pool.lifetime_stats["preparations"], 2)
                self.assertEqual(actual.stats["manifest_cache_hits"], 1)
                self.assertGreater(pool.lifetime_stats["decode_hits"], 0)
        self.assertEqual(
            {p.pid for p in multiprocessing.active_children()}, before
        )
        self.assertFalse(pool.manifests)


class ConfigurationCompatibilityTests(unittest.TestCase):
    def test_lazy_config_decoder_preserves_reference_relation_configuration(
        self,
    ):
        from dataclasses import asdict

        from v2.portfolio import PortfolioConfig
        from v2.qs import SIQSConfig
        from v2.tests.resource_stops import (
            historical_decode_siqs as _decode_siqs_config,
        )

        config = SIQSConfig()
        self.assertEqual(_decode_siqs_config(asdict(config)), config)
        self.assertEqual(
            PortfolioConfig(siqs=asdict(config)).siqs,
            config,
        )


class PackedBackingTests(unittest.TestCase):
    def test_deep_program_backing_keeps_stream_and_reuses_blocks(self):
        config = configuration(ecm_tiers=((3000, 9000, 3),))
        settings = WorkerSettings(
            packed_cache_bytes=2**20, worker_imports="ecm"
        )
        cached, ordinary = CurveEngine(settings), CurveEngine()
        plan = Manifest(
            2**89 - 1, 7919, "root", arithmetic_options(config, settings)
        )
        self.assertEqual(plan.options["ecm_program_bytes"], 0)
        for position in range(3):
            request = task(plan, position)
            left = cached.execute(request, cancelled=lambda: False)
            right = ordinary.execute(request, cancelled=lambda: False)
            for key in (
                "id",
                "factor",
                "reason",
                "phase",
                "sigma",
                "chain_chunks",
            ):
                self.assertEqual(left[key], right[key])
        self.assertGreater(cached.programs.hits, 0)
        self.assertLessEqual(
            cached.programs.used_bytes, settings.packed_cache_bytes
        )
        self.assertIsNone(cached.programs.chains)

    def test_existing_accepted_program_store_gets_additive_retention(self):
        config = configuration(ecm_program_bytes=65536)
        settings = WorkerSettings(packed_cache_bytes=2**20)
        engine = CurveEngine(settings)
        plan = Manifest(
            2**61 - 1, 7919, "root", arithmetic_options(config, settings)
        )
        result = engine.execute(task(plan, 0), cancelled=lambda: False)
        self.assertEqual(result["reason"], "complete")
        self.assertEqual(engine.programs.memory_bytes, 65536 + 2**20)
        self.assertEqual(plan.options["ecm_program_bytes"], 65536)
        self.assertEqual(engine.preparations, 1)

    def test_unfunded_scratch_refuses_before_dispatch(
        self,
    ):
        settings = WorkerSettings(packed_cache_bytes=4096)
        with ECMExecutor(settings=settings, force_workers=True) as owner:
            result = owner.run(35, config=configuration(), budget=allowance())
            self.assertEqual(result.reason, "curves_exhausted")
            self.assertTrue(owner.slots)
        base = WorkerSettings()
        extra = WorkerSettings(packed_cache_bytes=2**20)
        self.assertEqual(extra.owned_bytes - base.owned_bytes, 2 * 2**20)


def _burn_cpu(seconds):
    start = time.process_time()
    while time.process_time() - start < seconds:
        pass


def _receive(connection):
    if not connection.poll(10):
        raise TimeoutError("CPU accounting test phase timed out")
    return connection.recv()


def _cpu_leaf(channel):
    burner = threading.Thread(target=_burn_cpu, args=(0.12,))
    burner.start()
    burner.join()
    channel.send(time.process_time())
    if _receive(channel) != "exit":
        raise ValueError("expected leaf exit")
    _burn_cpu(0.04)


def _cpu_worker(channel):
    with channel:
        _burn_cpu(0.15)
        channel.send(time.process_time())
        if _receive(channel) != "second":
            raise ValueError("expected second call")
        ctx = multiprocessing.get_context("spawn")
        worker_start = time.process_time()
        monitor, peer = ctx.Pipe()
        leaf = ctx.Process(target=_cpu_leaf, args=(peer,))
        with monitor, peer:
            leaf.start()
            peer.close()
            try:
                leaf_live = _receive(monitor)
                _burn_cpu(0.06)
                channel.send((leaf.pid, leaf_live, time.process_time()))
                if _receive(channel) != "reap":
                    raise ValueError("expected leaf reap")
                monitor.send("exit")
                leaf.join(10)
                if leaf.exitcode != 0:
                    raise RuntimeError("leaf failed to exit")
            finally:
                if leaf.is_alive():
                    leaf.kill()
                    leaf.join(10)
                leaf.close()
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        channel.send((worker_start, usage.ru_utime + usage.ru_stime))
        if _receive(channel) != "exit":
            raise ValueError("expected worker exit")
        _burn_cpu(0.06)
        channel.send(time.process_time())


@unittest.skipUnless(sys.platform == "darwin", "requires macOS libproc")
class RetainedWorkerCpuTests(unittest.TestCase):
    def test_second_call_charges_live_and_reaped_descendant_once(self):
        ctx = multiprocessing.get_context("spawn")
        monitor, peer = ctx.Pipe()
        worker = ctx.Process(target=_cpu_worker, args=(peer,))
        # Isolate child charges from the observing test process's CPU.
        first = Budget(rss_bytes=None, cpu_clock=lambda: 0)
        with monitor, peer:
            worker.start()
            peer.close()
            try:
                first.track_process(worker)
                first_report = _receive(monitor)
                first.acknowledge(worker.pid, first_report)
                first.close()

                second = Budget(rss_bytes=None, cpu_clock=lambda: 0)
                second.track_process(worker)
                second.acknowledge(worker.pid, first_report)
                # Two clocks bracket pipe handoffs. Allow 30 ms for their CPU
                # and bookkeeping, well below either burn or a double charge.
                tolerance = 0.03
                self.assertAlmostEqual(second.cpu_used, 0, delta=tolerance)
                monitor.send("second")
                leaf_pid, leaf_live, worker_live = _receive(monitor)
                second.acknowledge(worker.pid, worker_live)
                self.assertIn(leaf_pid, second._live_identities)
                self.assertAlmostEqual(
                    second.cpu_used,
                    worker_live - first_report + leaf_live,
                    delta=tolerance,
                )

                monitor.send("reap")
                worker_start, leaf_final = _receive(monitor)
                second.publication()
                self.assertNotIn(leaf_pid, second._live_identities)
                self.assertAlmostEqual(
                    second.cpu_used,
                    worker_live - worker_start + leaf_final,
                    delta=tolerance,
                )
                monitor.send("exit")
                worker_final = _receive(monitor)
                worker.join(10)
                self.assertEqual(worker.exitcode, 0)
                second.close()
                self.assertAlmostEqual(
                    second.cpu_used,
                    worker_final - worker_start + leaf_final,
                    delta=tolerance,
                )
            finally:
                # Closing the pipe also releases a child stalled at a phase.
                monitor.close()
                worker.join(10)
                if worker.is_alive():
                    worker.kill()
                    worker.join(10)
                worker.close()
