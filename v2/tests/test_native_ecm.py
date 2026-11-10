"""Tiny exact oracles, refusal and optional native ECM integration."""

import copy
import json
import multiprocessing
import os
import signal
import threading
import time
import unittest
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import RecordingBudget
from v2.benchmarks.support.prac_oracle import (
    affine_multiply,
    matches,
    twist_point,
)
from v2.common.arithmetic import BackendUnavailableError
from v2.ecm import core
from v2.ecm.chains import ChainPlan
from v2.ecm.native_packet import (
    ABI,
    RESERVE,
    SOURCE,
    NativePacket,
    source_identity,
)
from v2.execution.budget import BudgetExhaustedError
from v2.execution.ecm_worker import CurveEngine
from v2.execution.ecm_worker_state import ActionBudget, ProspectiveActionBudget
from v2.execution.ecm_workers import arithmetic_options
from v2.execution.portfolio_config import close_execution_context
from v2.execution.resources import RemovedResourceLimitError
from v2.execution.stage_jobs import advance_job, new_job
from v2.execution.worker_contract import (
    Manifest,
    Search,
    WorkerSettings,
    digest,
)
from v2.portfolio import PortfolioConfig, _execution_context, factorize_bounded


def allowance(work=10**8, **changes):
    return RecordingBudget(
        work_limit=work, seconds=10, cpu_seconds=10, **changes
    )


def config(**changes):
    return replace(
        PortfolioConfig(
            trial_bound=5,
            rho_attempts=0,
            pm1_attempts=0,
            siqs=None,
            ecm_tiers=((20, 200, 2),),
            segment_size=16,
            max_input_bits=512,
        ),
        **changes,
    )


def spawned_smoke(channel):
    """Import/load is explicit in the spawned interpreter too."""
    packet = NativePacket()
    try:
        packet.available(1013, allowance())
        result = packet.ladder(13, (2, 1), 1013, 2, allowance())
        channel.send((result, packet.build["abi"]))
    finally:
        packet.close()
        channel.close()


class OptionalPacketTests(unittest.TestCase):
    def test_optional_kernel_preserves_automatic_routing_and_prac(self):
        for digits in (30, 39, 40, 41, 50, 60, 79, 80, 81, 100):
            n = 10 ** (digits - 1) + 1
            python = PortfolioConfig().for_input(n)
            native = PortfolioConfig(ecm_kernel="c-gmp").for_input(n)
            self.assertEqual(native.ecm_tiers, python.ecm_tiers)
            self.assertEqual(native.ecm_chain_mode, python.ecm_chain_mode)
            self.assertEqual(native.siqs, python.siqs)
            self.assertEqual(
                native.workspace_reserve, python.workspace_reserve + RESERVE
            )

    def test_retired_memory_caps_rejected_before_loading(self):
        with self.assertRaises(RemovedResourceLimitError):
            NativePacket(memory_bytes=RESERVE - 1)
        with self.assertRaises(RemovedResourceLimitError):
            config(memory_bytes=2**20, ecm_kernel="c-gmp")
        with self.assertRaises(ValueError):
            config(ecm_kernel="unknown")

    def test_lazy_missing_build_explicit_error_and_auto_fallback(self):
        with patch("v2.ecm.native_packet.SOURCE", SOURCE / "absent"):
            for mode in ("auto", "c-gmp"):
                packet = NativePacket(mode=mode)
                self.assertIsNone(packet.lib)
                try:
                    if mode == "auto":
                        self.assertFalse(packet.available(1013, allowance()))
                    else:
                        with self.assertRaises(BackendUnavailableError):
                            packet.available(1013, allowance())
                finally:
                    packet.close()

    def test_ineligible_even_and_wide_moduli(self):
        for n in (2, 100, 2**1024 + 1):
            packet = NativePacket(mode="auto")
            try:
                self.assertFalse(packet.available(n, allowance()))
                self.assertIsNone(packet.lib)
            finally:
                packet.close()

    def test_cheap_portfolio_does_not_load_library_and_keeps_resume_identity(
        self,
    ):
        selected = config(ecm_kernel="c-gmp")
        with patch.object(
            NativePacket, "_load", side_effect=AssertionError("eager load")
        ):
            run = factorize_bounded(2**8, config=selected, budget=allowance())
        self.assertEqual(run.result.reconstruct(), 2**8)
        self.assertTrue(run.result.complete)
        self.assertEqual(
            run.checkpoint["payload"]["native_packet"], source_identity()
        )
        resumed = factorize_bounded(
            2**8,
            config=selected,
            budget=allowance(),
            checkpoint=json.loads(json.dumps(run.checkpoint)),
        )
        self.assertGreaterEqual(resumed.work_used, run.work_used)
        with self.assertRaises(ValueError):
            factorize_bounded(
                2**8,
                config=config(),
                budget=allowance(),
                checkpoint=run.checkpoint,
            )
        broken = copy.deepcopy(run.checkpoint)
        broken["payload"]["native_packet"] = "wrong"
        broken["sha256"] = digest(broken["payload"])
        with self.assertRaises(ValueError):
            factorize_bounded(
                2**8, config=selected, budget=allowance(), checkpoint=broken
            )

    def test_kernel_preserves_b12_curve_stream_and_binds_resume(self):
        settings = WorkerSettings()
        python = Manifest(
            1009 * 1013, 7919, "root", arithmetic_options(config(), settings)
        )
        native = Manifest(
            1009 * 1013,
            7919,
            "root",
            arithmetic_options(config(ecm_kernel="c-gmp"), settings),
        )
        self.assertEqual(
            [python.assignment(i) for i in range(2)],
            [native.assignment(i) for i in range(2)],
        )
        snapshot = Search(native, settings, fixed_work=True).snapshot()
        with self.assertRaises(ValueError):
            Search.restore(snapshot, python, settings, fixed_work=True)
        with self.assertRaises(ValueError):
            Search.restore(
                snapshot, native, replace(settings, workers=2), fixed_work=True
            )


@unittest.skipUnless(
    (SOURCE / "build/manifest.json").exists(), "optional native build absent"
)
class CompiledPacketTests(unittest.TestCase):
    def setUp(self):
        self.packet = NativePacket()
        self.packet.available(1013, allowance())

    def tearDown(self):
        self.packet.close()

    @unittest.skipUnless(
        os.environ.get("FACTOR_B4A_LIFECYCLE_QA") == "1",
        "future Stage Q actual PyPy runtime lane",
    )
    def test_actual_cffi_foreign_call_releases_the_gil(self):
        import cffi

        # The prepared Factor ABI is immutable; libc needs its own binding.
        ffi = cffi.FFI()
        ffi.cdef("int usleep(unsigned int);")
        libc = ffi.dlopen(None)
        ready, start = threading.Event(), threading.Event()
        observed = []

        def contender():
            ready.set()
            start.wait(1)
            time.sleep(0.01)
            observed.append(time.monotonic())

        thread = threading.Thread(target=contender)
        thread.start()
        try:
            self.assertTrue(ready.wait(1))
            before = time.monotonic()
            start.set()
            libc.usleep(50000)
            returned = time.monotonic()
            thread.join(1)
            self.assertFalse(thread.is_alive())
            self.assertTrue(before + 0.005 < observed[0] < returned)
        finally:
            start.set()
            thread.join(1)
            ffi.dlclose(libc)

    @unittest.skipUnless(
        os.environ.get("FACTOR_B4A_LIFECYCLE_QA") == "1",
        "future Stage Q actual PyPy runtime lane",
    )
    def test_signal_interrupt_discards_private_packet_before_reuse(self):
        class SignalStopError(Exception):
            pass

        def stop(*_):
            raise SignalStopError()

        previous = signal.signal(signal.SIGALRM, stop)
        started = time.monotonic()
        try:
            signal.setitimer(signal.ITIMER_REAL, 0.001)
            with self.assertRaises(SignalStopError):
                self.packet.ladder(
                    2**4095 + 19, (2, 1), 2**1024 - 109, 2, allowance()
                )
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, previous)
        self.assertLess(time.monotonic() - started, 1)
        self.assertEqual(
            self.packet.lib.b4_step(self.packet.context, 5, 128), -1
        )
        self.assertEqual(
            self.packet.ladder(13, (2, 1), 1013, 2, allowance()),
            core.scalar_multiply(13, 2, 1, 1013, 2),
        )

    def test_exact_ladders_full_width_scalars_and_aliasing(self):
        for n in (
            3,
            1013,
            1009**2,
            2**63 - 1,
            2**64 + 13,
            10**100 + 267,
            2**1024 - 109,
        ):
            for point in ((0, 1), (1, 0), (n - 1, n - 1), (2 % n, 1)):
                for scalar in (0, 1, 2, 3, 13, 2**79 + 19):
                    expected = core.scalar_multiply(scalar, *point, n, 2 % n)
                    self.assertEqual(
                        self.packet.ladder(
                            scalar, point, n, 2 % n, allowance()
                        ),
                        expected,
                    )

    def test_independent_affine_oracle(self):
        n = 1013
        setup = core.setup_curve(n, 8)
        point, a, b = twist_point(setup.point, setup.a24, n)
        for scalar in (1, 2, 7, 31, 97):
            actual = self.packet.ladder(
                scalar, setup.point, n, setup.a24, allowance()
            )
            self.assertTrue(
                matches(actual, affine_multiply(scalar, point, n, a, b), n)
            )

    def test_reused_ladder_signed_bracket_and_scalar_constants(self):
        self.assertEqual(
            bool(self.packet.lib.b4_fused_expressions()),
            self.packet.build["fuse_ordinary_expressions"],
        )
        for n in (1013, 2**127 - 1, 2**267 - 1, 2**333 - 9, 2**1024 - 109):
            for a24 in (0, 1, n - 1):
                # X=n-1,Z=1 makes AA-BB negative, including a24=n-1.
                for point in ((n - 1, 1), (1, n - 1), (2, 1)):
                    for scalar in (2, 3, 13, 2**31 + 19):
                        self.assertEqual(
                            self.packet.ladder(
                                scalar, point, n, a24, allowance()
                            ),
                            core.scalar_multiply(scalar, *point, n, a24),
                        )

    def test_malformed_abi_canonical_bounds_and_program_references(self):
        ffi, lib, ctx = self.packet.ffi, self.packet.lib, self.packet.context
        owners, values = self.packet._values([1013, 2, 1, 2, 13])
        for owner in (*owners, *values):
            self.addCleanup(ffi.release, owner)
        self.assertEqual(
            lib.b4_ladder(ctx, ABI + 1, *(v[0] for v in values)), -1
        )
        malformed = ffi.new("unsigned char[]", b"\0\1")
        self.addCleanup(ffi.release, malformed)
        values[0].data, values[0].size = malformed, 2
        self.assertEqual(lib.b4_ladder(ctx, ABI, *(v[0] for v in values)), -1)
        owners, values = self.packet._values([1013, 2, 1, 2])
        for owner in (*owners, *values):
            self.addCleanup(ffi.release, owner)
        code = ffi.new("unsigned char[]", b"\1\1\1\xff")
        masks = ffi.new("unsigned char[]", b"\3\3")
        records = ffi.new("b4_record[]", 1)
        for owner in (code, masks, records):
            self.addCleanup(ffi.release, owner)
        records[0].code, records[0].masks = code, masks
        records[0].code_size, records[0].masks_size = 4, 2
        records[0].slots, records[0].output = 2, 1
        self.assertEqual(
            lib.b4_records(ctx, ABI, *(v[0] for v in values), records, 1), -1
        )
        self.assertTrue(owners)
        with self.assertRaises(ValueError):
            self.packet.ladder(2**4096, (2, 1), 1013, 2, allowance())
        with self.assertRaises(ValueError):
            self.packet.ladder(3, (0, 0), 1013, 2, allowance())

    def test_refused_copy_budget_and_cancellation_leave_input_unchanged(self):
        point = [2, 1]
        with self.assertRaises(BudgetExhaustedError):
            self.packet.ladder(13, point, 1013, 2, allowance(work=0))
        with self.assertRaises(BudgetExhaustedError):
            self.packet.ladder(
                2**2047 + 1, point, 1013, 2, allowance(cancelled=lambda: True)
            )
        self.assertEqual(point, [2, 1])
        self.assertEqual(
            self.packet.ladder(13, point, 1013, 2, allowance()),
            core.scalar_multiply(13, *point, 1013, 2),
        )

    def test_foreign_boundaries_force_amortized_worker_budget_polling(self):
        for budget_type in (ActionBudget, ProspectiveActionBudget):
            budget = budget_type(
                seconds=10,
                cpu_seconds=10,
                poll_every=64,
                cancelled=lambda: True,
            )
            budget.until_poll = 64
            with self.assertRaises(BudgetExhaustedError):
                self.packet.ladder(13, (2, 1), 1013, 2, budget)
            self.assertEqual(
                self.packet.lib.b4_step(self.packet.context, 5, 128), -1
            )

    def test_frontier_products_and_strict_saturated_recovery(self):
        plan = ChainPlan(2000, "python-int", allowance())
        powers = [
            [p, entry[0]] for p, entry in list(plan.entries.items())[:16]
        ]
        for n, point in (
            (1013, (2, 1)),
            (1009 * 1013, (1009, 1)),
            (3 * 5 * 7, (2, 1)),
            (1009**2, (1009, 1)),
        ):
            expected = plan.execute(powers, point, n, 2, allowance())
            actual = plan.execute(
                powers, point, n, 2, allowance(), native_packet=self.packet
            )
            self.assertEqual(actual, expected)
            if actual[1] is not None:
                self.assertTrue(1 < actual[1] < n and n % actual[1] == 0)
        entries = [(p, *plan.entries[p]) for p, _ in powers]
        point, product = (2, 1), 1
        for _, _, action, _, _ in entries:
            point, guard = action.run(point, 1013, 2)
            product = product * guard % 1013
            if not product:
                break
        actual, guard, completed = self.packet.records(
            entries, (2, 1), 1013, 2, allowance()
        )
        self.assertEqual((actual, guard), (point, product))
        self.assertTrue(1 <= completed <= 16)

    def test_stage_jobs_preserve_canonical_boundaries_and_pay_copy_work(self):
        for kernel in ("python", "c-gmp"):
            selected = config(ecm_kernel=kernel)
            budget = allowance()
            context, programs = _execution_context({}, selected, budget)
            job = new_job("ecm", 1009 * 1013, 7919, 20, 200)
            try:
                for _ in range(500):
                    advance_job(job, budget, programs or context, selected)
                    if job["done"]:
                        break
                self.assertTrue(job["done"])
            finally:
                close_execution_context(programs or context)
            if kernel == "python":
                expected, work = copy.deepcopy(job), budget.used
            else:
                self.assertEqual(job, expected)
                self.assertGreater(budget.used, work)

    def test_production_prac_chunk_retains_program_identity(self):
        n = 174094069460753795728688384823652336103409078814789169060711
        outcomes = []
        for kernel in ("python", "c-gmp"):
            selected = config(
                ecm_kernel=kernel,
                ecm_tiers=((2000, 200000, 8),),
                ecm_chain_mode="auto",
            )
            budget = allowance()
            context, programs = _execution_context({}, selected, budget)
            job = new_job("ecm", n, 7919, 2000, 200000)
            job.update(phase="stage_one", value=[2, 1], a24=2)
            try:
                for _ in range(10):
                    advance_job(job, budget, programs or context, selected)
                    if job.get("chain_chunks"):
                        break
                self.assertEqual(job["chain_chunks"], 1)
                outcomes.append(job)
            finally:
                close_execution_context(programs or context)
        self.assertEqual(outcomes[0], outcomes[1])

    def test_interrupted_private_packet_leaves_serializable_job_unchanged(
        self,
    ):
        selected = config(ecm_kernel="c-gmp")
        budget = allowance()
        context, programs = _execution_context({}, selected, budget)
        packet = (programs or context).native_packet
        packet.available(1013, budget)
        job = new_job("ecm", 1013, 7919, 20, 200)
        job.update(phase="stage_one", value=[2, 1], a24=2)
        try:
            for _ in range(10):
                if len(job["powers"]) == 8:
                    break
                advance_job(job, budget, programs or context, selected)
            saved = json.loads(json.dumps(job))
            used = budget.used
            with patch.object(
                packet,
                "_finish",
                side_effect=BudgetExhaustedError("cancelled"),
            ):
                with self.assertRaises(BudgetExhaustedError):
                    advance_job(job, budget, programs or context, selected)
            self.assertEqual(job, saved)
            self.assertGreater(budget.used, used)
            self.assertEqual(packet.lib.b4_step(packet.context, 5, 128), -1)
            advance_job(job, budget, programs or context, selected)
            self.assertNotEqual(job, saved)
        finally:
            close_execution_context(programs or context)

    def test_corrupt_program_scalar_or_frontier_identity_is_rejected(self):
        plan = ChainPlan(2000, "python-int", allowance())
        prime, entry = next(iter(plan.entries.items()))
        power, action, unit, reservation = entry
        for changed in (
            SimpleNamespace(
                record=replace(action.record, scalar=power + 1),
                masks=action.masks,
            ),
            SimpleNamespace(
                record=action.record, masks=b"\4" + action.masks[1:]
            ),
        ):
            with self.assertRaises(ValueError):
                self.packet.records(
                    [(prime, power, changed, unit, reservation)],
                    (2, 1),
                    1013,
                    2,
                    allowance(),
                )

    def test_repeated_context_use_closes_and_rejects_active_ownership(self):
        self.packet._lock.acquire()
        try:
            with self.assertRaises(RuntimeError):
                self.packet.close()
        finally:
            self.packet._lock.release()
        for _ in range(3):
            self.packet.ladder(13, (2, 1), 1013, 2, allowance())
        self.packet.close()
        with self.assertRaises(RuntimeError):
            self.packet.available(1013, allowance())

    def test_spawned_worker_loads_the_optional_build(self):
        context = multiprocessing.get_context("spawn")
        parent, child = context.Pipe(duplex=False)
        process = context.Process(target=spawned_smoke, args=(child,))
        process.start()
        child.close()
        try:
            self.assertTrue(parent.poll(5))
            self.assertEqual(
                parent.recv(), (core.scalar_multiply(13, 2, 1, 1013, 2), ABI)
            )
            process.join(5)
            self.assertEqual(process.exitcode, 0)
        finally:
            if process.is_alive():
                process.kill()
                process.join(1)
            parent.close()
            process.close()

    def test_worker_curve_engine_native_result_matches_python(self):
        signatures = []
        for kernel in ("python", "c-gmp"):
            selected = config(ecm_kernel=kernel)
            manifest = Manifest(
                1009 * 1013,
                7919,
                "root",
                arithmetic_options(selected, WorkerSettings()),
            )
            engine = CurveEngine()
            try:
                message = dict(
                    assignment=manifest.assignment(0),
                    options=manifest.options,
                    schedule=manifest.metadata,
                    token=1,
                    attempt=1,
                    work=10**8,
                    cpu=5,
                    deadline=time.monotonic() + 5,
                )
                result = engine.execute(message, cancelled=lambda: False)
                signatures.append(
                    tuple(
                        result[key]
                        for key in (
                            "id",
                            "factor",
                            "sigma",
                            "phase",
                            "chain_chunks",
                        )
                    )
                )
            finally:
                engine.close()
        self.assertEqual(signatures[0], signatures[1])


if __name__ == "__main__":
    unittest.main()
