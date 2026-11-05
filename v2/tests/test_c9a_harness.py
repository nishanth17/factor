"""Finite C9a planning, durable failures and independently checked replay."""

import json
import subprocess
import sys
import tempfile
import time
import unittest
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.qs.c9a.accounting import StudyLedger
from v2.benchmarks.qs.c9a.analysis import observer_summary, summarize
from v2.benchmarks.qs.c9a.configuration import admission_envelope
from v2.benchmarks.qs.c9a.gates import outcomes
from v2.benchmarks.qs.c9a.inputs import specification
from v2.benchmarks.qs.c9a.journal import (
    Journal,
    RecorderExhaustedError,
    atom_record,
    observe,
    records,
)
from v2.benchmarks.qs.c9a.launcher import (
    OUTPUT_BYTES,
    collect_output,
    validate_result,
)
from v2.benchmarks.qs.c9a.plan import build, cell, observer_cells, validate
from v2.benchmarks.qs.c9a.replay import replay
from v2.benchmarks.qs.c10a04.accounting import atomic
from v2.qs import SIQSConfig
from v2.tests.test_c9a_yield import collector, parity_witness, policy
from v2.tests.test_qs_dlp import admit


def resource_limits():
    return dict(
        A=dict(wall=100, cpu=100, work=10000, warm=20, disk_bytes=20 * 2**20),
        total=dict(wall=100, cpu=100, work=10000),
        disk_bytes=32 * 2**20,
    )


def worker_cell():
    return dict(
        sample=0,
        fixture="oracle",
        digits=5,
        n=10403,
        seed=7,
        arm="C1",
        seconds=5,
        margin=2,
        warm_seconds=4,
        warm_work_limit=100,
        work_limit=500,
        observer_work_limit=0,
    )


class PlanningTests(unittest.TestCase):
    def test_numerical_plan_counts_domains_and_no_full_80(self):
        fixtures, resolved = [], {}
        for digits, uneven, residue, name in specification():
            n = 10 ** (digits - 1) + residue
            fixtures.append(dict(id=name, n=n, digits=digits, uneven=uneven))
            resolved[name] = dict(bound=100000, multiplier=1)
        plan = build(dict(fixtures=fixtures), resolved)
        self.assertEqual(validate(plan)["A"], 40)
        self.assertEqual(len(observer_cells(plan)), 18)
        with self.assertRaises(ValueError):
            cell(
                dict(id="bad", n=10**79 + 7, digits=80),
                dict(bound=600000),
                "F",
                "C1",
            )
        for bound in (30000, 100000, 189020, 245000, 400000, 600000):
            limits = admission_envelope(bound)
            self.assertGreaterEqual(
                limits["large_product_bound"], 115200000000
            )
            self.assertGreaterEqual(
                limits["small_prime_bound"] ** 2
                + 2 * limits["small_prime_bound"]
                + 1,
                115200000000,
            )
            self.assertEqual(limits["candidate_bound"], 115200000000)

    def test_parent_reconstructs_with_independent_proofs(self):
        row = dict(n=10403, factors=[101, 103], uneven=False)
        selected = dict(worker_cell(), rss_bytes=100000)
        result = dict(
            n=10403,
            seed=7,
            arm="C1",
            divisor=101,
            factors=[101, 103],
            remaining=[],
            timing_valid=True,
            rss_bytes=1000,
            seconds=1,
            cpu_seconds=1,
        )
        self.assertTrue(validate_result(result, row, selected)["complete"])
        warning = validate_result(
            dict(result, timing_valid=False), row, selected
        )
        self.assertTrue(warning["resource_valid"])
        self.assertFalse(warning["timing_valid"])
        bad = dict(result, factors=[101, 107])
        with self.assertRaises(ArithmeticError):
            validate_result(bad, row, selected)


class AccountingTests(unittest.TestCase):
    def test_controller_restarts_reserve_pool_before_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            ledger = StudyLedger(directory, "oracle", resource_limits())
            with patch.dict(
                "v2.benchmarks.qs.c9a.accounting.CONTROLLER_SECONDS", A=15
            ):
                ledger.call_controller("A", 10, 100, lambda: dict(work=40))
                with self.assertRaises(RuntimeError):
                    ledger.call_controller(
                        "A", 10, 100, lambda: self.fail("executed")
                    )
            self.assertEqual(ledger.used(), 10)
            self.assertEqual(ledger.used(field="work"), 40)

    def test_work_settlement_and_missing_work_preserve_charge(self):
        with tempfile.TemporaryDirectory() as directory:
            ledger = StudyLedger(directory, "oracle", resource_limits())
            value = dict(
                measured_work=27, warm_work=13, observer_work=0, work=40
            )
            self.assertEqual(
                ledger.call_cell("A", "one", worker_cell(), lambda: value),
                value,
            )
            self.assertEqual(ledger.used(field="work"), 40)
            self.assertEqual(ledger.used(), 11)
            with self.assertRaises(ValueError):
                ledger.call_cell("A", "missing", worker_cell(), lambda: {})
            self.assertEqual(ledger.used(field="work"), 640)
            with self.assertRaises(RuntimeError):
                StudyLedger(directory, "oracle", resource_limits())

    def test_warm_and_output_reservations_precede_worker(self):
        with tempfile.TemporaryDirectory() as directory:
            limits = resource_limits()
            limits["A"]["warm"] = 3
            ledger = StudyLedger(directory, "oracle", limits)
            with self.assertRaises(RuntimeError):
                ledger.call_cell(
                    "A",
                    "refused",
                    worker_cell(),
                    lambda: self.fail("executed"),
                )
            self.assertEqual(ledger.data["calls"], [])
            limits["A"]["warm"] = 20
            limits["A"]["disk_bytes"] = 1024
            with self.assertRaises(RuntimeError):
                ledger.call_cell(
                    "A", "full", worker_cell(), lambda: self.fail("executed")
                )

    def test_receipt_before_ledger_recovers_without_reexecution(self):
        with tempfile.TemporaryDirectory() as directory:
            ledger = StudyLedger(directory, "oracle", resource_limits())
            original = ledger.write
            calls = 0

            def interrupted_write():
                nonlocal calls
                calls += 1
                if calls > 1:
                    raise OSError("injected ledger publication interruption")
                return original()

            with patch.object(ledger, "write", side_effect=interrupted_write):
                with self.assertRaises(OSError):
                    ledger.call_cell(
                        "A",
                        "one",
                        worker_cell(),
                        lambda: dict(
                            work=40,
                            measured_work=30,
                            warm_work=10,
                            observer_work=0,
                        ),
                    )
            recovered = StudyLedger(directory, "oracle", resource_limits())
            self.assertEqual(recovered.data["calls"][0]["status"], "complete")
            self.assertEqual(recovered.used(field="work"), 40)

    def test_lost_worker_is_censored_and_never_selectively_retried(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            identity = dict(worker_cell(), cap=5)
            atomic(
                directory / "ledger.json",
                dict(
                    identity="oracle",
                    limits=resource_limits(),
                    calls=[
                        dict(
                            index=0,
                            phase="A",
                            receipt="lost.json",
                            identity=identity,
                            status="reserved",
                            wall=11,
                            cpu=11,
                            work=600,
                        )
                    ],
                ),
            )
            ledger = StudyLedger(directory, "oracle", resource_limits())
            result = json.loads((directory / "lost.json").read_text())[
                "result"
            ]
            self.assertFalse(result["complete"])
            self.assertEqual(result["reason"], "interrupted")
            self.assertEqual(ledger.used(field="work"), 600)
            view = outcomes(ledger, "A")
            self.assertEqual(view[0]["fixture"], "oracle")
            self.assertFalse(view[0]["timing_valid"])
            self.assertEqual(
                summarize(view, samples=100)["classes"][0]["missing_pairs"], 1
            )


class JournalTests(unittest.TestCase):
    def test_live_toy_graph_journal_keeps_publication_lifetimes(self):
        import v2.qs as qs

        original = qs.SieveCollector._admit
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "toy"
            writer = Journal(path, {})
            with observe(qs, writer):
                target = collector(
                    cycle_policy="short-v1", record_policy="packed-v1"
                )
                for position, pair in (
                    (-5, (3, 11)),
                    (-4, (5, 11)),
                    (1, (3, 5)),
                ):
                    admit(target, position, pair)
            writer.close(dict(reason="toy_qa"), complete=True)
            events = list(records(path))
            graphs = [r["data"] for r in events if r["kind"] == "graph"]
            self.assertEqual(len(graphs), 3)
            self.assertEqual(len(graphs[-1]["pinned"]), 2)
        self.assertIs(qs.SieveCollector._admit, original)

    def test_recorder_exhaustion_is_distinct_and_framed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "stream"
            journal = Journal(path, {}, byte_limit=8192)
            with self.assertRaises(RecorderExhaustedError):
                journal.emit("atom_begin", {"padding": "x" * 8000})
            journal.close({"reason": "recorder_limit"}, complete=False)
            self.assertTrue(journal.exhausted)
            self.assertEqual(
                list(records(path, require_complete=False))[-1]["kind"],
                "footer",
            )
            with self.assertRaises(ValueError):
                list(records(path))
            self.assertFalse(issubclass(RecorderExhaustedError, MemoryError))

    def test_order_hash_and_pending_publication_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "stream"
            journal = Journal(path, {})
            journal.emit("atom_begin", {})
            journal.close({}, complete=True)
            with self.assertRaises(ValueError):
                list(records(path))
            lines = path.read_bytes().splitlines(keepends=True)
            path.write_bytes(b"".join(reversed(lines)))
            with self.assertRaises(ValueError):
                list(records(path))

    def test_exact_replay_of_equal_parity_factor_witness(self):
        job = parity_witness()
        config = SIQSConfig(
            base_bound=40,
            half_width=200,
            max_half_width=200,
            solve_policy="collection",
            collector=policy(),
        )
        selected = dict(
            fixture="oracle",
            digits=50,
            n=10403,
            seed=7,
            configuration=asdict(config),
            resolved=dict(bound=40, multiplier=1),
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "stream"
            writer = Journal(path, dict(format="c9a-atoms-v1", cell=selected))
            polynomial = job.collector.polynomial
            writer.emit(
                "window",
                dict(
                    polynomial=[polynomial.a, polynomial.b, 1],
                    lo=0,
                    hi=200,
                    memory_bytes=64 * 2**20,
                ),
            )
            for atom in job.collector._atoms.values():
                begin = writer.emit("atom_begin", atom_record(atom))
                writer.emit("atom_end", dict(begin=begin))
            begin = writer.emit(
                "solve_begin", dict(final=True, next_position=200, hi=200)
            )
            writer.emit("solve_end", dict(begin=begin))
            writer.close(
                dict(reason="factor_found", divisor=101), complete=True
            )
            result = replay(path, "C1", seconds=5, work_limit=10**9)
            self.assertIn(result["divisor"], (101, 103))
            self.assertGreater(result["work"], 0)


class UncertaintyTests(unittest.TestCase):
    def test_seeds_are_nested_and_invalid_times_are_explicit(self):
        results = []
        for modulus in range(3):
            for seed in (7, 11):
                for arm, seconds in (("C0", 2), ("C1", 1)):
                    results.append(
                        dict(
                            fixture=str(modulus),
                            seed=seed,
                            arm=arm,
                            digits=60,
                            uneven=True,
                            complete=True,
                            seconds=seconds,
                            capped_seconds=seconds,
                            timing_valid=True,
                        )
                    )
        result = summarize(results, samples=100)["classes"][0]
        self.assertEqual(result["independent_moduli"], 3)
        self.assertGreater(result["zero_event_regression_upper95"], 0.05)
        results[0]["timing_valid"] = False
        changed = summarize(results, samples=100)["classes"][0]
        self.assertEqual(changed["invalid_timing_pairs"], 0)
        self.assertEqual(changed["legacy_canary_warning_pairs"], 1)
        self.assertEqual(changed["complete_cost_moduli"], 3)
        historical = summarize(results, samples=100, legacy_timing=True)
        self.assertEqual(historical["classes"][0]["invalid_timing_pairs"], 1)
        self.assertEqual(historical["classes"][0]["complete_cost_moduli"], 2)

        results[0]["resource_valid"] = False
        invalid = summarize(results, samples=100)["classes"][0]
        self.assertEqual(invalid["invalid_timing_pairs"], 1)
        self.assertEqual(invalid["complete_cost_moduli"], 2)

    def test_observer_keeps_canary_warnings_and_rejects_recorder_failures(
        self,
    ):
        results = []
        for _ in range(9):
            for observation, seconds in (("clean", 1), ("capture", 2)):
                results.append(
                    dict(
                        observation=observation,
                        seconds=seconds,
                        cpu_seconds=seconds,
                        divisor=101,
                        reason="factor_found",
                        next_position=17,
                        measured_work=100,
                        stats={},
                        timing_valid=False,
                        resource_valid=True,
                        recorder=dict(exhausted=False, bytes=200),
                    )
                )
        result = observer_summary(results)
        self.assertEqual(result["overhead_ratio"], 2)
        self.assertEqual(result["legacy_canary_warning_pairs"], 9)
        results[1]["resource_valid"] = False
        with self.assertRaises(ValueError):
            observer_summary(results)
        results[1]["resource_valid"] = True
        results[1]["recorder"]["exhausted"] = True
        with self.assertRaises(ValueError):
            observer_summary(results)


class WorkerFailureTests(unittest.TestCase):
    def test_rewarm_policy_preserves_the_original_failed_screen_and_result(
        self,
    ):
        from v2.benchmarks.qs.c9a import worker

        config = SIQSConfig(collector=policy(), memory_bytes=64 * 2**20)
        selected = dict(
            worker_cell(),
            source="candidate",
            configuration=asdict(config),
            warm_seconds=30,
            end_reserve=10,
            warm_work_limit=10**6,
            work_limit=10**9,
            max_blocks=1,
            rss_bytes=3 * 2**30,
            timing_policy="retire-rewarm-v1",
        )
        checked = dict(work=10, signature="checked", wall=[1] * 9, cpu=[1] * 9)
        drifted = dict(checked, wall=[2] * 9, cpu=[2] * 9)
        result = SimpleNamespace(
            divisor=101,
            reason="factor_found",
            stats={"blocks": 1},
            next_position=17,
        )
        job = SimpleNamespace(
            multiplier=3,
            base=SimpleNamespace(bound=40),
            engine=None,
            run=lambda **kwargs: result,
        )
        clock = 0

        def tick():
            nonlocal clock
            clock += 0.01
            return clock

        with (
            patch.object(
                worker,
                "time",
                SimpleNamespace(monotonic=tick, process_time=tick),
            ),
            patch.object(worker, "warm_once", return_value=("checked", 1)),
            patch.object(
                worker, "sentinel", side_effect=[checked, drifted, checked]
            ),
            patch.object(worker.candidate, "SIQSJob", return_value=job),
            patch.object(worker.gc, "collect") as collect,
        ):
            output = worker.run(selected)

        self.assertFalse(output["legacy_drift"]["valid"])
        self.assertTrue(output["timing_valid"])
        self.assertEqual(output["factors"], [101, 103])
        self.assertEqual(output["next_position"], 17)
        self.assertEqual(output["actual_multiplier"], 3)
        self.assertEqual(output["actual_base_bound"], 40)
        self.assertGreaterEqual(output["retirement"]["rewarm_wall"], 3)
        self.assertGreaterEqual(output["retirement"]["rewarm_cpu"], 3)
        self.assertGreater(output["warm_work"], 30)
        collect.assert_called_once()

    def test_small_warm_fixture_is_checked_in_both_runtime_paths(self):
        from v2.benchmarks.qs.c9a import worker
        from v2.benchmarks.qs.c9a.control import control_modules

        template = asdict(
            SIQSConfig(
                collector=policy(
                    preparation_policy="target-v1",
                    cycle_policy="short-v1",
                    duplicate_policy="pair-v1",
                    record_policy="packed-v1",
                )
            )
        )
        for qs, budgets, selected in (
            (worker.candidate, worker.candidate_budget, template),
            (*control_modules(), asdict(SIQSConfig())),
        ):
            first = worker.warm_once(
                qs, budgets, selected, time.monotonic() + 5
            )
            second = worker.warm_once(
                qs, budgets, selected, time.monotonic() + 5
            )
            self.assertEqual(first, second)
            self.assertGreater(first[1], 0)

    def test_worker_output_exhaustion_preserves_bounded_transcripts(self):
        with tempfile.TemporaryDirectory() as directory:
            out_path, err_path = (
                Path(directory) / "out",
                Path(directory) / "err",
            )
            process = subprocess.Popen(
                [
                    sys.executable,
                    "-c",
                    "import sys; sys.stdout.write('x'*300000)",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                start_new_session=True,
            )
            with out_path.open("wb") as out, err_path.open("wb") as err:
                reason = collect_output(
                    process, out, err, time.monotonic() + 5, 2**30
                )
            self.assertEqual(reason, "output_limit")
            self.assertEqual(out_path.stat().st_size, OUTPUT_BYTES)
            self.assertLessEqual(err_path.stat().st_size, OUTPUT_BYTES)

    def test_start_and_end_work_exhaustion_are_real_worker_failures(self):
        from v2.benchmarks.qs.c9a import worker

        config = SIQSConfig(
            base_bound=40,
            half_width=32,
            max_half_width=32,
            collector=policy(),
            memory_bytes=64 * 2**20,
        )
        selected = dict(
            worker_cell(),
            source="candidate",
            configuration=asdict(config),
            warm_seconds=20,
            end_reserve=10,
            warm_work_limit=100,
            work_limit=10**9,
            max_blocks=1,
            rss_bytes=3 * 2**30,
        )
        good = dict(work=10, signature="checked", wall=[1] * 9, cpu=[1] * 9)
        bad = dict(good, work=1000)
        for sentinels, message in (
            ([bad], "end sentinel reservation would"),
            ([good, bad], "end sentinel reservation exhausted"),
        ):
            count = 0

            def tick():
                nonlocal count
                count += 1
                return count

            with (
                patch.object(
                    worker,
                    "time",
                    SimpleNamespace(monotonic=tick, process_time=tick),
                ),
                patch.object(worker, "warm_once", return_value=("checked", 1)),
                patch.object(worker, "sentinel", side_effect=sentinels),
            ):
                with self.assertRaisesRegex(RuntimeError, message):
                    worker.run(selected)
