"""Independent clock, concurrency and search-state checks for live grants."""

import io
import json
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from contextlib import ExitStack, redirect_stdout
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.qs.c10a01 import live_limits, pilots_live, plan_v4
from v2.benchmarks.qs.c10a01.runner import state_digest
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.qs import SIQSJob
from v2.tests.test_c10a01_v4 import config_v4


def call(name="candidate", seconds=2):
    return dict(
        id=name,
        digits=60,
        seconds=seconds,
        reservation_seconds=seconds + 5,
        work=plan_v4.call_work(60, seconds),
        phase="measured",
    )


class LiveGrantTests(unittest.TestCase):
    def test_worker_module_entrypoint_is_importable(self):
        result = subprocess.run(
            [sys.executable, "-m", pilots_live.MODULE, "worker", "--help"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--expect-revision", result.stdout)

    def test_pending_worker_ack_and_monotone_cpu_grant(self):
        with tempfile.TemporaryDirectory() as output:
            wall, cpu = [100.0], [10.0]
            with patch.object(time, "monotonic", side_effect=lambda: wall[0]):
                with patch.object(
                    time, "process_time", side_effect=lambda: cpu[0]
                ):
                    live_limits.initialize(output, "source", [call()])
                    budget = Budget(_wall_start=100, _cpu_start=10)
                    live_limits.LiveBudget(
                        output, "source", "candidate", budget
                    )
                    wall[0], cpu[0] = 100.5, 11.9
                    live_limits.extend(
                        output, "source", ["candidate"], 4, 0, "cpu"
                    )
                    cpu[0] = 12.1
                    budget.consume(1)
                    self.assertEqual(budget.cpu_seconds, 4)
                    self.assertEqual(budget.used, 1)
                    cpu[0] = 13.99
                    budget.consume(1)
                    with self.assertRaises(ValueError):
                        live_limits.extend(
                            output, "source", ["candidate"], 6, 1, "closing"
                        )
                    cpu[0] = 14.1
                    with self.assertRaises(BudgetExhaustedError):
                        budget.consume(1)
                    self.assertEqual(budget.reason, "cpu_limit")

    def test_live_worker_orchestration_receipts_and_close(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            row = dict(call(), fixture="tiny", tier=1, arm="nominee")
            freeze = dict(calls=[row])
            path = output / "freeze.json"
            path.write_text(json.dumps(freeze))
            identity = pilots_live.digest(path)
            live_limits.initialize(output, identity, [row])
            samples = [dict(wall=1, cpu=1, digest="same", work=1)] * 9
            sentinel = dict(samples=samples, work=9)

            def measured(*args, **kwargs):
                budget = Budget()
                live_limits.LiveBudget(output, identity, "candidate", budget)
                budget.consume(7)
                return dict(
                    work=7,
                    seconds=1,
                    cpu_seconds=1,
                    reason="factor_found",
                    complete=True,
                    factors=[101, 103],
                    cofactor=103,
                    divisor=101,
                    unresolved=[],
                    stats={},
                )

            patches = [
                patch.object(pilots_live, "OUTPUT", output),
                patch.object(pilots_live, "FREEZE", path),
                patch.object(pilots_live, "source_check", return_value=freeze),
                patch.object(pilots_live, "measured", side_effect=measured),
                patch.object(pilots_live.original, "host_admission"),
                patch.object(
                    pilots_live.original, "sentinel", return_value=sentinel
                ),
                patch.object(
                    pilots_live.execute_v4, "identity", return_value={}
                ),
                patch.object(
                    pilots_live.execute_v4, "runtime_modules", return_value=()
                ),
                patch.object(pilots_live.execute_v4, "OUTPUT", output),
                patch.object(plan_v4, "RSS", dict(plan_v4.RSS)),
            ]
            with ExitStack() as stack:
                for manager in patches:
                    stack.enter_context(manager)
                with redirect_stdout(io.StringIO()):
                    pilots_live.worker(0)
            state = live_limits.read(output)
            self.assertEqual(state["calls"]["candidate"]["status"], "complete")
            for reference in state["calls"]["candidate"]["receipts"].values():
                self.assertEqual(
                    pilots_live.digest(output / reference["name"]),
                    reference["sha256"],
                )

    def test_extension_crosses_old_deadline_without_resetting_budget(self):
        with tempfile.TemporaryDirectory() as output:
            wall, cpu = [100.0], [10.0]
            with patch.object(time, "monotonic", side_effect=lambda: wall[0]):
                with patch.object(
                    time, "process_time", side_effect=lambda: cpu[0]
                ):
                    live_limits.initialize(output, "source", [call()])
                    budget = Budget(
                        seconds=2,
                        cpu_seconds=2,
                        work_limit=plan_v4.call_work(60, 2),
                        _wall_start=100,
                        _cpu_start=10,
                    )
                    live_limits.LiveBudget(
                        output, "source", "candidate", budget
                    )
                    budget.consume(123)
                    wall[0], cpu[0] = 101.9, 11.9
                    live_limits.extend(
                        output, "source", ["candidate"], 4, 0, "test"
                    )
                    pending = live_limits.read(output)
                    self.assertEqual(
                        pending["calls"]["candidate"]["acknowledged"], 0
                    )
                    self.assertEqual(
                        live_limits.watchdog_deadline(pending, "candidate", 0),
                        184,
                    )
                    wall[0], cpu[0] = 102.1, 12.1
                    budget.consume(5)
                    self.assertEqual(budget.used, 128)
                    self.assertEqual(
                        (budget._wall_start, budget._cpu_start), (100, 10)
                    )
                    self.assertEqual(
                        (budget.seconds, budget.cpu_seconds), (4, 4)
                    )
                    state = live_limits.read(output)
                    self.assertEqual(
                        state["calls"]["candidate"]["acknowledged"], 1
                    )
                    self.assertEqual(
                        state["calls"]["candidate"]["used_at_ack"], 123
                    )
                    self.assertEqual(
                        state["limits"]["phases"]["measured"]["wall"], 9
                    )
                    wall[0] = 104.1
                    with self.assertRaises(BudgetExhaustedError):
                        budget.consume(1)
                    self.assertEqual(budget.used, 128)
                    with self.assertRaises(ValueError):
                        live_limits.extend(
                            output, "source", ["candidate"], 6, 1, "late"
                        )

    def test_cpu_limit_and_cancellation_remain_effective(self):
        with tempfile.TemporaryDirectory() as output:
            live_limits.initialize(output, "source", [call()])
            budget = Budget(seconds=2, cpu_seconds=2, cancelled=lambda: True)
            live_limits.LiveBudget(output, "source", "candidate", budget)
            with self.assertRaises(BudgetExhaustedError):
                budget.consume(1)
            self.assertEqual(budget.reason, "cancelled")
            self.assertEqual(budget.used, 0)

    def test_concurrent_extensions_have_one_cas_winner(self):
        with tempfile.TemporaryDirectory() as output:
            live_limits.initialize(output, "source", [call()])
            barrier, results = threading.Barrier(2), []

            def request(seconds):
                barrier.wait()
                try:
                    live_limits.extend(
                        output, "source", ["candidate"], seconds, 0, "race"
                    )
                    results.append("accepted")
                except ValueError:
                    results.append("stale")

            threads = [
                threading.Thread(target=request, args=(s,)) for s in (3, 4)
            ]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join(timeout=5)
                self.assertFalse(thread.is_alive())
            self.assertEqual(sorted(results), ["accepted", "stale"])
            self.assertEqual(live_limits.read(output)["revision"], 1)

    def test_refusal_is_atomic_and_interruption_keeps_latest_reservation(self):
        with tempfile.TemporaryDirectory() as output:
            live_limits.initialize(output, "source", [call(), call("control")])
            initial = (Path(output) / "live.json").read_bytes()
            initial_limits = live_limits.read(output)["limits"]["total"]
            for seconds, revision, identity in (
                (2, 0, "source"),
                (3, 1, "source"),
                (3, 0, "other"),
            ):
                with self.assertRaises(ValueError):
                    live_limits.extend(
                        output,
                        identity,
                        ["candidate"],
                        seconds,
                        revision,
                        "refuse",
                    )
                self.assertEqual(
                    (Path(output) / "live.json").read_bytes(), initial
                )
            live_limits.extend(
                output, "source", ["candidate", "control"], 4, 0, "both"
            )
            extended_limits = live_limits.read(output)["limits"]["total"]
            self.assertEqual(
                extended_limits["wall"], initial_limits["wall"] + 14
            )
            self.assertEqual(
                extended_limits["work"],
                initial_limits["work"] + 128_000_000_000,
            )
            live_limits.interrupt(output, "source", "candidate", "lost worker")
            state = live_limits.read(output)
            self.assertEqual(
                state["calls"]["candidate"]["reserved_interrupted"],
                dict(wall=9, cpu=9, work=plan_v4.call_work(60, 4)),
            )
            with self.assertRaises(ValueError):
                live_limits.extend(
                    output, "source", ["control"], 6, 1, "closed"
                )

    def test_closed_calls_and_corrupt_grants_are_refused(self):
        with tempfile.TemporaryDirectory() as output:
            live_limits.initialize(output, "source", [call()])
            budget = Budget()
            live_limits.LiveBudget(output, "source", "candidate", budget)
            live_limits.close_measurement(
                output,
                "source",
                "candidate",
                "result.json",
                "hash",
                dict(work=0, seconds=0, cpu_seconds=0, reason="factor_found"),
            )
            with self.assertRaises(ValueError):
                live_limits.extend(
                    output, "source", ["candidate"], 4, 0, "closed"
                )
            path = Path(output) / "live.json"
            state = json.loads(path.read_text())
            state["calls"]["candidate"]["work"] -= 1
            path.write_text(json.dumps(state))
            with self.assertRaises(ValueError):
                live_limits.read(output)

    def test_finite_journal_and_state_size(self):
        with tempfile.TemporaryDirectory() as output:
            live_limits.initialize(output, "source", [call()])
            budget = Budget()
            live_limits.LiveBudget(output, "source", "candidate", budget)
            before = (Path(output) / "live.json").read_bytes()
            with patch.object(
                live_limits, "CONTROL_HEADROOM", live_limits.STATE_BYTES
            ):
                with self.assertRaises(RuntimeError):
                    live_limits.extend(
                        output, "source", ["candidate"], 4, 0, "full"
                    )
            self.assertEqual((Path(output) / "live.json").read_bytes(), before)
            budget.consume(1)
            self.assertEqual(budget.used, 1)
        with tempfile.TemporaryDirectory() as output:
            with patch.object(live_limits, "MAX_REVISIONS", 1):
                live_limits.initialize(output, "source", [call()])
                live_limits.extend(
                    output, "source", ["candidate"], 3, 0, "one"
                )
                with self.assertRaises(RuntimeError):
                    live_limits.extend(
                        output, "source", ["candidate"], 4, 1, "two"
                    )
        with tempfile.TemporaryDirectory() as output:
            with patch.object(live_limits, "STATE_BYTES", 10):
                with self.assertRaises(RuntimeError):
                    live_limits.initialize(output, "source", [call()])
            self.assertFalse((Path(output) / "live.json").exists())

    def test_hook_preserves_exact_siqs_semantic_state_and_work(self):
        for multiplier in (1, 2, 3):
            with self.subTest(multiplier=multiplier):
                fingerprints = []
                for live in (False, True):
                    with tempfile.TemporaryDirectory() as output:
                        budget = Budget(
                            seconds=10,
                            cpu_seconds=10,
                            work_limit=plan_v4.call_work(60, 10),
                        )
                        job = SIQSJob(
                            10403,
                            seed=7,
                            config=config_v4(multiplier=multiplier),
                            budget=budget,
                        )
                        if live:
                            live_limits.initialize(
                                output, "source", [call(seconds=10)]
                            )
                            live_limits.LiveBudget(
                                output, "source", "candidate", budget
                            )
                        result = job.run(max_blocks=16)
                        self.assertEqual(
                            (result.divisor or 1) * result.cofactor, 10403
                        )
                        fingerprints.append(
                            (
                                result.divisor,
                                result.reason,
                                budget.used,
                                state_digest(job),
                            )
                        )
                self.assertEqual(fingerprints[0], fingerprints[1])

    def test_running_subprocess_keeps_pid_and_progress_past_old_cap(self):
        with tempfile.TemporaryDirectory() as output:
            live_limits.initialize(output, "source", [call()])
            code = """
import json,os,sys,time
from pathlib import Path
from v2.execution.budget import Budget
from v2.benchmarks.qs.c10a01 import live_limits,plan_v4
out=sys.argv[1]
b=Budget(seconds=2,cpu_seconds=2,work_limit=plan_v4.call_work(60,2))
origin=b._wall_start
live_limits.LiveBudget(out,'source','candidate',b)
Path(out,'ready').write_text(str(os.getpid()))
while b.wall_used<2.5:
 b.consume(1)
print(json.dumps(dict(pid=os.getpid(),origin_unchanged=origin==b._wall_start,used=b.used,seconds=b.seconds,wall=b.wall_used)))
"""
            process = subprocess.Popen(
                [sys.executable, "-c", code, output],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            try:
                deadline = time.monotonic() + 10
                while not (Path(output) / "ready").exists():
                    if (
                        process.poll() is not None
                        or time.monotonic() >= deadline
                    ):
                        self.fail("worker did not become ready")
                    time.sleep(0.01)
                pid = int((Path(output) / "ready").read_text())
                live_limits.extend(
                    output, "source", ["candidate"], 5, 0, "subprocess"
                )
                stdout, stderr = process.communicate(timeout=10)
                self.assertEqual(process.returncode, 0, stderr)
                result = json.loads(stdout)
                self.assertEqual(result["pid"], pid)
                self.assertTrue(result["origin_unchanged"])
                self.assertEqual(result["seconds"], 5)
                self.assertGreater(result["wall"], 2)
                self.assertGreater(result["used"], 0)
                state = live_limits.read(output)
                self.assertEqual(
                    state["calls"]["candidate"]["acknowledged"], 1
                )
                self.assertGreater(
                    result["used"], state["calls"]["candidate"]["used_at_ack"]
                )
            finally:
                if process.poll() is None:
                    process.kill()
                    process.wait()


if __name__ == "__main__":
    unittest.main()
