"""Unrun pilot isolation, native censoring and resource contracts."""

import fcntl
import json
import os
import subprocess
import sys
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.qs.c10a01 import pilots
from v2.benchmarks.qs.c10a01.configuration import configuration


class ImmediateLedger:
    def __init__(self):
        self.names = []

    def call(self, phase, name, seconds, work, identity, function):
        self.names.append((phase, name))
        return function()


class PilotIsolationTests(unittest.TestCase):
    def test_prospective_memory_suballocation_preserves_policies(self):
        plan = json.loads(pilots.PLAN.read_text())
        self.assertTrue(pilots.verify_plan(plan))
        old = json.loads(pilots.PLAN.with_name("pilots.json").read_text())
        for name, stage in plan["stages"].items():
            self.assertEqual(stage["limits"], old["stages"][name]["limits"])
            for call in stage["calls"]:
                expected = asdict(
                    configuration(
                        call["digits"],
                        call["arm"],
                        multiplier=call["multiplier"],
                    )
                )
                self.assertEqual(
                    expected.pop("root_update_policy"), "standard"
                )
                expected["memory_bytes"] -= pilots.PARENT_OWNED
                self.assertEqual(call["configuration"], expected)
        self.assertNotEqual(plan["identity"], old["identity"])

    def test_memory_censoring_does_not_run_post_sentinel_or_poison_next_cell(
        self,
    ):
        call = dict(
            digits=70,
            arm="control",
            expected_multiplier=1,
            seconds=7200,
            work=1,
        )
        before = dict(
            warmup=dict(
                last=dict(
                    multiplier=1,
                    stats=dict(capacity=dict(target_in_product_envelope=True)),
                )
            ),
            samples=[dict(wall=1, cpu=1)],
        )
        ledger = ImmediateLedger()
        with (
            patch.object(pilots, "sentinel", return_value=before),
            patch.object(
                pilots, "post_sentinel", return_value=before
            ) as after,
            patch.object(pilots, "complete_call") as complete,
        ):
            complete.return_value = dict(
                complete=False, reason="rss_limit", stats={"large": [1] * 1000}
            )
            first = pilots.execute_cell(ledger, {}, call, "first", {})
            after.assert_not_called()
            self.assertFalse(first["drift_pass"])
            self.assertIsNone(first["drift"])
            self.assertNotIn("stats", first["result"])
            self.assertEqual(first["result"]["reason"], "rss_limit")

            complete.return_value = dict(complete=False, reason="time_limit")
            second = pilots.execute_cell(ledger, {}, call, "second", {})
            after.assert_called_once()
            self.assertTrue(second["drift_pass"])
            self.assertEqual(second["result"]["reason"], "time_limit")
        self.assertNotIn(("warmup", "first-after"), ledger.names)
        self.assertIn(("measured", "second"), ledger.names)

    def test_failed_post_drift_is_censored_only_for_native_resource_failure(
        self,
    ):
        with patch.object(pilots, "sentinel", side_effect=MemoryError("cap")):
            result = pilots.post_sentinel({}, {"digits": 70})
        self.assertTrue(result["missing"])
        self.assertEqual(result["work"], 290000000000)
        with patch.object(pilots, "sentinel", side_effect=RuntimeError("bad")):
            with patch.object(pilots.runner, "rss", return_value=1):
                with self.assertRaisesRegex(RuntimeError, "bad"):
                    pilots.post_sentinel({}, {"digits": 70})

    def test_every_cell_executes_a_new_worker_and_failure_is_not_retried(self):
        call = dict(digits=70, seconds=7200, arm="control")
        row = dict(call=call, process_id=os.getpid() + 1)
        result = subprocess.CompletedProcess([], 0, json.dumps(row), "")
        with (
            patch.object(pilots.runner, "rss", return_value=1),
            patch.object(pilots.subprocess, "run", return_value=result) as run,
        ):
            for index in range(2):
                self.assertEqual(
                    pilots.launch_worker("C", index, call, 17), row
                )
            self.assertEqual(run.call_count, 2)
            for index, invocation in enumerate(run.call_args_list):
                argv = invocation.args[0]
                self.assertEqual(argv[0], pilots.sys.executable)
                self.assertIn("worker", argv)
                self.assertEqual(
                    argv[argv.index("--cell-index") + 1], str(index)
                )
                self.assertEqual(invocation.kwargs["pass_fds"], (17,))
                self.assertEqual(invocation.kwargs["timeout"], 7275)
            run.reset_mock()
            run.side_effect = subprocess.TimeoutExpired("worker", 7275)
            with self.assertRaises(subprocess.TimeoutExpired):
                pilots.launch_worker("C", 0, call, 17)
            run.assert_called_once()

    def test_parent_refusal_precedes_worker_launch(self):
        with (
            patch.object(
                pilots.runner, "rss", return_value=pilots.PARENT_RSS + 1
            ),
            patch.object(pilots.subprocess, "run") as run,
        ):
            with self.assertRaisesRegex(RuntimeError, "parent exceeds"):
                pilots.launch_worker("C", 0, {}, 17)
            run.assert_not_called()


class PilotWorkerTests(unittest.TestCase):
    def test_real_exec_lock_and_ledger_with_stubbed_numerical_calls(self):
        # Exercise exec, inherited flock and durable worker receipts without
        # running any warmup or complete-factor pilot.
        script = r"""
import json
import os
import sys
from pathlib import Path
from unittest.mock import patch
from v2.benchmarks.qs.c10a01 import pilots

output, lock_fd, index = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
plan = json.loads(pilots.PLAN.read_text())
call = plan["stages"]["B"]["calls"][index]
fixture = dict(id=call["fixture"], n=35)
real_stat = os.stat

def expected_lock(path, *args, **kwargs):
    if str(path) == "/private/tmp/factor-performance.lock":
        return os.fstat(lock_fd)
    return real_stat(path, *args, **kwargs)

before = dict(
    work=17,
    warmup=dict(last=dict(multiplier=call["expected_multiplier"],
        stats=dict(capacity=dict(target_in_product_envelope=True)))),
    samples=[dict(wall=1, cpu=1)],
)
result = dict(n=35, complete=False, cofactor=35, factors=[], certainty=[],
    unresolved=[35], work=23, seconds=1, cpu_seconds=1,
    reason="rss_limit" if index == 0 else "time_limit")
with (
    patch.object(pilots, "readiness", return_value=(plan,
        dict(fixtures=[fixture]), {})),
    patch.object(pilots, "output_path", return_value=output),
    patch.object(pilots, "plan_identity", return_value=dict(test=True)),
    patch.object(pilots.os, "stat", side_effect=expected_lock),
    patch.object(pilots, "sentinel", return_value=before),
    patch.object(pilots, "complete_call", return_value=result),
):
    pilots.worker("B", index, lock_fd)
"""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with (root / "performance.lock").open("a+") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                rows = []
                for index in range(2):
                    process = subprocess.run(
                        [
                            sys.executable,
                            "-c",
                            script,
                            str(root / "receipts"),
                            str(lock.fileno()),
                            str(index),
                        ],
                        cwd=pilots.runner.ROOT,
                        pass_fds=(lock.fileno(),),
                        capture_output=True,
                        text=True,
                        check=True,
                        timeout=30,
                    )
                    rows.append(json.loads(process.stdout))
            self.assertNotEqual(rows[0]["process_id"], rows[1]["process_id"])
            self.assertTrue(
                all(row["process_id"] != os.getpid() for row in rows)
            )
            self.assertIsNone(rows[0]["drift"])
            self.assertTrue(rows[1]["drift_pass"])
            for row in rows:
                self.assertEqual(
                    row["worker_rss_cap"], 5 * 2**30 - 256 * 2**20
                )
                self.assertEqual(row["parent_owned_reservation"], 32 * 2**20)
            ledger = json.loads(
                (root / "receipts" / "ledger.json").read_text()
            )
            self.assertEqual(len(ledger["calls"]), 5)
            self.assertTrue(
                all(row["status"] == "complete" for row in ledger["calls"])
            )
            self.assertEqual(sum(row["work"] for row in ledger["calls"]), 97)
