"""Opt-in SSS dispatch, checked full checkpoints and CLI reconstruction."""

import copy
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from v2 import utils
from v2.budget import Budget
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import SIQSConfig
from v2.qs.sss import SSSConfig, SSSJob

ROOT = Path(__file__).resolve().parents[2]


def allowance(work=10**9):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


def configuration(mode="sss"):
    return PortfolioConfig(
        trial_bound=2,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=(),
        memory_bytes=80 * 1024 * 1024,
        sss=SSSConfig(mode=mode, base_bound=400),
    )


def reseal(checkpoint, payload):
    changed = copy.deepcopy(checkpoint)
    changed["blob"] = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    )
    changed["sha256"] = hashlib.sha256(changed["blob"].encode()).hexdigest()
    return changed


class SSSDispatchTests(unittest.TestCase):
    def test_selected_dispatch_recursion_and_resource_caps(self):
        for mode in ("sss", "sssf"):
            config = configuration(mode)

            for n in (
                4001 * 4003,
                -4001 * 4003,
                (4001 * 4003) ** 2,
                2**4 * 4001 * 4003,
            ):
                run = factorize_bounded(
                    n, seed=31, config=config, budget=allowance()
                )

                self.assertTrue(run.result.complete)
                self.assertEqual(run.result.reconstruct(), n)
                self.assertTrue(run.result.proven)
                self.assertIn(mode, [e["stage"] for e in run.events])

            run = factorize_bounded(
                4001 * 4003, config=config, budget=allowance(0)
            )

            self.assertEqual(run.reason, "work_limit")
            self.assertEqual(run.result.remaining, (4001 * 4003,))

        self.assertIsNone(PortfolioConfig().sss)
        with self.assertRaises(MemoryError):
            replace(configuration(), memory_bytes=64 * 1024 * 1024)
        with self.assertRaises(ValueError):
            replace(configuration(), siqs=SIQSConfig())

    def test_fallback_after_ecm_and_checkpoint_resume(self):
        config = replace(configuration(), ecm_tiers=((5, 7, 1),))
        n = 4001 * 4003

        first = factorize_bounded(
            n, config=config, seed=31, budget=allowance(150000)
        )

        self.assertFalse(first.result.complete)
        self.assertEqual(first.result.reconstruct(), n)
        self.assertEqual(
            first.checkpoint["payload"]["state"]["current"]["stage"], "sss"
        )

        resumed = factorize_bounded(
            n, config=config, checkpoint=first.checkpoint, budget=allowance()
        )

        self.assertTrue(resumed.result.complete)
        self.assertGreater(resumed.work_used, first.work_used)
        stages = [event["stage"] for event in resumed.events]

        self.assertLess(stages.index("ecm"), stages.index("sss"))
        self.assertEqual(resumed.result.reconstruct(), n)

    def test_full_checkpoint_rebuilds_interrupted_assignment_and_solver(self):
        config = SSSConfig(base_bound=400)
        n = 4001 * 4003

        for stage in ("assignment", "solver", "extractor"):
            job = SSSJob(n, config=config, budget=allowance())

            def cancel():
                engine = job.pipeline
                if engine is None:
                    return False
                if stage == "assignment":
                    return (
                        engine.collector._assignment is not None
                        and engine.collector._cursor >= 1
                    )

                if stage == "solver":
                    return (
                        engine.solver is not None and engine.solver.xors >= 2
                    )
                return engine.extractor is not None

            job.budget.cancelled = cancel

            result = job.run()

            self.assertEqual(result.reason, "cancelled", stage)
            checkpoint = json.loads(json.dumps(job.checkpoint()))

            restored = SSSJob.from_checkpoint(checkpoint, budget=allowance())

            result = restored.run()

            self.assertTrue(utils.valid_divisor(result.divisor, n))
            self.assertEqual(result.divisor * result.cofactor, n)
            self.assertGreater(result.stats["work_used"], job.budget.used)

    def test_checkpoint_rejects_rehashed_arithmetic_and_progress(self):
        job = SSSJob(
            4001 * 4003, config=SSSConfig(base_bound=400), budget=allowance()
        )
        job.run(batch_limit=1)
        checkpoint = job.checkpoint()

        for field in ("atom", "base", "cursor", "stats"):
            payload = json.loads(checkpoint["blob"])
            if field == "atom":
                payload["store"]["atoms"][0][2] ^= 1
            elif field == "base":
                payload["base_identity"] = "wrong"
            elif field == "cursor":
                payload["engine"]["next_position"] = 1000000
            else:
                payload["engine"]["stats"]["stage_seconds"] = []

            with self.assertRaises(ValueError):
                SSSJob.from_checkpoint(
                    reseal(checkpoint, payload), budget=allowance()
                )

        with self.assertRaises(ValueError):
            SSSJob.from_checkpoint(checkpoint, budget=allowance(0))
        with self.assertRaises(ValueError):
            SSSJob.from_checkpoint(
                checkpoint,
                budget=allowance(),
                config=SSSConfig(base_bound=500),
            )

        small = SSSJob(
            4001 * 5003,
            config=SSSConfig(base_bound=400, checkpoint_bytes=4096),
            budget=allowance(),
        )
        small.run(batch_limit=1)
        with self.assertRaises(MemoryError):
            small.checkpoint()

        self.assertIsNotNone(small.run().divisor)
        with self.assertRaises(MemoryError):
            SSSJob(4001 * 5003, config=SSSConfig(memory_bytes=1)).checkpoint()

    def test_version_three_checkpoint_without_sss_remains_readable(self):
        config = PortfolioConfig(rho_attempts=0, pm1_attempts=0, ecm_tiers=())

        first = factorize_bounded(
            1009 * 1013, config=config, budget=allowance(0)
        )
        checkpoint = copy.deepcopy(first.checkpoint)
        checkpoint["payload"]["version"] = 3
        checkpoint["payload"]["config"].pop("sss")
        encoded = json.dumps(
            checkpoint["payload"], sort_keys=True, separators=(",", ":")
        )
        checkpoint["sha256"] = hashlib.sha256(encoded.encode()).hexdigest()

        resumed = factorize_bounded(
            1009 * 1013,
            config=config,
            checkpoint=checkpoint,
            budget=allowance(),
        )

        self.assertTrue(resumed.result.complete)
        self.assertEqual(resumed.checkpoint["payload"]["version"], 5)

    def test_cli_selection_full_output_and_checkpoint_resume(self):
        n = 100003 * 100019

        for mode in ("sss", "sssf"):
            command = [
                sys.executable,
                "-B",
                "-m",
                "v2.factor",
                str(n),
                "--method",
                mode,
                "--sss-base-bound",
                "400",
            ]
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "progress.json"

                first = subprocess.run(
                    command
                    + ["--work-limit", "100000", "--checkpoint", str(path)],
                    cwd=ROOT,
                    capture_output=True,
                    text=True,
                )

                self.assertEqual(first.returncode, 1, first.stderr)
                self.assertIn("unresolved", first.stdout)

                resumed = subprocess.run(
                    command + ["--resume", str(path)],
                    cwd=ROOT,
                    capture_output=True,
                    text=True,
                )

                self.assertEqual(resumed.returncode, 0, resumed.stderr)
                self.assertIn("100003^1 * 100019^1", resumed.stdout)


if __name__ == "__main__":
    unittest.main()
