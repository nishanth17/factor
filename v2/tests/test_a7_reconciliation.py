"""R5 forced-factor recovery, historical migration and cumulative resume."""

import importlib
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.infrastructure.performance.performance_audit import (
    materialize_baseline,
)
from v2.benchmarks.qs.phase_three.phase_three_sss_upstream import OutputCheck
from v2.execution.budget import BudgetExhaustedError
from v2.qs.checkpoint import _solver_digest
from v2.qs.linear_algebra import DependencySolver, FilteredMatrix
from v2.qs.sss import SSSConfig, SSSJob
from v2.tests.test_sss import search_collector, trial_residual
from v2.tests.test_sss_dispatch import ROOT, allowance, reseal


class ReconciliationTests(unittest.TestCase):
    def test_upstream_stdout_claims_require_an_exact_proper_split(self):
        for claim in ("1 | 77", "77 | 1", "7 | 13"):
            output = OutputCheck(77)
            with self.assertRaises(AssertionError):
                output.write("Proper factors found: " + claim)

        output = OutputCheck(77)
        output.write("SSS finished: 1000000 relations found")
        self.assertIsNone(output.divisor)
        output.write("Proper factors found: 7 | 11")
        self.assertEqual(output.divisor, 7)

    def test_forced_quotient_recovers_every_original_valuation(self):
        for mode in ("sss", "sssf"):
            collector = search_collector(mode=mode)
            quotients = dict(collector.assignment(0))

            result = collector.collect(0, 1)

            self.assertTrue(result.atoms)
            forced_seen = repeated_seen = False
            for atom in result.atoms:
                value = abs(atom.polynomial.value(atom.position))
                quotient = quotients[atom.position]
                self.assertEqual(value % quotient, 0)
                forced_seen |= value // quotient > 1
                remaining, expected = value, []
                for prime in collector.factor_base.primes:
                    exponent = 0
                    while remaining % prime == 0:
                        remaining //= prime
                        exponent += 1
                    if exponent:
                        expected.append((prime, exponent))
                        repeated_seen |= exponent > 1

                self.assertEqual(atom.exponents, tuple(expected))
                self.assertEqual(atom.residual, remaining)
                self.assertEqual(
                    trial_residual(value, collector.factor_base.primes),
                    remaining,
                )

            self.assertTrue(forced_seen and repeated_seen)

    def test_in_memory_budget_replacement_is_rejected_in_both_modes(self):
        for mode in ("sss", "sssf"):
            job = SSSJob(
                4001 * 4003, config=SSSConfig(mode=mode), budget=allowance(0)
            )
            self.assertEqual(job.run().reason, "work_limit")
            original = job.budget
            job.budget = allowance()
            with self.assertRaisesRegex(ValueError, "original Budget"):
                job.run()

            job.budget = original
            original.work_limit = 10**9
            result = job.run()
            self.assertEqual(result.divisor * result.cofactor, job.n)

    def test_restore_charges_reconstruction_and_retains_all_resources(self):
        for mode in ("sss", "sssf"):
            job = SSSJob(
                4001 * 4003,
                config=SSSConfig(mode=mode, base_bound=400),
                budget=allowance(),
            )
            job.run(batch_limit=1)
            checkpoint = job.checkpoint()
            resources = checkpoint["resources"]
            refused = allowance(resources["work_used"])

            with self.assertRaises(BudgetExhaustedError):
                SSSJob.from_checkpoint(checkpoint, budget=refused)
            self.assertEqual(refused.used, resources["work_used"])
            self.assertEqual(refused.reason, "work_limit")

            budget = allowance()
            restored = SSSJob.from_checkpoint(checkpoint, budget=budget)

            self.assertGreater(budget.used, resources["work_used"])
            self.assertGreaterEqual(budget.wall_used, resources["wall_used"])
            self.assertGreaterEqual(budget.cpu_used, resources["cpu_used"])
            second = restored.checkpoint()
            again = SSSJob.from_checkpoint(second, budget=allowance())
            self.assertGreater(again.budget.used, budget.used)
            result = again.run()
            self.assertEqual(result.divisor * result.cofactor, job.n)

    def test_terminal_evidence_round_trips_without_restarting_search(self):
        for mode in ("sss", "sssf"):
            for n, config, reason in (
                (10403, SSSConfig(mode=mode), "factor_found"),
                (
                    1009 * 1013,
                    SSSConfig(mode=mode, base_bound=200, search_rounds=1),
                    "search_exhausted",
                ),
                (
                    1009 * 1013,
                    SSSConfig(mode=mode, memory_bytes=3 * 2**20),
                    "memory_limit",
                ),
            ):
                job = SSSJob(n, config=config, budget=allowance())
                before = job.run()
                self.assertEqual(before.reason, reason)
                restored = SSSJob.from_checkpoint(
                    job.checkpoint(), budget=allowance()
                )
                after = restored.run()
                self.assertEqual(after.reason, before.reason)
                self.assertEqual(after.divisor, before.divisor)
                self.assertEqual(after.cofactor, before.cofactor)
                used = restored.budget.used
                self.assertEqual(restored.run().stats["work_used"], used)

    def test_real_version_one_and_version_two_migration(self):
        for version, filename in (
            (1, "p38_r3_baseline.json"),
            (2, "performance_followup_baseline.json"),
        ):
            with tempfile.TemporaryDirectory(prefix="a7-legacy-") as directory:
                materialize_baseline(
                    ROOT / "v2/benchmarks/inputs/baselines" / filename,
                    Path(directory),
                )
                package = types.ModuleType("_a7_legacy_" + str(version))
                package.__path__ = [str(Path(directory) / "v2")]
                sys.modules[package.__name__] = package
                old = importlib.import_module(package.__name__ + ".qs.sss")
                budgets = importlib.import_module(package.__name__ + ".budget")

                for mode in ("sss", "sssf"):
                    job = old.SSSJob(
                        100003 * 100019,
                        config=old.SSSConfig(mode=mode, base_bound=400),
                        budget=budgets.Budget(
                            work_limit=10**9, seconds=None, cpu_seconds=None
                        ),
                    )
                    job.run(batch_limit=1)
                    checkpoint = job.checkpoint()
                    self.assertEqual(checkpoint["version"], version)

                    restored = SSSJob.from_checkpoint(
                        checkpoint, budget=allowance()
                    )

                    self.assertEqual(restored.checkpoint()["version"], 3)
                    result = restored.run()
                    self.assertEqual(result.divisor * result.cofactor, job.n)

    def test_sss_checkpoint_hashes_wide_valid_lifted_masks(self):
        # Isolate packing's solver fingerprint from relation-store encoding.
        # A valid original-row dependency can exceed PyPy's decimal limit.
        count = 16000
        matrix = FilteredMatrix(
            original_rows=(1,) * count,
            rows=(1, 1),
            masks=(1 << (count - 2), 1 << (count - 1)),
            zero_dependencies=(),
            stats={},
            workspace_bytes=8 * 2**20,
        )
        solver = DependencySolver(matrix, budget=allowance())
        self.assertEqual(solver.run(), (3 << (count - 2),))
        job = SSSJob(4001 * 4003, budget=allowance())
        job._setup()
        job.pipeline.solver = solver

        checkpoint = job.checkpoint()

        prefix = json.loads(checkpoint["blob"])["engine"]["solver"]
        self.assertEqual(prefix["digest_encoding"], "hex-v1")
        self.assertEqual(
            prefix["digest"], _solver_digest(solver, encoding="hex-v1")
        )

    def test_sss_legacy_decimal_solver_fingerprint_stays_readable(self):
        for mode in ("sss", "sssf"):
            job = SSSJob(
                4001 * 4003,
                config=SSSConfig(mode=mode, base_bound=400),
                budget=allowance(),
            )
            with patch(
                "v2.qs.pipeline.DependencySolver.run",
                side_effect=BudgetExhaustedError("work_limit"),
            ):
                job.budget.reason = "work_limit"
                self.assertEqual(job.run().reason, "work_limit")
            checkpoint = job.checkpoint()
            payload = json.loads(checkpoint["blob"])
            prefix = payload["engine"]["solver"]
            self.assertEqual(prefix.pop("digest_encoding"), "hex-v1")
            prefix["digest"] = _solver_digest(job.pipeline.solver)

            restored = SSSJob.from_checkpoint(
                reseal(checkpoint, payload), budget=allowance()
            )

            result = restored.run()
            self.assertEqual(result.divisor * result.cofactor, job.n)
            prefix["digest_encoding"] = "unknown"
            with self.assertRaisesRegex(ValueError, "encoding"):
                SSSJob.from_checkpoint(
                    reseal(checkpoint, payload), budget=allowance()
                )

    def test_backend_identity_and_cumulative_deadlines_are_checked(self):
        job = SSSJob(4001 * 4003, budget=allowance())
        job.run(batch_limit=1)
        checkpoint = job.checkpoint()
        payload = json.loads(checkpoint["blob"])
        payload["backend"] = "unknown-build"
        with self.assertRaisesRegex(ValueError, "backend"):
            SSSJob.from_checkpoint(
                reseal(checkpoint, payload), budget=allowance()
            )

        for field, reason in (
            ("seconds", "wall_limit"),
            ("cpu_seconds", "cpu_limit"),
        ):
            budget = allowance()
            setattr(budget, field, 0)
            with self.assertRaises(BudgetExhaustedError):
                SSSJob.from_checkpoint(checkpoint, budget=budget)
            self.assertEqual(budget.reason, reason)
            self.assertEqual(budget.used, checkpoint["resources"]["work_used"])


if __name__ == "__main__":
    unittest.main()
