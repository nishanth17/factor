"""Wide dependency fingerprints and legacy SIQS checkpoint compatibility."""

import hashlib
import json
import unittest
from unittest.mock import patch

from v2.budget import BudgetExhaustedError
from v2.qs import SIQSJob
from v2.qs.checkpoint import _solver_digest
from v2.qs.linear_algebra import DependencySolver, FilteredMatrix
from v2.tests.test_siqs import allowance, configuration, mutate


class CheckpointMaskTests(unittest.TestCase):
    """Fingerprint binary masks without relying on decimal digit limits."""

    def test_wide_valid_dependency_digest_and_corruption(self):
        count = 16000
        rows = (1,) * count
        matrix = FilteredMatrix(
            original_rows=rows,
            rows=(1, 1),
            masks=(1 << (count - 2), 1 << (count - 1)),
            zero_dependencies=(),
            stats={},
            workspace_bytes=8 * 1024**2,
        )
        solver = DependencySolver(matrix, budget=allowance())

        dependencies = solver.run()

        self.assertEqual(dependencies, ((3 << (count - 2)),))
        digest = _solver_digest(solver, encoding="hex-v1")

        self.assertEqual(len(digest), 64)
        self.assertEqual(digest, _solver_digest(solver, encoding="hex-v1"))
        solver.xors += 1

        self.assertNotEqual(digest, _solver_digest(solver, encoding="hex-v1"))
        with self.assertRaises(ValueError):
            _solver_digest(solver, encoding="unknown")

    def test_streamed_hex_digest_matches_prior_canonical_encoding(self):
        matrix = FilteredMatrix(
            original_rows=(3, 5, 6),
            rows=(3, 5, 6),
            masks=(1, 2, 4),
            zero_dependencies=(),
            stats={},
            workspace_bytes=8192,
        )
        solver = DependencySolver(matrix, budget=allowance())
        solver.run()
        state = [
            solver.next_row,
            solver.xors,
            solver.pending,
            sorted(solver.pivots.items()),
            solver.dependencies,
        ]

        def encode(value):
            if type(value) is int:
                return hex(value)
            if value is None:
                return None
            return [encode(item) for item in value]

        expected = hashlib.sha256(
            json.dumps(
                encode(state), sort_keys=True, separators=(",", ":")
            ).encode()
        ).hexdigest()

        self.assertEqual(_solver_digest(solver, encoding="hex-v1"), expected)

    def test_legacy_decimal_solver_prefix_remains_readable(self):
        job = SIQSJob(4001 * 5003, config=configuration(), budget=allowance())
        with patch(
            "v2.qs.pipeline.DependencySolver.run",
            side_effect=BudgetExhaustedError("work_limit"),
        ):
            job.budget.reason = "work_limit"

            self.assertEqual(job.run().reason, "work_limit")

        checkpoint = job.checkpoint()
        legacy_digest = _solver_digest(job.engine.solver)

        def legacy(payload):
            prefix = payload["engine"]["solver"]

            self.assertEqual(prefix.pop("digest_encoding"), "hex-v1")
            prefix["digest"] = legacy_digest

        restored = SIQSJob.from_checkpoint(
            mutate(checkpoint, legacy), budget=allowance()
        )

        result = restored.run()

        self.assertEqual(result.reason, "factor_found")
        self.assertEqual({result.divisor, result.cofactor}, {4001, 5003})
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                mutate(
                    checkpoint,
                    lambda payload: payload["engine"]["solver"].update(
                        digest_encoding="unknown"
                    ),
                ),
                budget=allowance(),
            )


if __name__ == "__main__":
    unittest.main()
