"""Internal collection status preserves traversal and snapshot ownership."""

import gc
import json
import subprocess
import sys
import unittest
import weakref
from dataclasses import replace
from functools import partial
from unittest.mock import PropertyMock, patch

from v2.qs import (
    BoundedYieldSieveConfig,
    GrowingDoubleLargeSieveConfig,
    QSJob,
    SieveCollector,
    SieveConfig,
    build_factor_base,
    qs_polynomial,
)
from v2.qs.large_prime_checkpoint import pack_store
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_qs import unlimited_budget


class CollectionProgressTests(unittest.TestCase):
    """Compare independently owned stores through completion and refusal."""

    def make_worker(self, config):
        base = build_factor_base(104729, bound=40).factor_base
        return SieveCollector(
            qs_polynomial(base),
            base,
            config=config,
            budget=unlimited_budget(),
        )

    def assert_progress_equal(self, snapshot, progress):
        for name in (
            "divisor",
            "next_position",
            "reason",
            "stats",
            "workspace_bytes",
        ):
            self.assertEqual(getattr(snapshot, name), getattr(progress, name))

    def assert_store_equal(self, first, second):
        if first._graph is not None:
            self.assertEqual(
                json.loads(json.dumps(pack_store(first))),
                json.loads(json.dumps(pack_store(second))),
            )
        else:
            self.assertEqual(first._atoms, second._atoms)
            self.assertEqual(first._full, second._full)
            self.assertEqual(first._combined, second._combined)
            self.assertEqual(first._pending, second._pending)
        self.assertEqual(first.partial_ids, second.partial_ids)
        self.assertEqual(first._workspace, second._workspace)
        self.assertEqual(first.budget.used, second.budget.used)

    def test_signed_batches_preserve_stats_work_and_serialized_store(self):
        common = dict(residual_bound=500, block_width=17, division="full")
        configs = (
            SieveConfig(max_partials=2, **common),
            GrowingDoubleLargeSieveConfig(
                large_prime_bound=500, large_product_bound=500**2, **common
            ),
            BoundedYieldSieveConfig(
                large_prime_bound=500, large_product_bound=500**2, **common
            ),
        )
        for config in configs:
            with self.subTest(config=type(config).__name__):
                first, second = (self.make_worker(config) for _ in range(2))
                for lo, hi in ((-128, -31), (-31, 34), (34, 129)):
                    snapshot = first.collect(lo, hi)
                    progress = second._collect_progress(lo, hi)

                    self.assertEqual(snapshot.reason, "complete")
                    self.assert_progress_equal(snapshot, progress)
                    self.assert_store_equal(first, second)

                self.assertTrue(first._atoms)
                self.assertTrue(first._combined)

    def test_direct_divisor_and_empty_window_keep_progress(self):
        for n, lo, hi in ((101**2, 0, 1), (101 * 103, -1, 0)):
            base = build_factor_base(n, bound=40).factor_base
            first, second = (
                SieveCollector(
                    qs_polynomial(base),
                    base,
                    config=SieveConfig(residual_bound=500),
                    budget=unlimited_budget(),
                )
                for _ in range(2)
            )
            self.assert_progress_equal(
                first.collect(lo, lo), second._collect_progress(lo, lo)
            )
            snapshot = first.collect(lo, hi)
            progress = second._collect_progress(lo, hi)

            self.assertEqual(snapshot.divisor, 101)
            self.assertEqual(snapshot.reason, "factor_found")
            self.assertEqual(n % progress.divisor, 0)
            self.assert_progress_equal(snapshot, progress)
            self.assert_store_equal(first, second)

    def test_refusal_and_resume_preserve_cursor_state_and_charges(self):
        config = SieveConfig(residual_bound=97, block_width=17)
        cases = (
            ("work_limit", lambda: Budget(work_limit=200)),
            ("cancelled", lambda: Budget(cancelled=lambda: True)),
            ("wall_limit", lambda: Budget(seconds=0)),
            ("cpu_limit", lambda: Budget(cpu_seconds=0)),
            ("atom_limit", None),
            ("relation_limit", None),
            ("memory_limit", None),
        )
        for reason, budget in cases:
            with self.subTest(reason=reason):
                first, second = (self.make_worker(config) for _ in range(2))
                for worker in (first, second):
                    if budget is not None:
                        worker.budget = budget()
                    elif reason == "memory_limit":
                        worker.config = replace(config)
                    else:
                        worker.config = replace(
                            config,
                            **{
                                "max_atoms"
                                if reason == "atom_limit"
                                else "max_relations": 0
                            },
                        )

                snapshot = first.collect(-128, 129)
                progress = second._collect_progress(-128, 129)

                self.assertEqual(snapshot.reason, reason)
                self.assert_progress_equal(snapshot, progress)
                self.assert_store_equal(first, second)

                for worker in (first, second):
                    worker.config = config
                    worker.budget = Budget(
                        work_limit=200_000_000,
                        used=worker.budget.used,
                        seconds=None,
                        cpu_seconds=None,
                    )
                resumed = first.collect(snapshot.next_position, 129)
                progress = second._collect_progress(
                    progress.next_position, 129
                )

                self.assertEqual(resumed.reason, "complete")
                self.assert_progress_equal(resumed, progress)
                self.assert_store_equal(first, second)

    def test_default_factory_avoids_snapshot_expansion_in_pipeline(self):
        base = build_factor_base(104729, bound=100).factor_base
        job = QSJob(
            qs_polynomial(base),
            base,
            -64,
            65,
            batch_width=16,
            config=SieveConfig(residual_bound=1),
            collector_class=partial(SieveCollector, precomputed_roots=None),
            budget=unlimited_budget(),
        )
        with patch.object(
            SieveCollector,
            "partial_ids",
            new_callable=PropertyMock,
            side_effect=AssertionError("snapshot expansion on internal path"),
        ):
            result = job.run()

        self.assertEqual(result.reason, "window_exhausted")
        self.assertEqual(result.next_position, 65)
        self.assertTrue(job.collector._atoms)
        snapshot = job.collector.collect(65, 65)
        self.assertEqual(snapshot.atoms, tuple(job.collector._atoms.values()))
        self.assertEqual(snapshot.partial_ids, job.collector.partial_ids)

    def test_progress_does_not_pin_store_but_public_snapshot_does(self):
        # JIT frames may retain deleted temporaries. Check Python ownership in
        # a child with JIT off, without altering the suite's runtime settings.
        command = (
            "from v2.tests.test_qs_collection_progress import "
            "CollectionProgressTests; "
            "CollectionProgressTests()._check_result_ownership()"
        )
        process = subprocess.run(
            [sys.executable, "--jit", "off", "-B", "-c", command],
            capture_output=True,
            text=True,
            timeout=10,
        )

        self.assertEqual(
            process.returncode, 0, process.stdout + process.stderr
        )

    def _check_result_ownership(self):
        config = SieveConfig(max_partials=2, residual_bound=97)
        worker = self.make_worker(config)
        progress = worker._collect_progress(-128, 129)
        atom = weakref.ref(next(iter(worker._atoms.values())))
        del worker
        gc.collect()

        self.assertIsNone(atom())
        self.assertEqual(progress.reason, "complete")

        worker = self.make_worker(config)
        snapshot = worker.collect(-128, 0)
        atoms, partials = snapshot.atoms, snapshot.partial_ids
        worker._collect_progress(0, 129)
        atom = weakref.ref(snapshot.atoms[0])
        del worker
        gc.collect()

        self.assertIsNotNone(atom())
        self.assertEqual(snapshot.atoms, atoms)
        self.assertEqual(snapshot.partial_ids, partials)
