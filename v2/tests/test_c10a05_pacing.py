"""Independent arithmetic witnesses and charged collection-boundary replay."""

import json
import unittest
from collections import Counter
from dataclasses import replace
from math import gcd, isqrt
from unittest.mock import patch

from v2.execution.budget import BudgetExhaustedError
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs import (
    AtomicRelation,
    DependencySolver,
    DoubleLargeSieveConfig,
    GrowingDoubleLargeSieveConfig,
    QSJob,
    SieveConfig,
    SIQSConfig,
    SIQSJob,
    build_factor_base,
    extraction,
    filter_matrix,
    qs_polynomial,
    verify_atomic,
)
from v2.tests.test_qs_pipeline import dense_kernel, span
from v2.tests.test_siqs import allowance, mutate


def witness_job(*, paced=True, maximum=64):
    """139²-N and 251²-N have equal parity and give gcd 101."""
    base = build_factor_base(10403, bound=40).factor_base
    polynomial = qs_polynomial(base)
    job = QSJob(
        polynomial,
        base,
        200,
        201,
        budget=allowance(),
        config=SieveConfig(max_relations=maximum),
        row_excess=32,
        filter_row_growth=32,
        weight_two=False,
        collection_continues=paced,
    )
    atoms = (
        AtomicRelation(polynomial, 37, 1, ((2, 1), (7, 3), (13, 1))),
        AtomicRelation(polynomial, 149, 1, ((2, 1), (7, 1), (13, 1), (17, 2))),
    )
    for atom in atoms:
        verify_atomic(atom, base, budget=job.budget)
        assert job.collector._admit(atom, Counter()) is None
    return job


def terminal_job():
    config = SIQSConfig(
        base_bound=40,
        half_width=1,
        max_half_width=1,
        factor_count=1,
        family_count=1,
        pool_size=2,
        solve_policy="collection",
        row_excess=32,
        filter_row_growth=32,
        weight_two=False,
        collector=SieveConfig(max_relations=64),
    )
    job = SIQSJob(10403, config=config, budget=allowance())
    job._setup()
    job.engine = witness_job()
    job.engine.budget = job.budget
    job._stop_collection("families_exhausted")
    return job


class PacingTests(unittest.TestCase):
    def test_witness_and_kernel_are_independent_of_schedule(self):
        job = witness_job()
        x = 139 * 251
        product = (139**2 - 10403) * (251**2 - 10403)
        y = isqrt(product)
        self.assertEqual(y * y, product)
        self.assertEqual(gcd(x - y, 10403), 101)
        self.assertEqual(dense_kernel((14, 14)), {0, 3})
        for merge in (False, True):
            matrix = filter_matrix((14, 14), weight_two=merge)
            self.assertEqual(span(DependencySolver(matrix).run()), {0, 3})

        # A prefix may have useful dependencies despite surplus below 32.
        self.assertIsNone(job._solve(False))
        self.assertEqual(len(job.collector._rows), 2)
        self.assertEqual(job._solve(True), 101)

    def test_window_deferral_and_direct_qs_finalization(self):
        for paced in (True, False):
            job = witness_job(paced=paced)
            job.last_count = 1
            job.next_position = job.hi
            result = job.run()
            self.assertEqual(result.divisor, None if paced else 101)
            self.assertTrue(job.final_solve_done)
            if paced:
                before = job.budget.used
                self.assertEqual(job.run().reason, "window_exhausted")
                self.assertEqual(job.budget.used, before)
                self.assertEqual(job._solve(True), 101)

    def test_square_bypasses_cadence_and_checks_proper_divisor(self):
        base = build_factor_base(77, bound=5).factor_base
        poly = qs_polynomial(base)
        job = QSJob(
            poly,
            base,
            0,
            8,
            budget=allowance(),
            collection_continues=True,
            row_excess=32,
            filter_row_growth=32,
        )
        job.last_count = 0
        atom = AtomicRelation(poly, 0, 1, ((2, 2),))
        self.assertEqual(poly.u_value(0), 9)
        self.assertEqual(gcd(9 - 2, 77), 7)
        verify_atomic(atom, base, budget=job.budget)
        job.collector._admit(atom, Counter())
        self.assertEqual(job._solve(False), 7)

    def test_singletons_are_retained_for_a_later_partner(self):
        job = witness_job()
        second = job.collector._rows.pop()
        job.collector._full.pop()
        job.collector._atoms.pop(second.relation_id)
        self.assertIsNone(job._solve(False))
        self.assertEqual(job.stats["last_filter"]["singletons_removed"], 1)
        self.assertEqual(len(job.collector._rows), 1)
        job.collector._admit(second, Counter())
        self.assertEqual(job._solve(True), 101)

    def test_storage_cap_forces_the_retained_kernel(self):
        for reason in ("relation_limit", "atom_limit", "memory_limit"):
            job = witness_job()
            job.last_count = 1
            job.storage_reason = reason
            result = job.run()
            self.assertEqual(result.divisor, 101)
            self.assertTrue(job.storage_solve_done)

    def test_terminal_reasons_drain_once_even_below_excess(self):
        for reason in (
            "families_exhausted",
            "assignment_space_exhausted",
            "stalled_yield",
            "trivial_dependency_limit",
            "coefficient_limit",
        ):
            job = terminal_job()
            job.pending_final_reason = reason
            result = job.run()
            self.assertEqual(result.divisor, 101)
            self.assertIsNone(job.pending_final_reason)
            before = job.budget.used
            self.assertEqual(job.run().divisor, 101)
            self.assertEqual(job.budget.used, before)

    def test_refused_terminal_work_and_serialized_repeated_resume(self):
        for kind in ("work", "wall", "cpu", "cancel"):
            job = terminal_job()
            budget = job.budget
            if kind == "work":
                budget.work_limit = budget.used
            elif kind == "wall":
                budget.seconds = 0
            elif kind == "cpu":
                budget.cpu_seconds = 0
            else:
                budget.cancelled = lambda: True
            result = job.run()
            self.assertIsNone(result.divisor)
            self.assertIsNone(job.finished_reason)
            self.assertEqual(job.pending_final_reason, "families_exhausted")
            for _ in range(3):
                checkpoint = json.loads(json.dumps(job.checkpoint()))
                used = checkpoint["resources"]["work_used"]
                job = SIQSJob.from_checkpoint(checkpoint, budget=allowance())
                self.assertGreater(job.budget.used, used)
                self.assertGreaterEqual(
                    job.budget.prior_wall, checkpoint["resources"]["wall_used"]
                )
                self.assertEqual(job.run(max_blocks=0).reason, "paused")
            self.assertEqual(job.run().divisor, 101)

    def test_interrupted_stages_preserve_pending_terminal_work(self):
        for stage in ("preparation", "filtering", "elimination", "extraction"):
            job = terminal_job()
            original = job.engine._timed

            def interrupt(name, call):
                if name == stage:
                    job.budget.reason = "work_limit"
                    raise BudgetExhaustedError("work_limit")
                return original(name, call)

            with patch.object(job.engine, "_timed", interrupt):
                self.assertEqual(job.run().reason, "work_limit")
            checkpoint = json.loads(json.dumps(job.checkpoint()))
            restored = SIQSJob.from_checkpoint(checkpoint, budget=allowance())
            self.assertEqual(
                restored.pending_final_reason, "families_exhausted"
            )
            self.assertEqual(restored.run().divisor, 101)
            self.assertGreater(restored.budget.used, job.budget.used)

    def test_real_polynomial_pauses_and_schema_compatibility(self):
        config = SIQSConfig(
            base_bound=40,
            half_width=16,
            max_half_width=16,
            factor_count=2,
            family_count=4,
            pool_size=6,
            row_excess=32,
            filter_row_growth=32,
            solve_policy="collection",
        )
        direct = SIQSJob(10403, config=config, budget=allowance()).run()
        job = SIQSJob(10403, config=config, budget=allowance())
        for _ in range(100):
            result = job.run(max_blocks=1)
            if result.reason != "paused":
                break
            cp = json.loads(json.dumps(job.checkpoint()))
            self.assertEqual(cp["version"], 6)
            job = SIQSJob.from_checkpoint(cp, budget=allowance())
        self.assertEqual(result.divisor, direct.divisor)
        self.assertEqual((result.divisor or 1) * result.cofactor, 10403)
        legacy = SIQSJob(
            10403,
            config=replace(config, solve_policy="window"),
            budget=allowance(),
        )
        legacy.run(max_blocks=1)
        cp = legacy.checkpoint()
        self.assertEqual(cp["version"], 3)
        self.assertNotIn("solve_policy", json.loads(cp["blob"])["config"])
        self.assertEqual(
            SIQSJob.from_checkpoint(
                cp, budget=allowance()
            ).config.solve_policy,
            "window",
        )

    def test_corrupt_scheduling_state_and_changed_policy_are_rejected(self):
        job = terminal_job()
        checkpoint = job.checkpoint()
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                mutate(
                    checkpoint,
                    lambda p: p["engine"].update(zero_checked_count=3),
                ),
                budget=allowance(),
            )
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                checkpoint,
                budget=allowance(),
                allow_extension=True,
                config=replace(job.config, solve_policy="window"),
            )

    def test_partial_solver_and_trivial_extractor_resume(self):
        job = terminal_job()
        original_run = DependencySolver.run

        def partial(solver):
            solver.step()
            job.budget.reason = "work_limit"
            raise BudgetExhaustedError("work_limit")

        with patch.object(DependencySolver, "run", partial):
            self.assertEqual(job.run().reason, "work_limit")
        self.assertEqual(job.engine.solver.next_row, 1)
        restored = SIQSJob.from_checkpoint(
            job.checkpoint(), budget=allowance()
        )
        self.assertIs(DependencySolver.run, original_run)
        self.assertEqual(restored.run().divisor, 101)

        job = terminal_job()
        collector = job.engine.collector
        # Opposite U values make a known trivial square; the third atom
        # supplies the independent proper-factor witness from the first test.
        first, second = collector._rows
        collector._rows.pop()
        collector._full.pop()
        collector._atoms.pop(second.relation_id)
        negative = AtomicRelation(first.polynomial, -241, 1, first.exponents)
        verify_atomic(negative, collector.factor_base, budget=job.budget)
        collector._admit(negative, Counter())
        collector._admit(second, Counter())
        original = extraction.extract_dependency
        calls = []

        def interrupted(prepared, mask, **options):
            if calls:
                job.budget.reason = "work_limit"
                raise BudgetExhaustedError("work_limit")
            result = original(prepared, mask, **options)
            calls.append(result)
            return result

        with patch.object(extraction, "extract_dependency", interrupted):
            self.assertEqual(job.run().reason, "work_limit")
        self.assertIsNone(calls[0].divisor)
        self.assertEqual(job.engine.extractor.next_dependency, 1)
        cp = json.loads(json.dumps(job.checkpoint()))
        restored = SIQSJob.from_checkpoint(cp, budget=allowance())
        self.assertEqual(restored.engine.extractor.next_dependency, 1)
        self.assertEqual(restored.run().divisor, 101)

    def test_trivial_square_is_not_repeated_and_search_stops_drain(self):
        base = build_factor_base(97, bound=5).factor_base
        poly = qs_polynomial(base)
        job = QSJob(
            poly,
            base,
            0,
            1,
            budget=allowance(),
            collection_continues=True,
            filter_row_growth=32,
        )
        atom = AtomicRelation(poly, 39, 1, ((2, 8), (3, 2)))
        verify_atomic(atom, base, budget=job.budget)
        job.collector._admit(atom, Counter())
        self.assertIsNone(job._solve(False))
        self.assertEqual(job.stats["trivial_dependencies"], 1)
        before = job.budget.used
        self.assertIsNone(job._solve(False))
        self.assertEqual(job.budget.used, before)

        for stop in ("families", "stalled", "trivial"):
            job = terminal_job()
            job.pending_final_reason = None
            job.engine.last_count = 1
            job.engine.next_position = job.engine.hi
            if stop == "families":
                job.active = False
                job.family_index = job.config.family_count
                job.assignments = ()
            else:
                job.active = True
                job.poly_start_rows = 2
                if stop == "stalled":
                    job.stalled = job.config.max_stalled
                else:
                    job.engine.stats["trivial_dependencies"] = (
                        job.config.max_trivial
                    )
            self.assertEqual(job.run().divisor, 101)

    def test_new_schema_collector_types_and_legacy_portfolio_identity(self):
        for collector in (
            DoubleLargeSieveConfig(
                large_prime_bound=100, large_product_bound=10000
            ),
            GrowingDoubleLargeSieveConfig(
                large_prime_bound=100, large_product_bound=10000
            ),
        ):
            for policy in ("window", "collection"):
                config = SIQSConfig(solve_policy=policy, collector=collector)
                job = SIQSJob(4001 * 5003, config=config, budget=allowance())
                cp = json.loads(json.dumps(job.checkpoint()))
                restored = SIQSJob.from_checkpoint(cp, budget=allowance())
                self.assertEqual(restored.config, config)
                self.assertIs(type(restored.config.collector), type(collector))

        for policy in ("window", "collection"):
            config = PortfolioConfig(
                siqs=SIQSConfig(solve_policy=policy), memory_bytes=64 * 2**20
            )
            first = factorize_bounded(
                4001 * 5003, config=config, budget=allowance(0)
            )
            saved = first.checkpoint["payload"]["config"]["siqs"]
            self.assertEqual("solve_policy" in saved, policy == "collection")
            for explicit in (None, config):
                restored = factorize_bounded(
                    4001 * 5003,
                    checkpoint=first.checkpoint,
                    config=explicit,
                    budget=allowance(0),
                )
                self.assertEqual(restored.result.reconstruct(), 4001 * 5003)


if __name__ == "__main__":
    unittest.main()
