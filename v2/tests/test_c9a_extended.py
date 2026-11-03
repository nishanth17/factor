"""Independent arithmetic, kernel and fault oracles for the expanded arm."""

import copy
import json
import random
import tempfile
import unittest
from collections import Counter
from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from v2.execution.budget import BudgetExhaustedError
from v2.qs import (
    AtomicRelation,
    ExtendedYieldSieveConfig,
    SieveCollector,
    SIQSConfig,
    SIQSJob,
    build_factor_base,
    qs_polynomial,
    verify_atomic,
    verify_combined,
)
from v2.qs.cofactor_policy import SplitMemo, split_cofactor, split_residual
from v2.qs.cold_factors import SpillStore
from v2.qs.extended_checkpoint import pack_store, restore_store
from v2.qs.extraction import DependencyExtractor
from v2.qs.incremental_structure import IncrementalStructure
from v2.qs.indexed_relations import (
    IndexedFactors,
    decode_factors,
    encode_factors,
)
from v2.qs.lazy_relations import LazyCombinedRelation
from v2.qs.linear_algebra import DependencySolver, filter_matrix
from v2.qs.recovery_batch import BaseRecovery, PendingHints
from v2.qs.residual_kernel import ResidualKernel
from v2.qs.residual_relations import ResidualAtom
from v2.qs.yield_feedback import OutcomeTarget, ReportingController
from v2.tests.test_qs_dlp import allowance
from v2.tests.test_qs_pipeline import dense_kernel, span


def policy(**changes):
    values = dict(
        large_prime_bound=1000,
        large_product_bound=10**9,
        residual_bound=1000,
        memory_bytes=128 * 2**20,
        division="full",
        block_width=16,
        duplicate_policy="identity",
    )
    values.update(changes)
    return ExtendedYieldSieveConfig(**values)


def collector(**changes):
    base = build_factor_base(91, bound=3).factor_base
    return SieveCollector(
        qs_polynomial(base), base, config=policy(**changes), budget=allowance()
    )


def atom_at(worker, position, *, opaque=False):
    poly = worker.polynomial
    value = poly.a * poly.value(position)
    remaining = abs(value) // poly.square_coefficient**2
    factors = []
    for prime in worker.factor_base.primes:
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        if exponent:
            factors.append((prime, exponent))
    if remaining == 1:
        return AtomicRelation(
            poly, position, -1 if value < 0 else 1, tuple(factors)
        )
    if opaque:
        return ResidualAtom(
            poly,
            position,
            -1 if value < 0 else 1,
            tuple(factors),
            opaque=remaining,
        )
    primes, divisor = [], 2
    while divisor * divisor <= remaining:
        while remaining % divisor == 0:
            remaining //= divisor
            primes.append(divisor)
        divisor += 1
    if remaining > 1:
        primes.append(remaining)
    return ResidualAtom(
        poly,
        position,
        -1 if value < 0 else 1,
        tuple(factors),
        primes=tuple(primes),
    )


class ExpandedArithmeticTests(unittest.TestCase):
    def test_index_recovery_has_exact_tuple_interchange(self):
        worker = collector()
        base = worker.factor_base
        values = IndexedFactors(((1, 7),), base.primes)
        self.assertEqual(tuple(values), ((base.primes[0], 7),))
        encoded = encode_factors(values, base)
        self.assertEqual(
            tuple(decode_factors(encoded, base.primes, 1)), tuple(values)
        )
        with self.assertRaises(ValueError):
            IndexedFactors(((1, 1), (1, 2)), base.primes)

    def test_three_endpoints_and_opaque_contract_are_explicit(self):
        worker = collector()
        triple = atom_at(worker, 28)
        self.assertEqual(triple.large_primes, (3, 11, 41))
        with self.assertRaises(ValueError):
            verify_atomic(
                triple,
                worker.factor_base,
                residual_bound=1000,
                large_prime_bound=1000,
                large_product_bound=10**9,
            )
        for atom in (triple, atom_at(worker, 1, opaque=True)):
            self.assertTrue(
                verify_atomic(
                    atom,
                    worker.factor_base,
                    residual_bound=1000,
                    large_prime_bound=1000,
                    large_product_bound=10**9,
                    residual_arity=3,
                    opaque_pairs=True,
                )
            )

    def test_malformed_arity_is_rejected_before_product_arithmetic(self):
        with self.assertRaises(ValueError):
            policy(residual_arity=10**9)
        worker = collector()
        with self.assertRaises(ValueError):
            verify_atomic(
                atom_at(worker, 1, opaque=True),
                worker.factor_base,
                residual_arity=10**9,
                opaque_pairs=True,
            )

    def test_failed_split_memo_does_not_hide_a_new_proper_hint(self):
        worker = collector()
        memo = SplitMemo(4)
        memo.remember(
            (15, 1000, "portfolio-v2", 3, True), ((), "split_failed")
        )
        primes, _, _ = split_residual(
            15,
            1000,
            worker.budget,
            worker._checked,
            "portfolio-v2",
            hint=3,
            memo=memo,
        )
        self.assertEqual(primes, (3, 5))

    def test_lazy_collection_and_checkpoint_do_not_compute_roots(self):
        worker = collector()
        with patch.object(
            LazyCombinedRelation,
            "materialize",
            side_effect=AssertionError("eager arithmetic"),
        ):
            for position in (1, -21):
                self.assertIsNone(
                    worker._admit(atom_at(worker, position), Counter())
                )
            row = worker._rows[0]
            self.assertIsInstance(row, LazyCombinedRelation)
            worker._checked.row(
                row, worker._atoms, worker.budget, worker.config.memory_limit
            )
            payload = json.loads(json.dumps(pack_store(worker)))
        self.assertTrue(
            verify_combined(
                row, worker.factor_base, worker._atoms, budget=allowance()
            )
        )
        restored = collector()
        restore_store(payload, restored, restored.budget)
        self.assertEqual(pack_store(restored), payload)

    def test_lazy_arithmetic_refusal_leaves_no_partial_cache(self):
        worker = collector()
        for position in (1, -21):
            worker._admit(atom_at(worker, position), Counter())
        row = worker._rows[0]
        row.bind_budget(allowance(0))
        with self.assertRaises(BudgetExhaustedError):
            row.materialize()
        row.bind_budget(allowance())
        self.assertEqual(row.materialize(), row.materialize())

    def test_opaque_conversion_preserves_old_atoms_and_replays(self):
        worker = collector()
        worker._admit(atom_at(worker, 1, opaque=True), Counter())
        worker._admit(atom_at(worker, -21), Counter())
        self.assertEqual(len(worker._conversion_events), 1)
        self.assertTrue(worker._rows)
        for row in worker._rows:
            verify_combined(
                row, worker.factor_base, worker._atoms, budget=allowance()
            )
        payload = json.loads(json.dumps(pack_store(worker)))
        restored = collector()
        restore_store(payload, restored, restored.budget)
        self.assertEqual(pack_store(restored), payload)
        damaged = copy.deepcopy(payload)
        damaged["conversion_events"][0][2] = [3, 7]
        with self.assertRaises(ValueError):
            other = collector()
            restore_store(damaged, other, other.budget)

    def test_recompute_matches_retained_exact_factors(self):
        worker = collector(recompute_prime_limit=13)
        source = atom_at(worker, 1)
        worker._admit(source, Counter())
        stored = worker._atoms[source.relation_id]
        self.assertEqual(tuple(stored.exponents), source.exponents)
        verify_atomic(
            stored,
            worker.factor_base,
            residual_bound=1000,
            large_prime_bound=1000,
            large_product_bound=10**9,
            residual_arity=3,
            opaque_pairs=True,
        )

    def test_three_prime_split_and_memo_revalidate_proofs(self):
        worker = collector()
        memo = SplitMemo(2)
        factors, _, _ = split_residual(
            3 * 11 * 41,
            1000,
            allowance(),
            worker._checked,
            "portfolio-v2",
            memo=memo,
        )
        self.assertEqual(factors, (3, 11, 41))
        with patch(
            "v2.qs.cofactor_policy.split_cofactor",
            side_effect=AssertionError("repeat split"),
        ):
            self.assertEqual(
                split_residual(
                    3 * 11 * 41,
                    1000,
                    allowance(),
                    worker._checked,
                    "portfolio-v2",
                    memo=memo,
                )[0],
                factors,
            )

    def test_ecm_fallback_is_forced_after_both_failed_walks(self):
        worker = collector()
        with (
            patch("v2.qs.cofactor_policy.squfof", return_value=None),
            patch(
                "v2.qs.cofactor_policy.factorize_rho", return_value=None
            ) as rho,
            patch("v2.ecm.core.factorize_ecm", return_value=97) as ecm,
        ):
            factors, _, _ = split_cofactor(
                97 * 101, 1000, allowance(), worker._checked, "ecm-v2"
            )
        self.assertEqual(factors, (97, 101))
        self.assertEqual(rho.call_count, 2)
        self.assertEqual(ecm.call_count, 1)


class ExpandedKernelTests(unittest.TestCase):
    def test_incidence_kernel_matches_dense_oracle_for_mixed_arities(self):
        generator = random.Random(822)
        config = policy(cycle_reduction_choices=4)
        for _ in range(40):
            rows = tuple(generator.randrange(32) for _ in range(8))
            kernel, dependencies = ResidualKernel(), []
            for handle, bits in enumerate(rows):
                atom = SimpleNamespace(
                    residual_primes=tuple(
                        i + 2 for i in range(5) if bits & (1 << i)
                    ),
                    exponents=(),
                )
                plan = kernel.plan_atom(atom, handle, config, allowance())
                if plan.cycle is not None:
                    dependencies.append(sum(1 << i for i in plan.cycle))
                kernel.commit_atom(plan, atom, config.cycle_reduction_choices)
            self.assertEqual(span(dependencies), dense_kernel(rows))

    def test_sparse_high_weight_and_panel_solver_preserve_full_kernel(self):
        generator = random.Random(911)
        for _ in range(40):
            rows = tuple(generator.randrange(64) for _ in range(8))
            expected = dense_kernel(rows)
            for sparse in (False, True):
                matrix = filter_matrix(
                    rows,
                    weight_two=True,
                    merge_policy="fill-v2",
                    merge_weight=4,
                    merge_fill=2,
                    sparse_incidence=sparse,
                    panel_bits=3,
                    budget=allowance(),
                )
                for pivot in ("highest", "lowest"):
                    result = DependencySolver(
                        matrix, pivot=pivot, budget=allowance()
                    ).run()
                    self.assertEqual(span(result), expected)

    def test_incremental_singletons_revive_and_match_batch_kernel(self):
        structure = IncrementalStructure(2**20)
        self.assertEqual(structure.observe((1,), allowance()), 0)
        self.assertEqual(structure.observe((1, 1), allowance()), 3)
        rows = (1, 1, 6, 10, 12, 0)
        for count in range(2, len(rows) + 1):
            prefix = rows[:count]
            mask = structure.observe(prefix, allowance())
            filtered = filter_matrix(
                prefix, survivor_mask=mask, budget=allowance()
            )
            self.assertEqual(
                span(DependencySolver(filtered, budget=allowance()).run()),
                dense_kernel(prefix),
            )

    def test_incremental_refusal_falls_back_without_dropping_rows(self):
        structure = IncrementalStructure(0)
        self.assertIsNone(structure.observe((1, 1), allowance()))
        self.assertIsNone(structure.observe((1, 1, 0), allowance()))

    def test_panel_pause_replays_the_same_action_prefix(self):
        rows = (0b1101, 0b1010, 0b0111, 0b1110, 0b0011, 0b1001)
        matrix = filter_matrix(rows, panel_bits=3, budget=allowance())
        solver = DependencySolver(matrix, budget=allowance())
        for _ in range(3):
            solver.step()
        state = copy.deepcopy((solver.next_row, solver.pending, solver.pivots))
        solver.budget = allowance(0)
        with self.assertRaises(BudgetExhaustedError):
            solver.step()
        # Initial pending selection is permitted; no row/pivot action commits.
        self.assertEqual(
            (solver.next_row, solver.pivots), (state[0], state[2])
        )
        solver.budget = allowance()
        self.assertEqual(span(solver.run()), dense_kernel(rows))


class ExpandedResourceTests(unittest.TestCase):
    def test_schema13_resumes_pending_subset_solver_extractor_and_full_retry(
        self,
    ):
        for phase in ("solver", "extractor", "full_retry"):
            config = SIQSConfig(
                mode="qs",
                multiplier=1,
                base_bound=40,
                half_width=512,
                max_half_width=512,
                memory_bytes=128 * 2**20,
                solve_policy="collection",
                row_excess=0,
                filter_row_growth=1,
                collector=policy(),
            )
            job = SIQSJob(10403, config=config, budget=allowance())
            job._setup()
            self.assertTrue(job._activate())
            worker = job.engine.collector
            for position, factors in (
                (37, ((2, 1), (7, 3), (13, 1))),
                (149, ((2, 1), (7, 1), (13, 1), (17, 2))),
                (-241, ((2, 1), (7, 3), (13, 1))),
                (-353, ((2, 1), (7, 1), (13, 1), (17, 2))),
            ):
                atom = AtomicRelation(worker.polynomial, position, 1, factors)
                self.assertIsNone(worker._admit(atom, Counter()))
            job.engine.next_position = 150
            calls = []

            def refuse(state):
                calls.append(job.engine._selected_subset)
                if phase == "full_retry" and len(calls) == 1:
                    return None
                job.budget.reason = "work_limit"
                raise BudgetExhaustedError("work_limit")

            cls = (
                DependencySolver if phase == "solver" else DependencyExtractor
            )
            with patch.object(cls, "run", refuse):
                with self.assertRaises(BudgetExhaustedError):
                    job.engine._solve(False)
            self.assertEqual(
                calls, [True, False] if phase == "full_retry" else [True]
            )
            checkpoint = json.loads(json.dumps(job.checkpoint()))
            restored = SIQSJob.from_checkpoint(checkpoint, budget=allowance())
            self.assertIsNotNone(restored.engine.solver)
            divisor = restored.engine._solve(False)
            if divisor is not None:
                self.assertTrue(1 < divisor < 10403 and 10403 % divisor == 0)
            self.assertFalse(restored.engine._retry_full_view)

    def test_piecewise_refinement_never_rejects_a_reportable_norm(self):
        worker = collector(
            residual_bound=10,
            candidate_bound=15,
            refinement_width=4,
            score_policy="powers",
        )
        stats = Counter()
        worker._sieve(-30, -14, stats)
        for offset, position in enumerate(range(-30, -14)):
            value = abs(worker.polynomial.value(position))
            remaining = value
            for prime in worker.factor_base.primes:
                while remaining and remaining % prime == 0:
                    remaining //= prime
            if remaining <= worker._candidate_bound:
                self.assertTrue(worker._piecewise_passes(offset))

    def test_interrupted_batch_preserves_first_uncommitted_position(self):
        worker = collector(smooth_batch_size=0)
        with patch(
            "v2.qs.cofactor_policy.batch_divisors",
            side_effect=BudgetExhaustedError,
        ):
            result = worker.collect(1, 4)
        self.assertEqual(result.next_position, 1)
        self.assertFalse(worker._atoms)
        spent = worker.budget.used
        worker.budget = allowance()
        resumed = worker.collect(1, 4)
        self.assertGreater(spent, 0)
        self.assertGreater(resumed.next_position, 1)

    def test_spill_limits_corruption_and_explicit_lifetime(self):
        with tempfile.TemporaryDirectory() as directory:
            owner = SpillStore(directory, 6, 65536, allowance())
            reference = owner.append(b"abc")
            owner.append(b"def")
            with self.assertRaises(MemoryError):
                owner.append(b"g")
            manifest = owner.manifest()
            owner.close()
            restored = SpillStore(directory, 6, 65536, allowance())
            restored.restore_manifest(manifest)
            self.assertEqual(restored.read(reference), b"abc")
            path = Path(directory) / reference[0]
            path.write_bytes(b"xbcdef")
            with self.assertRaises(ValueError):
                restored.read(reference)
            restored.cleanup()
            self.assertFalse(path.exists())

    def test_spilled_collector_restores_without_duplicate_disk_charge(self):
        with tempfile.TemporaryDirectory() as directory:
            worker = collector(spill_directory=directory, spill_bytes=2**20)
            for position in (1, -21):
                worker._admit(atom_at(worker, position), Counter())
            payload = json.loads(json.dumps(pack_store(worker)))
            used = worker._arena.spill.used
            restored = collector(spill_directory=directory, spill_bytes=2**20)
            restore_store(payload, restored, restored.budget)
            self.assertEqual(restored._arena.spill.used, used)
            self.assertEqual(pack_store(restored), payload)
            worker._arena.spill.close()
            restored._arena.spill.cleanup()

    def test_spill_reads_follow_the_resumed_allowance(self):
        with tempfile.TemporaryDirectory() as directory:
            worker = collector(spill_directory=directory, spill_bytes=2**20)
            source = atom_at(worker, 1)
            worker._admit(source, Counter())
            worker.budget = allowance(0)
            with self.assertRaises(BudgetExhaustedError):
                tuple(worker._atoms[source.relation_id].exponents)
            worker.budget = allowance()
            self.assertEqual(
                tuple(worker._atoms[source.relation_id].exponents),
                source.exponents,
            )
            worker._arena.spill.cleanup()

    def test_cold_count_is_checked_before_reading_an_extent(self):
        from v2.qs.cold_factors import ColdFactors

        with tempfile.TemporaryDirectory() as directory:
            owner = SpillStore(directory, 2**20, 65536, allowance())
            with self.assertRaises(ValueError):
                ColdFactors(
                    owner, ("c9a-cold-test.bin", 0, 10000, "0" * 64), (2,), 1
                )

    def test_component_contractions_reuse_untouched_originals(self):
        cache = IncrementalStructure(16 * 2**20)
        rows = (0b11, 0b11, 0b110000, 0b110000)
        for suffix in ((), (0,), (0b10010,), (0b10010,)):
            rows += suffix
            previous = dict(cache.cores)
            matrix = filter_matrix(
                rows,
                weight_two=True,
                budget=allowance(),
                _incremental=cache,
                memory_bytes=16 * 2**20,
            )
            self.assertTrue(matrix.stats["incremental_core"])
            dependencies = DependencySolver(matrix, budget=allowance()).run()
            self.assertEqual(span(dependencies), dense_kernel(rows))
            if suffix == (0,):
                for root, value in previous.items():
                    self.assertIs(cache.cores[root], value)

    def test_reporting_restore_refusal_does_not_publish_progress(self):
        worker = collector()
        payload = pack_store(worker)
        payload["reporting"]["positions"] = 17
        original = ReportingController.restore

        def refuse_after_validation(controller, state):
            original(controller, state)
            raise BudgetExhaustedError()

        with patch.object(
            ReportingController, "restore", refuse_after_validation
        ):
            with self.assertRaises(BudgetExhaustedError):
                restore_store(payload, worker, worker.budget)
        self.assertEqual(worker._reporting.positions, 0)

    def test_base_batch_recovers_all_valuations_independently(self):
        worker = collector()
        helper = BaseRecovery(worker.factor_base, 2**20, allowance())
        values = (2**20 * 15, 2**7 * 3, 1, 2**20 * 15)
        self.assertEqual(helper.residuals(values, allowance()), (15, 3, 1, 15))

    def test_background_hints_match_serial_exact_gcd(self):
        values = (3 * 5, 3 * 7, 5 * 11, 13 * 17)
        results = PendingHints(values, allowance()).result()
        for value, hint in zip(values, results):
            if hint is not None:
                self.assertTrue(1 < hint < value and value % hint == 0)
        self.assertTrue(any(results))

    def test_outcome_feedback_and_reporting_state_are_bounded(self):
        target = OutcomeTarget()
        target.extracted(100, [SimpleNamespace(divisor=None)], 0, None, 5, 200)
        self.assertEqual(target.streak, 1)
        self.assertTrue(100 < target.next_target <= 125)
        report = ReportingController(100, 400, 32, 4)
        report.observe(32, 32, 0)
        self.assertEqual(report.bound, 100)
        report.observe(32, 32, 0)
        self.assertEqual(report.bound, 200)
        clone = ReportingController(100, 400, 32, 4)
        clone.restore(report.state())
        self.assertEqual(clone.state(), report.state())

    def test_schema13_roundtrip_and_policy_binding(self):
        config = SIQSConfig(
            base_bound=100,
            half_width=256,
            max_half_width=256,
            memory_bytes=128 * 2**20,
            solve_policy="collection",
            collector=policy(),
        )
        job = SIQSJob(
            10007 * 10009, seed=73, config=config, budget=allowance()
        )
        job.run(max_blocks=1)
        checkpoint = json.loads(json.dumps(job.checkpoint()))
        self.assertEqual(checkpoint["version"], 13)
        restored = SIQSJob.from_checkpoint(
            checkpoint, config=config, budget=allowance()
        )
        self.assertEqual(
            restored.engine.next_position, job.engine.next_position
        )
        self.assertEqual(
            pack_store(restored.engine.collector),
            pack_store(job.engine.collector),
        )
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                checkpoint,
                config=replace(
                    config,
                    collector=replace(config.collector, lazy_cycles=False),
                ),
                budget=allowance(),
            )

    def test_no_production_default_is_replaced(self):
        self.assertNotIn("extension_policy", asdict(SIQSConfig().collector))


if __name__ == "__main__":
    unittest.main()
