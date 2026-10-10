"""Independent coverage, exponent recovery, and bounded store failure tests."""

import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.qs import Polynomial, build_factor_base, qs_polynomial, verify_combined
from v2.qs.sieve_collector import SieveCollector, SieveConfig
from v2.tests.test_qs import reference_positions, unlimited_budget


def signature(result):
    """Consume retained atom payloads without depending on score internals."""
    return {
        atom.position: (atom.sign, atom.exponents, atom.residual)
        for atom in result.atoms
    }


def collector(polynomial, base, **options):
    """Use no eviction for coverage checks and no timer noise."""
    options.setdefault("residual_bound", 97)
    config = SieveConfig(
        max_partials=4096, max_relations=4096, max_atoms=4096, **options
    )
    return SieveCollector(
        polynomial, base, config=config, budget=unlimited_budget()
    )


class SieveCoverageTests(unittest.TestCase):
    """Independent generic factoring and modular enumeration references."""

    def test_all_modes_signed_windows_and_tails(self):
        """All buffer/mark/division modes retain identical verified atoms."""
        base = build_factor_base(10403, bound=100).factor_base

        for a, b in ((1, 102), (7, 1), (49, 8)):
            polynomial = Polynomial(10403, 1, a, b)
            expected = reference_positions(polynomial, base, -71, 80, 97)

            for backend in ("list", "bytearray", "array"):
                for marking in ("dense", "sparse", "bucket"):
                    for division in ("full", "roots", "bucket", "resieve"):
                        with self.subTest(
                            a=a,
                            backend=backend,
                            marking=marking,
                            division=division,
                        ):
                            run = collector(
                                polynomial,
                                base,
                                score_backend=backend,
                                marking=marking,
                                division=division,
                                block_width=17,
                                metadata_chunk=3,
                            ).collect(-71, 80)

                            self.assertEqual(run.reason, "complete")
                            self.assertEqual(signature(run), expected)
                            store = {a.relation_id: a for a in run.atoms}

                            for relation in run.combined_relations:
                                self.assertTrue(
                                    verify_combined(
                                        relation,
                                        base,
                                        store,
                                        budget=unlimited_budget(),
                                    )
                                )

    def test_cutoffs_blocks_and_high_prime_powers(self):
        """Omitted primes and high valuations cannot defeat safe scores."""
        for n, h in ((10403, 1), (10403, 9), (1022117, 3)):
            base = build_factor_base(n, multiplier=h, bound=100).factor_base
            polynomial = qs_polynomial(base)
            expected = reference_positions(polynomial, base, -111, 146, 97)

            for width in (1, 8, 64, 257):
                for cutoff in (0, 3, 8, 100):
                    run = collector(
                        polynomial,
                        base,
                        block_width=width,
                        small_prime_cutoff=cutoff,
                    ).collect(-111, 146)

                    self.assertEqual(signature(run), expected)

        # N is a square only for this arithmetic diagnostic. F(2)=2**40
        # when B=2**38-1; byte scores must saturate conservatively.
        b = 2**38 - 1
        polynomial = Polynomial(b * b, 1, 1, b)
        base = build_factor_base(b * b, bound=3).factor_base
        for backend in ("list", "bytearray", "array"):
            run = collector(polynomial, base, score_backend=backend).collect(
                2, 3
            )

            self.assertEqual(run.atoms[0].exponents, ((2, 40),))

    def test_saturation_large_threshold(self):
        """Clip both sides; a 300-bit exact power still becomes a candidate."""
        b = 2**298 - 1
        polynomial = Polynomial(b * b, 1, 1, b)
        base = build_factor_base(b * b, bound=3).factor_base

        run = collector(polynomial, base, score_backend="bytearray").collect(
            2, 3
        )

        self.assertEqual(run.stats["threshold_max"], 255)
        self.assertEqual(run.atoms[0].exponents, ((2, 300),))

    def test_adaptive_scores_and_candidate_resieving(self):
        """Bypassed zero thresholds and resieving retain all exact payloads."""
        base = build_factor_base(1022117, bound=100).factor_base
        polynomial = qs_polynomial(base)
        expected = reference_positions(polynomial, base, -128, 129, 500)

        for policy in ("adaptive", "conservative", "candidate"):
            for cutoff in (0, 8, 100):
                worker = collector(
                    polynomial,
                    base,
                    residual_bound=500,
                    division="resieve",
                    score_policy=policy,
                    small_prime_cutoff=cutoff,
                    block_width=17,
                )

                run = worker.collect(-128, 129)

                self.assertEqual(signature(run), expected)
                self.assertFalse(worker._resieved)
                if policy == "adaptive":
                    self.assertGreater(run.stats["skipped_score_blocks"], 0)

        worker = collector(polynomial, base, division="resieve")
        worker.budget = Budget(work_limit=0)

        stopped = worker.collect(-128, 129)

        self.assertEqual(stopped.next_position, -128)
        self.assertFalse(stopped.atoms)
        worker.budget = unlimited_budget()

        resumed = worker.collect(stopped.next_position, 129)
        expected_small = reference_positions(polynomial, base, -128, 129, 97)

        self.assertEqual(signature(resumed), expected_small)

    def test_refined_thresholds_against_independent_oracle(self):
        """Per-norm thresholds preserve payloads across bounded backends."""
        base = build_factor_base(10403, bound=40).factor_base
        polynomial = qs_polynomial(base)
        expected = reference_positions(polynomial, base, -128, 129, 97)

        for backend in ("list", "bytearray", "array"):
            for division in ("full", "roots", "bucket", "resieve"):
                for cutoff in (0, 8, 100):
                    result = collector(
                        polynomial,
                        base,
                        score_policy="candidate",
                        score_backend=backend,
                        division=division,
                        small_prime_cutoff=cutoff,
                    ).collect(-128, 129)

                    self.assertEqual(signature(result), expected)
                    if cutoff == 0:
                        self.assertGreater(
                            result.stats["refined_rejections"], 0
                        )

    def test_roots_against_independent_enumeration(self):
        """Hit exclusion includes singular A, 2 and degenerate roots."""
        base = build_factor_base(10403, multiplier=9, bound=100).factor_base

        for polynomial in (qs_polynomial(base), Polynomial(10403, 9, 9, 0)):
            worker = collector(polynomial, base)

            for roots in worker._roots:
                expected = {
                    x
                    for x in range(roots.prime)
                    if (
                        polynomial.a * x * x
                        + 2 * polynomial.b * x
                        + polynomial.c
                    )
                    % roots.prime
                    == 0
                }
                actual = (
                    set(range(roots.prime))
                    if roots.all_positions
                    else (set(roots.roots))
                )

                self.assertEqual(actual, expected)

    def test_resieving_powers_singular_roots_and_scratch_refusal(self):
        """Recovery handles high valuations and atomic scratch caps."""
        b = 2**298 - 1
        polynomial = Polynomial(b * b, 1, 1, b)
        base = build_factor_base(b * b, bound=3).factor_base

        for policy in ("adaptive", "candidate"):
            run = collector(
                polynomial,
                base,
                division="resieve",
                score_backend="bytearray",
                score_policy=policy,
            ).collect(2, 3)

            self.assertEqual(run.atoms[0].exponents, ((2, 300),))

        base = build_factor_base(10403, multiplier=9, bound=100).factor_base
        polynomial = Polynomial(10403, 9, 9, 0)
        expected = reference_positions(polynomial, base, -71, 80, 97)
        worker = collector(polynomial, base, division="resieve")

        result = worker.collect(-71, 80)

        self.assertEqual(signature(result), expected)
        worker = collector(polynomial, base, division="resieve")
        # Fault-inject a nearly occupied reservation before candidate scratch.
        worker._workspace = worker.config.memory_bytes - 1

        result = worker.collect(-71, 80)

        self.assertEqual(result.reason, "memory_limit")
        self.assertEqual(result.next_position, -71)
        self.assertFalse(result.atoms)

    def test_window_larger_than_working_block_limit(self):
        """A full window spans reusable blocks and retains its exact tail."""
        base = build_factor_base(1022117, bound=100).factor_base
        polynomial = qs_polynomial(base)
        expected = reference_positions(polynomial, base, 0, 4097, 500)

        run = collector(
            polynomial,
            base,
            block_width=256,
            residual_bound=500,
            memory_bytes=32 * 1024 * 1024,
        ).collect(0, 4097)

        self.assertEqual(run.reason, "complete")
        self.assertEqual(run.stats["blocks"], 17)
        self.assertEqual(run.next_position, 4097)
        self.assertEqual(signature(run), expected)

    def test_zero_empty_and_factor(self):
        """Zeros are checked before division; proper GCD stops explicitly."""
        base = build_factor_base(101**2, bound=40).factor_base
        polynomial = qs_polynomial(base)
        worker = collector(polynomial, base)

        empty = worker.collect(-3, -3)

        self.assertEqual(empty.stats["scanned"], 0)

        zero = worker.collect(0, 1)

        self.assertEqual(zero.divisor, 101)
        self.assertEqual(zero.next_position, 1)
        self.assertEqual(zero.reason, "factor_found")
        base = build_factor_base(101 * 103, bound=40).factor_base
        worker = collector(qs_polynomial(base), base, residual_bound=500)

        run = worker.collect(-1, 0)

        self.assertEqual(run.divisor, 101)

    def test_composite_residual_and_lossy_threshold(self):
        """Extra thresholds lose yield only in explicitly lossy mode."""
        base = build_factor_base(10403, bound=100).factor_base
        polynomial = qs_polynomial(base)

        safe = collector(polynomial, base).collect(-128, 129)

        lossy = collector(polynomial, base, threshold_extra=1000).collect(
            -128, 129
        )

        self.assertGreater(len(safe.atoms), len(lossy.atoms))
        small_base = build_factor_base(9471, bound=3).factor_base

        composite = collector(
            Polynomial(9471, 1, 1, 100), small_base, residual_bound=1000
        ).collect(0, 1)

        self.assertGreater(composite.stats["composite_residuals"], 0)


class SieveStoreTests(unittest.TestCase):
    """FIFO eviction, provenance retention, cancellation and atomic refusal."""

    def setUp(self):
        """Use enough partial matches to force eviction and combinations."""
        self.base = build_factor_base(10403, bound=40).factor_base
        self.polynomial = qs_polynomial(self.base)

    def test_eviction_never_deletes_combined_provenance(self):
        """Tiny FIFO caps evict unmatched IDs and preserve every pinned ID."""
        worker = SieveCollector(
            self.polynomial,
            self.base,
            config=SieveConfig(max_partials=2, residual_bound=97),
            budget=unlimited_budget(),
        )

        run = worker.collect(-128, 129)

        self.assertEqual(run.reason, "complete")
        self.assertGreater(run.stats["evictions"], 0)
        self.assertGreater(len(run.combined_relations), 0)
        self.assertLessEqual(len(run.partial_ids), 2)
        store = {atom.relation_id: atom for atom in run.atoms}
        for relation in run.combined_relations:
            verify_combined(
                relation, self.base, store, budget=unlimited_budget()
            )
        other = SieveCollector(
            self.polynomial,
            self.base,
            config=SieveConfig(max_partials=2, residual_bound=97),
            budget=unlimited_budget(),
        ).collect(-128, 129)

        self.assertEqual(run.atoms, other.atoms)
        self.assertEqual(run.stats, other.stats)

    def test_caps_refuse_before_publication(self):
        """Atom/relation caps leave the first refused position resumable."""
        for config, reason in (
            (SieveConfig(max_atoms=0), "atom_limit"),
            (SieveConfig(max_relations=0, residual_bound=1), "relation_limit"),
        ):
            worker = SieveCollector(
                self.polynomial,
                self.base,
                config=config,
                budget=unlimited_budget(),
            )

            run = worker.collect(-128, 129)

            self.assertEqual(run.reason, reason)
            self.assertFalse(run.atoms)
            worker.config = replace(config, max_atoms=2048, max_relations=1024)

            resumed = worker.collect(run.next_position, 129)

            expected = collector(
                self.polynomial,
                self.base,
                residual_bound=config.residual_bound,
            ).collect(run.next_position, 129)

            self.assertEqual(resumed.atoms, expected.atoms)

    def test_memory_refusal_and_zero_partial_cap(self):
        """Check storage before mutation and count dropped partials."""
        worker = collector(self.polynomial, self.base)
        worker.config = replace(worker.config, memory_bytes=worker._workspace)

        run = worker.collect(-128, 129)

        self.assertEqual(run.reason, "memory_limit")
        self.assertFalse(run.atoms)

        run = SieveCollector(
            self.polynomial,
            self.base,
            config=SieveConfig(max_partials=0, residual_bound=97),
            budget=unlimited_budget(),
        ).collect(-128, 129)

        self.assertGreater(run.stats["dropped_partials"], 0)
        self.assertTrue(all(atom.residual == 1 for atom in run.atoms))

    def test_budget_refusal_combination_and_resume(self):
        """Fault injection at combination leaves partner intact on retry."""
        worker = collector(self.polynomial, self.base)

        expected = collector(self.polynomial, self.base).collect(-128, 129)
        worker.budget.reason = "work_limit"
        with patch(
            "v2.qs.sieve_collector.combine_relations",
            side_effect=BudgetExhaustedError,
        ):
            run = worker.collect(-128, 129)

        self.assertEqual(run.reason, "work_limit")
        self.assertTrue(run.partial_ids)
        worker.budget = unlimited_budget()

        resumed = worker.collect(run.next_position, 129)

        self.assertEqual(resumed.atoms, expected.atoms)
        self.assertEqual(
            resumed.combined_relations, expected.combined_relations
        )

    def test_budget_stops_at_first_uncommitted_position(self):
        """Zero work/cancellation/deadlines do not publish a partial block."""
        for budget, reason in (
            (Budget(work_limit=0), "work_limit"),
            (Budget(cancelled=lambda: True), "cancelled"),
            (Budget(seconds=0), "wall_limit"),
            (Budget(cpu_seconds=0), "cpu_limit"),
        ):
            worker = collector(self.polynomial, self.base)
            worker.budget = budget

            run = worker.collect(-128, 129)

            self.assertEqual(run.reason, reason)
            self.assertEqual(run.next_position, -128)
            self.assertFalse(run.atoms)

    def test_verifier_failure_does_not_publish_atom(self):
        """Rejected arithmetic leaves the accepted store untouched."""
        worker = collector(self.polynomial, self.base)
        with patch(
            "v2.qs.sieve_collector.verify_atomic",
            side_effect=ValueError("corrupt exponent"),
        ):
            with self.assertRaises(ValueError):
                worker.collect(-128, 129)

        self.assertFalse(worker._atoms)
        self.assertFalse(worker._pending)
        self.assertFalse(worker._combined)

    def test_duplicate_and_invalid_inputs(self):
        """Rescans cannot admit or combine a pinned atom twice."""
        worker = collector(self.polynomial, self.base)

        first = worker.collect(-128, 129)

        second = worker.collect(-128, 129)

        self.assertEqual(first.atoms, second.atoms)
        self.assertEqual(first.combined_relations, second.combined_relations)
        self.assertGreater(second.stats["duplicates"], 0)
        for options in (
            {"block_width": 0},
            {"score_backend": "float"},
            {"residual_bound": 2**64},
            {"metadata_chunk": 1025},
        ):
            with self.assertRaises(ValueError):
                SieveConfig(**options)

        with self.assertRaises(MemoryError):
            SieveCollector(
                self.polynomial, self.base, config=SieveConfig(memory_bytes=0)
            )
        with self.assertRaises(ValueError):
            worker.collect(0, 1_000_001)
        with self.assertRaises(ValueError):
            SieveCollector(Polynomial(10403, 1, 101, 0), self.base)
        with patch(
            "v2.qs.sieve_collector.polynomial_roots",
            side_effect=ValueError("inverse failure"),
        ):
            with self.assertRaises(ValueError):
                collector(self.polynomial, self.base)
