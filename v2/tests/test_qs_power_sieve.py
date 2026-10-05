"""Exhaustive p-adic scoring and independent candidate-coverage checks."""

import contextlib
import hashlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from v2.budget import Budget, BudgetExhaustedError
from v2.prime_sieve import prime_sieve
from v2.qs import Polynomial, build_factor_base, verify_atomic
from v2.qs.extraction import DependencyExtractor, prepare_relations
from v2.qs.factor_base import FactorBase, FactorBaseEntry, modular_square_roots
from v2.qs.linear_algebra import DependencySolver, filter_matrix
from v2.qs.polynomial import PolynomialRoots
from v2.qs.power_sieve import prime_power_roots
from v2.qs.relations import AtomicRelation, verify_combined
from v2.qs.sieve_collector import SieveCollector, SieveConfig
from v2.tests.test_qs import reference_positions, unlimited_budget
from v2.tests.test_qs_sieve import signature


class PrimePowerSieveTests(unittest.TestCase):
    def test_lifted_roots_against_exhaustive_modular_enumeration(self):
        for n in (9, 25, 101, 10403):
            polynomial = Polynomial(n, 1, 1, 0)

            for prime in (2, 3, 5, 7, 11):
                values = tuple(
                    x for x in range(prime) if polynomial.value(x) % prime == 0
                )
                roots = PolynomialRoots(prime, values, False)
                marks = list(
                    prime_power_roots(
                        polynomial, roots, 2048, 1, unlimited_budget()
                    )
                )

                for modulus, lifted, weight in marks:
                    self.assertEqual(weight, 1)
                    expected = {
                        x
                        for x in range(modulus)
                        if polynomial.value(x) % modulus == 0
                    }

                    self.assertEqual(set(lifted), expected)

    def test_window_pruning_matches_direct_valuations(self):
        for n in (9, 25, 101, 10403, 3**12):
            polynomial = Polynomial(n, 1, 1, 0)

            for prime in (2, 3, 5, 7, 11):
                roots = PolynomialRoots(
                    prime,
                    tuple(
                        x
                        for x in range(prime)
                        if polynomial.value(x) % prime == 0
                    ),
                    False,
                )

                for lo, hi in ((-105, -92), (-7, 9), (47, 48), (513, 540)):
                    maximum = max(
                        abs(polynomial.value(x)) for x in range(lo, hi)
                    )
                    marks = list(
                        prime_power_roots(
                            polynomial,
                            roots,
                            maximum,
                            1,
                            unlimited_budget(),
                            lo=lo,
                            hi=hi,
                        )
                    )

                    for position in range(lo, hi):
                        value = abs(polynomial.value(position))
                        if not value:
                            continue
                        exponent = 0
                        while value % prime == 0:
                            exponent += 1
                            value //= prime
                        score = sum(
                            weight
                            for modulus, residues, weight in marks
                            if position % modulus in residues
                        )

                        self.assertGreaterEqual(score, exponent)
                        if not any(weight > 1 for _, _, weight in marks):
                            self.assertEqual(score, exponent)

        empty = PolynomialRoots(11, (3,), False)
        budget = unlimited_budget()

        self.assertEqual(
            list(
                prime_power_roots(
                    Polynomial(9, 1, 1, 0),
                    empty,
                    10**30,
                    1,
                    budget,
                    lo=0,
                    hi=1,
                )
            ),
            [],
        )
        self.assertLess(budget.used, 10)

    def test_singular_branch_cap_keeps_every_valuation_covered(self):
        polynomial = Polynomial(3**12, 1, 1, 0)
        roots = PolynomialRoots(3, (0,), False)
        maximum = 3**12
        marks = list(
            prime_power_roots(
                polynomial, roots, maximum, 2, unlimited_budget()
            )
        )

        self.assertTrue(any(weight > 2 for _, _, weight in marks))
        for position in range(-512, 513):
            value = abs(polynomial.value(position))
            if value == 0 or value > maximum:
                continue
            exponent = 0
            while value % 3 == 0:
                value //= 3
                exponent += 1
            score = sum(
                weight
                for modulus, residues, weight in marks
                if position % modulus in residues
            )

            self.assertGreaterEqual(score, 2 * exponent)

    def test_power_scoring_preserves_independent_signed_window_coverage(self):
        base = build_factor_base(10403, bound=100).factor_base

        for a, b in ((1, 102), (7, 1), (49, 8)):
            polynomial = Polynomial(10403, 1, a, b)
            expected = reference_positions(polynomial, base, -71, 80, 97)

            for backend in ("list", "array", "bytearray"):
                for division in ("full", "roots", "bucket", "resieve"):
                    for cutoff in (0, 5):
                        config = SieveConfig(
                            score_policy="powers",
                            score_backend=backend,
                            division=division,
                            block_width=17,
                            small_prime_cutoff=cutoff,
                            residual_bound=97,
                            max_atoms=4096,
                            max_partials=4096,
                            max_relations=4096,
                            memory_bytes=64 * 1024 * 1024,
                        )

                        result = SieveCollector(
                            polynomial,
                            base,
                            config=config,
                            budget=unlimited_budget(),
                        ).collect(-71, 80)

                        self.assertEqual(result.reason, "complete")
                        self.assertEqual(signature(result), expected)

    def test_sparse_matches_and_temporary_memory_remain_bounded(self):
        from dataclasses import replace

        small = build_factor_base(10403, bound=40).factor_base
        polynomial = Polynomial(10403, 1, 1, 102)

        initial = SieveCollector(
            polynomial,
            small,
            config=SieveConfig(residual_bound=97),
            budget=unlimited_budget(),
        ).collect(-128, 129)
        store = {a.relation_id: a for a in initial.atoms}
        pair = tuple(store[i] for i in initial.combined_relations[0].atom_ids)
        entries = list(small.entries)

        for prime in prime_sieve(100000):
            if prime < 50000:
                continue
            roots = modular_square_roots(
                10403, prime, budget=unlimited_budget()
            )
            if roots:
                entries.append(FactorBaseEntry(prime, roots))
            if len(entries) >= 1200:
                break

        base = FactorBase(10403, 1, 100000, tuple(entries))
        worker = SieveCollector(
            polynomial,
            base,
            config=SieveConfig(
                residual_bound=97, memory_bytes=32 * 1024 * 1024
            ),
            budget=unlimited_budget(),
        )
        stats = dict(duplicates=0, admitted_atoms=0, matches=0)

        for index in range(200):
            shift = 1000 * index
            translated = tuple(
                AtomicRelation(
                    Polynomial(10403, 1, 1, a.polynomial.b + shift),
                    a.position - shift,
                    a.sign,
                    a.exponents,
                    a.residual,
                )
                for a in pair
            )

            self.assertIsNone(worker._admit(translated[0], stats))
            if index == 199:
                prior = dict(worker._atoms)
                config = worker.config
                worker.config = replace(config, memory_bytes=worker._workspace)

                self.assertEqual(
                    worker._admit(translated[1], stats), "memory_limit"
                )
                self.assertEqual(worker._atoms, prior)
                worker.config = config

            self.assertIsNone(worker._admit(translated[1], stats))

        self.assertEqual(len(worker._combined), 200)
        self.assertLess(
            worker._workspace + worker._scratch_peak_bytes, 32 * 1024 * 1024
        )
        for relation in worker._combined:
            self.assertTrue(
                verify_combined(
                    relation,
                    base,
                    worker._atoms,
                    budget=unlimited_budget(),
                    memory_bytes=32 * 1024 * 1024,
                )
            )

    def test_larger_verified_store_reaches_extraction(self):
        base = build_factor_base(10403, bound=100).factor_base
        atoms = tuple(
            AtomicRelation(Polynomial(10403, 1, 1, 102 + i), -i, 1, ())
            for i in range(5000)
        )
        store = {atom.relation_id: atom for atom in atoms}
        prepared = prepare_relations(
            atoms,
            base,
            store,
            budget=unlimited_budget(),
            memory_bytes=128 * 1024 * 1024,
        )
        matrix = filter_matrix(
            prepared.rows,
            budget=unlimited_budget(),
            memory_bytes=128 * 1024 * 1024,
        )

        dependencies = DependencySolver(
            matrix, budget=unlimited_budget()
        ).run()

        result = DependencyExtractor(
            prepared, dependencies, budget=unlimited_budget()
        ).run()

        self.assertIn(result, (101, 103))
        with self.assertRaises(TypeError):
            base._columns[2] = 99

    def test_large_single_prime_residual_keeps_exact_identity(self):
        base = build_factor_base(4001, bound=100).factor_base
        polynomial = Polynomial(4001, 1, 1, 1002)
        config = SieveConfig(
            score_policy="powers",
            residual_bound=10**10,
            memory_bytes=64 * 1024 * 1024,
        )

        result = SieveCollector(
            polynomial, base, config=config, budget=unlimited_budget()
        ).collect(0, 1)

        self.assertEqual(len(result.atoms), 1)
        atom = result.atoms[0]

        self.assertEqual(atom.residual, 1000003)
        self.assertTrue(all(atom.residual % d for d in range(2, 1001)))
        self.assertTrue(verify_atomic(atom, base, residual_bound=10**10))

    def test_historic_collector_control_keeps_certified_coverage(self):
        from v2.benchmarks.qs_snapshot import load_qs_arm

        arm = load_qs_arm("_m31_power_loader", ("qs.sieve_collector",))
        base = arm.qs.factor_base.build_factor_base(
            10403, bound=100
        ).factor_base
        polynomial = arm.qs.polynomial.qs_polynomial(base)
        config = arm.qs.sieve_collector.SieveConfig(
            score_policy="powers",
            residual_bound=97,
            memory_bytes=64 * 1024 * 1024,
            max_atoms=4096,
            max_partials=4096,
            max_relations=4096,
        )

        result = arm.qs.sieve_collector.SieveCollector(
            polynomial, base, config=config, budget=unlimited_budget()
        ).collect(-71, 80)

        self.assertEqual(
            signature(result),
            reference_positions(polynomial, base, -71, 80, 97),
        )

    def test_benchmark_checkpoint_survives_disk_and_cumulative_resume(self):
        from v2.benchmarks.phase_three_large import configuration, run_one

        fixture = dict(
            id="checkpoint_control",
            kind="balanced",
            digits=5,
            n=10403,
            factors=[101, 103],
        )
        config = configuration(10403, 100, 32, 1000)

        self.assertEqual(config.checkpoint_bytes, 16 * 1024 * 1024)
        from dataclasses import replace

        with self.assertRaises(ValueError):
            replace(config, checkpoint_bytes=16 * 1024 * 1024 + 1)
        with tempfile.TemporaryDirectory() as directory:
            checkpoints = Path(directory) / "checkpoints"
            with contextlib.redirect_stdout(io.StringIO()):
                paused = run_one(
                    fixture,
                    7,
                    config,
                    "siqs",
                    seconds=0,
                    checkpoint_dir=checkpoints,
                )
                metadata = paused["checkpoint"]
                encoded = (Path(directory) / metadata["path"]).read_bytes()

                self.assertEqual(
                    hashlib.sha256(encoded).hexdigest(), metadata["sha256"]
                )
                resumed = run_one(
                    fixture,
                    7,
                    config,
                    "siqs",
                    seconds=3,
                    checkpoint=json.loads(encoded),
                    checkpoint_dir=checkpoints,
                )

            self.assertFalse(paused["complete"])
            self.assertTrue(resumed["complete"])
            self.assertTrue(resumed["resumed"])
            self.assertGreater(
                resumed["total_wall_used"], paused["total_wall_used"]
            )
            self.assertEqual(resumed["factors"], [101, 103])

    def test_lifting_refusal_is_charged_and_large_residuals_are_exact(self):
        polynomial = Polynomial(10403, 1, 1, 0)
        roots = PolynomialRoots(2, (1,), False)
        marks = prime_power_roots(
            polynomial,
            roots,
            1000,
            1,
            Budget(work_limit=0, seconds=None, cpu_seconds=None),
        )
        next(marks)
        with self.assertRaises(BudgetExhaustedError):
            next(marks)
        SieveConfig(residual_bound=10**12, max_atoms=65536)
        with self.assertRaises(ValueError):
            SieveConfig(residual_bound=10**12 + 1)
        with self.assertRaises(ValueError):
            SieveConfig(max_atoms=65537)


if __name__ == "__main__":
    unittest.main()
