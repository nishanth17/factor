"""Independent n-1 proof checks and generation failure boundaries."""

import copy
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from v2.benchmarks.ecm.c3 import round2_prime_proofs as proofs
from v2.execution.budget import Budget, BudgetExhaustedError


class PrimeProofTests(unittest.TestCase):
    def test_integer_checker_accepts_full_and_partial_proofs(self):
        certificates = {
            "2": dict(kind="trial"),
            "17": dict(kind="trial"),
            "65537": dict(kind="lucas", factors=[[2, 16]], witness=3),
            "103": dict(kind="pocklington", q=17, witness=2),
        }
        with patch.object(
            proofs.utils, "is_prime", side_effect=AssertionError
        ):
            proofs.verify_proofs(certificates)

    def test_checker_rejects_changed_factors_witnesses_and_kinds(self):
        original = {
            "2": dict(kind="trial"),
            "65537": dict(kind="lucas", factors=[[2, 16]], witness=3),
        }
        for field, value in (
            ("factors", [[2, 15]]),
            ("factors", [[2, 8], [2, 8]]),
            ("witness", 65536),
            ("kind", "probable"),
        ):
            changed = copy.deepcopy(original)
            changed["65537"][field] = value
            with self.assertRaises(ValueError):
                proofs.verify_proofs(changed)

    def test_generation_proves_instead_of_trusting_terminal_labels(self):
        source = proofs.UniformPrimeSource(17)
        source.prove(65537)
        self.assertEqual(source.certificates["65537"]["kind"], "lucas")
        proofs.verify_proofs(source.certificates)
        self.assertGreater(source.budget.used, 0)

    def test_incomplete_proof_aborts_instead_of_resampling(self):
        source = proofs.UniformPrimeSource(17)
        result = SimpleNamespace(complete=False)
        with patch.object(
            proofs,
            "factorize_bounded",
            return_value=SimpleNamespace(result=result),
        ):
            with self.assertRaisesRegex(RuntimeError, "did not complete"):
                source.prove(65537)
        self.assertNotIn("65537", source.certificates)

    def test_sampling_is_deterministic_and_respects_shared_allowance(self):
        left = proofs.UniformPrimeSource(29)
        right = proofs.UniformPrimeSource(29)
        self.assertEqual(left.prime(4), right.prime(4))
        proofs.verify_proofs(left.certificates)
        source = proofs.UniformPrimeSource(17, budget=Budget(work_limit=0))
        with self.assertRaises(BudgetExhaustedError):
            source.prove(17)

    def test_recursive_proof_rechecks_storage_before_parent(self):
        source = proofs.UniformPrimeSource(17)
        with patch.object(proofs, "MAX_NODES", 1):
            with self.assertRaisesRegex(RuntimeError, "node allowance"):
                source.prove(65537)
        self.assertEqual(set(source.certificates), {"2"})

    def test_fermat_pseudoprime_is_not_accepted_as_a_proven_prime(self):
        certificates = {
            "2": dict(kind="trial"),
            "5": dict(kind="trial"),
            "17": dict(kind="trial"),
            "341": dict(
                kind="lucas", factors=[[2, 2], [5, 1], [17, 1]], witness=2
            ),
        }
        self.assertEqual(pow(2, 340, 341), 1)
        with self.assertRaisesRegex(ValueError, "invalid n-1 prime proof"):
            proofs.verify_proofs(certificates)

    def test_multiple_proofs_keep_one_cumulative_parent_allowance(self):
        source = proofs.UniformPrimeSource(17)
        source.prove(65537)
        previous = source.budget.used
        source.prove(131071)
        self.assertGreater(source.budget.used, previous)
        proofs.verify_proofs(source.certificates)

    def test_child_work_is_charged_before_parent_deadline_check(self):
        # Explicit starts avoid default factories bound to original clocks.
        parent = Budget(
            work_limit=100, seconds=1, cpu_seconds=None, _wall_start=100
        )
        source = proofs.UniformPrimeSource(17, budget=parent)
        clock = [100]

        def finish(*args, budget, **kwargs):
            budget.consume(7)
            clock[0] = 102
            return SimpleNamespace()

        with (
            patch("v2.execution.budget.time.monotonic", lambda: clock[0]),
            patch.object(proofs, "factorize_bounded", side_effect=finish),
        ):
            with self.assertRaises(BudgetExhaustedError):
                source._factor_predecessor(65537)
        self.assertEqual(parent.used, 7)
