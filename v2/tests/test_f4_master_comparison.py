"""Failure accounting independent of ECM arithmetic or factor discovery."""

import unittest
from types import SimpleNamespace

from v2.benchmarks.ecm.b12_f4.selection_child import (
    run_accounted,
    validate_reconstruction,
)


class FailureAccountingTests(unittest.TestCase):
    def test_result_serialization_refusal_keeps_prior_and_current_charges(
        self,
    ):
        budget = SimpleNamespace(used=0, prior_cpu=0.0)
        totals = dict(work=100, children=7.0)

        def refusing(*, budget):
            budget.used = 23
            budget.prior_cpu = 11.5
            raise MemoryError("checkpoint refused after charged work")

        with self.assertRaises(MemoryError):
            run_accounted(refusing, budget, totals)
        self.assertEqual(totals, dict(work=123, children=18.5))

    def test_success_and_other_failure_each_charge_once(self):
        budget = SimpleNamespace(used=3, prior_cpu=2.0)
        totals = dict(work=0, children=0.0)
        self.assertEqual(
            run_accounted(lambda **kwargs: 17, budget, totals), 17
        )

        def failing(**kwargs):
            raise ValueError("malformed mathematical result")

        with self.assertRaises(ValueError):
            run_accounted(failing, budget, totals)
        self.assertEqual(totals, dict(work=6, children=4.0))

    def test_refusal_preserves_input_without_claiming_completion(self):
        row = dict(
            n=35,
            complete=False,
            factorization_complete=False,
            factors=[],
            remaining=[35],
            reconstruct=35,
            endpoint="complete",
        )
        validate_reconstruction(row)
        for field, invalid in (("remaining", []), ("complete", "false")):
            damaged = dict(row, **{field: invalid})
            with self.assertRaises(ValueError):
                validate_reconstruction(damaged)


if __name__ == "__main__":
    unittest.main()
