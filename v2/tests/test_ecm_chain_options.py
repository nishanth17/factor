"""Independent records and finite recovery for optional native chains."""

import unittest
from unittest.mock import patch

from v2 import arithmetic, ecm, ecm_chains, prac, utils
from v2.benchmarks import c6_cf, c6_fast
from v2.benchmarks.prac_oracle import (
    affine_multiply,
    historical_points,
    matches,
)
from v2.budget import Budget, BudgetExhaustedError
from v2.ecm_chain_options import (
    ChainPlan,
    ChainPlans,
    identity,
    verify_progress,
)
from v2.schedules import SieveContext


def allowance(work=10**8):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


class OptionalChainTests(unittest.TestCase):
    def test_records_and_independent_affine_oracle(self):
        catalogs = {
            "lucas": c6_fast.load_catalog()["lucas"],
            "cf": {
                scalar: row.certified
                for scalar, row in c6_cf.load_catalog().items()
            },
        }
        for family, catalog in catalogs.items():
            plan = ChainPlan(2000, "python-int", allowance(), family)
            self.assertLessEqual(plan.owned_bytes, ecm_chains.PLAN_BYTES)
            for _, (power, action, unit, _) in plan.entries.items():
                self.assertEqual(
                    action.record.code, catalog[power].record.code
                )
                self.assertEqual(action.masks, catalog[power].masks)
                self.assertEqual(unit.record.scalar, unit.strict.record.scalar)
                for modulus, curve_a, point in historical_points():
                    a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
                    expected = affine_multiply(power, point, modulus, curve_a)
                    value, divisor, _ = plan.execute(
                        [(unit.record.scalar, power)],
                        (point[0], 1),
                        modulus,
                        a24,
                        allowance(),
                    )
                    if value is not None:
                        self.assertTrue(matches(value, expected, modulus))
                    self.assertTrue(
                        divisor is None
                        or utils.valid_divisor(divisor, modulus)
                    )

    def test_refusal_precedes_arithmetic_and_saturation_keeps_divisor(self):
        for family in ("lucas", "cf"):
            plan = ChainPlan(2000, "python-int", allowance(), family)
            power, action, _, _ = plan.entries[2]
            with patch.object(action, "run", return_value=((1, 1), 35)) as run:
                with self.assertRaises(BudgetExhaustedError):
                    plan.execute([[2, power]], (2, 1), 35, 2, allowance(1))
                run.assert_not_called()
                with patch.object(
                    action, "strict", side_effect=prac.NonunitPointError(5)
                ):
                    value, divisor, replayed = plan.execute(
                        [[2, power]], (2, 1), 35, 2, allowance()
                    )
                self.assertEqual((value, divisor, replayed), (None, 5, True))

    def test_finite_cache_failure_cancellation_and_eviction(self):
        for family in ("lucas", "cf"):
            store = ChainPlans(
                ecm_chains.MIN_MEMORY_BYTES, "python-int", (1999, 2000), family
            )
            budget = allowance()
            first = store.get(2000, budget)
            before = budget.used
            self.assertIs(first, store.get(2000, budget))
            self.assertEqual(budget.used, before + 1)
            store.get(1999, budget)
            self.assertEqual(store.evictions, 1)
            self.assertLessEqual(store.used_bytes, store.memory_bytes)
            with self.assertRaises(BudgetExhaustedError):
                store.get(2000, Budget(cancelled=lambda: True))
            failed = ChainPlans(
                ecm_chains.MIN_MEMORY_BYTES, "python-int", (2000,), family
            )
            with self.assertRaises(BudgetExhaustedError):
                failed.get(2000, allowance(1))
            self.assertFalse(failed.plans)
            self.assertEqual(failed.used_bytes, ecm_chains.SCRATCH_BYTES)

    def test_optional_identity_and_prefix_do_not_change_legacy(self):
        for family in ("lucas", "cf"):
            n = 1000003
            setup = ecm.setup_curve(n, 17)
            job = dict(
                phase="stage_two_setup",
                b1=2000,
                n=n,
                value=list(setup.point),
                chain_chunks=19,
                chain_identity=identity(2000, "python-int", family),
            )
            verifier = SieveContext(2501)
            verify_progress(job, "python-int", verifier, family)
            other = "lucas" if family == "cf" else "cf"
            with self.assertRaises(ValueError):
                verify_progress(job, "python-int", verifier, other)
            self.assertIn("/batch16/optional-", job["chain_identity"])
        self.assertEqual(
            ecm_chains.identity(2000, "python-int"),
            ecm_chains.CHAIN_VERSION
            + "/2000/prac-reduced/"
            + "python-int/batch16",
        )

    def test_optional_gmp_records_keep_exact_coordinates(self):
        try:
            backend = arithmetic.get_backend("gmpy2-mpz")
        except arithmetic.BackendUnavailableError:
            self.skipTest("optional GMP backend unavailable")
        n = backend.integer(1000003)
        setup = ecm.setup_curve(n, 17)
        for family in ("lucas", "cf"):
            plan = ChainPlan(2000, backend.name, allowance(), family)
            for prime, (power, _, _, _) in list(plan.entries.items())[:16]:
                value, divisor, _ = plan.execute(
                    [[prime, power]], setup.point, n, setup.a24, allowance()
                )
                if value is not None:
                    expected = ecm.scalar_multiply(
                        power, *setup.point, n, setup.a24
                    )
                    self.assertEqual(
                        (value[0] * expected[1] - value[1] * expected[0]) % n,
                        0,
                    )
                self.assertTrue(
                    divisor is None or utils.valid_divisor(divisor, n)
                )


if __name__ == "__main__":
    unittest.main()
