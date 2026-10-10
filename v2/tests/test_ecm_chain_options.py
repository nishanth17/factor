"""Independent records and finite recovery for optional native chains."""

import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.benchmarks.ecm.c6 import c6_cf, c6_fast
from v2.benchmarks.support.prac_oracle import (
    affine_multiply,
    historical_points,
    matches,
)
from v2.common import arithmetic, utils
from v2.ecm import chains as ecm_chains
from v2.ecm import core as ecm
from v2.ecm import prac
from v2.ecm.chain_options import (
    ChainPlan,
    ChainPlans,
    default_options,
    identity,
    verify_progress,
)
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.schedules import SieveContext
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.tests.test_ecm_chains import configuration


def allowance(work=10**8):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


class OptionalChainTests(unittest.TestCase):
    def test_finite_default_decision_and_fallbacks(self):
        base = PortfolioConfig(memory_bytes=16 * 1024**2, ecm_chain_mode="off")
        options = default_options(base)
        candidate = replace(base, **options)
        automatic = replace(base, ecm_chain_mode="auto")
        self.assertEqual(automatic, candidate)
        self.assertEqual(candidate.ecm_chain_mode, "reuse")
        self.assertEqual(
            candidate.ecm_chain_bytes, ecm_chains.MIN_MEMORY_BYTES
        )
        self.assertEqual(candidate.ecm_program_bytes, 512 * 1024)
        self.assertGreaterEqual(
            candidate.memory_bytes - candidate.workspace_reserve, 8192
        )
        fresh = PortfolioConfig()
        self.assertEqual(fresh.ecm_chain_mode, "reuse")
        self.assertEqual(fresh.memory_bytes, 16 * 1024**2)
        self.assertEqual(
            PortfolioConfig(memory_bytes=8 * 1024**2).ecm_chain_mode, "off"
        )
        self.assertEqual(
            PortfolioConfig(ecm_tiers=()).memory_bytes, 8 * 1024**2
        )
        for changes in (
            dict(memory_bytes=8 * 1024**2),
            dict(chunk_size=32),
            dict(ecm_tiers=((2001, 2500, 32),)),
            dict(ecm_tiers=((2000, 2500, 7),)),
            dict(ecm_tiers=()),
        ):
            unsupported = replace(base, **changes)
            self.assertEqual(
                default_options(unsupported)["ecm_chain_mode"], "off"
            )
        try:
            gmp = replace(base, backend="gmpy2-mpz")
        except arithmetic.BackendUnavailableError:
            return
        self.assertEqual(default_options(gmp)["ecm_chain_mode"], "off")
        selected = replace(gmp, ecm_chain_family="lucas")
        self.assertEqual(default_options(selected)["ecm_chain_mode"], "reuse")

    def test_optional_portfolio_resume_pins_family_and_reconstructs(self):
        n = 6120168563605791616423380424731852610871
        for family in ("lucas", "cf"):
            config = configuration(ecm_chain_family=family)
            full = factorize_bounded(
                n, seed=19, config=config, budget=allowance()
            )
            partial = factorize_bounded(
                n, seed=19, config=config, budget=allowance(250000)
            )
            self.assertEqual(partial.checkpoint["payload"]["version"], 11)
            restored = factorize_bounded(
                n,
                config=config,
                checkpoint=partial.checkpoint,
                budget=allowance(),
            )
            self.assertEqual(restored.result, full.result)
            self.assertEqual(restored.result.reconstruct(), n)
            other = replace(
                config, ecm_chain_family="cf" if family == "lucas" else "lucas"
            )
            with self.assertRaises(ValueError):
                factorize_bounded(
                    n,
                    config=other,
                    checkpoint=partial.checkpoint,
                    budget=allowance(),
                )

    def test_implicit_resume_keeps_saved_caps_and_legacy_executor(self):
        # Exact terminal input keeps this test about serialized defaults;
        # arithmetic-bearing legacy/resume coverage is in test_ecm_chains.
        for config in (
            PortfolioConfig(memory_bytes=8 * 1024**2, ecm_chain_mode="off"),
            replace(
                PortfolioConfig(memory_bytes=16 * 1024**2),
                ecm_chain_mode="reuse",
                ecm_program_bytes=512 * 1024,
                ecm_chain_bytes=ecm_chains.MIN_MEMORY_BYTES,
            ),
        ):
            first = factorize_bounded(101, config=config, budget=allowance(0))
            restored = factorize_bounded(
                101, checkpoint=first.checkpoint, budget=allowance()
            )
            self.assertEqual(restored.result.reconstruct(), 101)
            saved = restored.checkpoint["payload"]["config"]
            self.assertEqual(saved["memory_bytes"], config.memory_bytes)
            self.assertEqual(
                saved.get("ecm_chain_mode", "off"), config.ecm_chain_mode
            )

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
