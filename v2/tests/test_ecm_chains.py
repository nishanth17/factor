"""Production scalar/coverage, budget, cache and legacy-resume gates."""

import copy
import unittest
from dataclasses import asdict, replace
from math import gcd, prod
from unittest.mock import patch

from v2 import arithmetic, ecm, prac, utils
from v2.benchmarks import c6_b4, c6_fast
from v2.benchmarks.p41_campaign import PYTHON_BACKEND
from v2.benchmarks.prac_oracle import (
    affine_multiply,
    historical_points,
    matches,
)
from v2.budget import Budget, BudgetExhaustedError
from v2.ecm_chain_records import Record, verify_frontier
from v2.ecm_chains import (
    MIN_MEMORY_BYTES,
    SCRATCH_BYTES,
    ChainPlan,
    ChainPlans,
    identity,
    supports_modulus,
    verify_progress,
)
from v2.ecm_programs import ECMPrograms
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.schedules import SieveContext
from v2.stage_jobs import advance_job, new_job
from v2.tests.test_phase_two import reseal


def allowance(work=10**8):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


def configuration(**changes):
    return replace(
        PortfolioConfig(
            trial_bound=5,
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=((2000, 2500, 8),),
            segment_size=128,
            max_input_bits=512,
            memory_bytes=32 * 1024**2,
            ecm_program_bytes=262144,
            ecm_chain_mode="reuse",
            ecm_chain_bytes=MIN_MEMORY_BYTES,
        ),
        **changes,
    )


class ProductionChainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.native = ChainPlan(2000, "python-int", allowance())
        cls.catalog = c6_fast.load_catalog()
        cls.reference = c6_b4.build_program(
            "prac", "reduced", PYTHON_BACKEND, 16
        )

    def test_all_records_retain_independent_scalar_and_factor_coverage(self):
        for prime, (power, action, unit, _) in self.native.entries.items():
            self.assertEqual(action.record.scalar, power)
            self.assertEqual(unit.record.scalar, prime)
            certified = self.catalog["prac"][power]
            self.assertEqual(action.record.code, certified.record.code)
            self.assertEqual(action.masks, certified.masks)
            self.assertTrue(verify_frontier(action.record, action.masks))
            broken = Record(
                power + 1,
                action.record.code,
                action.record.output,
                action.record.slots,
            )
            with self.assertRaises(ValueError):
                verify_frontier(broken, action.masks)

    def test_every_record_matches_independent_affine_and_c6_kernels(self):
        for prime, (power, action, _, _) in self.native.entries.items():
            for modulus, curve_a, point in historical_points():
                a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
                actual, _ = action.run((point[0], 1), modulus, a24)
                # A coverage-unit aggregate certifies equality; exceptional
                # outputs must take the same finite strict recovery as C6.
                expected = affine_multiply(power, point, modulus, curve_a)
                guarded = self.native.execute(
                    [(prime, power)], (point[0], 1), modulus, a24, allowance()
                )
                if guarded[0] is not None:
                    self.assertTrue(matches(guarded[0], expected, modulus))
                reference_action = next(
                    row[2] for row in self.reference.entries if row[0] == prime
                )
                self.assertEqual(
                    actual,
                    reference_action.run((point[0], 1), modulus, a24)[0],
                )

    def test_atomic_refusal_and_saturated_mixed_factor_recovery(self):
        powers = [
            [prime, row[0]]
            for prime, row in list(self.native.entries.items())[:16]
        ]
        setup = ecm.setup_curve(1000003, 17)
        original = list(setup.point)
        budget = allowance(1)
        with self.assertRaises(BudgetExhaustedError):
            self.native.execute(powers, original, 1000003, setup.a24, budget)
        self.assertEqual(list(setup.point), original)
        self.assertEqual(budget.used, 0)

        # Product = n cannot hide proper factors on different old coordinates.
        action = self.native.entries[2][1]
        with (
            patch.object(action, "run", return_value=((1, 1), 35)),
            patch.object(
                action, "strict", side_effect=prac.NonunitPointError(5)
            ),
        ):
            value, divisor, replayed = self.native.execute(
                [[2, self.native.entries[2][0]]], (2, 1), 35, 2, allowance()
            )
        self.assertIsNone(value)
        self.assertEqual(divisor, 5)
        self.assertTrue(replayed)

    def test_catalog_corruption_missing_schedule_and_caps(self):
        with patch("v2.ecm_chains.CATALOG_SHA256", "0" * 64):
            budget = allowance()
            with self.assertRaisesRegex(ValueError, "catalog identity"):
                ChainPlan(2000, "python-int", budget)
            self.assertGreater(budget.used, 0)
        with self.assertRaises(ValueError):
            self.native.execute([[2, 2]], (2, 1), 101, 2, allowance())
        with self.assertRaises(MemoryError):
            configuration(ecm_chain_bytes=4096)
        with self.assertRaises(ValueError):
            configuration(ecm_chain_mode="off")

    def test_cache_hit_eviction_cancellation_and_failed_build(self):
        store = ChainPlans(MIN_MEMORY_BYTES, "python-int", (1999, 2000))
        budget = allowance()
        first = store.get(2000, budget)
        preparation = budget.used
        self.assertIs(first, store.get(2000, budget))
        self.assertEqual(budget.used, preparation + 1)
        store.get(1999, budget)
        self.assertEqual(store.evictions, 1)
        self.assertLessEqual(store.used_bytes, store.memory_bytes)
        store.get(2000, budget)
        self.assertEqual(store.evictions, 2)
        self.assertIsNone(store.get(2001, allowance()))
        with self.assertRaises(BudgetExhaustedError):
            store.get(2000, Budget(cancelled=lambda: True))
        failed = ChainPlans(MIN_MEMORY_BYTES, "python-int", (2000,))
        with self.assertRaises(BudgetExhaustedError):
            failed.get(2000, allowance(202462))
        self.assertFalse(failed.plans)
        self.assertEqual(failed.used_bytes, SCRATCH_BYTES)

    def test_actual_jobs_preserve_factors_replay_and_fallback(self):
        config = configuration()
        for n, seed in (
            (1009 * 1013, 6),
            (1000000000039 * 1000000000061, 17),
            (33554520197234177, 2046841451),
            (6120168563605791616423380424731852610871, 17),
        ):
            context = SieveContext(
                config.max_hi, segment_size=config.segment_size
            )
            store = ECMPrograms(context, memory_bytes=config.ecm_program_bytes)
            store.chains = ChainPlans(MIN_MEMORY_BYTES, "python-int", (2000,))
            job, budget = new_job("ecm", n, seed, 2000, 2500), allowance()
            for _ in range(10000):
                advance_job(job, budget, store, config)
                if job["done"]:
                    break
            else:
                self.fail("production chain job failed to terminate")
            if job["factor"] is not None:
                self.assertTrue(utils.valid_divisor(job["factor"], n))
            if "chain_identity" in job:
                self.assertEqual(
                    job["chain_identity"], identity(2000, "python-int")
                )

        fallback = configuration(ecm_tiers=((2001, 2500, 8),))
        n = 1000000000039 * 1000000000061
        candidate = factorize_bounded(
            n, seed=7, config=fallback, budget=allowance()
        )
        control = factorize_bounded(
            n,
            seed=7,
            config=replace(fallback, ecm_chain_mode="off", ecm_chain_bytes=0),
            budget=allowance(),
        )
        self.assertEqual(candidate.result, control.result)
        self.assertEqual(
            candidate.checkpoint["payload"]["work_used"],
            control.checkpoint["payload"]["work_used"],
        )

    def test_pause_rebuild_cumulative_allowance_identity_and_legacy(self):
        config, n = configuration(), 6120168563605791616423380424731852610871
        full = factorize_bounded(n, seed=19, config=config, budget=allowance())
        self.assertEqual(full.result.reconstruct(), n)
        partial = factorize_bounded(
            n, seed=19, config=config, budget=allowance(250000)
        )
        self.assertEqual(partial.reason, "work_limit")
        self.assertEqual(partial.checkpoint["payload"]["version"], 10)
        resumed = factorize_bounded(
            n, checkpoint=partial.checkpoint, config=config, budget=allowance()
        )
        self.assertEqual(resumed.result.reconstruct(), n)
        self.assertEqual(resumed.result, full.result)
        self.assertGreater(
            resumed.checkpoint["payload"]["work_used"],
            full.checkpoint["payload"]["work_used"],
        )
        corrupt = copy.deepcopy(partial.checkpoint)
        corrupt["payload"]["chains"] += "wrong"
        reseal(corrupt)
        with self.assertRaisesRegex(ValueError, "incompatible"):
            factorize_bounded(
                n, checkpoint=corrupt, config=config, budget=allowance()
            )
        for field, value in (
            ("value", [0, 0]),
            ("chain_chunks", 999),
            ("chain_identity", "wrong"),
        ):
            corrupt = copy.deepcopy(partial.checkpoint)
            corrupt["payload"]["state"]["current"]["job"][field] = value
            reseal(corrupt)
            with self.assertRaises(ValueError):
                factorize_bounded(
                    n, checkpoint=corrupt, config=config, budget=allowance()
                )
        off = replace(config, ecm_chain_mode="off", ecm_chain_bytes=0)
        old = factorize_bounded(
            n, seed=19, config=off, budget=allowance(10000)
        )
        self.assertNotIn("chains", old.checkpoint["payload"])
        self.assertNotIn("ecm_chain_mode", old.checkpoint["payload"]["config"])
        self.assertLess(old.checkpoint["payload"]["version"], 10)
        restored = factorize_bounded(
            n, checkpoint=old.checkpoint, config=off, budget=allowance()
        )
        self.assertEqual(restored.result.reconstruct(), n)

    def test_separate_gmp_tuple_lucas_and_native_policy(self):
        try:
            backend = arithmetic.get_backend("gmpy2-mpz")
        except arithmetic.BackendUnavailableError:
            self.skipTest("optional GMP backend unavailable")
        plan = ChainPlan(2000, backend.name, allowance())
        prime, (power, action, _, _) = next(iter(plan.entries.items()))
        self.assertEqual(
            action.record.code, self.catalog["lucas"][power].record.code
        )
        n = backend.integer(1000003)
        setup = ecm.setup_curve(n, 17)
        point, divisor, _ = plan.execute(
            [[prime, power]], setup.point, n, setup.a24, allowance()
        )
        if point is not None:
            expected = ecm.scalar_multiply(power, *setup.point, n, setup.a24)
            self.assertEqual(
                (point[0] * expected[1] - point[1] * expected[0]) % n, 0
            )
            self.assertEqual(gcd(int(point[1]), int(n)), 1)
        self.assertTrue(divisor is None or utils.valid_divisor(divisor, n))
        self.assertEqual(
            prod(
                [
                    f.value**f.exponent
                    for f in factorize_bounded(
                        1009 * 1013,
                        seed=7,
                        config=configuration(backend=backend.name),
                        budget=allowance(),
                    ).result.factors
                ]
            ),
            1009 * 1013,
        )

    def test_committed_mainline_bidirectional_ladder_resume(self):
        from v2.benchmarks.b3_production import baseline

        old = baseline()
        config = configuration(ecm_chain_mode="off", ecm_chain_bytes=0)
        options = asdict(config)
        options.pop("ecm_chain_mode")
        options.pop("ecm_chain_bytes")
        old_config = old.PortfolioConfig(**options)
        n = 1000000000039 * 1000000000061
        first = old.factorize_bounded(
            n,
            seed=19,
            config=old_config,
            budget=old.Budget(
                work_limit=10000, seconds=None, cpu_seconds=None
            ),
        )
        restored = factorize_bounded(
            n, config=config, checkpoint=first.checkpoint, budget=allowance()
        )
        self.assertEqual(restored.result.reconstruct(), n)
        second = factorize_bounded(
            n, seed=19, config=config, budget=allowance(10000)
        )
        self.assertEqual(
            {
                k: v
                for k, v in first.checkpoint["payload"]["state"].items()
                if k != "stage_seconds"
            },
            {
                k: v
                for k, v in second.checkpoint["payload"]["state"].items()
                if k != "stage_seconds"
            },
        )
        old_restored = old.factorize_bounded(
            n,
            config=old_config,
            checkpoint=second.checkpoint,
            budget=old.Budget(
                work_limit=10**8, seconds=None, cpu_seconds=None
            ),
        )
        self.assertEqual(old_restored.result.reconstruct(), n)
        self.assertEqual(restored.work_used, old_restored.work_used)

    def test_supported_size_band_uses_exact_boundaries(self):
        self.assertFalse(supports_modulus(10**39 - 1))
        self.assertTrue(supports_modulus(10**39))
        self.assertTrue(supports_modulus(10**80 - 1))
        self.assertFalse(supports_modulus(10**80))
        config = configuration()
        n = 1000000000039 * 1000000000061
        with patch(
            "v2.ecm_chains.ChainPlans.get",
            side_effect=AssertionError("unsupported modulus prepared chains"),
        ):
            candidate = factorize_bounded(
                n, seed=19, config=config, budget=allowance()
            )
        control = factorize_bounded(
            n,
            seed=19,
            config=replace(config, ecm_chain_mode="off", ecm_chain_bytes=0),
            budget=allowance(),
        )
        self.assertEqual(candidate.result, control.result)
        self.assertEqual(candidate.work_used, control.work_used)

    def test_durable_replay_progress_rejects_corrupt_recovery(self):
        config = configuration()
        n = 6120168563605791616423380424731852610871
        context = SieveContext(config.max_hi, segment_size=config.segment_size)
        store = ECMPrograms(context, memory_bytes=config.ecm_program_bytes)
        store.chains = ChainPlans(MIN_MEMORY_BYTES, "python-int", (2000,))
        job = new_job("ecm", n, 17, 2000, 2500)
        with patch.object(
            ChainPlan, "execute", return_value=(None, None, True)
        ):
            for _ in range(100):
                advance_job(job, allowance(), store, config)
                if job["phase"] == "replay":
                    break
            else:
                self.fail("saved-chunk replay was not reached")
        verify_progress(job, "python-int", context)
        advance_job(job, allowance(), store, config)
        self.assertFalse(job["done"])
        verify_progress(job, "python-int", context)
        for field, value in (
            ("replay_index", 17),
            ("replay_power", 3),
            ("replay_value", [1, 0]),
        ):
            corrupt = copy.deepcopy(job)
            corrupt[field] = value
            with self.assertRaises(ValueError):
                verify_progress(corrupt, "python-int", context)
