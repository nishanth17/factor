"""Independent wheel-cell coverage, affine products and resumable execution."""

import copy
import json
import unittest
from dataclasses import replace
from math import gcd
from unittest.mock import patch

from v2 import arithmetic, ecm_paired, ecm_wheel, portfolio, utils
from v2.budget import BudgetExhaustedError
from v2.stage_jobs import advance_job, new_job, prime_cursor
from v2.tests.test_ecm_paired import continuation, store
from v2.tests.test_ecm_programs import allowance, configuration, primes_between
from v2.tests.test_phase_one import affine_add
from v2.tests.test_phase_two import reseal


def settings(**changes):
    return configuration(ecm_pair_wheel=30, **changes)


class WheelExecutionTests(unittest.TestCase):
    def test_mixed_factor_recovery_and_plus_only_tail_survive_resume(self):
        # Independent affine order certificates on A=6: x=29 mod 151
        # has order 19 and x=55 mod 163 has order 41. Their shared center
        # is 30, so one offset-11 product saturates over the CRT modulus.
        for modulus, point, order in ((151, (29, 87), 19), (163, (55, 8), 41)):
            value = None
            for index in range(1, order + 1):
                value = affine_add(value, point, modulus, 6)
                self.assertEqual(value is None, index == order)
        n = 151 * 163
        x = 29 + 151 * ((55 - 29) * pow(151, -1, 163) % 163)
        backends = ["python-int"]
        try:
            arithmetic.get_backend("gmpy2-mpz")
            backends.append("gmpy2-mpz")
        except arithmetic.BackendUnavailableError:
            pass
        for backend in backends:
            integer = arithmetic.get_backend(backend).integer
            for b1, expected in ((18, 151), (30, 163)):
                config = settings(
                    backend=backend, ecm_tiers=((b1, 41, 1),), gcd_batch=1
                )
                context = store(config)
                job = continuation(
                    integer(n),
                    (integer(x), integer(1)),
                    b1=b1,
                    b2=41,
                    a24=integer(2),
                )
                phases = set()
                while not job["done"]:
                    resumed = arithmetic.canonical(copy.deepcopy(job))
                    from v2.stage_jobs import promote_job

                    promote_job(resumed, arithmetic.get_backend(backend))
                    advance_job(job, allowance(), context, config)
                    advance_job(resumed, allowance(), store(config), config)
                    self.assertEqual(job, resumed)
                    ecm_paired.verify_progress(
                        arithmetic.canonical(job), config, context.context
                    )
                    phases.add(job["phase"])
                self.assertIn("pair_scalar_replay", phases)
                self.assertEqual(job["factor"], expected)
                self.assertTrue(utils.valid_divisor(job["factor"], n))

    def test_admission_and_legacy_identity(self):
        for wheel in (0, 1, 3, 34):
            with self.assertRaises(ValueError):
                replace(settings(), ecm_pair_wheel=wheel)
        with self.assertRaises(TypeError):
            replace(settings(), ecm_pair_wheel=True)
        with self.assertRaises(ValueError):
            settings(ecm_pair_distance=8)
        with self.assertRaises(ValueError):
            settings(ecm_program_bytes=0)
        with self.assertRaises(MemoryError):
            settings(memory_bytes=65536)
        for distance, version in ((None, 5), (8, 7)):
            config = configuration(ecm_pair_distance=distance)
            run = portfolio.factorize_bounded(
                1000000000039 * 1000000000061,
                config=config,
                budget=allowance(10),
            )
            self.assertEqual(run.checkpoint["payload"]["version"], version)
            self.assertNotIn(
                "ecm_pair_wheel", run.checkpoint["payload"]["config"]
            )

    def test_every_prime_once_and_no_pair_split_at_segment_boundaries(self):
        for wheel in (2, 6, 30, 210):
            for width in (wheel // 2, wheel // 2 + 7, 233):
                for b1, b2 in ((2, 251), (13, 223), (50, 251), (211, 211)):
                    config = replace(
                        settings(),
                        ecm_pair_wheel=wheel,
                        segment_size=width,
                        ecm_tiers=((b1, b2, 1),),
                    )
                    context = store(config)
                    cursor = prime_cursor(b1 + 1, b2 + 1)
                    records = []
                    while (
                        ecm_wheel.peek_prime(
                            cursor, context, allowance(), wheel
                        )
                        is not None
                    ):
                        records.extend(
                            context.coverage(
                                cursor,
                                b1=b1,
                                b2=b2,
                                distance=wheel // 2,
                                wheel=wheel,
                                budget=allowance(),
                            )
                        )
                        cursor["index"] = len(cursor["values"])
                    primes = primes_between(b1 + 1, b2 + 1)
                    self.assertEqual(
                        sorted(p for r in records for p in r[2:] if p), primes
                    )
                    # Independent exhaustive pair oracle: use the affine
                    # symmetry equation and wheel grid, not the compiler.
                    expected = {
                        (p, q)
                        for p in primes
                        for q in primes
                        if p < q
                        and (p + q) % (2 * wheel) == 0
                        and q - p < wheel
                        and gcd(p, wheel) == 1
                    }
                    actual = {(r[2], r[3]) for r in records if r[2] and r[3]}
                    self.assertEqual(actual, expected)
                    self.assertLessEqual(
                        context.used_bytes, context.memory_bytes
                    )

    def test_sparse_points_and_products_against_affine_oracle(self):
        modulus, curve_a, point = 1009, 6, (3, 293)
        affine = [None]
        for _ in range(420):
            affine.append(affine_add(affine[-1], point, modulus, curve_a))
        for wheel in (2, 6, 30, 210):
            config = replace(
                settings(),
                ecm_pair_wheel=wheel,
                segment_size=max(16, wheel // 2),
                ecm_tiers=((2, 251, 1),),
            )
            context = store(config)
            job = continuation(modulus, (3, 1), b1=2, b2=251)
            covered = []
            while not job["done"]:
                start = len(job["terms"])
                advance_job(job, allowance(), context, config)
                ecm_paired.verify_progress(job, config, context.context)
                for term, record in zip(
                    job["terms"][start:], job["term_records"][start:]
                ):
                    center, offset, minus, plus = record
                    covered.extend(p for p in (minus, plus) if p)
                    if offset:
                        baby = job["baby"][job["baby_offsets"].index(offset)]
                        giant = job["giant"]
                        if affine[center] and affine[offset]:
                            self.assertEqual(
                                giant[0],
                                affine[center][0] * giant[1] % modulus,
                            )
                            self.assertEqual(
                                baby[0], affine[offset][0] * baby[1] % modulus
                            )
                            self.assertEqual(
                                term,
                                (affine[center][0] - affine[offset][0])
                                * giant[1]
                                * baby[1]
                                % modulus,
                            )
                    else:
                        self.assertEqual(term == 0, affine[center] is None)
            self.assertEqual(sorted(covered), primes_between(3, 252))
            self.assertEqual(
                job["baby_offsets"],
                [d for d in range(1, wheel // 2 + 1) if gcd(d, wheel) == 1],
            )

    def test_every_action_resume_refusal_cancellation_and_regeneration(self):
        config = settings(ecm_program_bytes=4096 + 256 * 16)
        context = store(config)
        n = 1000000000039 * 1000000000061
        job = new_job("ecm", n, 17, 50, 1000)
        while not job["done"]:
            before = copy.deepcopy(job)
            cancelled = allowance()
            cancelled.cancelled = lambda: True
            for refused_budget in (allowance(0), cancelled):
                refused = copy.deepcopy(job)
                try:
                    advance_job(refused, refused_budget, context, config)
                except BudgetExhaustedError:
                    self.assertEqual(refused, before)
            advance_job(job, allowance(), context, config)
            advance_job(before, allowance(), store(config), config)
            self.assertEqual(job, before)
            if job["phase"].startswith("pair_"):
                ecm_paired.verify_progress(job, config, context.context)
        self.assertTrue(context.unretained)
        self.assertFalse(context.coverage_blocks)
        self.assertTrue(
            job["factor"] is None or utils.valid_divisor(job["factor"], n)
        )

    def test_gmp_actions_and_canonical_campaign_resume(self):
        backends = ["python-int"]
        try:
            arithmetic.get_backend("gmpy2-mpz")
            backends.append("gmpy2-mpz")
        except arithmetic.BackendUnavailableError:
            pass
        n = 1000000000039 * 1000000000061
        configs = [settings(backend=backend) for backend in backends]
        contexts = [store(config) for config in configs]
        jobs = [
            new_job(
                "ecm",
                arithmetic.get_backend(c.backend).integer(n),
                17,
                50,
                1000,
            )
            for c in configs
        ]
        budgets = [allowance() for _ in configs]
        while not jobs[0]["done"]:
            for job, budget, context, config in zip(
                jobs, budgets, contexts, configs
            ):
                advance_job(job, budget, context, config)
            self.assertEqual(jobs[0], arithmetic.canonical(jobs[-1]))
            self.assertEqual(budgets[0].used, budgets[-1].used)
        for config in configs:
            whole = portfolio.factorize_bounded(
                n, config=config, budget=allowance()
            )
            for phase in ("stage_one", "pair_baby", "pair_terms"):
                advance = portfolio.advance_job

                def pause(job, budget, context, options):
                    advance(job, budget, context, options)
                    if job["phase"] == phase:
                        budget.cancelled = lambda: True

                with patch.object(portfolio, "advance_job", pause):
                    first = portfolio.factorize_bounded(
                        n, config=config, budget=allowance()
                    )
                self.assertEqual(first.reason, "cancelled")
                self.assertEqual(first.checkpoint["payload"]["version"], 8)
                resumed = portfolio.factorize_bounded(
                    n,
                    config=config,
                    budget=allowance(),
                    checkpoint=json.loads(json.dumps(first.checkpoint)),
                )
                self.assertEqual(resumed.result, whole.result)
                self.assertEqual(resumed.reason, whole.reason)
                self.assertEqual(
                    [e["seed"] for e in resumed.events],
                    [e["seed"] for e in whole.events],
                )
                self.assertGreaterEqual(resumed.work_used, whole.work_used)

    def test_resealed_wheel_metadata_corruption_is_rejected(self):
        config = settings()
        n = 1000000000039 * 1000000000061
        advance = portfolio.advance_job

        def pause(job, budget, context, options):
            advance(job, budget, context, options)
            if job.get("pair_records") and job["terms"]:
                budget.cancelled = lambda: True

        with patch.object(portfolio, "advance_job", pause):
            run = portfolio.factorize_bounded(
                n, config=config, budget=allowance()
            )
        for key in ("baby_offsets", "baby_next", "pair_records", "center"):
            bad = copy.deepcopy(run.checkpoint)
            job = bad["payload"]["state"]["current"]["job"]
            if key == "baby_offsets":
                job[key][0] = True
            elif key == "pair_records":
                job[key][0][2] += 2
            else:
                job[key] += 1
            with self.assertRaises(ValueError):
                portfolio.factorize_bounded(
                    n,
                    config=config,
                    budget=allowance(),
                    checkpoint=reseal(bad),
                )
        for changed in (
            replace(config, ecm_pair_wheel=6),
            replace(config, ecm_tiers=((51, 1000, 4),)),
        ):
            with self.assertRaises(ValueError):
                portfolio.factorize_bounded(
                    n,
                    config=changed,
                    budget=allowance(),
                    checkpoint=run.checkpoint,
                )

    def test_finite_work_extension_completes_and_reconstructs(self):
        n, config = 25013 * 25031, settings()
        run = portfolio.factorize_bounded(
            n, config=config, budget=allowance(1000)
        )
        for limit in range(2000, 100001, 1000):
            if run.reason in ("complete", "exhausted"):
                break
            run = portfolio.factorize_bounded(
                n,
                config=config,
                budget=allowance(limit),
                checkpoint=run.checkpoint,
            )
        self.assertTrue(run.result.complete)
        self.assertEqual(run.result.reconstruct(), n)
