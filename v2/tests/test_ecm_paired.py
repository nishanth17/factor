"""Independent coverage, affine, saturation and paired resume gates."""

import copy
import json
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2 import portfolio
from v2.common import arithmetic, utils
from v2.ecm.programs import ECMPrograms
from v2.execution.budget import BudgetExhaustedError
from v2.execution.schedules import SieveContext
from v2.execution.stage_jobs import advance_job, new_job, prime_cursor
from v2.tests.test_ecm_programs import allowance, configuration, primes_between
from v2.tests.test_phase_one import affine_add
from v2.tests.test_phase_two import reseal


def settings(**changes):
    return configuration(ecm_pair_distance=8, **changes)


def store(config):
    return ECMPrograms(
        SieveContext(config.max_hi, segment_size=config.segment_size),
        memory_bytes=config.ecm_program_bytes,
    )


def continuation(n, point, *, b1=13, b2=100, a24=2):
    job = new_job("ecm", n, 0, b1, b2)
    job.update(
        value=list(point),
        a24=a24,
        phase="stage_two_setup",
        cursor=prime_cursor(b1 + 1, b2 + 1),
    )
    return job


class PairedExecutionTests(unittest.TestCase):
    def test_admission_and_default_compatibility(self):
        for distance in (-2, 1, 3, 24):
            with self.assertRaises(ValueError):
                replace(
                    settings(),
                    ecm_pair_distance=distance,
                    ecm_tiers=((29, 1000, 1),),
                )
        with self.assertRaises(TypeError):
            replace(settings(), ecm_pair_distance=True)
        with self.assertRaises(ValueError):
            replace(settings(), ecm_program_bytes=0)
        with self.assertRaises(MemoryError):
            replace(settings(), memory_bytes=65536)
        self.assertIsNone(portfolio.PortfolioConfig().ecm_pair_distance)

    def test_exact_coverage_and_affine_products_including_tails(self):
        modulus, curve_a, point = 1009, 6, (3, 293)
        affine = [None]
        for _ in range(260):
            affine.append(affine_add(affine[-1], point, modulus, curve_a))

        for b1, b2, distance in ((2, 101, 0), (13, 127, 4), (29, 251, 12)):
            for width in (1, 7, 32):
                config = replace(
                    settings(),
                    ecm_pair_distance=distance,
                    ecm_tiers=((b1, b2, 1),),
                    segment_size=width,
                )
                context = store(config)
                job = continuation(modulus, (3, 1), b1=b1, b2=b2)
                covered = []
                for _ in range(10000):
                    before = len(job["terms"])
                    advance_job(job, allowance(), context, config)
                    for term, record in zip(
                        job["terms"][before:],
                        job.get("term_records", [])[before:],
                    ):
                        center, offset, minus, plus = record
                        covered.extend(p for p in (minus, plus) if p)
                        if offset and affine[center] and affine[offset]:
                            giant, baby = (
                                job["giant"],
                                job["baby"][offset // 2],
                            )
                            self.assertEqual(
                                term,
                                (affine[center][0] - affine[offset][0])
                                * giant[1]
                                * baby[1]
                                % modulus,
                            )
                            self.assertNotEqual(
                                giant[1] * baby[1] % modulus, 0
                            )
                        elif not offset:
                            self.assertEqual(term == 0, affine[center] is None)
                    if job["done"]:
                        break
                else:
                    self.fail("paired continuation did not terminate")
                self.assertEqual(
                    sorted(covered), primes_between(b1 + 1, b2 + 1)
                )
                self.assertEqual(len(covered), len(set(covered)))

    def test_single_pair_mixed_factor_saturation_recovers_both_signs(self):
        # Independent affine orders: (29,87) mod 151 has order 19;
        # (13,153) mod 367 has order 23, both on A=6. CRT forms Q.
        for modulus, point, order in (
            (151, (29, 87), 19),
            (367, (13, 153), 23),
        ):
            value = None
            for index in range(1, order + 1):
                value = affine_add(value, point, modulus, 6)
                self.assertEqual(value is None, index == order)
        n = 151 * 367
        x = 29 + 151 * ((13 - 29) * pow(151, -1, 367) % 367)
        config = replace(
            settings(),
            ecm_pair_distance=4,
            ecm_tiers=((13, 23, 1),),
            gcd_batch=1,
        )
        context = store(config)
        job = continuation(n, (x, 1), b2=23)
        phases = set()
        for _ in range(1000):
            before = copy.deepcopy(job)
            advance_job(job, allowance(), context, config)
            # Rebuilt program stores never own or replace curve points.
            advance_job(before, allowance(), store(config), config)
            self.assertEqual(job, before)
            phases.add(job["phase"])
            if job["done"]:
                break
        self.assertIn("pair_scalar_replay", phases)
        self.assertEqual(job["factor"], 151)
        self.assertTrue(utils.valid_divisor(job["factor"], n))

        # A segment tail can retain only the plus certificate of the same
        # saturated relation. Its missing side must not suppress recovery.
        tail = continuation(n, (x, 1), b2=23)
        tail.update(
            phase="pair_term_replay",
            terms=[0],
            product=0,
            term_records=[[21, 2, 0, 23]],
        )
        for _ in range(10):
            advance_job(tail, allowance(), store(config), config)
            if tail["done"]:
                break
        self.assertEqual(tail["factor"], 367)

    def test_every_action_atomic_and_store_rebuilt_without_point_sharing(self):
        config = settings()
        context = store(config)
        n = 1000000000039 * 1000000000061
        job = new_job("ecm", n, 17, 50, 1000)
        for _ in range(10000):
            before = copy.deepcopy(job)
            refused = copy.deepcopy(job)
            try:
                advance_job(refused, allowance(0), context, config)
            except BudgetExhaustedError:
                self.assertEqual(refused, before)
            else:
                # Empty stage transitions and completion need no arithmetic.
                self.assertTrue(
                    refused["done"] or refused["phase"] == "stage_two_setup"
                )
            advance_job(job, allowance(), context, config)
            advance_job(before, allowance(), store(config), config)
            self.assertEqual(job, before)
            if job["done"]:
                break
        self.assertTrue(job["done"])
        self.assertGreater(context.coverage_misses, 0)
        self.assertLessEqual(context.used_bytes, context.memory_bytes)
        other = new_job("ecm", n, 19, 50, 1000)
        while not other["done"]:
            advance_job(other, allowance(), context, config)
        self.assertGreater(context.coverage_hits, 0)
        self.assertIsNot(job["baby"], other["baby"])

    def test_program_regeneration_and_cancelled_cached_read(self):
        config = settings(ecm_program_bytes=4096 + 256 * 16)
        context = store(config)
        cursor = prime_cursor(51, 1001)
        values = context.program_segment(51, 83, allowance())
        cursor.update(left=51, next=83, values=values)
        first = context.coverage(
            cursor, b1=50, b2=1000, distance=8, budget=allowance()
        )
        second = context.coverage(
            cursor, b1=50, b2=1000, distance=8, budget=allowance()
        )
        self.assertEqual(first, second)
        self.assertFalse(context.coverage_blocks)
        self.assertEqual(context.coverage_misses, 2)
        budget = allowance()
        budget.cancelled = lambda: True
        with self.assertRaisesRegex(BudgetExhaustedError, "cancelled"):
            context.coverage(cursor, b1=50, b2=1000, distance=8, budget=budget)
        self.assertEqual(context.coverage_misses, 2)

    def test_partial_compilation_and_cached_cancellation_do_not_publish(self):
        config = settings()
        context = store(config)
        values = context.program_segment(51, 83, allowance())
        cursor = dict(left=51, next=83, hi=1001, values=values, index=0)
        used = context.used_bytes
        with self.assertRaises(BudgetExhaustedError):
            context.coverage(
                cursor,
                b1=50,
                b2=1000,
                distance=8,
                budget=allowance(2 * len(values)),
            )
        self.assertEqual(context.used_bytes, used)
        self.assertEqual(context.coverage_misses, 0)
        self.assertFalse(context.coverage_blocks)
        context.coverage(
            cursor, b1=50, b2=1000, distance=8, budget=allowance()
        )
        cancelled = allowance()
        cancelled.cancelled = lambda: True
        with self.assertRaises(BudgetExhaustedError):
            context.coverage(
                cursor, b1=50, b2=1000, distance=8, budget=cancelled
            )
        self.assertEqual(context.coverage_hits, 0)
        self.assertEqual(cancelled.used, 0)
        with self.assertRaisesRegex(ValueError, "one prime segment"):
            context.coverage(
                dict(cursor, next=1000),
                b1=50,
                b2=1000,
                distance=8,
                budget=allowance(),
            )


class PairedResumeTests(unittest.TestCase):
    def test_supported_gmp_matches_every_paired_action(self):
        try:
            backend = arithmetic.get_backend("gmpy2-mpz")
        except arithmetic.BackendUnavailableError as error:
            self.skipTest(str(error))
        n = 1000000000039 * 1000000000061
        configs = [settings(), settings(backend="gmpy2-mpz")]
        contexts = [store(config) for config in configs]
        jobs = [
            new_job("ecm", value, 17, 50, 1000)
            for value in (n, backend.integer(n))
        ]
        budgets = [allowance(), allowance()]
        for _ in range(10000):
            for job, budget, context, config in zip(
                jobs, budgets, contexts, configs
            ):
                advance_job(job, budget, context, config)
            self.assertEqual(jobs[0], arithmetic.canonical(jobs[1]))
            self.assertEqual(budgets[0].used, budgets[1].used)
            if jobs[0]["done"]:
                break
        self.assertTrue(jobs[0]["done"])

    def test_each_phase_checkpoint_cumulative_campaign_and_backend(self):
        backends = ["python-int"]
        try:
            arithmetic.get_backend("gmpy2-mpz")
            backends.append("gmpy2-mpz")
        except arithmetic.BackendUnavailableError:
            pass
        n = 1000000000039 * 1000000000061
        for backend in backends:
            config = settings(backend=backend)
            whole = portfolio.factorize_bounded(
                n, config=config, budget=allowance()
            )
            self.assertEqual(whole.reason, "exhausted")
            self.assertEqual(len(whole.events), 4)
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
                self.assertEqual(first.checkpoint["payload"]["version"], 7)
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
                self.assertEqual(resumed.result.reconstruct(), n)

    def test_resealed_coverage_table_and_identity_corruption_rejected(self):
        config = settings()
        n = 1000000000039 * 1000000000061
        advance = portfolio.advance_job

        def pause(job, budget, context, options):
            advance(job, budget, context, options)
            if job.get("pair_records") and job["terms"]:
                budget.cancelled = lambda: True

        with patch.object(portfolio, "advance_job", pause):
            first = portfolio.factorize_bounded(
                n, config=config, budget=allowance()
            )
        for corrupt_key in ("record", "index", "table", "product", "center"):
            bad = copy.deepcopy(first.checkpoint)
            job = bad["payload"]["state"]["current"]["job"]
            if corrupt_key == "record":
                job["pair_records"][0][2] += 2
            elif corrupt_key == "index":
                job["pair_index"] = 10**6
            elif corrupt_key == "table":
                job["baby"].append([1, 1])
            elif corrupt_key == "product":
                job["product"] += 1
            else:
                job["center"] += 1
            with self.assertRaises(ValueError):
                portfolio.factorize_bounded(
                    n,
                    config=config,
                    budget=allowance(),
                    checkpoint=reseal(bad),
                )
        for modified in (
            replace(config, ecm_pair_distance=4),
            replace(config, ecm_tiers=((51, 1000, 4),)),
            replace(config, ecm_tiers=((50, 1000, 5),)),
        ):
            with self.assertRaises(ValueError):
                portfolio.factorize_bounded(
                    n,
                    config=modified,
                    budget=allowance(),
                    checkpoint=first.checkpoint,
                )

    def test_complete_factoring_and_finite_budget_extension(self):
        n = 25013 * 25031
        config = settings()
        run = portfolio.factorize_bounded(
            n, config=config, budget=allowance(1000)
        )
        self.assertEqual(run.result.reconstruct(), n)
        for limit in range(2000, 100001, 1000):
            if run.reason in ("complete", "exhausted"):
                break
            run = portfolio.factorize_bounded(
                n,
                config=config,
                budget=allowance(limit),
                checkpoint=run.checkpoint,
            )
        whole = portfolio.factorize_bounded(
            n, config=config, budget=allowance()
        )
        self.assertEqual(run.result, whole.result)
        self.assertTrue(run.result.complete)
        self.assertEqual(run.result.reconstruct(), n)
