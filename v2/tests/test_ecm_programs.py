"""Independent schedule, coverage, accounting and resume gates for P5.2 A3."""

import copy
import io
import unittest
from contextlib import redirect_stdout
from dataclasses import FrozenInstanceError, replace
from math import gcd, isqrt, prod
from struct import pack
from unittest.mock import patch

from v2.common import utils
from v2.ecm import core as ecm
from v2.ecm.programs import ECMPrograms, ProgramBlock, pair_coverage
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.schedules import ScheduleCache, SieveContext
from v2.execution.stage_jobs import (
    advance_job,
    new_job,
    peek_prime,
    prime_cursor,
)
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.tests.test_phase_one import affine_add
from v2.tests.test_phase_two import reseal


def allowance(work=10**8):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


def primes_between(lo, hi):
    """Trial division supplies an oracle independent of segmented schedules."""
    return [
        value
        for value in range(max(2, lo), hi)
        if all(value % divisor for divisor in range(2, isqrt(value) + 1))
    ]


def configuration(**changes):
    return replace(
        PortfolioConfig(
            trial_bound=5,
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=((50, 1000, 4),),
            segment_size=16,
            max_input_bits=128,
            ecm_program_bytes=131072,
        ),
        **changes,
    )


class ProgramTests(unittest.TestCase):
    def test_exact_powers_and_immutable_blocks(self):
        for bound in (2, 3, 8, 16, 50, 243):
            for width in (1, 7, 16):
                store = ECMPrograms(
                    SieveContext(bound + 1, segment_size=width),
                    memory_bytes=131072,
                )
                powers, primes = [], []
                for lo in range(2, bound + 1, 2 * width):
                    hi = min(bound + 1, lo + 2 * width)
                    values = store.program_segment(
                        lo, hi, allowance(), bound=bound
                    )
                    primes.extend(values)
                    powers.extend(
                        store.power_values(lo, hi, bound, 0, len(values))
                    )

                expected = 1
                for value in range(1, bound + 1):
                    expected = expected // gcd(expected, value) * value
                self.assertEqual(primes, primes_between(2, bound + 1))
                self.assertEqual(prod(powers), expected)
                block = store.last_block
                with self.assertRaises(FrozenInstanceError):
                    block.lo = 0
                self.assertIsInstance(block.primes, bytes)
                self.assertIsInstance(block.powers, bytes)

    def test_refusal_does_not_publish_or_advance(self):
        context = SieveContext(100, segment_size=16)
        store = ECMPrograms(context, memory_bytes=131072)
        cursor = prime_cursor(2, 51)
        original = copy.deepcopy(cursor)
        budget = allowance(16 + len(context.base_primes))

        with self.assertRaises(BudgetExhaustedError):
            peek_prime(cursor, store, budget, power_bound=50)

        self.assertEqual(cursor, original)
        self.assertFalse(store.blocks)
        self.assertIsNone(store.last_block)
        self.assertEqual(store.misses, 0)
        self.assertGreater(budget.used, 0)
        budget.work_limit += 100
        self.assertEqual(peek_prime(cursor, store, budget, power_bound=50), 2)

    def test_caps_reuse_and_regeneration(self):
        context = SieveContext(500, segment_size=16)
        minimum = 4096 + 256 * context.segment_size
        store = ECMPrograms(context, memory_bytes=minimum + 700)
        costs = []
        for _ in range(2):
            budget = allowance()
            result = []
            for lo in range(2, 500, 32):
                result.extend(
                    store.program_segment(lo, min(lo + 32, 500), budget)
                )
            costs.append(budget.used)

            self.assertEqual(result, primes_between(2, 500))
            self.assertLessEqual(store.used_bytes, store.memory_bytes)
        self.assertEqual(len(store.blocks), 1)
        self.assertGreater(store.unretained, 0)
        self.assertEqual(store.hits, 1)
        self.assertLess(costs[1], costs[0])

    def test_schedule_cache_composition_and_packed_limits(self):
        context = ScheduleCache(SieveContext(100, segment_size=16))
        store = ECMPrograms(context, memory_bytes=131072)
        self.assertEqual(
            store.program_segment(2, 34, allowance()), primes_between(2, 34)
        )
        self.assertEqual(
            store.program_segment(2, 34, allowance()), primes_between(2, 34)
        )
        self.assertEqual(store.hits, 1)
        with self.assertRaises(ValueError):
            store.program_segment(2**64, 2**64 + 1, allowance())
        with self.assertRaises(MemoryError):
            ECMPrograms(context, memory_bytes=4096)
        with self.assertRaises(MemoryError):
            configuration(ecm_program_bytes=4096)
        with self.assertRaises(ValueError):
            configuration(ecm_tiers=((50, 2**64, 1),))
        with self.assertRaises(TypeError):
            configuration(ecm_program_bytes=True)

    def test_wide_words_and_cancelled_cached_reads(self):
        lo = 2**32 - 16
        context = SieveContext(lo + 64, segment_size=16)
        store = ECMPrograms(context, memory_bytes=131072)
        actual = store.program_segment(lo, lo + 32, allowance())
        self.assertEqual(actual, primes_between(lo, lo + 32))
        self.assertTrue(any(value >= 2**32 for value in actual))
        before = (store.hits, store.misses, store.used_bytes)
        cancelled = Budget(cancelled=lambda: True)

        with self.assertRaisesRegex(BudgetExhaustedError, "cancelled"):
            store.program_segment(lo, lo + 32, cancelled)

        self.assertEqual((store.hits, store.misses, store.used_bytes), before)
        self.assertEqual(cancelled.used, 0)

    def test_each_candidate_action_matches_streamed_control(self):
        for b1, b2 in ((2, 100), (10, 200), (50, 1000)):
            for seed, n in (
                (0, 1000003),
                (6, 1009 * 1013),
                (17, 1000000000039 * 1000000000061),
            ):
                config = configuration(ecm_tiers=((b1, b2, 3),))
                context = SieveContext(config.max_hi, segment_size=16)
                store = ECMPrograms(context, memory_bytes=131072)
                budgets = [allowance(), allowance()]
                jobs = [new_job("ecm", n, seed, b1, b2) for _ in range(2)]

                for _ in range(10000):
                    advance_job(jobs[0], budgets[0], context, config)
                    advance_job(jobs[1], budgets[1], store, config)
                    self.assertEqual(jobs[0], jobs[1])
                    if jobs[0]["done"]:
                        break
                else:
                    self.fail("finite candidate did not terminate")
                divisor = jobs[0]["factor"]
                if divisor is not None:
                    self.assertTrue(utils.valid_divisor(divisor, n))


class CoverageTests(unittest.TestCase):
    def test_independent_coverage_initialization_exceptions_and_tails(self):
        for b1, b2, distance in (
            (2, 101, 0),
            (11, 127, 4),
            (29, 503, 12),
            (50, 997, 24),
        ):
            store = ECMPrograms(
                SieveContext(b2 + 1, segment_size=16), memory_bytes=131072
            )
            covered = []
            paired = 0
            for lo in range(b1 + 1, b2 + 1, 32):
                hi = min(lo + 32, b2 + 1)
                store.program_segment(lo, hi, allowance())
                coverage = pair_coverage(
                    store.last_block,
                    b1=b1,
                    b2=b2,
                    distance=distance,
                    memory_bytes=16384,
                    budget=allowance(),
                )
                for center, offset, minus, plus in coverage.records():
                    if offset:
                        self.assertGreater(center - 2 * distance, 0)
                        self.assertEqual(offset % 2, 0)
                        self.assertLessEqual(offset, distance)
                        self.assertTrue(not minus or minus == center - offset)
                        self.assertTrue(not plus or plus == center + offset)
                    else:
                        self.assertEqual(minus, center)
                        self.assertEqual(plus, 0)
                    covered.extend(value for value in (minus, plus) if value)
                    paired += bool(minus and plus)

            self.assertEqual(sorted(covered), primes_between(b1 + 1, b2 + 1))
            self.assertEqual(len(covered), len(set(covered)))
            if distance:
                self.assertGreater(paired, 0)

    def test_projective_cross_products_against_affine_oracle(self):
        modulus, curve_a = 1009, 6
        point = (3, 293)
        self.assertEqual(
            point[1] ** 2 % modulus,
            (point[0] ** 3 + curve_a * point[0] ** 2 + point[0]) % modulus,
        )
        points = [None]
        for _ in range(100):
            points.append(affine_add(points[-1], point, modulus, curve_a))
        a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus

        for center in (11, 19, 27, 35):
            for offset in (2, 4):
                giant = ecm.scalar_multiply(center, 3, 1, modulus, a24)
                baby = ecm.scalar_multiply(offset, 3, 1, modulus, a24)
                if points[center] is None or points[offset] is None:
                    continue
                self.assertNotEqual(giant[1], 0)
                self.assertNotEqual(baby[1], 0)
                cross = (giant[0] * baby[1] - baby[0] * giant[1]) % modulus
                affine_difference = (
                    points[center][0] - points[offset][0]
                ) % modulus
                self.assertEqual(
                    cross, affine_difference * giant[1] * baby[1] % modulus
                )
                self.assertEqual(
                    cross == 0, points[center][0] == points[offset][0]
                )

    def test_invalid_initialization_caps_and_duplicate_certificates(self):
        block = ProgramBlock(
            12, 30, None, b"".join(pack("<Q", p) for p in (13, 17, 19)), b""
        )
        for distance in (1, 3, 6):
            with self.assertRaises(ValueError):
                pair_coverage(
                    block,
                    b1=11,
                    b2=29,
                    distance=distance,
                    memory_bytes=16384,
                    budget=allowance(),
                )
        with self.assertRaises(MemoryError):
            pair_coverage(
                block,
                b1=11,
                b2=29,
                distance=4,
                memory_bytes=4096,
                budget=allowance(),
            )
        duplicate = replace(block, primes=pack("<QQ", 13, 13))
        with self.assertRaises(ValueError):
            pair_coverage(
                duplicate,
                b1=11,
                b2=29,
                distance=0,
                memory_bytes=16384,
                budget=allowance(),
            )


class ResumeTests(unittest.TestCase):
    def test_quiet_reconstruction_and_cumulative_resume(self):
        config = configuration()
        n = 25013 * 25031
        with redirect_stdout(io.StringIO()) as output:
            whole = factorize_bounded(
                n, seed=7, config=config, budget=allowance()
            )
            run = factorize_bounded(
                n, seed=7, config=config, budget=allowance(1000)
            )
            self.assertEqual(run.checkpoint["payload"]["version"], 5)
            self.assertGreater(run.work_used, 0)
            for limit in range(2000, 100001, 1000):
                if run.reason in ("complete", "exhausted"):
                    break
                previous = run.work_used
                run = factorize_bounded(
                    n,
                    config=config,
                    budget=allowance(limit),
                    checkpoint=run.checkpoint,
                )
                self.assertGreaterEqual(run.work_used, previous)
                self.assertEqual(run.result.reconstruct(), n)
            else:
                self.fail("resumed campaign did not terminate")

        self.assertEqual(output.getvalue(), "")
        self.assertEqual(run.result, whole.result)
        self.assertEqual(run.reason, whole.reason)
        self.assertEqual(
            [(e["seed"], e["outcome"]) for e in run.events],
            [(e["seed"], e["outcome"]) for e in whole.events],
        )
        self.assertGreaterEqual(run.work_used, whole.work_used)

    def test_rebuilt_store_handles_every_arithmetic_phase(self):
        config = configuration()
        context = SieveContext(config.max_hi, segment_size=16)
        original_store = ECMPrograms(context, memory_bytes=131072)
        original = new_job("ecm", 1000000000039 * 1000000000061, 17, 50, 1000)
        seen = set()
        for _ in range(10000):
            if original["done"]:
                break
            before = copy.deepcopy(original)
            advance_job(original, allowance(), original_store, config)
            rebuilt = ECMPrograms(context, memory_bytes=131072)
            advance_job(before, allowance(), rebuilt, config)

            self.assertEqual(before, original)
            seen.add(original["phase"])
        self.assertTrue(
            {"stage_one", "baby_steps", "giant_setup", "stage_two"} <= seen
        )

    def test_saturated_product_replay_survives_rebuilt_programs(self):
        config = configuration()
        context = SieveContext(config.max_hi, segment_size=16)
        job = new_job("ecm", 35, 0, 50, 1000)
        job.update(phase="stage_two", terms=[5, 7], product=0)
        job["cursor"] = prime_cursor(1001, 1001)

        advance_job(
            job, allowance(), ECMPrograms(context, memory_bytes=131072), config
        )
        self.assertEqual(job["phase"], "term_replay")
        resumed = copy.deepcopy(job)
        advance_job(
            resumed,
            allowance(),
            ECMPrograms(context, memory_bytes=131072),
            config,
        )
        self.assertEqual(resumed["factor"], 5)
        self.assertTrue(resumed["done"])

    def test_old_schema_and_corrupt_identity(self):
        config = configuration(ecm_program_bytes=0)
        run = factorize_bounded(
            25013 * 25031, config=config, budget=allowance(1000)
        )
        self.assertEqual(run.checkpoint["payload"]["version"], 4)
        self.assertNotIn(
            "ecm_program_bytes", run.checkpoint["payload"]["config"]
        )
        resumed = factorize_bounded(
            run.result.original,
            config=config,
            budget=allowance(),
            checkpoint=run.checkpoint,
        )
        self.assertEqual(resumed.result.reconstruct(), run.result.original)
        with self.assertRaises(ValueError):
            factorize_bounded(
                run.result.original,
                config=configuration(),
                budget=allowance(),
                checkpoint=run.checkpoint,
            )

        run = factorize_bounded(
            25013 * 25031, config=configuration(), budget=allowance(1000)
        )
        corrupt = copy.deepcopy(run.checkpoint)
        corrupt["payload"]["schedule"] = "wrong-program-version"
        with self.assertRaises(ValueError):
            factorize_bounded(
                run.result.original,
                config=configuration(),
                budget=allowance(),
                checkpoint=reseal(corrupt),
            )


class CampaignContractTests(unittest.TestCase):
    def test_resealed_bound_and_endpoint_changes_are_rejected(self):
        config = configuration()
        n = 1000000000039 * 1000000000061
        run = factorize_bounded(n, config=config, budget=allowance(3000))
        job = run.checkpoint["payload"]["state"]["current"]["job"]
        self.assertIsNotNone(job)

        for key, value in (("b1", 51), ("b2", 1001)):
            corrupt = copy.deepcopy(run.checkpoint)
            corrupt["payload"]["state"]["current"]["job"][key] = value
            with self.assertRaisesRegex(ValueError, "bounds disagree"):
                factorize_bounded(
                    n,
                    config=config,
                    budget=allowance(),
                    checkpoint=reseal(corrupt),
                )
        corrupt = copy.deepcopy(run.checkpoint)
        cursor = corrupt["payload"]["state"]["current"]["job"]["cursor"]
        cursor["hi"] += 1
        with self.assertRaisesRegex(ValueError, "wrong endpoint"):
            factorize_bounded(
                n,
                config=config,
                budget=allowance(),
                checkpoint=reseal(corrupt),
            )

    def test_finite_large_tiers_and_predeclared_curve_extension(self):
        config = configuration(
            ecm_tiers=((50000, 5000000, 10000),),
            max_input_bits=329,
            segment_size=1024,
            memory_bytes=16 * 2**20,
            ecm_program_bytes=8 * 2**20,
        )
        self.assertLess(config.workspace_reserve, config.memory_bytes)
        with self.assertRaises(MemoryError):
            replace(config, max_input_bits=4096)

        # More work resumes a declared campaign; changing its identity refuses.
        config = configuration(ecm_tiers=((50, 1000, 2),))
        n = 1000000000039 * 1000000000061
        first = factorize_bounded(n, config=config, budget=allowance(1000))
        resumed = factorize_bounded(
            n, config=config, budget=allowance(), checkpoint=first.checkpoint
        )
        self.assertEqual(resumed.result.reconstruct(), n)
        self.assertEqual(len(resumed.events), 2)
        self.assertEqual(len({event["seed"] for event in resumed.events}), 2)
        with self.assertRaises(ValueError):
            factorize_bounded(
                n,
                config=replace(config, ecm_tiers=((50, 1000, 3),)),
                budget=allowance(),
                checkpoint=first.checkpoint,
            )

    def test_prime_two_shares_the_smallest_resumable_segment(self):
        from v2 import portfolio

        config = configuration(segment_size=1)
        n = 1000000000039 * 1000000000061
        original = portfolio.advance_job

        def pause_first_buffer(job, budget, context, settings):
            original(job, budget, context, settings)
            if job["cursor"] and job["cursor"]["values"] == [2, 3]:
                budget.cancelled = lambda: True

        with patch.object(portfolio, "advance_job", pause_first_buffer):
            first = factorize_bounded(n, config=config, budget=allowance())

        self.assertEqual(first.reason, "cancelled")
        job = first.checkpoint["payload"]["state"]["current"]["job"]
        self.assertEqual(job["cursor"]["values"], [2, 3])
        resumed = factorize_bounded(
            n, config=config, budget=allowance(), checkpoint=first.checkpoint
        )
        self.assertEqual(resumed.result.reconstruct(), n)
        self.assertEqual(resumed.reason, "exhausted")

    def test_versioned_runner_inputs_and_small_schedule_oracle(self):
        from v2.benchmarks.ecm.p52 import p52_a3

        corpus = p52_a3.load_corpus()
        control, _ = p52_a3.load_control()
        self.assertEqual(len(corpus["fixtures"]), 10)
        self.assertNotIn(
            "ecm_program_bytes", control.PortfolioConfig.__dataclass_fields__
        )
        self.assertLess(
            p52_a3.options("small", "regenerated")["ecm_program_bytes"],
            p52_a3.options("small", "programs")["ecm_program_bytes"],
        )
        oracle = p52_a3.reference_schedule(8, 31)
        self.assertEqual(oracle["counts"], [4, 7])
