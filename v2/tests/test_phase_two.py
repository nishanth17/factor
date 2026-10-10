"""Independent oracles and interruption tests for the bounded portfolio."""

import copy
import hashlib
import json
import random
import unittest
from dataclasses import replace
from math import gcd, isqrt, prod
from unittest.mock import patch

from v2 import ecm, utils
from v2.budget import Budget, BudgetExhaustedError
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.preprocessing import fermat_step, integer_root, strip_twos
from v2.schedules import ScheduleCache, SieveContext, iter_primes, prime_powers
from v2.stage_jobs import advance_job, new_job


def reference_primes(lo, hi):
    """Independent exact trial division for short test intervals."""
    return [n for n in range(max(2, lo), hi) if utils.is_prime_bf(n)]


def small_config(**changes):
    """Small finite candidates make saturation/resume tests inexpensive."""
    config = PortfolioConfig(
        pm1_gap_mode="cached",
        pm1_chunk_size=None,
        trial_bound=5,
        rho_attempts=2,
        rho_evaluations=500,
        pm1_b1=10,
        pm1_b2=200,
        ecm_tiers=((50, 1000, 3),),
        segment_size=16,
        max_input_bits=512,
    )
    return replace(config, **changes)


def finish_job(job, config, work=1_000_000):
    """Drive a candidate under a finite guard; retain its final state."""
    context = SieveContext(config.max_hi, segment_size=config.segment_size)
    budget = Budget(work_limit=work, seconds=None, cpu_seconds=None)
    for _ in range(100_000):
        if job["done"]:
            return job, budget.used
        advance_job(job, budget, context, config)
    raise AssertionError("candidate failed to terminate")


def reseal(checkpoint):
    """Reseal corruption to exercise semantic checkpoint validation."""
    encoded = json.dumps(
        checkpoint["payload"], sort_keys=True, separators=(",", ":")
    )
    checkpoint["sha256"] = hashlib.sha256(encoded.encode()).hexdigest()
    return checkpoint


class BudgetTests(unittest.TestCase):
    """Limits must refuse actions before their state can be mutated."""

    def test_exact_charge_and_no_overshoot(self):
        """An exhausted reservation preserves the already-consumed count."""
        budget = Budget(work_limit=3, seconds=None, cpu_seconds=None)

        budget.consume(3)
        with self.assertRaises(BudgetExhaustedError):
            budget.consume()

        self.assertEqual(budget.used, 3)

    def test_expired_and_cancelled(self):
        """Wall, CPU, and cancellation checks also apply to zero-cost steps."""
        for budget, reason in (
            (Budget(seconds=0), "wall_limit"),
            (Budget(cpu_seconds=0), "cpu_limit"),
            (Budget(cancelled=lambda: True), "cancelled"),
        ):
            with self.assertRaisesRegex(BudgetExhaustedError, reason):
                budget.consume(0)

            self.assertEqual(budget.used, 0)

    def test_invalid_limits(self):
        """NaN, infinity, negative values, and Boolean work are invalid."""
        for value in (-1, float("nan"), float("inf"), True):
            with self.assertRaises(ValueError):
                Budget(seconds=value)
        with self.assertRaises(TypeError):
            Budget(work_limit=True)


class ScheduleTests(unittest.TestCase):
    """Streaming values, private ownership, cap boundaries, and restarts."""

    def test_exhaust_small_intervals_and_segments(self):
        """Both strike implementations agree with trial division."""
        for rolling in (False, True):
            for width in (1, 2, 7, 31):
                context = SieveContext(
                    503, segment_size=width, rolling=rolling
                )

                for lo in range(-2, 120, 3):
                    for span in (0, 1, 2, 15, 60):
                        hi = lo + span

                        self.assertEqual(
                            list(context.primes(lo, hi)),
                            reference_primes(lo, hi),
                            (rolling, width, lo, hi),
                        )

    def test_high_values_and_prime_squares(self):
        """Packed base primes never truncate output above 2**32."""
        context = SieveContext(2**32 + 200, rolling=True, segment_size=11)
        self.assertEqual(
            list(context.primes(2**32 - 50, 2**32 + 100)),
            reference_primes(2**32 - 50, 2**32 + 100),
        )
        self.assertEqual(list(iter_primes(3700, 3722)), [3701, 3709, 3719])

    def test_context_ownership_and_limits(self):
        """Abandoned streams must close; caps fail before allocating."""
        context = SieveContext(1000)
        first = context.primes(2, 100)

        self.assertEqual(next(first), 2)
        with self.assertRaises(RuntimeError):
            next(context.primes(2, 100))
        first.close()

        self.assertEqual(list(context.primes(2, 5)), [2, 3])
        with self.assertRaises(ValueError):
            list(context.primes(2, 1001))
        with self.assertRaises(MemoryError):
            SieveContext(10**12, memory_bytes=8192)

    def test_prime_power_schedule(self):
        """The schedule product equals an independent LCM construction."""
        for bound in (2, 3, 10, 243):
            scalar = prod(power for _, power in prime_powers(bound))
            expected = 1
            for value in range(1, bound + 1):
                expected = expected // gcd(expected, value) * value
            self.assertEqual(scalar, expected)

    def test_cache_caps_keys_and_abandoned_consumers(self):
        """Keep representations distinct and avoid caching partial streams."""
        context = SieveContext(10000, segment_size=11)
        cache = ScheduleCache(context, cache_bytes=4096, max_entries=2)

        self.assertEqual(list(cache.values(2, 30)), reference_primes(2, 30))
        self.assertEqual(list(cache.values(2, 30)), reference_primes(2, 30))
        self.assertEqual(cache.hits, 1)
        gaps = list(cache.values(2, 30, kind="gaps"))
        previous, reconstructed = 2, []
        for gap in gaps:
            previous += gap
            reconstructed.append(previous)

        self.assertEqual(reconstructed, reference_primes(2, 30))
        self.assertEqual(
            list(cache.values(2, 11, kind="powers", bound=10)), [8, 9, 5, 7]
        )
        self.assertLessEqual(cache.used_bytes, cache.cache_bytes)
        stream = cache.values(5000, 10000)
        next(stream)
        stream.close()

        self.assertFalse(context.active)
        self.assertFalse(cache.active)
        self.assertNotIn(
            (5000, 10000, "half-open", "primes", None), cache.entries
        )
        self.assertLessEqual(len(cache.entries), 2)

    def test_experimental_sieve_residues_and_restarts(self):
        """Unpromoted variants still require complete independent sequences."""
        from v2.benchmarks.sieve_candidates import (
            integer_bitset,
            presieved,
            wheel_thirty,
        )

        for function in (integer_bitset, presieved, wheel_thirty):
            for lo in range(-2, 65):
                for span in (0, 1, 2, 30, 103):
                    self.assertEqual(
                        function(lo, lo + span),
                        reference_primes(lo, lo + span),
                    )

            self.assertEqual(function(3700, 3722), [3701, 3709, 3719])


class PreprocessingTests(unittest.TestCase):
    """Integer roots distinguish powers and their immediate neighbors."""

    def test_roots_and_twos(self):
        """Use inequalities rather than the implementation's iteration."""
        generator = random.Random(20261003)

        for bits in (8, 64, 256, 1024):
            for exponent in (2, 3, 5, 7, 11, 31):
                n = generator.getrandbits(bits)
                root = integer_root(n, exponent)

                self.assertLessEqual(root**exponent, n)
                self.assertGreater((root + 1) ** exponent, n)

        for exponent in (0, 1, 400, 1000):
            self.assertEqual(strip_twos(15 << exponent), (15, exponent))

    def test_fermat_bound_and_exact_factor(self):
        """One exact close-factor step succeeds; wrong starts are rejected."""
        n = 1009 * 1013

        self.assertEqual(fermat_step(n, 1011), 1009)
        with self.assertRaises(ValueError):
            fermat_step(n, isqrt(n))


class CandidateTests(unittest.TestCase):
    """Verify chunk exponent action and candidate pause/resume paths."""

    def test_m12_atomic_boundary_state_and_work_are_unchanged(self):
        """Compare every committed action with the verified frozen control."""
        from v2.benchmarks.snapshot_loader import load_stage_jobs

        baseline = load_stage_jobs()

        for chunk in (1, 16):
            config = small_config(chunk_size=chunk)

            for kind, n, b1, b2 in (
                ("rho", 35, 0, 0),
                ("rho", 25013 * 25031, 0, 0),
                ("pm1", 13 * 19, 10, 200),
                ("pm1", 607 * 1019, 10, 200),
                ("ecm", 1009 * 1013, 50, 1000),
                ("ecm", 1000000000039 * 1000000000061, 50, 1000),
            ):
                jobs = [new_job(kind, n, 0, b1, b2) for _ in range(2)]
                contexts = [
                    SieveContext(config.max_hi, segment_size=16)
                    for _ in range(2)
                ]
                budgets = [
                    Budget(seconds=None, cpu_seconds=None) for _ in range(2)
                ]

                while not jobs[0]["done"]:
                    for action, job, budget, context in zip(
                        (advance_job, baseline.advance_job),
                        jobs,
                        budgets,
                        contexts,
                    ):
                        action(job, budget, context, config)

                    self.assertEqual(jobs[0], jobs[1])
                    self.assertEqual(budgets[0].used, budgets[1].used)

                divisor = jobs[0]["factor"]
                if divisor is not None:
                    self.assertTrue(utils.valid_divisor(divisor, n))

    def test_chunk_actions_against_independent_exponent(self):
        """Over a prime field, chunking preserves the LCM exponent action."""
        for kind in ("pm1", "ecm"):
            for length in (1, 3, 16):
                config = small_config(chunk_size=length)
                job = new_job(kind, 1000003, 6, 10, 10)
                context = SieveContext(config.max_hi, segment_size=16)
                budget = Budget(seconds=None, cpu_seconds=None)
                advance_job(job, budget, context, config)
                start = job["value"]
                while job["phase"] == "stage_one" and not job["done"]:
                    advance_job(job, budget, context, config)
                scalar = 2520  # Independent lcm(1,...,10).
                if kind == "pm1":
                    expected = pow(start, scalar, job["n"])

                    self.assertEqual(job["value"], expected)
                else:
                    expected = ecm.scalar_multiply(
                        scalar, *start, job["n"], job["a24"]
                    )

                    self.assertEqual(
                        (job["value"][0] * expected[1]) % job["n"],
                        (expected[0] * job["value"][1]) % job["n"],
                    )

    def test_every_candidate_action_can_resume(self):
        """JSON-roundtripping between actions preserves output and work."""
        config = small_config()

        for kind, n, b1, b2 in (
            ("rho", 25013 * 25031, 0, 0),
            ("pm1", 607 * 1019, 10, 200),
            ("ecm", 1009 * 1013, 50, 1000),
        ):
            original = new_job(kind, n, 0, b1, b2)
            expected, used = finish_job(copy.deepcopy(original), config)
            context = SieveContext(config.max_hi, segment_size=16)
            budget = Budget(seconds=None, cpu_seconds=None)
            job = original
            while not job["done"]:
                advance_job(job, budget, context, config)
                job = json.loads(json.dumps(job))

            self.assertEqual(job, expected)
            self.assertEqual(budget.used, used)
            if job["factor"] is not None:
                self.assertTrue(utils.valid_divisor(job["factor"], n))

    def test_saturated_chunk_replay(self):
        """Mixed smooth factors expose a proper factor during finer replay."""
        config = small_config(chunk_size=16)
        job, _ = finish_job(new_job("pm1", 13 * 19, 0, 10, 200), config)
        self.assertTrue(utils.valid_divisor(job["factor"], 13 * 19))

    def test_ready_chunks_check_before_generating_more_primes(self):
        """A ready factor must not require the budget for another segment."""
        config = small_config(chunk_size=1, gcd_batch=1)
        first = new_job("pm1", 35, 0, 10, 200)
        first.update(phase="stage_one", value=2, powers=[[2, 8]])
        budget = Budget(work_limit=5, seconds=None, cpu_seconds=None)
        advance_job(first, budget, None, config)

        self.assertEqual(first["factor"], 5)
        second = new_job("pm1", 35, 0, 10, 200)
        second.update(phase="stage_two", terms=[5], product=5)
        budget = Budget(work_limit=1, seconds=None, cpu_seconds=None)
        advance_job(second, budget, None, config)

        self.assertEqual(second["factor"], 5)

    def test_ecm_early_middle_tail_saturation_is_bounded(self):
        """Replay saturation at early, middle, and tail chunks."""
        from v2 import stage_jobs

        original = stage_jobs._apply

        for position in (1, 2, 4):
            calls = []

            def injected(job, value, scalar):
                calls.append(scalar)
                if len(calls) == position:
                    return [0, 0]
                return original(job, value, scalar)

            with patch.object(stage_jobs, "_apply", injected):
                job, work = finish_job(
                    new_job("ecm", 1000003, 6, 10, 200),
                    small_config(chunk_size=1),
                )

            self.assertTrue(job["done"])
            self.assertIsNone(job["factor"])
            self.assertLess(work, 10000)


class PortfolioTests(unittest.TestCase):
    """Reconstruction, arbitrary-size powers, shared budgets, and metadata."""

    def test_complete_small_sweep(self):
        """Verify terminal factors with an independent trial oracle."""
        config = small_config(trial_bound=100)

        for n in range(1, 501):
            run = factorize_bounded(n, config=config)

            self.assertTrue(run.result.proven, n)
            self.assertEqual(run.result.reconstruct(), n)
            self.assertTrue(
                all(
                    utils.is_prime_bf(factor.value)
                    for factor in run.result.factors
                )
            )

    def test_powers_signs_and_certainty(self):
        """Higher powers preserve multiplicities and probable-prime labels."""
        config = small_config(rho_attempts=0, pm1_attempts=0, ecm_tiers=())

        for n, expected in (
            (-(2**100) * 1009**7, [(2, 100), (1009, 7)]),
            ((2**127 - 1) ** 3, [(2**127 - 1, 3)]),
        ):
            run = factorize_bounded(n, config=config)

            self.assertEqual(run.result.reconstruct(), n)
            self.assertTrue(run.result.complete)
            self.assertEqual(
                [
                    (factor.value, factor.exponent)
                    for factor in run.result.factors
                ],
                expected,
            )
            self.assertEqual(run.result.proven, n < 0)

    def test_shared_budget_and_resume_identity(self):
        """Small resumed grants consume the same work as uninterrupted work."""
        config = small_config()
        n = -(2**9) * 1009 * 1013 * 1019

        full = factorize_bounded(n, seed=91, config=config)
        checkpoint = None

        for allowance in range(0, full.work_used + 1000, 97):
            run = factorize_bounded(
                n,
                seed=91,
                config=config,
                checkpoint=checkpoint,
                budget=Budget(
                    work_limit=allowance, seconds=None, cpu_seconds=None
                ),
            )

            self.assertLessEqual(run.work_used, allowance)
            self.assertEqual(run.result.reconstruct(), n)
            checkpoint = json.loads(json.dumps(run.checkpoint))
            if run.reason == "complete":
                break

        self.assertEqual(run.result, full.result)
        self.assertEqual(run.work_used, full.work_used)
        self.assertEqual(run.events, full.events)

    def test_expiry_cancellation_and_corruption(self):
        """Expired work preserves inputs and incompatible snapshots reject."""
        config = small_config()
        n = 1009 * 1013

        run = factorize_bounded(n, config=config, budget=Budget(seconds=0))

        self.assertEqual(run.reason, "wall_limit")
        self.assertEqual(run.result.remaining, (n,))

        cancelled = factorize_bounded(
            n, config=config, budget=Budget(cancelled=lambda: True)
        )

        self.assertEqual(cancelled.reason, "cancelled")
        checkpoint = copy.deepcopy(run.checkpoint)
        checkpoint["payload"]["state"]["original"] += 1
        with self.assertRaisesRegex(ValueError, "checksum"):
            factorize_bounded(n, config=config, checkpoint=checkpoint)
        with self.assertRaisesRegex(ValueError, "incompatible"):
            factorize_bounded(
                n,
                config=replace(config, chunk_size=3),
                checkpoint=run.checkpoint,
            )

    def test_explicit_exhaustion_and_trace_cap(self):
        """Keep failures and dropped diagnostics explicit and bounded."""
        config = small_config(
            rho_evaluations=1, pm1_attempts=0, ecm_tiers=(), trace_limit=1
        )

        run = factorize_bounded(1009 * 1013, config=config)

        self.assertEqual(run.reason, "exhausted")
        self.assertEqual(run.result.remaining, (1009 * 1013,))
        self.assertEqual(len(run.events), 1)
        self.assertEqual(run.dropped_events, 1)

    def test_factor_one_resume_preserves_full_assignment(self):
        """Preserve children, RNG, and work when stopping at a split."""
        config = small_config()
        n = 1009 * 1013 * 1019

        full = factorize_bounded(n, seed=81, config=config)

        split = factorize_bounded(
            n, seed=81, config=config, stop_after_split=True
        )

        self.assertEqual(split.reason, "factor_found")
        self.assertFalse(split.result.complete)

        resumed = factorize_bounded(
            n, config=config, checkpoint=split.checkpoint
        )

        self.assertEqual(resumed.result, full.result)
        self.assertEqual(resumed.work_used, full.work_used)
        self.assertEqual(resumed.events, full.events)

    def test_semantic_checkpoint_corruption_is_rejected(self):
        """Validate resealed cached claims, witnesses, and prime buffers."""
        config = small_config(trial_bound=100)
        n = 1009 * 1013

        run = factorize_bounded(n, config=config, budget=Budget(work_limit=0))
        checkpoint = copy.deepcopy(run.checkpoint)
        checkpoint["payload"]["state"]["classifications"][str(n)] = (
            utils.Primality.PROVEN.value
        )
        with self.assertRaisesRegex(ValueError, "cached primality"):
            factorize_bounded(n, config=config, checkpoint=reseal(checkpoint))
        run = factorize_bounded(
            n, config=config, budget=Budget(work_limit=200)
        )
        checkpoint = copy.deepcopy(run.checkpoint)
        cursor = checkpoint["payload"]["state"]["current"]["cursor"]

        self.assertTrue(cursor["values"])
        cursor["values"][0] = 9
        with self.assertRaisesRegex(ValueError, "prime values"):
            factorize_bounded(n, config=config, checkpoint=reseal(checkpoint))
        prime = 2147483647

        run = factorize_bounded(
            prime, config=config, budget=Budget(work_limit=100)
        )
        checkpoint = copy.deepcopy(run.checkpoint)
        witness = checkpoint["payload"]["state"]["current"]["prime_job"]
        witness["bases"] = [2]
        with self.assertRaisesRegex(ValueError, "primality progress"):
            factorize_bounded(
                prime, config=config, checkpoint=reseal(checkpoint)
            )


if __name__ == "__main__":
    unittest.main()
