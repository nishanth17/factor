"""Independent p-1 default, accounting and old-checkpoint oracles."""

import copy
import hashlib
import json
import unittest
from dataclasses import replace
from math import gcd

from v2 import portfolio, stage_jobs
from v2.benchmarks.a6_production import load_control
from v2.budget import Budget, BudgetExhaustedError
from v2.pm1_gaps import verify_powers
from v2.schedules import SieveContext


def allowance(work=1000000):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


def configuration(**changes):
    return portfolio.PortfolioConfig(
        trial_bound=5,
        rho_attempts=0,
        ecm_tiers=(),
        pm1_b1=7,
        pm1_b2=71,
        segment_size=16,
        gcd_batch=4,
        **changes,
    )


class ProductionPM1Tests(unittest.TestCase):
    def test_all_stage_two_relations_against_direct_powers(self):
        # Start from an arbitrary unit, independently of stage one. Exercise
        # small even gaps, the initial odd exponent and an oversized fallback.
        n, value = 1009 * 1013, 7
        primes = [11, 13, 17, 19, 23, 29, 31, 167, 173]
        for mode in ("cached", "recurrence"):
            config = configuration(pm1_gap_mode=mode)
            job = stage_jobs.new_job("pm1", n, 0, 7, 173)
            job.update(
                phase="stage_two",
                value=value,
                previous_prime=0,
                stage_two_value=1,
                cursor={
                    "left": 8,
                    "next": 174,
                    "hi": 174,
                    "values": primes,
                    "index": 0,
                },
            )
            context, ledger = SieveContext(174, segment_size=100), allowance()
            config = replace(config, gcd_batch=64)
            stage_jobs.advance_job(job, ledger, context, config)
            self.assertEqual(
                job["terms"], [(pow(value, q, n) - 1) % n for q in primes]
            )
            self.assertEqual(
                job["product"], __import__("math").prod(job["terms"]) % n
            )
            growth = len(job.get("even_powers", []))
            gaps = [primes[0]] + [b - a for a, b in zip(primes, primes[1:])]
            self.assertEqual(
                ledger.used, sum(g.bit_length() + 1 for g in gaps) + growth
            )
            verify_powers(job, ledger, mode)
            for k, power in enumerate(job.get("even_powers", []), 1):
                self.assertEqual(power, pow(value, 2 * k, n))

    def test_atomic_growth_refusal_and_cancellation(self):
        config = configuration()
        job = stage_jobs.new_job("pm1", 1009 * 1013, 0, 7, 71)
        job.update(
            phase="stage_two",
            value=7,
            previous_prime=11,
            stage_two_value=pow(7, 11, job["n"]),
            cursor={
                "left": 12,
                "next": 72,
                "hi": 72,
                "values": [13, 17, 19, 23],
                "index": 0,
            },
        )
        original = copy.deepcopy(job)
        for ledger in (allowance(0), Budget(cancelled=lambda: True)):
            with self.assertRaises(BudgetExhaustedError):
                stage_jobs.advance_job(job, ledger, SieveContext(72), config)
            self.assertEqual(job, original)

    def test_finite_saturation_nonunits_and_bound_cases(self):
        for n in (3 * 5, 7 * 11, 17 * 23, 29 * 47, 31 * 43, 101 * 107):
            for mode in ("cached", "recurrence"):
                for chunk in (1, 16, 64):
                    config = configuration(
                        pm1_gap_mode=mode, pm1_chunk_size=chunk
                    )
                    job = stage_jobs.new_job("pm1", n, 0, 7, 71)
                    context, ledger = (
                        SieveContext(72, segment_size=16),
                        allowance(),
                    )
                    for _ in range(1000):
                        if job["done"]:
                            break
                        stage_jobs.advance_job(job, ledger, context, config)
                    self.assertTrue(job["done"])
                    if job["factor"] is not None:
                        self.assertEqual(n % job["factor"], 0)
                        self.assertTrue(1 < job["factor"] < n)
        self.assertEqual(gcd(pow(2, 420, 23) - 1, 23), 1)
        self.assertEqual(gcd(pow(pow(2, 420, 23), 11, 23) - 1, 23), 23)

    def test_old_execution_matches_frozen_jobs_and_snapshot(self):
        old = load_control("portfolio")
        options = dict(
            trial_bound=5,
            rho_attempts=0,
            ecm_tiers=(),
            pm1_b1=7,
            pm1_b2=71,
            segment_size=16,
            gcd_batch=4,
        )
        oldconfig = old.PortfolioConfig(**options)
        legacy = portfolio.PortfolioConfig(
            **options, pm1_gap_mode="cached", pm1_chunk_size=None
        )
        n = 1009 * 1013
        oldrun = old.factorize_bounded(
            n, seed=7, config=oldconfig, budget=allowance(700)
        )
        newrun = portfolio.factorize_bounded(
            n, seed=7, config=legacy, budget=allowance(700)
        )
        for key in ("version", "config", "rng", "work_used", "schedule"):
            self.assertEqual(
                oldrun.checkpoint["payload"][key],
                newrun.checkpoint["payload"][key],
            )
        resumed = portfolio.factorize_bounded(
            n,
            config=legacy,
            budget=allowance(),
            checkpoint=json.loads(json.dumps(oldrun.checkpoint)),
        )
        whole = old.factorize_bounded(
            n, seed=7, config=oldconfig, budget=allowance()
        )
        self.assertEqual(
            (resumed.result, resumed.work_used),
            (whole.result, whole.work_used),
        )
        with self.assertRaisesRegex(ValueError, "incompatible"):
            portfolio.factorize_bounded(
                n,
                config=configuration(),
                budget=allowance(),
                checkpoint=oldrun.checkpoint,
            )

    def test_new_schema_resume_and_table_corruption(self):
        n, config = 1009 * 1013, configuration()
        whole = portfolio.factorize_bounded(
            n, seed=7, config=config, budget=allowance()
        )
        checkpoint = None
        for work in range(0, whole.work_used, 1):
            paused = portfolio.factorize_bounded(
                n, seed=7, config=config, budget=allowance(work)
            )
            current = paused.checkpoint["payload"]["state"]["current"]
            if current and (current.get("job") or {}).get("even_powers"):
                checkpoint = paused.checkpoint
                break
        self.assertIsNotNone(checkpoint)
        self.assertEqual(checkpoint["payload"]["version"], 9)
        resumed = portfolio.factorize_bounded(
            n,
            config=config,
            budget=allowance(),
            checkpoint=json.loads(json.dumps(checkpoint)),
        )
        self.assertEqual(resumed.result, whole.result)
        self.assertGreater(resumed.work_used, whole.work_used)
        same = portfolio.factorize_bounded(
            n,
            config=config,
            budget=allowance(checkpoint["payload"]["work_used"]),
            checkpoint=checkpoint,
        )
        self.assertEqual(same.reason, "work_limit")
        damaged = copy.deepcopy(checkpoint)
        damaged["payload"]["state"]["current"]["job"]["even_powers"][0] += 1
        encoded = portfolio._canonical(damaged["payload"])
        damaged["sha256"] = hashlib.sha256(encoded.encode()).hexdigest()
        with self.assertRaisesRegex(ValueError, "corrupt"):
            portfolio.factorize_bounded(
                n, config=config, budget=allowance(), checkpoint=damaged
            )

    def test_implicit_defaults_and_legacy_snapshot_inference(self):
        config = portfolio.PortfolioConfig()
        self.assertEqual(
            (config.pm1_chunk_size, config.pm1_gap_mode), (64, "recurrence")
        )
        self.assertEqual(config.chunk_size, 16)
        n = 1009 * 1013
        fresh = portfolio.factorize_bounded(n, budget=allowance(0))
        self.assertEqual(fresh.checkpoint["payload"]["version"], 9)
        resumed = portfolio.factorize_bounded(
            n, budget=allowance(0), checkpoint=fresh.checkpoint
        )
        self.assertEqual(resumed.reason, "work_limit")

        old = load_control("portfolio")
        paused = old.factorize_bounded(n, budget=allowance(0))
        resumed = portfolio.factorize_bounded(
            n, budget=allowance(0), checkpoint=paused.checkpoint
        )
        self.assertEqual(resumed.checkpoint["payload"]["version"], 4)
        self.assertNotIn(
            "pm1_gap_mode", resumed.checkpoint["payload"]["config"]
        )
        self.assertEqual(resumed.work_used, paused.work_used)

    def test_implicit_bounded_campaign_default_and_legacy_resume(self):
        from v2.pm1_bounded import factorize_pm1_bounded
        from v2.tests.test_a6_followup import frozen_control

        n = 1009 * 1013
        fresh = factorize_pm1_bounded(n, budget=allowance(), max_actions=0)
        saved = fresh.checkpoint["payload"]
        self.assertEqual(saved["execution"], "pm1-tuning-v1")
        self.assertEqual(
            (saved["config"]["chunk_size"], saved["config"]["gap_mode"]),
            (64, "recurrence"),
        )
        resumed = factorize_pm1_bounded(
            n, budget=allowance(), checkpoint=fresh.checkpoint, max_actions=0
        )
        self.assertEqual(resumed.result, fresh.result)

        old = frozen_control()
        paused = old.factorize_pm1_bounded(
            n, budget=allowance(), max_actions=0
        )
        resumed = factorize_pm1_bounded(
            n, budget=allowance(), checkpoint=paused.checkpoint, max_actions=0
        )
        self.assertEqual(
            resumed.checkpoint["payload"]["execution"], "pm1-campaign-v1"
        )
        self.assertEqual(
            resumed.checkpoint["payload"]["config"]["chunk_size"], 16
        )

    def test_malformed_implicit_resume_is_rejected(self):
        from v2.pm1_bounded import factorize_pm1_bounded

        for checkpoint in ([1], {"payload": []}, {"payload": {}}):
            for function in (
                portfolio.factorize_bounded,
                factorize_pm1_bounded,
            ):
                with self.assertRaises(ValueError):
                    function(
                        1009 * 1013, budget=allowance(), checkpoint=checkpoint
                    )
