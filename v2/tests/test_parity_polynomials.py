"""Independent CRT, integer norm and valuation oracles for parity SIQS."""

import itertools
import json
import unittest
from dataclasses import asdict, replace
from math import prod
from unittest.mock import patch

from v2.benchmarks.qs.c10.oracles import DiscriminantPolynomial
from v2.execution.budget import BudgetExhaustedError
from v2.portfolio import _decode_siqs_config
from v2.qs import Polynomial, SIQSJob, build_factor_base
from v2.qs.coordinate_family import CoordinateFamily
from v2.qs.families import PolynomialFamily
from v2.qs.multiplier import PARITY_MULTIPLIERS, select_multiplier
from v2.qs.power_sieve import prime_power_roots
from v2.qs.sieve_collector import SieveCollector, SieveConfig
from v2.qs.siqs import make_family
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_qs import reference_positions, unlimited_budget
from v2.tests.test_qs_sieve import signature
from v2.tests.test_siqs import allowance, configuration, mutate


def parity_configuration(**options):
    values = dict(polynomial_policy="parity2", sign_policy="swap_high")
    values.update(options)
    return configuration(**values)


class ParityPolynomialTests(unittest.TestCase):
    def test_cartesian_crt_norm_mapping_roots_and_high_valuations(self):
        for multiplier in (3, 7):
            base = build_factor_base(
                10403, multiplier=multiplier, bound=97
            ).factor_base
            primes = tuple(
                p for p in base.primes if p > 2 and base.n_prime % p
            )[:3]
            a = prod(primes)
            choices = [
                tuple(r for r in range(p) if r * r % p == base.n_prime % p)
                for p in primes
            ]
            expected = set()
            for residues in itertools.product(*choices):
                b = (
                    sum(
                        r * (a // p) * pow(a // p, -1, p)
                        for r, p in zip(residues, primes)
                    )
                    % a
                )
                if b % 2 == 0:
                    b += a
                expected.add(min(b, -b % (2 * a)))

            for policy in ("canonical", "swap_high"):
                family = make_family(
                    base,
                    primes,
                    config=parity_configuration(sign_policy=policy),
                    budget=allowance(),
                )
                actual = set()
                for _ in range(family.count):
                    step = family.next()
                    polynomial = step.polynomial
                    b = polynomial.b
                    self.assertEqual(polynomial.a, 2 * a)
                    self.assertEqual(b % 2, 1)
                    self.assertTrue(-a <= b < a)
                    actual.add(min(b % (2 * a), -b % (2 * a)))
                    oracle = DiscriminantPolynomial(base.n_prime, a, b)
                    b0 = (b + a // 2) % a - a // 2
                    delta = (b - b0) // a
                    standard = Polynomial(10403, multiplier, a, b0)

                    for x in range(-19, 22):
                        value = oracle.value(x)
                        self.assertEqual(polynomial.value(x), 2 * value)
                        self.assertEqual(
                            standard.value(2 * x + delta), 4 * value
                        )
                        self.assertEqual(
                            polynomial.u_value(x),
                            standard.u_value(2 * x + delta),
                        )
                        self.assertEqual(
                            polynomial.u_value(x) ** 2 - base.n_prime,
                            4 * a * value,
                        )

                    for roots in step.roots:
                        p = roots.prime
                        brute = tuple(
                            x for x in range(p) if polynomial.value(x) % p == 0
                        )
                        represented = (
                            tuple(range(p))
                            if roots.all_positions
                            else roots.roots
                        )
                        self.assertEqual(represented, brute)
                        maximum = 256 if p == 2 else p**2
                        for modulus, residues, weight in prime_power_roots(
                            polynomial, roots, maximum, 1, allowance()
                        ):
                            expected_roots = tuple(
                                x
                                for x in range(modulus)
                                if 2 * oracle.value(x) % modulus == 0
                            )
                            self.assertEqual(
                                tuple(sorted(residues)), expected_roots
                            )
                            self.assertEqual(weight, 1)
                self.assertIsNone(family.next())
                self.assertEqual(actual, expected)

    def test_exceptional_multiplier_fallback_and_coupled_scores(self):
        for n in (101, 103, 107, 109):
            choice = select_multiplier(
                n,
                candidates=PARITY_MULTIPLIERS,
                policy="parity2",
                budget=allowance(),
            )
            self.assertEqual(
                choice.multiplier,
                max(choice.scores, key=lambda item: (item[1], -item[0]))[0],
            )
            for k in (1, 2, 3, 5, 7):
                base = build_factor_base(n, multiplier=k, bound=97).factor_base
                primes = tuple(
                    p for p in base.primes if p > 2 and base.n_prime % p
                )[:2]
                family = make_family(
                    base,
                    primes,
                    config=parity_configuration(),
                    budget=allowance(),
                )
                self.assertEqual(family.parity, n * k % 4 == 1)
                if not family.parity:
                    literal = PolynomialFamily(
                        base, primes, budget=allowance()
                    )
                    self.assertEqual(family.next(), literal.next())

        choice = select_multiplier(15, candidates=(2, 6), policy="parity2")
        self.assertEqual(choice.divisor, 3)
        for candidates in ((2,), (4,), (9,), (2, 2)):
            with self.assertRaises(ValueError):
                select_multiplier(10403, candidates=candidates)
        with self.assertRaises(ValueError):
            select_multiplier(10403, candidates=(4,), policy="parity2")
        with self.assertRaises(ValueError):
            configuration(multiplier_policy="parity2")
        with self.assertRaises(ValueError):
            parity_configuration(mode="mpqs")

    def test_exact_division_queue_presieve_plan_and_refusal(self):
        base = build_factor_base(10403, multiplier=3, bound=97).factor_base
        family = make_family(
            base,
            tuple((p for p in base.primes if p > 2 and base.n_prime % p))[:3],
            config=parity_configuration(),
            budget=allowance(),
        )
        polynomial = family.next().polynomial
        expected = reference_positions(polynomial, base, -71, 80, 97)
        for marking, queue, plan in (
            ("sparse", 0, 0),
            ("sparse", 1, 1),
            ("sparse", 2**20, 2**20),
            ("density+presieve+scan_fused", 2**20, 2**20),
        ):
            cfg = SieveConfig(
                score_policy="powers",
                marking=marking,
                division="bucket",
                block_width=17,
                residual_bound=97,
                sparse_queue_bytes=queue,
                power_plan_bytes=plan,
                max_partials=4096,
                preparation_cache_bytes=min(2 * 2**20, 2048 * 4096),
                target_rows=4096,
            )
            collector = SieveCollector(
                polynomial, base, config=cfg, budget=unlimited_budget()
            )
            result = collector.collect(-71, 80)
            self.assertEqual(signature(result), expected)
            self.assertLessEqual(collector._workspace, cfg.memory_limit)

        collector = SieveCollector(
            polynomial, base, config=cfg, budget=allowance()
        )
        before = collector._workspace
        with patch(
            "v2.qs.sparse_queue.build",
            side_effect=BudgetExhaustedError("work"),
        ):
            result = collector.collect(-71, 80)
        self.assertEqual(result.next_position, -71)
        self.assertEqual(collector._workspace, before)
        self.assertEqual(signature(collector.collect(-71, 80)), expected)

    def test_family_checkpoint_and_atomic_refusal(self):
        base = build_factor_base(10403, multiplier=3, bound=97).factor_base
        for policy in ("canonical", "swap_high"):
            cfg = parity_configuration(sign_policy=policy)
            family = make_family(
                base,
                tuple((p for p in base.primes if p > 2 and base.n_prime % p))[
                    :3
                ],
                config=cfg,
                budget=allowance(),
            )
            family.next()
            saved = family.checkpoint()
            resources = saved["payload"]["resources"]
            budget = Budget(
                work_limit=10**9,
                used=resources["work_used"],
                prior_wall=resources["wall_used"],
                prior_cpu=resources["cpu_used"],
                seconds=None,
                cpu_seconds=None,
            )
            kind = (
                CoordinateFamily if policy == "swap_high" else PolynomialFamily
            )
            restored = kind.from_checkpoint(
                base, saved, budget=budget, memory_bytes=2**24
            )
            self.assertEqual(restored.next(), family.next())
            previous = family.current
            family.budget = allowance(0)
            with self.assertRaises(BudgetExhaustedError):
                family.next()
            self.assertEqual(family.current, previous)
            self.assertEqual(family.next_index, 2)

    def test_two_adic_density_and_singular_roots_independently(self):
        # Enumerate U residues rather than replaying the production scorer.
        modulus = 2**12
        for d in (17, 21, 19, 23, 26, 30):
            standard = (
                sum(
                    sum((u * u - d) % (2**j) == 0 for j in range(1, 13))
                    for u in range(modulus)
                )
                / modulus
            )
            expected = {1: 2, 5: 1, 3: 0.5, 7: 0.5, 2: 0.5, 6: 0.5}[d % 8]
            self.assertLessEqual(abs(standard - expected), 4 / modulus)
            if d % 4 == 1:
                parity = sum(
                    sum((u * u - d) % (2**j) == 0 for j in range(1, 13)) - 1
                    for u in range(1, modulus, 2)
                ) / (modulus // 2)
                self.assertLessEqual(
                    abs(parity - (3 if d % 8 == 1 else 1)), 8 / modulus
                )
        # A multiplier prime has one root, but cannot lift to p squared.
        for p in (3, 5, 7):
            d = p * 101
            self.assertEqual(sum((u * u - d) % p == 0 for u in range(p)), 1)
            self.assertFalse(
                any((u * u - d) % (p * p) == 0 for u in range(p * p))
            )

    def test_cancelled_pending_activation_replays_without_skipping(self):
        cfg = parity_configuration(multiplier=3)
        job = SIQSJob(4001 * 5003, config=cfg, budget=allowance())
        job._setup()
        original = job._next_step
        calls = []

        def cancel_after_step():
            step = original()
            calls.append(step)
            job.budget.cancelled = lambda: True
            return step

        with patch.object(job, "_next_step", side_effect=cancel_after_step):
            result = job.run(max_blocks=1)
        self.assertEqual(result.reason, "cancelled")
        self.assertEqual(len(calls), 1)
        self.assertIsNotNone(job.pending_step)
        prior = job.budget.used
        job.budget.cancelled = None
        saved = job.checkpoint()
        restored = SIQSJob.from_checkpoint(saved, budget=allowance())
        self.assertGreaterEqual(restored.budget.used, prior)
        expected = job.run()
        actual = restored.run()
        self.assertEqual(
            (actual.divisor, actual.cofactor),
            (expected.divisor, expected.cofactor),
        )
        self.assertEqual(actual.divisor * actual.cofactor, job.n)

    def test_job_extraction_charged_resume_and_policy_identity(self):
        for multiplier in (1, 2, 3, 7):
            cfg = parity_configuration(multiplier=multiplier)
            job = SIQSJob(4001 * 5003, config=cfg, seed=73, budget=allowance())
            job.run(max_blocks=1)
            saved = json.loads(json.dumps(job.checkpoint()))
            self.assertEqual(saved["version"], 8)
            restored = SIQSJob.from_checkpoint(saved, budget=allowance())
            self.assertEqual(restored.config, cfg)
            self.assertGreaterEqual(restored.budget.used, job.budget.used)
            left, right = job.run(), restored.run()
            self.assertEqual(
                (left.divisor, left.cofactor), (right.divisor, right.cofactor)
            )
            self.assertEqual(left.reason, "factor_found")
            self.assertEqual(left.divisor * left.cofactor, job.n)
            self.assertEqual(_decode_siqs_config(asdict(cfg)), cfg)
            with self.assertRaises(ValueError):
                SIQSJob.from_checkpoint(
                    saved,
                    budget=allowance(),
                    config=replace(cfg, polynomial_policy="standard"),
                    allow_extension=True,
                )

        legacy = SIQSJob(
            4001 * 5003, config=configuration(), budget=allowance()
        ).checkpoint()
        self.assertLess(legacy["version"], 8)
        self.assertNotIn(
            "polynomial_policy", json.loads(legacy["blob"])["config"]
        )
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                mutate(
                    legacy,
                    lambda p: p["config"].update(polynomial_policy="parity2"),
                ),
                budget=allowance(),
            )


if __name__ == "__main__":
    unittest.main()
