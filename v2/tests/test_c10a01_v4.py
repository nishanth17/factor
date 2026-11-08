"""Independent population, CRT, root-cache and v4 selection/resume checks."""

import itertools
import json
import unittest
from dataclasses import asdict, replace
from fractions import Fraction
from math import prod
from unittest.mock import patch

from v2.benchmarks.qs.c10a01.configuration_v3 import configuration as v3_config
from v2.benchmarks.qs.c10a01.configuration_v4 import configuration as v4_config
from v2.execution.budget import BudgetExhaustedError
from v2.qs import Polynomial, SIQSJob, build_factor_base
from v2.qs.assignment_stream import AssignmentStream
from v2.qs.coordinate_family import CoordinateFamily
from v2.qs.families import PolynomialFamily, _checksum
from v2.qs.multiplier import select_coupled_multiplier
from v2.qs.norm_multiplier import (
    NormPolicyRefusedError,
    preview_polynomials,
    preview_positions,
    select_norm_multiplier,
)
from v2.qs.product_window import feasible_core_window
from v2.qs.root_updates import root_update_workspace
from v2.qs.sign_schedule import SignSchedule, schedule_for
from v2.qs.siqs import make_family
from v2.tests.test_c10a01_coupled import coupled_config
from v2.tests.test_qs_families import resumed_budget
from v2.tests.test_siqs import allowance, mutate


def config_v4(**options):
    values = dict(
        root_update_policy="batch",
        assignment_policy="windowed_pair",
        multiplier_policy="parity2_norm",
    )
    values.update(options)
    return coupled_config(**values)


def direct_polynomial(base, primes, mask, parity):
    """Successive CRT congruence merging, independent of term/Gray updates."""
    b, modulus = 0, 1
    for index, prime in enumerate(primes):
        root = next(
            e.square_roots[0] for e in base.entries if e.prime == prime
        )
        residue = -root if mask & (1 << index) else root
        b += modulus * ((residue - b) * pow(modulus, -1, prime) % prime)
        modulus *= prime
    if parity:
        b += modulus if b % 2 == 0 else 0
        modulus *= 2
    b = (b + modulus // 2) % modulus - modulus // 2
    return Polynomial(base.n, base.multiplier, modulus, b)


class GrayTableTests(unittest.TestCase):
    def test_direct_crt_brute_roots_signed_coordinates_and_even_fallback(self):
        saw_wrap = False
        for k in (1, 2, 3, 6, 7, 10, 15, 21):
            base = build_factor_base(10403, multiplier=k, bound=97).factor_base
            primes = tuple(
                p for p in base.primes if p != 2 and base.n_prime % p
            )[:4]
            parity = base.n_prime % 4 == 1
            schedules = (
                schedule_for(4, "control"),
                schedule_for(4, "swap_high"),
                SignSchedule(4, (3, 1, 2), 10),
            )
            for schedule in schedules:
                options = dict(
                    schedule=schedule, parity=parity, memory_bytes=8 * 2**20
                )
                families = [
                    CoordinateFamily(
                        base,
                        primes,
                        budget=allowance(),
                        root_update_policy=policy,
                        **options,
                    )
                    for policy in ("standard", "cached", "batch")
                ]
                self.assertEqual(
                    families[1].root_updates, families[2].root_updates
                )
                for index in range(8):
                    steps = [family.next() for family in families]
                    expected = direct_polynomial(
                        base, primes, schedule.mask(index), parity
                    )
                    for step in steps:
                        self.assertEqual(step.polynomial, expected)
                        self.assertEqual(step.roots, steps[0].roots)
                    for entry, roots in zip(base.entries, steps[2].roots):
                        actual = (
                            tuple(range(entry.prime))
                            if roots.all_positions
                            else roots.roots
                        )
                        self.assertEqual(
                            actual,
                            tuple(
                                x
                                for x in range(entry.prime)
                                if expected.value(x) % entry.prime == 0
                            ),
                        )
                    for x in (-65, -1, 0, 1, 65):
                        self.assertEqual(
                            expected.u_value(x) ** 2 - base.n_prime,
                            expected.a * expected.value(x),
                        )
                    if index:
                        old_gray = (index - 1) ^ ((index - 1) >> 1)
                        bit = (
                            (index ^ (index >> 1)) ^ old_gray
                        ).bit_length() - 1
                        sign = -1 if old_gray & (1 << bit) else 1
                        saw_wrap |= (
                            steps[2].delta_b
                            != -2 * sign * families[2].terms[bit + 1]
                        )
                self.assertTrue(all(f.next() is None for f in families))
        self.assertTrue(saw_wrap)

    def test_quota_memory_cancel_and_both_checkpoint_types(self):
        base = build_factor_base(10403, multiplier=3, bound=97).factor_base
        primes = tuple(p for p in base.primes if p != 2 and base.n_prime % p)[
            :4
        ]
        standard = PolynomialFamily(
            base, primes, parity=True, budget=allowance()
        )
        for policy in ("cached", "batch"):
            for quota in (1, 2, 3, 5, 8):
                family = PolynomialFamily(
                    base,
                    primes,
                    parity=True,
                    root_update_policy=policy,
                    root_update_limit=quota,
                    budget=allowance(),
                )
                extra = root_update_workspace(len(base.entries), 4, quota)
                self.assertEqual(
                    family.workspace_bytes, standard.workspace_bytes + extra
                )
                self.assertEqual(
                    len(family.root_updates), (quota - 1).bit_length()
                )
                with self.assertRaises(TypeError):
                    PolynomialFamily(
                        base,
                        primes,
                        parity=True,
                        root_update_policy=policy,
                        root_update_limit=quota,
                        memory_bytes=family.workspace_bytes - 1,
                        budget=allowance(),
                    )
                family.next()
                saved = family.checkpoint()
                restored = PolynomialFamily.from_checkpoint(
                    base, saved, budget=resumed_budget(saved)
                )
                self.assertGreater(
                    restored.budget.used,
                    saved["payload"]["resources"]["work_used"],
                )
                self.assertEqual(restored.root_updates, family.root_updates)
                for _ in range(quota):
                    self.assertEqual(family.next(), restored.next())
            family = CoordinateFamily(
                base,
                primes,
                parity=True,
                schedule=SignSchedule(4, (3, 1, 2), 10),
                budget=allowance(),
                root_update_policy=policy,
                root_update_limit=5,
            )
            family.next()
            family.next()
            saved = family.checkpoint()
            restored = CoordinateFamily.from_checkpoint(
                base,
                saved,
                budget=resumed_budget(saved),
                memory_bytes=8 * 2**20,
            )
            self.assertEqual(family.next(), restored.next())
            old = (family.next_index, family.current)
            family.budget.cancelled = lambda: True
            with self.assertRaises(BudgetExhaustedError):
                family.next()
            self.assertEqual((family.next_index, family.current), old)
            family.budget.cancelled = None
            self.assertEqual(family.next(), restored.next())
            bad = json.loads(json.dumps(saved))
            bad["payload"]["schema"] = "c10a01-coordinate-family-v2"
            bad["sha256"] = _checksum(bad["payload"])
            with self.assertRaises(TypeError):
                CoordinateFamily.from_checkpoint(
                    base, bad, budget=allowance(), memory_bytes=8 * 2**20
                )
        ledger = allowance()
        ledger.cancelled = lambda: ledger.used > standard.budget.used + 10
        with self.assertRaises(BudgetExhaustedError):
            PolynomialFamily(
                base,
                primes,
                parity=True,
                root_update_policy="batch",
                budget=ledger,
            )


class ProductPopulationTests(unittest.TestCase):
    def test_every_small_window_and_cartesian_pair_population(self):
        core = (3, 11, 23, 37, 53, 71)
        lower, upper = (5, 17, 31, 59), (7, 19, 41, 67)
        for count in (1, 2, 3):
            for target in range(10, 200000, 997):

                def key(start):
                    stop = start + 4
                    pool = core[start:stop]
                    products = [
                        prod(c) for c in itertools.combinations(pool, count)
                    ]
                    misses = [
                        max(
                            0,
                            target - c * max(lower) * max(upper),
                            c * min(lower) * min(upper) - target,
                        )
                        for c in products
                    ]
                    return (
                        max(misses),
                        abs(pool[1] ** (count + 2) - target),
                        start,
                    )

                expected = min(range(3), key=key)
                actual, miss = feasible_core_window(
                    core, lower, upper, count, 4, target, allowance()
                )
                stop = expected + 4
                self.assertEqual(actual, core[expected:stop])
                self.assertEqual(miss, key(expected)[0])
        base = build_factor_base(4001 * 5003, bound=300).factor_base
        options = dict(
            factor_count=4,
            pool_size=6,
            seed=73,
            policy="windowed_pair",
            memory_bytes=8 * 2**20,
        )
        stream = AssignmentStream(
            base, 64, family_count=10000, budget=allowance(), **options
        )
        self.assertFalse(set(stream.pool) & set(stream.lower_flyers))
        self.assertFalse(set(stream.pool) & set(stream.upper_flyers))
        self.assertFalse(set(stream.lower_flyers) & set(stream.upper_flyers))
        expected = set()
        for core in itertools.combinations(stream.pool, 2):
            pair = min(
                itertools.product(stream.lower_flyers, stream.upper_flyers),
                key=lambda pair: (
                    abs(prod(core) * prod(pair) - stream.target),
                    *pair,
                ),
            )
            expected.add(tuple(sorted(core + pair)))
        actual = [stream[i] for i in range(len(stream))]
        self.assertEqual(set(actual), expected)
        self.assertEqual(len(actual), len(set(actual)))
        short = AssignmentStream(
            base, 64, family_count=3, budget=allowance(), **options
        )
        self.assertEqual(short.identity, stream.identity)
        self.assertEqual(actual[:3], [short[i] for i in range(3)])
        used = stream.budget.used
        stream.budget.work_limit = used
        with self.assertRaises(BudgetExhaustedError):
            stream[0]
        self.assertEqual(stream.budget.used, used)
        stream.budget.work_limit += 10000000
        self.assertEqual(stream[0], actual[0])
        with self.assertRaises(IndexError):
            stream[len(stream)]


class NormPreviewTests(unittest.TestCase):
    def test_preview_trials_and_independent_rational_ranking(self):
        n = 4001 * 5003
        cfg = config_v4()
        result = select_norm_multiplier(
            n, config=cfg, seed=73, budget=allowance()
        )
        density = dict(select_coupled_multiplier(n, budget=allowance()).scores)
        scores = []
        for receipt in result.previews:
            k = receipt["multiplier"]
            base = build_factor_base(
                n, multiplier=k, bound=cfg.base_bound, budget=allowance()
            ).factor_base
            stream = AssignmentStream(
                base,
                cfg.half_width * (2 if base.n_prime % 4 == 1 else 1),
                factor_count=cfg.factor_count,
                family_count=4,
                pool_size=cfg.pool_size,
                seed=73,
                policy=cfg.assignment_policy,
                budget=allowance(),
            )
            logs, corrections, samples = 0, Fraction(0), 0
            for index in range(4):
                primes = stream[index]
                family = make_family(
                    base, primes, config=cfg, budget=allowance()
                )
                expected = [
                    family.next().polynomial for _ in range(cfg.gray_limit)
                ]
                ordinals = sorted(
                    {(cfg.gray_limit - 1) * j // 3 for j in range(4)}
                )
                self.assertEqual(
                    preview_polynomials(base, primes, cfg, allowance()),
                    tuple(expected[j] for j in ordinals),
                )
                for j in ordinals:
                    poly = direct_polynomial(
                        base,
                        primes,
                        schedule_for(len(primes), "swap_high").mask(j),
                        base.n_prime % 4 == 1,
                    )
                    for x in sorted(
                        {
                            -cfg.half_width + 2 * cfg.half_width * i // 16
                            for i in range(17)
                        }
                    ):
                        value = abs(
                            poly.a * x * x
                            + 2 * poly.b * x
                            + (poly.b * poly.b - n * k) // poly.a
                        )
                        logs += (max(1, value) ** 256).bit_length() - 1
                        samples += 1
                corrections += sum(
                    (
                        Fraction((p**256).bit_length() - 1, p - 1)
                        for p in primes
                        if p < 1988
                    ),
                    Fraction(0),
                )
            score = (
                density[k]
                + Fraction((k**256).bit_length() - 1, 2)
                - Fraction(logs, samples)
                - corrections / 4
            )
            self.assertEqual(
                score,
                Fraction(
                    receipt["score_numerator"], receipt["score_denominator"]
                ),
            )
            scores.append((k, score))
        self.assertEqual(
            result.multiplier, max(scores, key=lambda x: (x[1], -x[0]))[0]
        )
        self.assertEqual(len(result.previews), 4)
        self.assertEqual(preview_positions(1), (-1, 0, 1))
        for k in (2, 3, 7, 10):
            job = SIQSJob(
                n, config=replace(cfg, multiplier=k), budget=allowance()
            )
            job._setup()
            primes = job.assignments[0]
            family = make_family(
                job.base, primes, config=cfg, budget=allowance()
            )
            self.assertEqual(
                preview_polynomials(job.base, primes, cfg, allowance()),
                tuple(family.next().polynomial for _ in range(4)),
            )

    def test_refusal_cancellation_and_exceptional_gcd(self):
        cfg = config_v4()
        with self.assertRaises(TypeError):
            select_norm_multiplier(
                4001 * 5003,
                config=cfg,
                seed=7,
                budget=allowance(),
                memory_bytes=1,
            )
        with patch(
            "v2.qs.norm_multiplier._preview",
            side_effect=NormPolicyRefusedError("fixture"),
        ) as preview:
            with self.assertRaises(NormPolicyRefusedError):
                select_norm_multiplier(
                    4001 * 5003,
                    config=cfg,
                    seed=7,
                    budget=allowance(),
                    memory_bytes=32 * 2**20,
                )
            self.assertEqual(preview.call_count, 4)
        result = select_norm_multiplier(
            101 * 4001, config=cfg, seed=7, budget=allowance()
        )
        self.assertEqual(result.divisor, 101)
        job = SIQSJob(4001 * 5003, config=cfg, budget=allowance())
        job.budget.cancelled = lambda: job.budget.used > 10000000
        self.assertEqual(job.run().reason, "cancelled")
        self.assertIsNone(job.multiplier)
        charged = job.budget.used
        job.budget.cancelled = None
        job._setup()
        self.assertGreater(job.budget.used, charged)


class V4IntegrationTests(unittest.TestCase):
    def test_literal_controls_and_each_opt_in_checkpoint_gate(self):
        for digits in (60, 62, 65, 70):
            self.assertEqual(
                v4_config(digits, "control"), v3_config(digits, "control")
            )
            new, old = (
                v4_config(digits, "coupled"),
                v3_config(digits, "coupled"),
            )
            self.assertEqual(new.memory_bytes, old.memory_bytes)
            expected = asdict(old)
            expected.update(
                root_update_policy="batch",
                assignment_policy="windowed_pair",
                multiplier_policy="parity2_norm",
            )
            self.assertEqual(asdict(new), expected)
        old = SIQSJob(
            4001 * 5003,
            config=coupled_config(multiplier=3),
            budget=allowance(),
        )
        saved = old.checkpoint()
        self.assertEqual(saved["version"], 9)
        self.assertNotIn(
            "root_update_policy", json.loads(saved["blob"])["config"]
        )
        for field, value in (
            ("root_update_policy", "batch"),
            ("assignment_policy", "windowed_pair"),
            ("multiplier_policy", "parity2_norm"),
        ):
            bad = mutate(saved, lambda p: p["config"].update({field: value}))
            with self.assertRaisesRegex(ValueError, "legacy checkpoint"):
                SIQSJob.from_checkpoint(bad, budget=allowance())

    def test_proper_extraction_resume_and_policy_substitution(self):
        for k in (0, 2, 3, 7):
            cfg = config_v4(multiplier=k)
            job = SIQSJob(4001 * 5003, config=cfg, seed=73, budget=allowance())
            job._setup()
            self.assertTrue(job._activate())
            saved = json.loads(json.dumps(job.checkpoint()))
            self.assertEqual(saved["version"], 10)
            restored = SIQSJob.from_checkpoint(
                saved, config=cfg, budget=allowance()
            )
            self.assertGreater(restored.budget.used, job.budget.used)
            left, right = job.run(), restored.run()
            self.assertEqual(
                (left.reason, left.divisor), (right.reason, right.divisor)
            )
            self.assertEqual(left.reason, "factor_found")
            self.assertEqual(left.divisor * left.cofactor, job.n)
            for field, value in (
                ("root_update_policy", "cached"),
                ("assignment_policy", "paired_flyer"),
                ("multiplier_policy", "parity2_wide"),
            ):
                with self.assertRaises(ValueError):
                    SIQSJob.from_checkpoint(
                        saved,
                        config=replace(cfg, **{field: value}),
                        allow_extension=True,
                        budget=allowance(),
                    )


class V4BoundaryTests(unittest.TestCase):
    def test_pending_transition_and_failed_private_table_are_retryable(self):
        cfg = config_v4(multiplier=3)
        job = SIQSJob(4001 * 5003, config=cfg, budget=allowance())
        job._setup()
        original = job._next_step

        def cancel_after_step():
            step = original()
            job.budget.cancelled = lambda: True
            return step

        with patch.object(job, "_next_step", side_effect=cancel_after_step):
            self.assertEqual(job.run().reason, "cancelled")
        self.assertIsNotNone(job.pending_step)
        saved = job.checkpoint()
        restored = SIQSJob.from_checkpoint(saved, budget=allowance())
        job.budget.cancelled = None
        left, right = job.run(), restored.run()
        self.assertEqual(
            (left.reason, left.divisor), (right.reason, right.divisor)
        )
        self.assertEqual(left.divisor * left.cofactor, job.n)

        pristine = SIQSJob(job.n, config=cfg, budget=allowance())
        pristine._setup()
        with patch(
            "v2.qs.root_updates.build_root_updates",
            side_effect=MemoryError("table fixture"),
        ):
            with self.assertRaises(MemoryError):
                pristine._next_step()
        self.assertIsNone(pristine.family)
        self.assertEqual(pristine.family_index, 0)
        self.assertIsNotNone(pristine._next_step())

    def test_canonical_standard_form_root_only_checkpoint(self):
        cfg = config_v4(
            multiplier=2,
            polynomial_policy="standard",
            multiplier_policy="standard",
            assignment_policy="flyer",
            sign_policy="canonical",
        )
        cfg = replace(cfg, collector=replace(cfg.collector, marking="sparse"))
        job = SIQSJob(4001 * 5003, config=cfg, budget=allowance())
        job._setup()
        self.assertTrue(job._activate())
        saved = job.checkpoint()
        self.assertEqual(saved["version"], 10)
        restored = SIQSJob.from_checkpoint(
            saved, config=cfg, budget=allowance()
        )
        self.assertEqual(restored.family.root_update_policy, "batch")
        self.assertEqual(restored.family.root_updates, job.family.root_updates)
        left, right = job.run(), restored.run()
        self.assertEqual(
            (left.reason, left.divisor), (right.reason, right.divisor)
        )
        self.assertEqual(left.divisor * left.cofactor, job.n)

    def test_norm_preview_grant_independence_and_invalid_configurations(self):
        cfg = config_v4()
        choices = [
            select_norm_multiplier(
                4001 * 5003,
                config=replace(cfg, family_count=count),
                seed=7,
                budget=allowance(),
            )
            for count in (1, 64)
        ]
        self.assertEqual(choices[0], choices[1])
        for changes in (
            dict(root_update_policy="unknown"),
            dict(root_update_policy="batch", mode="mpqs"),
            dict(assignment_policy="reference"),
            dict(factor_count=2),
        ):
            with self.assertRaises(ValueError):
                replace(cfg, **changes)
        tiny = replace(cfg, base_bound=3, multiplier=1)
        result = SIQSJob(4001 * 5003, config=tiny, budget=allowance()).run()
        self.assertEqual(result.reason, "configuration_refusal")
        self.assertIsNone(result.divisor)
        self.assertEqual(result.cofactor, 4001 * 5003)
