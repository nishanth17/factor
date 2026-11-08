"""Independent oracles and resume contracts for the combined C10a-01 arm."""

import itertools
import json
import unittest
from collections import defaultdict
from dataclasses import asdict, replace
from fractions import Fraction
from math import isqrt, prod
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.benchmarks.qs.c10a01.configuration import configuration as old_config
from v2.benchmarks.qs.c10a01.configuration_v3 import (
    configuration as study_config,
)
from v2.benchmarks.qs.c10a01.inputs import CORPUS
from v2.common.prime_sieve import small_sieve
from v2.execution.budget import BudgetExhaustedError
from v2.portfolio import _decode_siqs_config
from v2.qs import Polynomial, SIQSJob, build_factor_base
from v2.qs.assignment_stream import (
    AssignmentStream,
    nearest_factor_pair,
    paired_workspace,
)
from v2.qs.families import PolynomialFamily
from v2.qs.multiplier import (
    WIDE_MULTIPLIER_WORKSPACE,
    WIDE_MULTIPLIERS,
    select_coupled_multiplier,
)
from v2.qs.sieve_collector import SieveCollector, SieveConfig
from v2.qs.siqs import make_family
from v2.tests.resource_stops import InterruptBudget as Budget
from v2.tests.test_qs import reference_positions
from v2.tests.test_qs_sieve import signature
from v2.tests.test_siqs import allowance, configuration, mutate


def direct_scores(n):
    """Direct target Euler tests and Fraction arithmetic; no symbol reuse."""
    primes = tuple(
        p for p in range(2, 1988) if all(p % d for d in range(2, isqrt(p) + 1))
    )
    weights = {p: (p**256).bit_length() - 1 for p in primes}
    result = []
    for k in WIDE_MULTIPLIERS:
        residue = n * k % 8
        two = {1: 3, 5: 1}.get(residue, Fraction(1, 2))
        score = Fraction(-((k**256).bit_length() - 1), 2) + 256 * two
        for p in primes[1:]:
            value = n * k % p
            if value == 0:
                score += Fraction(weights[p], p)
            elif pow(value, (p - 1) // 2, p) == 1:
                score += Fraction(2 * weights[p], p - 1)
        result.append((k, score.numerator // score.denominator))
    return tuple(result)


def coupled_config(**options):
    values = dict(
        polynomial_policy="parity2",
        multiplier_policy="parity2_wide",
        assignment_policy="paired_flyer",
        multiplier=0,
        factor_count=3,
        pool_size=8,
        family_count=64,
        half_width=64,
        sign_policy="swap_high",
        collector=SieveConfig(
            marking="parity_two",
            score_policy="powers",
            division="bucket",
            block_width=17,
            sparse_queue_bytes=2**20,
            power_plan_bytes=2**18,
            residual_bound=500,
            max_partials=4096,
            preparation_cache_bytes=min(2 * 2**20, 2048 * 4096),
            target_rows=4096,
        ),
    )
    values.update(options)
    return configuration(**values)


class CoupledMultiplierTests(unittest.TestCase):
    def test_direct_integer_and_fraction_oracle_all_residue_classes(self):
        expected = tuple(
            k
            for k in range(1, 98)
            if all(k % (d * d) for d in range(2, isqrt(k) + 1))
        )
        self.assertEqual(WIDE_MULTIPLIERS, expected)
        self.assertEqual(len(expected), 61)
        self.assertEqual(sum(k % 2 == 0 for k in expected), 20)
        self.assertEqual(len(small_sieve(1988)), 300)
        primes = [
            p
            for p in range(4000, 4100)
            if all(p % d for d in range(2, isqrt(p) + 1))
        ]
        for residue in (1, 3, 5, 7):
            q = next(p for p in primes if p % 8 == residue)
            n = 4001 * q
            choice = select_coupled_multiplier(n, budget=allowance())
            oracle = direct_scores(n)
            self.assertIsNone(choice.divisor)
            self.assertEqual(choice.scores, oracle)
            self.assertEqual(
                choice.multiplier, max(oracle, key=lambda x: (x[1], -x[0]))[0]
            )

    def test_symbol_reuse_work_memory_and_exceptional_gcd(self):
        with patch("v2.qs.multiplier.pow", wraps=pow) as power:
            choice = select_coupled_multiplier(4001 * 5003)
        self.assertIsNone(choice.divisor)
        self.assertEqual(power.call_count, 299 * 26)
        self.assertLess(power.call_count, 299 * 61)
        for n in (101 * 4001, 3 * 4001, 1009 * 1013):
            split = select_coupled_multiplier(n)
            self.assertTrue(1 < split.divisor < n)
            self.assertEqual(n % split.divisor, 0)
        self.assertIsNone(select_coupled_multiplier(101).divisor)
        ledger = allowance()
        with self.assertRaises(TypeError):
            select_coupled_multiplier(
                4001 * 5003,
                budget=ledger,
                memory_bytes=WIDE_MULTIPLIER_WORKSPACE - 1,
            )
        self.assertEqual(ledger.used, 0)
        with self.assertRaises(BudgetExhaustedError):
            select_coupled_multiplier(4001 * 5003, budget=allowance(0))
        ledger.cancelled = lambda: ledger.used > 100000
        with self.assertRaises(BudgetExhaustedError):
            select_coupled_multiplier(4001 * 5003, budget=ledger)
        self.assertGreater(ledger.used, 0)
        self.assertEqual(select_coupled_multiplier(4001 * 5003), choice)

    def test_finite_default_charge_bound_at_large_input_width(self):
        radical = prod(small_sieve(1988))
        n = radical * 2 ** (4096 - radical.bit_length()) + 1
        self.assertLessEqual(n.bit_length(), 4096)
        ledger = Budget(work_limit=32_000_000, seconds=None, cpu_seconds=None)
        choice = select_coupled_multiplier(n, budget=ledger)
        self.assertEqual(len(choice.scores), 61)
        self.assertLess(ledger.used, ledger.work_limit)


class PairedAssignmentTests(unittest.TestCase):
    def test_monotone_minimum_against_cartesian_products_and_ties(self):
        for lower, upper in (
            ((3, 5, 7), (17, 19, 23)),
            ((2,), (11, 13)),
            ((3, 7), (19,)),
        ):
            for core in (1, 11, 1237):
                for target in itertools.chain(
                    range(1, 3000, 13),
                    (core * p * q for p in lower for q in upper),
                ):
                    actual = nearest_factor_pair(
                        lower, upper, core, target, allowance()
                    )
                    expected = min(
                        itertools.product(lower, upper),
                        key=lambda pair: (
                            abs(core * prod(pair) - target),
                            *pair,
                        ),
                    )
                    self.assertEqual(actual, expected)
        with self.assertRaises(ValueError):
            nearest_factor_pair((), (7,), 1, 100, allowance())

    def test_population_uniqueness_prefix_exhaustion_and_charged_retry(self):
        base = build_factor_base(
            4001 * 5003, multiplier=3, bound=300
        ).factor_base
        options = dict(
            factor_count=4,
            family_count=1000,
            pool_size=6,
            seed=17,
            policy="paired_flyer",
            memory_bytes=2**24,
        )
        stream = adapt_resource_configuration(
            AssignmentStream, base, 32, budget=allowance(), **options
        )
        short = adapt_resource_configuration(
            AssignmentStream,
            base,
            32,
            budget=allowance(),
            **dict(options, family_count=3),
        )
        values = [stream[i] for i in range(len(stream))]
        self.assertEqual(len(values), 15)
        self.assertEqual(len(set(map(prod, values))), 15)
        self.assertEqual(values[:3], [short[i] for i in range(3)])
        self.assertEqual(stream.identity, short.identity)
        cores = {tuple(p for p in row if p in stream.pool) for row in values}
        self.assertEqual(
            cores, set(itertools.combinations(sorted(stream.pool), 2))
        )
        for row in values:
            self.assertEqual(len(set(row)), 4)
            self.assertTrue(all(p != 2 and base.n_prime % p for p in row))
            core = prod(p for p in row if p in stream.pool)
            expected = min(
                itertools.product(stream.lower_flyers, stream.upper_flyers),
                key=lambda pair: (
                    abs(core * prod(pair) - stream.target),
                    *pair,
                ),
            )
            self.assertEqual(
                tuple(p for p in row if p not in stream.pool), expected
            )
        with self.assertRaises(IndexError):
            stream[len(stream)]
        before = dict(stream.__dict__)
        stream.budget = Budget(
            work_limit=10**9, seconds=None, cpu_seconds=None
        )
        stream.budget.cancelled = lambda: stream.budget.used > 40
        with self.assertRaises(BudgetExhaustedError):
            stream[2]
        stream.budget = allowance()
        self.assertEqual(stream[2], values[2])
        for key in before:
            if key != "budget":
                self.assertEqual(stream.__dict__[key], before[key])
        ledger = allowance()
        with self.assertRaises(TypeError):
            adapt_resource_configuration(
                AssignmentStream,
                base,
                32,
                budget=ledger,
                **dict(options, memory_bytes=stream.workspace_bytes - 1),
            )
        self.assertEqual(ledger.used, 0)
        self.assertGreaterEqual(
            stream.workspace_bytes, paired_workspace(len(base.entries), 4)
        )
        with self.assertRaises(ValueError):
            adapt_resource_configuration(
                AssignmentStream,
                base,
                32,
                budget=allowance(),
                **dict(options, factor_count=2),
            )


class ParityInitializationTests(unittest.TestCase):
    def compare_buffers(self, polynomial, base, lo, hi, **options):
        values = dict(
            block_width=17,
            score_policy="powers",
            division="bucket",
            memory_bytes=64 * 2**20,
            residual_bound=97,
            max_atoms=4096,
            max_relations=4096,
            max_partials=4096,
        )
        values.update(options)
        marking = values.pop("marking", "sparse")
        opt_marking = (
            "parity_two" if marking == "sparse" else marking + "+parity_two"
        )
        collectors = [
            SieveCollector(
                polynomial,
                base,
                budget=allowance(),
                config=adapt_resource_configuration(
                    SieveConfig, marking=kind, **values
                ),
            )
            for kind in (marking, opt_marking)
        ]
        for collector in collectors:
            if values.get("power_plan_bytes") or values.get(
                "sparse_queue_bytes"
            ):
                collector._set_plan_interval(lo, hi)
        seeded = 0
        for left in range(lo, hi, values["block_width"]):
            right = min(hi, left + values["block_width"])
            rows = []
            for collector in collectors:
                stats = defaultdict(int, threshold_min=2**63)
                threshold = collector._sieve(left, right, stats)
                rows.append(
                    (
                        threshold,
                        tuple(collector._scores[: right - left]),
                        tuple(collector._hits[: right - left]),
                        stats["root_hits"],
                        stats["power_hits"],
                    )
                )
                if collector is collectors[1]:
                    seeded += stats["parity_two_seeded_hits"]
            self.assertEqual(*rows)
        results = [collector.collect(lo, hi) for collector in collectors]
        self.assertEqual(signature(results[0]), signature(results[1]))
        if results[1].divisor is None:
            self.assertEqual(
                signature(results[1]),
                reference_positions(polynomial, base, lo, hi, 97),
            )
        else:
            self.assertEqual(results[0].divisor, results[1].divisor)
            self.assertTrue(1 < results[1].divisor < polynomial.n)
            self.assertEqual(polynomial.n % results[1].divisor, 0)
        return seeded

    def test_complete_buffers_roots_refinement_queues_plans_and_fallbacks(
        self,
    ):
        for multiplier in (3, 7):
            base = build_factor_base(
                10403, multiplier=multiplier, bound=97
            ).factor_base
            primes = tuple(
                p for p in base.primes if p != 2 and base.n_prime % p
            )[:3]
            polynomial = (
                PolynomialFamily(base, primes, parity=True, budget=allowance())
                .next()
                .polynomial
            )
            for options in (
                {},
                {"score_policy": "fixed"},
                {"power_plan_bytes": 2**20},
                {"sparse_queue_bytes": 2**20},
                {"power_plan_bytes": 1, "sparse_queue_bytes": 1},
            ):
                self.assertGreater(
                    self.compare_buffers(polynomial, base, -71, 80, **options),
                    0,
                )
            for options in (
                {"score_backend": "bytearray"},
                {"score_backend": "array"},
                {"score_policy": "candidate"},
                {"small_prime_cutoff": 3},
                {"marking": "density+presieve"},
                {"marking": "presieve_powers"},
            ):
                self.assertEqual(
                    self.compare_buffers(polynomial, base, -31, 40, **options),
                    0,
                )
            for width in (1, 2):
                self.compare_buffers(
                    polynomial,
                    base,
                    -5,
                    7,
                    block_width=width,
                    sparse_queue_bytes=2**20,
                )
        base = build_factor_base(9, bound=3).factor_base
        polynomial = Polynomial(9, 1, 2, 3)
        self.compare_buffers(polynomial, base, 0, 1, block_width=1)
        self.compare_buffers(
            polynomial, base, -3, 4, block_width=1, power_plan_bytes=2**20
        )

    def test_high_two_valuation_retains_every_higher_level(self):
        n = 2**41 + 1
        base = build_factor_base(n, bound=3).factor_base
        polynomial = Polynomial(n, 1, 2, 1)
        self.assertEqual(polynomial.value(0), -(2**40))
        for options in ({}, {"power_plan_bytes": 2**20}):
            self.assertGreater(
                self.compare_buffers(polynomial, base, -5, 6, **options), 0
            )

    def test_cancelled_initialization_leaves_position_retryable(self):
        base = build_factor_base(10403, multiplier=3, bound=97).factor_base
        primes = tuple(p for p in base.primes if p > 2 and base.n_prime % p)[
            :3
        ]
        polynomial = (
            PolynomialFamily(base, primes, parity=True, budget=allowance())
            .next()
            .polynomial
        )
        collector = SieveCollector(
            polynomial,
            base,
            budget=allowance(),
            config=SieveConfig(
                marking="parity_two",
                score_policy="powers",
                block_width=17,
                residual_bound=97,
            ),
        )
        used = collector.budget.used
        collector.budget.cancelled = lambda: collector.budget.used > used + 10
        result = collector.collect(-7, 10)
        self.assertEqual(result.next_position, -7)
        self.assertEqual(result.reason, "cancelled")
        collector.budget.cancelled = None
        result = collector.collect(-7, 10)
        self.assertEqual(
            signature(result),
            reference_positions(polynomial, base, -7, 10, 97),
        )


class CoupledIntegrationTests(unittest.TestCase):
    def test_configurations_preserve_literal_control_and_matched_limits(self):
        for digits in (60, 62, 65, 70):
            control, coupled = (
                study_config(digits, arm) for arm in ("control", "coupled")
            )
            self.assertEqual(control, old_config(digits, "control"))
            expected = asdict(old_config(digits, "parity2"))
            expected["assignment_policy"] = "paired_flyer"
            expected["multiplier_policy"] = "parity2_wide"
            expected["collector"]["marking"] = "parity_two"
            self.assertEqual(asdict(coupled), expected)
            self.assertEqual(coupled.memory_bytes, control.memory_bytes)
            self.assertEqual(_decode_siqs_config(asdict(coupled)), coupled)

    def test_four_size_structural_setup_and_exact_norms_without_sieving(self):
        fixtures = json.loads(CORPUS.read_text())["fixtures"]
        for fixture in fixtures:
            digits = len(str(fixture["n"]))
            job = SIQSJob(
                fixture["n"],
                config=study_config(digits, "coupled"),
                budget=allowance(2_000_000_000),
            )
            job._setup()
            self.assertIsNone(job.divisor)
            self.assertEqual(len(job.base.entries), job.config.base_count)
            self.assertTrue(
                job.stats["capacity"]["paired_target_in_completion_envelope"]
            )
            for index in (0, 1, 7):
                primes = job.assignments[index]
                self.assertEqual(len(primes), job.config.factor_count)
                family = make_family(
                    job.base, primes, config=job.config, budget=job.budget
                )
                polynomial = family.next().polynomial
                for x in (-65536, -1, 0, 1, 65536):
                    self.assertEqual(
                        polynomial.u_value(x) ** 2 - job.base.n_prime,
                        polynomial.a * polynomial.value(x),
                    )
            self.assertIsNone(job.engine)
            self.assertEqual(job.stats["blocks"], 0)

    def test_each_new_selector_requires_version_nine(self):
        old = coupled_config(
            multiplier=3,
            multiplier_policy="parity2",
            assignment_policy="flyer",
        )
        old = replace(old, collector=replace(old.collector, marking="sparse"))
        job = SIQSJob(4001 * 5003, config=old, budget=allowance())
        job.run(max_blocks=1)
        saved = job.checkpoint()
        self.assertEqual(saved["version"], 8)
        changes = (
            lambda p: p["config"].update(multiplier_policy="parity2_wide"),
            lambda p: p["config"].update(assignment_policy="paired_flyer"),
            lambda p: p["config"]["collector"].update(marking="parity_two"),
        )
        for change in changes:
            with self.assertRaisesRegex(ValueError, "legacy checkpoint"):
                SIQSJob.from_checkpoint(
                    mutate(saved, change), budget=allowance()
                )

    def test_complete_factor_resume_pending_and_policy_substitution(self):
        for multiplier in (0, 2, 3, 7):
            cfg = coupled_config(multiplier=multiplier)
            job = SIQSJob(4001 * 5003, config=cfg, seed=73, budget=allowance())
            paused = job.run(max_blocks=1)
            self.assertEqual(paused.reason, "paused")
            saved = json.loads(json.dumps(job.checkpoint()))
            self.assertEqual(saved["version"], 9)
            restored = SIQSJob.from_checkpoint(
                saved, budget=allowance(), config=cfg
            )
            self.assertGreaterEqual(restored.budget.used, job.budget.used)
            left, right = job.run(), restored.run()
            self.assertEqual(
                (left.reason, left.divisor), (right.reason, right.divisor)
            )
            self.assertEqual(left.reason, "factor_found")
            self.assertEqual(left.divisor * left.cofactor, job.n)
            for changed in (
                replace(cfg, multiplier_policy="parity2"),
                replace(cfg, assignment_policy="flyer"),
                replace(
                    cfg, collector=replace(cfg.collector, marking="sparse")
                ),
            ):
                with self.assertRaises(ValueError):
                    SIQSJob.from_checkpoint(
                        saved,
                        config=changed,
                        budget=allowance(),
                        allow_extension=True,
                    )
            bad = mutate(saved, lambda p: p.update(version=8))
            bad["version"] = 8
            with self.assertRaisesRegex(ValueError, "legacy checkpoint"):
                SIQSJob.from_checkpoint(bad, budget=allowance())
        job = SIQSJob(
            4001 * 5003,
            config=coupled_config(multiplier=3),
            budget=allowance(),
        )
        job._setup()
        original = job._next_step

        def cancel_after_step():
            step = original()
            job.budget.cancelled = lambda: True
            return step

        with patch.object(job, "_next_step", side_effect=cancel_after_step):
            self.assertEqual(job.run(max_blocks=1).reason, "cancelled")
        self.assertIsNotNone(job.pending_step)
        saved = job.checkpoint()
        job.budget.cancelled = None
        restored = SIQSJob.from_checkpoint(saved, budget=allowance())
        left, right = job.run(), restored.run()
        self.assertEqual(
            (left.divisor, left.cofactor), (right.divisor, right.cofactor)
        )
        self.assertEqual(left.divisor * left.cofactor, job.n)
