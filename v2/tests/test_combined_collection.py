"""Independent CRT/coverage oracles and policy-bound production resume."""

import hashlib
import itertools
import json
import unittest
from dataclasses import asdict, replace
from math import prod
from pathlib import Path
from unittest.mock import patch

from v2.execution.budget import BudgetExhaustedError
from v2.portfolio import _decode_siqs_config
from v2.qs import Polynomial, SIQSJob, build_factor_base
from v2.qs.coordinate_family import CoordinateFamily
from v2.qs.defaults import default_siqs_config
from v2.qs.sieve_collector import SieveCollector, SieveConfig
from v2.qs.sign_schedule import schedule_for
from v2.tests.qs_default_controls import previous_siqs_config
from v2.tests.test_qs import reference_positions, unlimited_budget
from v2.tests.test_qs_sieve import signature
from v2.tests.test_siqs import allowance, configuration


class CombinedCollectionTests(unittest.TestCase):
    def test_production_matches_retained_measured_combined_prefixes(self):
        path = Path(__file__).resolve().parents[1] / (
            "benchmarks/inputs/controls/c10a07_combined_adoption_checks.json"
        )
        for fixture in json.loads(path.read_text())["fixtures"]:
            job = SIQSJob(
                fixture["n"],
                config=previous_siqs_config(fixture["n"]),
                seed=fixture["seed"],
                budget=allowance(10**13),
            )
            result = job.run(max_blocks=fixture["blocks"])
            self.assertEqual(result.reason, "paused")
            collector = job.engine.collector
            value = dict(
                family=job.family_index,
                polynomials=job.stats["polynomials"],
                poly=collector.polynomial.identity,
                position=job.engine.next_position,
                roots=[asdict(root) for root in collector._roots],
                atoms=sorted(
                    (key, atom.sign, atom.exponents, atom.residual_product)
                    for key, atom in collector._atoms.items()
                ),
                rows=[asdict(row) for row in collector._rows],
                partials=sorted(collector.partial_ids),
                divisor=job.divisor,
            )
            digest = hashlib.sha256(
                json.dumps(value, sort_keys=True).encode()
            ).hexdigest()
            self.assertEqual(digest, fixture["arithmetic_signature"])
            for name, expected in fixture["counters"].items():
                self.assertEqual(
                    result.stats["engine"]["collector_counts"][name], expected
                )
            self.assertGreaterEqual(result.stats["workspace_bytes"], 0)

    def test_swapped_family_matches_independent_cartesian_crt(self):
        base = build_factor_base(1009 * 1013, bound=100).factor_base
        primes = tuple(p for p in base.primes if p != 2 and base.n_prime % p)[
            :4
        ]
        a = prod(primes)
        local = [
            min(x for x in range(p) if x * x % p == base.n_prime % p)
            for p in primes
        ]
        expected = set()
        for signs in itertools.product((1, -1), repeat=len(primes) - 1):
            roots = [local[0]] + [s * r for s, r in zip(signs, local[1:])]
            b = (
                sum(
                    r * (a // p) * pow(a // p, -1, p)
                    for r, p in zip(roots, primes)
                )
                % a
            )
            expected.add(min(b, -b % a))
        family = CoordinateFamily(
            base,
            primes,
            schedule=schedule_for(4, "swap_high"),
            budget=unlimited_budget(),
        )
        actual = set()
        for _ in expected:
            step = family.next()
            b = step.polynomial.b % a
            actual.add(min(b, -b % a))
            for roots in step.roots:
                brute = tuple(
                    x
                    for x in range(roots.prime)
                    if step.polynomial.value(x) % roots.prime == 0
                )
                self.assertEqual(tuple(sorted(roots.roots)), brute)
        self.assertEqual(actual, expected)

    def test_queue_matches_exact_division_on_signed_windows_and_tails(self):
        base = build_factor_base(10403, bound=100).factor_base
        for a, b in ((1, 102), (7, 1), (49, 8)):
            poly = Polynomial(10403, 1, a, b)
            expected = reference_positions(poly, base, -71, 80, 97)
            for cap in (0, 1, 2**20):
                cfg = SieveConfig(
                    score_policy="powers",
                    division="bucket",
                    block_width=17,
                    residual_bound=97,
                    sparse_queue_bytes=cap,
                    max_partials=4096,
                    preparation_cache_bytes=min(2 * 2**20, 2048 * 4096),
                    target_rows=4096,
                )
                collector = SieveCollector(
                    poly, base, config=cfg, budget=unlimited_budget()
                )
                result = collector.collect(-71, 80)
                self.assertEqual(signature(result), expected)
                self.assertGreaterEqual(collector._workspace, 0)
                if cap == 1:
                    self.assertGreater(result.stats["queue_refusals"], 0)

    def test_interrupted_queue_build_releases_reservation_and_replays(self):
        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 1, 102)
        cfg = SieveConfig(
            score_policy="powers",
            division="bucket",
            block_width=17,
            sparse_queue_bytes=2**20,
            residual_bound=97,
        )
        collector = SieveCollector(
            poly, base, config=cfg, budget=unlimited_budget()
        )
        before = collector._workspace
        with patch(
            "v2.qs.sparse_queue.build",
            side_effect=BudgetExhaustedError("work"),
        ):
            result = collector.collect(-71, 80)
        self.assertEqual(result.next_position, -71)
        self.assertEqual(collector._workspace, before)
        self.assertIsNone(collector._sparse_queue_key)
        result = collector.collect(-71, 80)
        self.assertEqual(
            signature(result), reference_positions(poly, base, -71, 80, 97)
        )

    def test_combined_job_checkpoint_preserves_policy_and_resources(self):
        cfg = configuration(
            sign_policy="swap_high",
            collector=SieveConfig(
                score_policy="powers",
                division="bucket",
                sparse_queue_bytes=2**20,
                residual_bound=500,
                max_partials=4096,
                preparation_cache_bytes=min(2 * 2**20, 2048 * 4096),
                target_rows=4096,
            ),
        )
        job = SIQSJob(4001 * 5003, config=cfg, seed=73, budget=allowance())
        job.run(max_blocks=1)
        saved = json.loads(json.dumps(job.checkpoint()))
        self.assertEqual(saved["version"], 7)
        restored = SIQSJob.from_checkpoint(saved, budget=allowance())
        self.assertEqual(restored.config, cfg)
        self.assertGreaterEqual(restored.budget.used, job.budget.used)
        left, right = job.run(), restored.run()
        self.assertEqual(
            (left.divisor, left.cofactor), (right.divisor, right.cofactor)
        )
        self.assertEqual(left.divisor * left.cofactor, job.n)
        with self.assertRaises(ValueError):
            SIQSJob.from_checkpoint(
                saved,
                budget=allowance(),
                config=replace(cfg, sign_policy="canonical"),
            )
        self.assertEqual(_decode_siqs_config(asdict(cfg)), cfg)

    def test_legacy_checkpoint_stays_canonical_and_factory_reaches_all_sizes(
        self,
    ):
        for digits in (20, 40, 45, 50, 51, 55, 60, 61, 65, 80):
            cfg = default_siqs_config(10 ** (digits - 1) + 123)
            self.assertEqual(cfg.sign_policy, "swap_high")
            self.assertEqual(_decode_siqs_config(asdict(cfg)), cfg)
        cfg = configuration()
        job = SIQSJob(4001 * 5003, config=cfg, budget=allowance())
        saved = job.checkpoint()
        payload = json.loads(saved["blob"])
        self.assertNotIn("sign_policy", payload["config"])
        self.assertNotIn("sparse_queue_bytes", payload["config"]["collector"])
        restored = SIQSJob.from_checkpoint(saved, budget=allowance())
        self.assertEqual(restored.config.sign_policy, "canonical")
        self.assertEqual(restored.config.collector.sparse_queue_bytes, 0)


if __name__ == "__main__":
    unittest.main()
