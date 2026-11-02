"""Independent coverage of the selected growing-DLP fixed-score path."""

import json
import unittest
from dataclasses import replace
from unittest.mock import patch

from v2.benchmarks.qs.a1b import CONTROL
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.portfolio import _decode_siqs_config
from v2.qs import Polynomial, SieveCollector, SIQSJob, build_factor_base
from v2.tests.test_qs import reference_factors


def allowance():
    return Budget(work_limit=10**11, seconds=15, cpu_seconds=15)


def selected(digits=55):
    data = json.loads(CONTROL.read_text())
    return _decode_siqs_config(data["configurations"][str(digits)]["fixed"])


def oracle(poly, base, position, config):
    """Recover exponents and the allowed SLP/DLP domain by trial division."""
    value = poly.a * poly.value(position)
    if value == 0:
        return None
    factors = reference_factors(value)
    inside = tuple((p, e) for p, e in factors if p in base.primes)
    outside = tuple(
        p for p, e in factors if p not in base.primes for _ in range(e)
    )
    admissible = (
        not outside
        or (len(outside) == 1 and outside[0] <= config.residual_bound)
        or (
            len(outside) == 2
            and outside[-1] <= config.large_prime_bound
            and outside[0] * outside[1] <= config.large_product_bound
        )
    )
    if not admissible:
        return None
    return (
        -1 if value < 0 else 1,
        inside,
        outside[0] if len(outside) == 1 else 1,
        outside if len(outside) == 2 else (),
    )


class A1bContractTests(unittest.TestCase):
    def test_selected_growing_dlp_windows_exceptions_and_recovery(self):
        saw_pair = False
        for digits in (50, 55, 60):
            cfg = replace(
                selected(digits).collector,
                memory_bytes=selected(digits).memory_bytes,
            )
            for multiplier in (1, 3):
                base = build_factor_base(
                    10403, bound=100, multiplier=multiplier
                ).factor_base
                coefficients = [(1, 102)]
                for entry in base.entries:
                    if entry.prime != 2 and len(coefficients) < 4:
                        coefficients.append(
                            (entry.prime, entry.square_roots[0])
                        )
                for a, b in coefficients:
                    poly = Polynomial(10403, multiplier, a, b)
                    worker = SieveCollector(
                        poly, base, config=cfg, budget=allowance()
                    )
                    expected = {
                        x: oracle(poly, base, x, cfg) for x in range(-71, 80)
                    }
                    stats = dict.fromkeys(
                        (
                            "blocks",
                            "root_hits",
                            "power_hits",
                            "plan_builds",
                            "plan_replays",
                            "plan_refusals",
                            "slice_bytes",
                            "slice_allocations",
                            "threshold_max",
                        ),
                        0,
                    )
                    stats["threshold_min"] = 2**63
                    with patch(
                        "v2.qs.sparse_queue.build",
                        side_effect=AssertionError("fixed queued"),
                    ):
                        worker._set_plan_interval(-71, 80)
                        threshold = worker._sieve(-71, 80, stats)
                        for x, truth in expected.items():
                            if truth is not None:
                                self.assertGreaterEqual(
                                    worker._scores[x + 71], threshold
                                )
                                self.assertTrue(
                                    worker._candidate_passes(
                                        poly.value(x), x + 71
                                    )
                                )
                        result = worker.collect(-71, 80)
                    self.assertIn(result.reason, ("complete", "factor_found"))
                    self.assertNotIn("queue_builds", result.stats)
                    for atom in result.atoms:
                        actual = (
                            atom.sign,
                            atom.exponents,
                            atom.residual,
                            atom.large_primes,
                        )
                        self.assertEqual(actual, expected[atom.position])
                        saw_pair |= bool(atom.large_primes)
                    if result.divisor is not None:
                        self.assertTrue(1 < result.divisor < base.n)
                        self.assertEqual(base.n % result.divisor, 0)
        self.assertTrue(saw_pair)

    def test_fixed_growing_high_valuation_clipping_and_zero(self):
        b = 2**298 - 1
        poly = Polynomial(b * b, 1, 1, b)
        base = build_factor_base(b * b, bound=3).factor_base
        for backend in ("list", "array", "bytearray"):
            cfg = replace(
                selected().collector,
                score_backend=backend,
                memory_bytes=selected().memory_bytes,
            )
            worker = SieveCollector(poly, base, config=cfg, budget=allowance())
            result = worker.collect(2, 3)
            self.assertEqual(result.reason, "complete")
            self.assertEqual(result.atoms[0].exponents, ((2, 300),))
            if backend != "list":
                clip = 255 if backend == "bytearray" else 2**32 - 1
                worker._add((0,), 2**40, {})
                self.assertEqual(worker._scores[0], clip)
                self.assertEqual(worker._clip_threshold(2**40), clip)
            zero = worker.collect(0, 1)
            self.assertEqual(zero.reason, "factor_found")
            self.assertEqual(zero.divisor, b)
            self.assertEqual(zero.divisor * (base.n // zero.divisor), base.n)

    def test_fixed_growing_refusal_cancellation_and_exact_resume(self):
        base = build_factor_base(10403, bound=100).factor_base
        poly = Polynomial(10403, 1, 49, 8)
        cfg = replace(
            selected().collector, memory_bytes=selected().memory_bytes
        )
        worker = SieveCollector(poly, base, config=cfg, budget=allowance())
        worker.budget.work_limit = worker.budget.used
        refused = worker.collect(-71, 80)
        self.assertEqual(refused.reason, "work_limit")
        self.assertEqual(refused.next_position, -71)
        worker.budget = allowance()
        worker.budget.cancelled = lambda: True
        cancelled = worker.collect(-71, 80)
        self.assertEqual(cancelled.reason, "cancelled")
        self.assertEqual(cancelled.next_position, -71)
        worker.budget = allowance()
        resumed = worker.collect(-71, 80)
        clean = SieveCollector(
            poly, base, config=cfg, budget=allowance()
        ).collect(-71, 80)
        self.assertEqual(resumed, clean)
        with self.assertRaises(MemoryError):
            SieveCollector(
                poly,
                base,
                config=replace(cfg, memory_bytes=1),
                budget=allowance(),
            )

    def test_literal_selected_job_native_outcome_and_charged_restore(self):
        data = json.loads(CONTROL.read_text())
        n = next(f["n"] for f in data["fixtures"] if f["digits"] == 55)
        for arm in ("control", "fixed"):
            cfg = _decode_siqs_config(data["configurations"]["55"][arm])
            job = SIQSJob(n, seed=7, config=cfg, budget=allowance())
            first = job.run(max_blocks=1)
            self.assertEqual(first.reason, "paused")
            self.assertIsNone(first.divisor)
            self.assertEqual(first.cofactor, n)
            saved = json.loads(json.dumps(job.checkpoint()))
            restored = SIQSJob.from_checkpoint(saved, budget=allowance())
            self.assertEqual(restored.config, cfg)
            self.assertGreater(
                restored.budget.used, saved["resources"]["work_used"]
            )
            for change in (
                dict(score_policy="powers" if arm == "fixed" else "fixed"),
                dict(score_backend="array"),
            ):
                with self.assertRaises(ValueError):
                    SIQSJob.from_checkpoint(
                        saved,
                        budget=allowance(),
                        config=replace(
                            cfg, collector=replace(cfg.collector, **change)
                        ),
                    )
            stopped = allowance()
            stopped.cancelled = lambda: True
            with self.assertRaises(BudgetExhaustedError):
                SIQSJob.from_checkpoint(saved, budget=stopped)
            left, right = job.run(max_blocks=1), restored.run(max_blocks=1)
            self.assertEqual(
                (left.divisor, left.cofactor, left.next_position),
                (right.divisor, right.cofactor, right.next_position),
            )
            self.assertEqual(
                job.engine.collector._atoms, restored.engine.collector._atoms
            )


if __name__ == "__main__":
    unittest.main()
