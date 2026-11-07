"""Independent CRT peers, paid-prefix refusal and cached-merge resume."""

import unittest
from dataclasses import replace

from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.work_budget import PollingBudget
from v2.qs.coefficient_cache import CoefficientCache
from v2.qs.current_parallel import CurrentParallelSIQSJob
from v2.qs.worker_packet import polynomial_at
from v2.tests.test_current_qs_workers import (
    N,
    allowance,
    configuration,
    prepared,
    settings,
)


class CoefficientCacheTests(unittest.TestCase):
    def test_every_gray_sign_and_parity_matches_uncached_crt(self):
        for multiplier in (1, 2, 3, 5):
            for sign in ("canonical", "swap_high"):
                job = CurrentParallelSIQSJob(
                    N,
                    seed=29,
                    config=replace(
                        configuration(multiplier=multiplier), sign_policy=sign
                    ),
                    workers=settings(),
                    budget=allowance(),
                )
                job._setup()
                if job.source.base is None:
                    continue
                base, primes = job.source.base, job.source.assignments[0]
                cache = CoefficientCache(base, primes, job.config)
                for ordinal in range(job.config.gray_limit):
                    left, right = allowance(), allowance()
                    expected = polynomial_at(
                        base, primes, ordinal, job.config, left
                    )
                    actual = polynomial_at(
                        base,
                        primes,
                        ordinal,
                        job.config,
                        right,
                        coefficient_cache=cache,
                    )
                    self.assertEqual(actual, expected)
                    self.assertEqual(left.used, right.used)

    def test_cache_hits_keep_poll_mutation_and_refusal_boundaries(self):
        job = prepared()
        base, primes = job.source.base, job.source.assignments[0]
        cache = CoefficientCache(base, primes, job.config)
        polynomial_at(
            base, primes, 0, job.config, allowance(), coefficient_cache=cache
        )
        for limit in (0, 100, 1000, 10**8):
            results = []
            for cached in (False, True):
                ledger = Budget(
                    work_limit=limit, seconds=None, cpu_seconds=None
                )
                polls = []

                def shorten(shared):
                    polls.append(shared.used)
                    if len(polls) == 2:
                        shared.work_limit = shared.used

                budget = PollingBudget(ledger, interval=1, poll=shorten)
                try:
                    value = polynomial_at(
                        base,
                        primes,
                        1,
                        job.config,
                        budget,
                        coefficient_cache=cache if cached else None,
                    )
                except BudgetExhaustedError:
                    value = None
                results.append((value, ledger.used, ledger.reason, polls))
            self.assertEqual(*results)

    def test_template_cannot_cross_base_or_configuration_ownership(self):
        job = prepared()
        base, primes = job.source.base, job.source.assignments[0]
        cache = CoefficientCache(base, primes, job.config)
        with self.assertRaises(ValueError):
            polynomial_at(
                base,
                primes,
                0,
                replace(job.config),
                allowance(),
                coefficient_cache=cache,
            )
        self.assertNotIn("coefficient_cache", settings().record())
        with self.assertRaises(ValueError):
            replace(settings(), coefficient_cache=1)

    def test_cached_batch_merge_resume_preserves_exact_atoms_and_identity(
        self,
    ):
        reference = prepared(True)
        reference.run(fixed_work=True)
        options = replace(
            settings(),
            coefficient_cache=True,
            packed_transport=True,
            lease_budget=True,
            reuse_state=True,
            verification="batch",
        )
        job = prepared(True, worker_settings=options)
        job.run(max_assignments=1, fixed_work=True)
        saved = job.checkpoint()
        restored = CurrentParallelSIQSJob.from_checkpoint(
            saved, budget=allowance()
        )
        restored.run(fixed_work=True)
        self.assertEqual(
            dict(restored.engine.collector._atoms),
            dict(reference.engine.collector._atoms),
        )
        with self.assertRaises(ValueError):
            CurrentParallelSIQSJob.from_checkpoint(
                saved,
                budget=allowance(),
                workers=replace(options, coefficient_cache=False),
            )


if __name__ == "__main__":
    unittest.main()
