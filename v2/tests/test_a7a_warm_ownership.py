"""Warm the measured entry point and reconcile exceptional child charges."""

import itertools
import unittest
from unittest.mock import patch

from v2.benchmarks.infrastructure.resource_control import Budget
from v2.benchmarks.qs.a7a.cell import complete_call, execute
from v2.benchmarks.qs.a7a.combined import cell


class WarmOwnershipTests(unittest.TestCase):
    def test_literal_endpoint_override_uses_checked_collector_configuration(
        self,
    ):
        spec = cell(
            "S0",
            policy="Y0",
            endpoint="complete",
            lifetime="cold",
            samples=1,
            tier_a=False,
            observe=False,
            collector_overrides={
                "small_prime_bound": 2400000,
                "large_prime_bound": 9600000,
                "large_product_bound": 1728000000000000000,
            },
        )
        number = 10**79 + 7
        partial = dict(
            complete=False, factors=[], remaining=[number], reason="bounded"
        )

        with patch(
            "v2.benchmarks.qs.a7a.cell.complete", return_value=partial
        ) as complete:
            result = execute(dict(n=number, cell=spec))

        observed = complete.call_args.args[1].collector
        self.assertEqual(observed.small_prime_bound, 2400000)
        self.assertEqual(observed.large_prime_bound, 9600000)
        self.assertEqual(
            result["configuration"]["collector"]["large_product_bound"],
            1728000000000000000,
        )
        self.assertEqual(observed.candidate_bound, 115200000000)

        spec["collector_overrides"]["small_prime_bound"] = 9600001
        with patch("v2.benchmarks.qs.a7a.cell.complete") as complete:
            with self.assertRaises(ValueError):
                execute(dict(n=number, cell=spec))
        complete.assert_not_called()

    def test_complete_warmup_never_binds_a_bare_prefix_pool(self):
        spec = cell(
            "M0",
            policy="G0",
            endpoint="complete",
            lifetime="warm",
            samples=1,
            tier_a=False,
            observe=False,
            wall=120,
            cpu=1200,
            warmup_seconds=10,
        )
        calls = []

        def full_entry(n, config, workers, budget, pool, endpoint, **kwargs):
            calls.append((workers, pool, endpoint))
            budget.used += 5
            return dict(
                complete=False, factors=[], remaining=[n], reason="bounded"
            )

        with (
            patch("v2.benchmarks.qs.a7a.cell.prefix") as prefix,
            patch(
                "v2.benchmarks.qs.a7a.cell.complete", side_effect=full_entry
            ),
            patch(
                "v2.benchmarks.qs.a7a.cell.time.monotonic",
                side_effect=itertools.count(100, 1),
            ),
        ):
            result = execute(dict(n=10403, cell=spec))
        prefix.assert_not_called()
        self.assertGreaterEqual(len(calls), 2)
        self.assertTrue(all(pool is calls[0][1] for _, pool, _ in calls))
        self.assertTrue(
            all(endpoint == "complete" for _, _, endpoint in calls)
        )
        self.assertGreaterEqual(result["warm_seconds"], 10)
        self.assertEqual(result["warm_endpoint"], "complete")
        self.assertEqual(result["warmup_target"], 10)
        first, last = result["samples"][0]["measurement_interval"]
        self.assertLessEqual(first, last)

    def test_unfunded_or_unadmitted_warmup_is_rejected_before_configuration(
        self,
    ):
        for duration, wall in ((True, 120), (4, 120), (10, 20), (20, 30)):
            spec = cell(
                "S0", lifetime="warm", warmup_seconds=duration, wall=wall
            )
            with patch(
                "v2.benchmarks.qs.a7a.cell.configuration"
            ) as configuration:
                with self.assertRaises(ValueError):
                    execute(dict(n=10403, cell=spec))
            configuration.assert_not_called()

    def test_exception_retains_work_and_publishes_child_cpu_once(self):
        budget = Budget(
            work_limit=100, seconds=30, cpu_seconds=30, used=10, prior_cpu=1
        )

        def fail(n, config, workers, child, pool, endpoint, **kwargs):
            self.assertIs(child, budget)
            self.assertIsNotNone(child.cancelled)
            child.consume(5)
            child.prior_cpu += 7
            raise ValueError("owned failure")

        with patch("v2.benchmarks.qs.a7a.cell.complete", side_effect=fail):
            with self.assertRaisesRegex(ValueError, "owned failure"):
                complete_call(
                    10403,
                    None,
                    None,
                    budget,
                    None,
                    dict(seed=7, endpoint="complete", memory=2**26),
                    seconds_cap=3,
                )
        self.assertEqual(budget.used, 15)
        self.assertEqual(budget.prior_cpu, 8)


if __name__ == "__main__":
    unittest.main()
