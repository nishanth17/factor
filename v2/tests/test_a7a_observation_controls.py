"""Memory observation must not force per-relation production counters."""

import unittest
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from v2.benchmarks.qs.a7a import cell
from v2.benchmarks.qs.a7a.combined import cell as make_cell


class ObservationControlTests(unittest.TestCase):
    def test_rss_stays_enabled_when_detailed_worker_counters_are_disabled(
        self,
    ):
        for explicit in (False, True):
            spec = make_cell(
                "C1",
                observe=True,
                endpoint="prefix",
                lifetime="cold",
                samples=1,
                tasks=1,
                tier_a=False,
                features={"assignment_cpu": 5},
            )
            if explicit:
                spec["worker_observe"] = False
                spec["rss_stop_bytes"] = 18 * 2**30
                spec["memory"] = 24 * 2**30
            observer = MagicMock()
            observer.error = None
            observer.exceeded = False
            observer.rows = [dict(bytes=123)]
            observer.backend = "test"
            pool = SimpleNamespace(close=lambda: None)

            with (
                patch("v2.qs.parallel.CollectionPool", return_value=pool),
                patch(
                    "v2.benchmarks.qs.a7a.observer.OwnedRSS",
                    return_value=observer,
                ) as rss,
                patch.object(
                    cell,
                    "prefix",
                    return_value=dict(complete=False, reason="paused"),
                ) as collect,
            ):
                result = cell.execute(dict(n=15, cell=spec))

            self.assertEqual(result["rss_samples"], [dict(bytes=123)])
            self.assertEqual(
                result["worker_configuration"]["observe"], not explicit
            )
            self.assertEqual(collect.call_args.args[2].observe, not explicit)
            self.assertEqual(collect.call_args.args[2].assignment_cpu, 5)
            rss.assert_called_once()
            self.assertEqual(
                rss.call_args.kwargs["limit"],
                spec.get("rss_stop_bytes", spec["memory"]),
            )
            observer.__exit__.assert_called_once()

    def test_guard_requires_an_active_observer_and_a_finite_lower_cap(self):
        for observed, limit in ((False, 1), (True, 0), (True, 25 * 2**30)):
            spec = make_cell(
                "S0", observe=observed, rss_stop_bytes=limit, memory=24 * 2**30
            )
            with self.assertRaisesRegex(ValueError, "live RSS guard"):
                cell.execute(dict(n=91, cell=spec))

    def test_cooperative_cancel_preserves_a_valid_resumable_checkpoint(self):
        from v2.execution.budget import Budget
        from v2.portfolio import PortfolioConfig, factorize_bounded

        n = 91
        config = cell.configuration(n, "G0")
        budget = Budget(work_limit=10**8, seconds=2, cpu_seconds=2)
        budget.cancelled = lambda: True
        spec = dict(
            seed=7,
            endpoint="complete",
            memory=24 * 2**30,
            preserve_checkpoint=True,
        )

        result = cell.complete_call(n, config, None, budget, None, spec)

        self.assertEqual(result["reason"], "cancelled")
        self.assertEqual(result["remaining"], [n])
        resumed = factorize_bounded(
            n,
            seed=7,
            config=PortfolioConfig(
                siqs=config,
                rho_attempts=0,
                pm1_attempts=0,
                ecm_tiers=(),
                memory_bytes=spec["memory"],
            ),
            budget=Budget(work_limit=10**8, seconds=2, cpu_seconds=2),
            checkpoint=result["checkpoint"],
        )
        self.assertTrue(resumed.result.complete)
