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
            observer.__exit__.assert_called_once()
