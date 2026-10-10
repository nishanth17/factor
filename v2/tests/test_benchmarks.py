"""Independent checks for corpus validation and benchmark accounting."""

import io
import json
import unittest
from contextlib import redirect_stdout
from unittest.mock import MagicMock, patch

from v2.benchmarks.infrastructure.parallel import parallel_candidates
from v2.benchmarks.suites import phase_two


class BenchmarkTests(unittest.TestCase):
    """Reject incorrect answers even during untimed workload warmup."""

    def test_worker_cpu_requires_every_worker(self):
        """Duplicate PIDs cannot silently omit another worker's CPU use."""
        executor = MagicMock()
        executor.submit.return_value.result.return_value = (17, 1.0)
        with self.assertRaisesRegex(RuntimeError, "every worker"):
            parallel_candidates._worker_cpu(executor, 2)

    def test_worker_cpu_preserves_process_readings(self):
        """Use process totals rather than only candidate-body thread time."""
        futures = [MagicMock(), MagicMock()]
        for future, reading in zip(futures, ((17, 1.5), (23, 2.75))):
            future.result.return_value = reading
        executor = MagicMock()
        executor.submit.side_effect = futures

        self.assertEqual(
            parallel_candidates._worker_cpu(executor, 2),
            {17: 1.5, 23: 2.75},
        )

    def test_simultaneous_timeouts_do_not_interrupt_cleanup(self):
        """Inject a second timer expiry while the first is being disarmed."""
        handlers = {}

        def install(signum, handler):
            handlers[signum] = handler

        def disarm(timer, seconds):
            if not seconds:
                handlers[phase_two.signal.SIGPROF](
                    phase_two.signal.SIGPROF, None
                )

        def expire(*args, **kwargs):
            handlers[phase_two.signal.SIGALRM](phase_two.signal.SIGALRM, None)

        request = {
            "config": {},
            "n": 15,
            "seed": 7,
            "engine": "bounded",
            "mode": "complete",
            "seconds": 0.05,
            "cpu_seconds": 0.05,
            "work_limit": 100,
            "rss_cap_bytes": 2**40,
        }
        with (
            patch.object(phase_two.signal, "signal", side_effect=install),
            patch.object(phase_two.signal, "setitimer", side_effect=disarm),
            patch.object(phase_two, "factorize_bounded", side_effect=expire),
        ):
            row = phase_two._sample(request)

        self.assertEqual(row["reason"], "wall_limit")
        self.assertEqual(row["remaining"], [15])
        self.assertTrue(row["reconstructs"])
        self.assertFalse(row["success"])

    def test_warmup_checks_hidden_oracle(self):
        """A reconstructible composite terminal factor fails the oracle."""
        request = {
            "warmup": 3,
            "requests": [{"n": 15}],
            "fixtures": [{"n": 15, "factors": [[3, 1], [5, 1]]}],
        }
        invalid = {
            "reconstructs": True,
            "factors": [[15, 1]],
            "complete": True,
        }
        with (
            patch.object(phase_two.sys, "stdin", [json.dumps(request)]),
            patch.object(phase_two, "_sample", return_value=invalid),
            self.assertRaisesRegex(AssertionError, "terminal factor"),
        ):
            phase_two._worker()

    def test_cold_sample_preserves_requested_mode(self):
        """Cold factor-one requests include startup and shutdown costs."""
        fixture = {"id": "close", "band": "close_small", "factors": []}
        request = {"engine": "bounded", "mode": "factor_one", "seed": 7}
        row = {
            "cpu_seconds": 0.0,
            "reconstructs": True,
            "complete": False,
            "factors": [],
        }
        with patch.object(phase_two, "Worker") as factory:
            factory.return_value.call.return_value = row
            measured = phase_two._cold_sample(request, fixture)

            factory.return_value.call.assert_called_once_with(request)
            factory.return_value.close.assert_called_once_with()

        self.assertEqual(measured["mode"], "factor_one")
        self.assertEqual(measured["temperature"], "cold")
        self.assertGreaterEqual(measured["total_seconds"], 0)
        self.assertGreaterEqual(measured["cpu_seconds"], 0)

    def test_summary_includes_timeouts_and_failures(self):
        """Fast successful samples cannot hide censored slow outcomes."""
        rows = [
            {
                "engine": "bounded",
                "mode": "complete",
                "temperature": "warm",
                "band": "random_small",
                "total_seconds": seconds,
                "success": reason == "complete",
                "reason": reason,
                "cpu_seconds": seconds,
                "peak_rss_bytes": 1024,
            }
            for seconds, reason in (
                (1, "complete"),
                (9, "wall_limit"),
                (10, "work_limit"),
            )
        ]
        summary = phase_two.summarize(rows)[0]

        self.assertEqual(summary["median_seconds"], 9)
        self.assertEqual(summary["p95_seconds"], 10)
        self.assertEqual(summary["completion_fraction"], 1 / 3)
        self.assertEqual(summary["censored_timeouts"], 1)
        self.assertEqual(summary["work_exhaustions"], 1)

    def test_warmup_returns_validated_counts(self):
        """The worker consumes complete and factor-one warmup requests."""
        fixture = {"n": 15, "factors": [[3, 1], [5, 1]]}
        requests = [{"mode": mode} for mode in ("complete", "factor_one")]
        payload = {
            "warmup": 3,
            "requests": requests,
            "fixtures": [fixture, fixture],
        }
        answer = {
            "reconstructs": True,
            "factors": [[3, 1], [5, 1]],
            "complete": True,
        }
        output = io.StringIO()
        with (
            patch.object(phase_two.sys, "stdin", [json.dumps(payload)]),
            patch.object(phase_two, "_sample", return_value=answer) as sample,
            patch.object(
                phase_two.time, "perf_counter", side_effect=[0, 0, 4, 4]
            ),
            redirect_stdout(output),
        ):
            phase_two._worker()

        self.assertEqual(sample.call_count, 2)
        self.assertEqual(
            json.loads(output.getvalue()), {"seconds": 4, "calls": 2}
        )


if __name__ == "__main__":
    unittest.main()
