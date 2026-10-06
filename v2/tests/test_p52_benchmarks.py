"""Independent workload proofs and exclusive benchmark-window controls."""

import json
import subprocess
import unittest
from unittest.mock import Mock, patch

from v2 import ecm
from v2.benchmarks import p52_a3, p52_realistic
from v2.budget import Budget
from v2.portfolio import PortfolioConfig


class RealisticBenchmarkTests(unittest.TestCase):
    def test_frozen_control_ignores_active_helpers_and_pauses(self):
        control, stages = p52_a3.load_control()
        expected = stages.ecm.scalar_multiply(97, 3, 1, 1009, 2)
        with patch.object(
            ecm, "scalar_multiply", side_effect=AssertionError("active helper")
        ):
            self.assertEqual(
                stages.ecm.scalar_multiply(97, 3, 1, 1009, 2), expected
            )
        self.assertIsNot(control.Budget, Budget)

        n = 25013 * 25031
        config = control.PortfolioConfig(
            trial_bound=5,
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=((50, 1000, 2),),
            segment_size=16,
            max_input_bits=329,
            memory_bytes=16 * 2**20,
        )
        run = control.factorize_bounded(
            n,
            config=config,
            budget=control.Budget(work_limit=100, seconds=5, cpu_seconds=5),
        )

        self.assertEqual(run.reason, "work_limit")
        self.assertEqual(run.result.reconstruct(), n)

    def test_certified_strata_and_finite_deep_campaign(self):
        corpus = p52_realistic.load_corpus()
        strata = {(f["digits"], f["shape"]) for f in corpus["fixtures"]}
        self.assertEqual(
            strata,
            {
                (digits, shape)
                for digits in (30, 40, 60, 80)
                for shape in ("small10", "target", "balanced")
            },
        )
        config = PortfolioConfig(**p52_realistic.options("deep", "programs"))

        self.assertEqual(config.ecm_tiers, ((11000, 1900000, 256),))
        self.assertLess(config.workspace_reserve, config.memory_bytes)
        self.assertLessEqual(
            max(f["n"].bit_length() for f in corpus["fixtures"]),
            config.max_input_bits,
        )

    def test_missing_factor_proof_is_rejected(self):
        corpus = p52_realistic.load_corpus()
        corpus["fixtures"][0]["factors"][0][0] = 4
        with patch.object(
            type(p52_realistic.CORPUS),
            "read_text",
            return_value=json.dumps(corpus),
        ):
            with self.assertRaisesRegex(ValueError, "lacks.*certificate"):
                p52_realistic.load_corpus()

    def test_monitor_tracks_process_ownership(self):
        listing = "\n".join(
            (
                "100 1 /opt/pypy3 -m v2.benchmarks.p52_realistic",
                "101 100 /opt/pypy3 -m v2.benchmarks.p52_realistic --worker",
                "200 1 /opt/python -m v2.benchmarks.p43_sizes",
                "201 1 /opt/pypy3 -m unittest discover",
                "202 1 /opt/python -m http.server",
                "203 1 /opt/pypy3 -m v2.benchmarks.p52_realistic",
                "malformed row",
            )
        )
        competing = p52_realistic.competing_processes(listing, {100, 101})

        self.assertEqual([p["pid"] for p in competing], [200, 201, 203])

    def test_overlap_stops_only_the_owned_worker(self):
        process = Mock(pid=101)
        process.communicate.side_effect = [
            subprocess.TimeoutExpired("worker", 1),
            ("", "stopped"),
        ]
        with (
            patch.object(
                p52_realistic.subprocess, "Popen", return_value=process
            ),
            patch.object(
                p52_realistic,
                "check_quiet",
                side_effect=[None, RuntimeError("benchmark/test overlap")],
            ),
        ):
            with self.assertRaisesRegex(RuntimeError, "overlap"):
                p52_realistic.launch("middle", "baseline", 3, 9, quiet=True)

        process.terminate.assert_called_once_with()
        self.assertEqual(process.communicate.call_count, 2)
