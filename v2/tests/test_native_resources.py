"""Native root observation equivalence, fallback and refusing-poll lifetime."""

import errno
import os
import subprocess
import sys
import unittest
from types import SimpleNamespace
from unittest.mock import Mock

from v2.ecm.native_packet import SOURCE, NativePacket
from v2.ecm.native_resources import (
    NativeRootProvider,
    close_resource_checks,
    poll_native_budget,
)
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.resources import (
    MacResources,
    ProcessSample,
    ResourceProviderError,
)
from v2.tests.test_native_correctness import deterministic_budget


@unittest.skipUnless(
    (SOURCE / "build/manifest.json").exists(), "scalar native build required"
)
class NativeResourceTests(unittest.TestCase):
    def packet(self):
        packet = NativePacket()
        packet.available(1013, deterministic_budget())
        self.addCleanup(packet.close)
        return packet

    def fake(self, status=0, rss=10):
        packet = self.packet()

        def observe(pid, output):
            self.assertEqual(pid, os.getpid())
            output.birth, output.rss = 7, rss
            output.user, output.system = 1_234_567_890, 2_000_000_000
            return status

        owner = SimpleNamespace(
            ffi=packet.ffi,
            lib=SimpleNamespace(b4_observe_root=Mock(side_effect=observe)),
        )
        return owner

    def test_exact_fields_and_expected_birth(self):
        owner = self.fake()
        provider = NativeRootProvider(owner)
        self.addCleanup(provider.close)

        sample = provider.sample({os.getpid(): None})

        self.assertEqual(
            sample, [ProcessSample(os.getpid(), 7, 10, 3.23456789)]
        )
        self.assertEqual(provider.sample({os.getpid(): 8}), [])
        self.assertEqual(provider.sample({os.getpid(): 7}), sample)

    def test_child_and_multiple_root_paths_use_original_tree_provider(self):
        for status, roots in (
            (1, {os.getpid(): None}),
            (0, {os.getpid(): None, 123456: None}),
        ):
            owner = self.fake(status)
            provider = NativeRootProvider(owner)
            self.addCleanup(provider.close)
            fallback = SimpleNamespace(sample=Mock(return_value=["tree"]))
            provider.fallback = fallback

            self.assertEqual(provider.sample(roots), ["tree"])

            fallback.sample.assert_called_once_with(roots)

    def test_missing_process_and_observation_failure_remain_distinct(self):
        for status in (-errno.ESRCH, -errno.EPERM):
            provider = NativeRootProvider(self.fake(status))
            self.addCleanup(provider.close)
            if status == -errno.ESRCH:
                self.assertEqual(provider.sample({os.getpid(): None}), [])
            else:
                with self.assertRaises(ResourceProviderError):
                    provider.sample({os.getpid(): None})

    def test_rss_refusal_restores_original_provider(self):
        original = object.__new__(MacResources)
        original.sample = Mock(
            return_value=[ProcessSample(os.getpid(), 7, 10, 0)]
        )
        budget = Budget(
            provider=original,
            rss_bytes=100,
            wall_clock=lambda: 0,
            cpu_clock=lambda: 0,
        )
        owner = self.fake(rss=100)
        self.addCleanup(close_resource_checks, owner)

        with self.assertRaisesRegex(BudgetExhaustedError, "rss_limit"):
            poll_native_budget(budget, owner)

        self.assertIs(budget.provider, original)
        self.assertEqual(budget.observed_rss, 100)
        self.assertIsNone(owner._resource_provider.fallback)

    def test_custom_provider_and_cancellation_keep_original_semantics(self):
        budget = deterministic_budget()
        original = budget.provider
        owner = self.fake()

        poll_native_budget(budget, owner)

        self.assertIs(budget.provider, original)
        owner.lib.b4_observe_root.assert_not_called()
        budget.cancelled = lambda: True
        with self.assertRaisesRegex(BudgetExhaustedError, "cancelled"):
            poll_native_budget(budget, owner)
        self.assertIs(budget.provider, original)

    def test_actual_os_fields_match_original_usage_units(self):
        packet = self.packet()
        original = MacResources()
        provider = NativeRootProvider(packet)
        self.addCleanup(provider.close)
        provider.fallback = SimpleNamespace(
            sample=Mock(
                side_effect=AssertionError("childless fast path fell back")
            )
        )
        roots = {os.getpid(): None}

        before = original.sample(roots)[0]
        actual = provider.sample(roots)[0]
        after = original.sample(roots)[0]

        self.assertEqual(actual.birth, before.birth)
        self.assertEqual(actual.birth, after.birth)
        self.assertGreater(actual.rss, 0)
        self.assertLessEqual(before.cpu, actual.cpu)
        self.assertLessEqual(actual.cpu, after.cpu)
        provider.fallback.sample.assert_not_called()

    def test_actual_live_child_is_included_by_fallback(self):
        packet = self.packet()
        provider = NativeRootProvider(packet)
        self.addCleanup(provider.close)
        original = MacResources()
        original.sample = Mock(wraps=original.sample)
        provider.fallback = original
        child = subprocess.Popen(
            [sys.executable, "-B", "-c", "import time; time.sleep(5)"]
        )
        try:
            samples = provider.sample({os.getpid(): None})

            self.assertIn(child.pid, {sample.pid for sample in samples})
            original.sample.assert_called_once()
        finally:
            child.terminate()
            child.wait(timeout=2)
