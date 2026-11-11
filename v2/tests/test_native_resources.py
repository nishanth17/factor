"""Native root observation equivalence, fallback and refusing-poll lifetime."""

import errno
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

from v2.ecm.native_curve import BlockStore
from v2.ecm.native_packet import SOURCE, NativePacket
from v2.ecm.native_resources import (
    NativeRootProvider,
    close_resource_checks,
    poll_native_budget,
)
from v2.execution.budget import Budget, BudgetExhaustedError
from v2.execution.resources import (
    MacResources,
    ProcessObservationDeniedError,
    ProcessSample,
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

    @staticmethod
    def resources():
        provider = object.__new__(MacResources)
        provider._cpu_numer = 1
        provider._cpu_denom = 1_000_000_000
        return provider

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

    @unittest.skipUnless(
        os.environ.get("FACTOR_B4A_LIFECYCLE_QA") == "1",
        "compiled libproc proof requires the locked native QA window",
    )
    def test_checked_child_query_distinguishes_empty_and_denied(self):
        manifest = json.loads((SOURCE / "build/manifest.json").read_text())
        probe = (
            Path(__file__).resolve().parents[1]
            / "benchmarks/ecm/b4a/resource_query_probe.c"
        )
        prefix = Path(manifest["gmp_prefix"])
        with tempfile.TemporaryDirectory(
            prefix="b4a-resource-proof-"
        ) as folder:
            binary = Path(folder) / "probe"
            command = [
                manifest["command"][0],
                "-std=c11",
                "-O2",
                "-Wall",
                "-Wextra",
                "-Werror",
                "-isysroot",
                manifest["sysroot"],
                '-DB4_PACKET_SOURCE="' + str(SOURCE / "packet.c") + '"',
                str(probe),
                "-I" + str(prefix / "include"),
                "-L" + str(prefix / "lib"),
                "-lgmp",
                "-o",
                str(binary),
            ]
            subprocess.run(
                command, check=True, capture_output=True, timeout=20
            )
            subprocess.run([str(binary)], check=True, timeout=3)

    @unittest.skipUnless(
        os.environ.get("FACTOR_B4A_LIFECYCLE_QA") == "1",
        "compiled pthread proof requires the locked native QA window",
    )
    def test_entry_observer_progresses_without_pausing_workers(self):
        manifest = json.loads((SOURCE / "build/manifest.json").read_text())
        probe = (
            Path(__file__).resolve().parents[1]
            / "benchmarks/ecm/b4a/notification_probe.c"
        )
        prefix = Path(manifest["gmp_prefix"])
        with tempfile.TemporaryDirectory(prefix="b4a-notify-proof-") as folder:
            binary = Path(folder) / "probe"
            command = [
                manifest["command"][0],
                "-std=c11",
                "-O2",
                "-Wall",
                "-Wextra",
                "-Werror",
                "-pthread",
                "-DB4_PTHREAD",
                "-isysroot",
                manifest["sysroot"],
                '-DB4_PACKET_SOURCE="' + str(SOURCE / "packet.c") + '"',
                '-DB4_TEAM_SOURCE="' + str(SOURCE / "team.c") + '"',
                str(probe),
                "-I" + str(prefix / "include"),
                "-L" + str(prefix / "lib"),
                "-lgmp",
                "-o",
                str(binary),
            ]
            subprocess.run(
                command, check=True, capture_output=True, timeout=20
            )
            subprocess.run([str(binary)], check=True, timeout=3)

    def test_exact_fields_and_expected_birth(self):
        owner = self.fake()
        provider = NativeRootProvider(owner, self.resources())
        self.addCleanup(provider.close)

        sample = provider.sample({os.getpid(): None})

        self.assertEqual(
            sample, [ProcessSample(os.getpid(), 7, 10, 3.23456789)]
        )
        self.assertEqual(provider.sample({os.getpid(): 8}), [])
        self.assertEqual(provider.sample({os.getpid(): 7}), sample)

    def test_native_ticks_use_the_original_provider_timebase(self):
        original = self.resources()
        original._cpu_numer = 125
        original._cpu_denom = 3_000_000_000
        provider = NativeRootProvider(self.fake(), original)
        self.addCleanup(provider.close)

        sample = provider.sample({os.getpid(): None})[0]

        self.assertAlmostEqual(sample.cpu, 3.23456789 * 125 / 3)

    def test_child_and_multiple_root_paths_use_original_tree_provider(self):
        for status, roots in (
            (1, {os.getpid(): None}),
            (0, {os.getpid(): None, 123456: None}),
        ):
            owner = self.fake(status)
            provider = NativeRootProvider(owner, self.resources())
            self.addCleanup(provider.close)
            fallback = SimpleNamespace(sample=Mock(return_value=["tree"]))
            provider.fallback = fallback

            self.assertEqual(provider.sample(roots), ["tree"])

            fallback.sample.assert_called_once_with(roots)

    def test_missing_process_and_observation_failure_remain_distinct(self):
        for status in (-errno.ESRCH, -errno.EPERM, -errno.EACCES):
            provider = NativeRootProvider(self.fake(status), self.resources())
            self.addCleanup(provider.close)
            if status == -errno.ESRCH:
                self.assertEqual(provider.sample({os.getpid(): None}), [])
            else:
                with self.assertRaises(ProcessObservationDeniedError):
                    provider.sample({os.getpid(): None})

    def test_rss_refusal_restores_original_provider(self):
        original = self.resources()
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

    def block_store(self, rss=10):
        packet = self.packet()
        owner = self.fake(rss=rss)
        store = BlockStore()
        store.bind(
            SimpleNamespace(
                ffi=packet.ffi,
                lib=SimpleNamespace(
                    b4_observe_root=owner.lib.b4_observe_root,
                    b4_block_create=packet.lib.b4_block_create,
                    b4_block_destroy=packet.lib.b4_block_destroy,
                    b4_block_prepare=packet.lib.b4_block_prepare,
                ),
            )
        )
        self.addCleanup(store.close)
        return store

    def test_registration_reuses_observer_and_releases_it_with_store(self):
        original = self.resources()
        original.sample = Mock(
            return_value=[ProcessSample(os.getpid(), 7, 10, 0)]
        )
        budget = Budget(
            provider=original,
            rss_bytes=10000,
            wall_clock=lambda: 0,
            cpu_clock=lambda: 0,
        )
        store = self.block_store()

        first = store.register(
            "first", (3, 5, 7), 3, 11, 31, None, budget, b1=31, b2=997
        )
        provider = store._resource_provider
        cached = store.register(
            "first", (), 3, 11, 31, None, budget, b1=31, b2=997
        )
        store.register(
            "second", (11, 13), 11, 17, 31, None, budget, b1=31, b2=997
        )

        self.assertIs(first, cached)
        self.assertIs(store._resource_provider, provider)
        self.assertEqual(store.lib.b4_observe_root.call_count, 2)
        self.assertIs(budget.provider, original)
        self.assertIsNone(provider.fallback)
        store.close()
        self.assertIsNone(store._resource_provider)
        self.assertEqual(store.blocks, {})
        store.close()

    def test_registration_refusal_keeps_checked_block_and_restores_provider(
        self,
    ):
        original = self.resources()
        original.sample = Mock(
            return_value=[ProcessSample(os.getpid(), 7, 10, 0)]
        )
        budget = Budget(
            provider=original,
            rss_bytes=100,
            wall_clock=lambda: 0,
            cpu_clock=lambda: 0,
        )
        store = self.block_store(rss=100)

        with self.assertRaisesRegex(BudgetExhaustedError, "rss_limit"):
            store.register(
                "first", (3, 5), 3, 7, 31, None, budget, b1=31, b2=997
            )

        self.assertIn("first", store.blocks)
        self.assertEqual(budget.observed_rss, 100)
        self.assertIs(budget.provider, original)
        self.assertIsNone(store._resource_provider.fallback)
        store.close()
        self.assertIsNone(store._resource_provider)

    def test_actual_os_fields_match_original_usage_units(self):
        if os.environ.get("FACTOR_RESOURCE_TEST_CHILD") != "1":
            # Full-suite multiprocessing helpers may remain live. Prove the
            # childless fast path in its own interpreter rather than weakening
            # that assertion or rejecting the correct full-tree fallback.
            test = (
                "v2.tests.test_native_resources.NativeResourceTests."
                "test_actual_os_fields_match_original_usage_units"
            )
            script = "\n".join(
                (
                    "from pathlib import Path",
                    "import unittest",
                    "from v2.ecm import native_packet",
                    "native_packet.SOURCE = Path(" + repr(str(SOURCE)) + ")",
                    "suite = unittest.defaultTestLoader.loadTestsFromName("
                    + repr(test)
                    + ")",
                    "result = unittest.TextTestRunner().run(suite)",
                    "assert result.testsRun == 1 and not result.skipped",
                    "raise SystemExit(not result.wasSuccessful())",
                )
            )
            result = subprocess.run(
                [sys.executable, "-B", "-c", script],
                cwd=Path(__file__).resolve().parents[2],
                env=dict(os.environ, FACTOR_RESOURCE_TEST_CHILD="1"),
                capture_output=True,
                text=True,
                timeout=10,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            return

        packet = self.packet()
        original = MacResources()
        provider = NativeRootProvider(packet, original)
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
        original = MacResources()
        provider = NativeRootProvider(packet, original)
        self.addCleanup(provider.close)
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
