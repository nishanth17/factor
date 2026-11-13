"""Deterministic native deadline sampling, including slow arithmetic gaps."""

import unittest

from v2.ecm import core
from v2.ecm.native_packet import SOURCE, NativePacket
from v2.tests.test_native_correctness import deterministic_budget


@unittest.skipUnless(
    (SOURCE / "build/manifest.json").exists(), "scalar native build required"
)
class NativePollingTests(unittest.TestCase):
    def packet(self, walls, cpus):
        packet = NativePacket()
        self.addCleanup(packet.close)
        packet.prepare_ladder(255, (2, 1), 1013, 2, deterministic_budget())
        self.addCleanup(packet.end)
        ffi, lib = packet.ffi, packet.lib
        wall = ffi.new("unsigned long long[]", walls)
        cpu = ffi.new("unsigned long long[]", cpus)
        self.assertEqual(
            lib.b4_test_clocks(
                packet.context, wall, len(walls), cpu, len(cpus)
            ),
            0,
        )
        return packet

    def step(self, packet, wall=10**9, cpu=10**9, actions=128):
        return packet.lib.b4_slice(
            packet.context, wall, cpu, actions, packet.ffi.NULL
        )

    def reads(self, packet):
        return tuple(
            packet.lib.b4_test_clock_reads(packet.context, kind)
            for kind in (0, 1)
        )

    def test_cpu_reads_amortize_and_canonical_arithmetic_is_unchanged(self):
        packet = self.packet([100], [10])
        self.assertEqual(self.step(packet), 1)
        self.assertEqual(self.reads(packet), (7, 1))
        result = packet.checked(packet._export(), 1013)
        self.assertEqual(result, core.scalar_multiply(255, 2, 1, 1013, 2))

    def test_slow_operation_checks_wall_at_the_very_next_boundary(self):
        packet = self.packet([100, 10**9], [10])
        self.assertEqual(self.step(packet), 0)
        self.assertEqual(self.reads(packet), (2, 1))
        self.assertEqual(self.step(packet), 0)
        self.assertEqual(self.reads(packet), (3, 1))

    def test_cpu_sample_after_one_ms_and_again_on_the_next_slice(self):
        packet = self.packet([100, 101, 1000100], [10, 20])
        self.assertEqual(self.step(packet, cpu=20), 0)
        self.assertEqual(self.reads(packet), (3, 2))
        self.assertEqual(self.step(packet, cpu=20), 0)
        self.assertEqual(self.reads(packet), (4, 3))

    def test_wall_cpu_and_clock_failures_precede_arithmetic(self):
        for walls, cpus, wall_end, cpu_end, expected in (
            ([100], [10], 100, 20, 0),
            ([100], [20], 200, 20, 0),
            ([0], [10], 200, 20, -1),
            ([100], [0], 200, 20, -1),
        ):
            packet = self.packet(walls, cpus)
            self.assertEqual(self.step(packet, wall_end, cpu_end), expected)
            self.assertNotEqual(
                packet.lib.b4_result(packet.context, packet.output), 0
            )

    def test_action_cap_and_backward_clock_force_new_cpu_sample(self):
        packet = self.packet([100, 99], [10])
        self.assertEqual(self.step(packet, actions=2), 0)
        self.assertEqual(self.reads(packet), (2, 2))
        self.assertEqual(self.step(packet), 1)

    def test_test_hooks_reject_invalid_script_sizes(self):
        packet = self.packet([100], [10])
        self.assertEqual(
            packet.lib.b4_test_clocks(
                packet.context, packet.ffi.NULL, 17, packet.ffi.NULL, 0
            ),
            -1,
        )
