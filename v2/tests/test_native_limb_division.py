"""Signed/aliased modular products at every experimental width transition."""

import random
import unittest

from v2.ecm.native_packet import SOURCE, NativePacket
from v2.tests.test_native_correctness import deterministic_budget


@unittest.skipUnless(
    (SOURCE / "build/manifest.json").exists(), "scalar native build required"
)
class NativeLimbDivisionTests(unittest.TestCase):
    def test_exact_signed_aliases_boundaries_and_fallback(self):
        rng = random.Random(7311)
        widths = (191, 192, 193, 255, 256, 257, 319, 320, 321, 383, 384, 385)
        moduli = [3, (1 << 1024) - 109]
        for bits in widths:
            moduli.extend(((1 << bits) - 1, (1 << (bits - 1)) + 3))
        for n in moduli:
            packet = NativePacket()
            try:
                packet.prepare_ladder(1, (2, 1), n, 2, deterministic_budget())
                ffi, lib = packet.ffi, packet.lib
                enabled = bool(lib.b4_limb_division())
                self.assertEqual(enabled, packet.manifest["limb_division"])
                values = [0, 1, n - 1, n, n + 1, 2 * n, -1, 1 - n]
                values.extend(rng.randrange(-2 * n, 2 * n) for _ in range(12))
                for a in values:
                    for b in (a, rng.choice(values)):
                        for alias in (0, 1, 2, 3) if a == b else (0, 1, 2):
                            records, owners = [], []
                            for value in (a, b):
                                raw = abs(value).to_bytes(
                                    max(1, (abs(value).bit_length() + 7) // 8),
                                    "big",
                                )
                                data = ffi.new("unsigned char[]", raw)
                                record = ffi.new("b4_value *")
                                record.data, record.size = data, len(raw)
                                records.append(record)
                                owners.append(data)
                            output = ffi.new("b4_atom *")
                            try:
                                path = lib.b4_test_product(
                                    packet.context,
                                    records[0][0],
                                    records[1][0],
                                    int(a < 0),
                                    int(b < 0),
                                    alias,
                                    output,
                                )
                                width = (n.bit_length() + 63) // 64
                                fast = (
                                    enabled
                                    and 4 <= width <= 6
                                    and abs(a).bit_length() <= 64 * width
                                    and abs(b).bit_length() <= 64 * width
                                )
                                self.assertEqual(path, int(fast))
                                actual = int.from_bytes(
                                    bytes(
                                        ffi.buffer(output.data, output.size)
                                    ),
                                    "big",
                                )
                                self.assertEqual(actual, a * b % n)
                            finally:
                                for item in (output, *records, *owners):
                                    ffi.release(item)
            finally:
                packet.end()
                packet.close()
