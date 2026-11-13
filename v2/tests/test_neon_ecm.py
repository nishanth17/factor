"""Exact, untimed checks for the optional two-curve ECM NEON prototype."""

import os
import random
import re
import unittest
from unittest.mock import patch

from v2.common.arithmetic import BackendUnavailableError
from v2.ecm.neon_probe import ARMS, SOURCE, WIDTHS, NEONProbe


class NEONAdmissionTests(unittest.TestCase):
    def test_only_selected_full_widths_are_supported(self):
        for width in (64, 199, 201, 265, 267, 332, 334, True):
            with self.assertRaises(ValueError):
                NEONProbe(width, 2**200 - 1)
        for modulus in (0, 2**199 - 1, 2**200 - 2, True):
            with self.assertRaises(ValueError):
                NEONProbe(200, modulus)

    def test_missing_build_is_explicit(self):
        with patch("v2.ecm.neon_probe.SOURCE", SOURCE / "absent"):
            with self.assertRaises(BackendUnavailableError):
                NEONProbe(200, 2**200 - 1)


@unittest.skipUnless(
    (SOURCE / "build/neon/manifest.json").exists()
    and os.environ.get("FACTOR_B4A_NEON_QA") == "1",
    "optional ECM NEON prototype is not built",
)
class CompiledNEONTests(unittest.TestCase):
    def test_boundaries_random_products_squares_and_complete_paths(self):
        generator = random.Random(2026110407)
        for width in WIDTHS:
            maximum = 2**width - 1
            probe = NEONProbe(width, maximum)
            try:
                boundaries = [0, 1, 2, maximum, maximum - 1]
                for bits in (26, 52, 64, 130, width - 1):
                    boundaries.extend((2**bits - 1, 2**bits))
                pairs = [(a, maximum - a, maximum - a, a) for a in boundaries]
                pairs.append((maximum, maximum, maximum, maximum))
                pairs += [
                    tuple(generator.getrandbits(width) for _ in range(4))
                    for _ in range(24)
                ]
                # Correctness admission must never enter timed batches/warmup.
                with patch.object(
                    probe, "batch", side_effect=AssertionError("timing in QA")
                ):
                    for values in pairs:
                        probe.set(0, values)
                        for operation in ("multiply", "square"):
                            products = tuple(
                                values[2 * lane]
                                * values[2 * lane + (operation == "multiply")]
                                for lane in range(2)
                            )
                            for arm in ARMS:
                                expected = (
                                    products
                                    if arm.startswith("arithmetic-")
                                    else tuple(
                                        value % maximum for value in products
                                    )
                                )
                                self.assertEqual(
                                    probe.once(arm, operation), expected
                                )
            finally:
                probe.close()
                probe.close()
            self.assertIsNone(probe.context)

    def test_single_lane_tail_ignores_absent_descriptors(self):
        for width in WIDTHS:
            probe = NEONProbe(width, 2**width - 1)
            try:
                values = (2**width - 1, 2**width - 2)
                probe.set(0, values, active=1)
                # Poison the unused descriptors. C must not inspect them.
                for index in (2, 3):
                    probe.values[index].data = probe.ffi.NULL
                    probe.values[index].size = 2**32 - 1
                self.assertEqual(
                    probe.lib.b4n_set(probe.context, 0, probe.values, 1), 0
                )
                for arm in ARMS:
                    for operation in ("multiply", "square"):
                        product = values[0] * (
                            values[1] if operation == "multiply" else values[0]
                        )
                        if not arm.startswith("arithmetic-"):
                            product %= probe.modulus
                        self.assertEqual(
                            probe.once(arm, operation), (product, 0)
                        )
            finally:
                probe.close()

    def test_malformed_canonical_input_invalidates_prior_publication(self):
        probe = NEONProbe(200, 2**200 - 1)
        try:
            probe.set(0, (1, 2, 3, 4))
            probe.once(ARMS[2], "multiply")
            probe.ffi.memmove(probe.values[0].data, b"\0\1", 2)
            probe.values[0].size = 2
            self.assertEqual(
                probe.lib.b4n_set(probe.context, 0, probe.values, 2), -1
            )
            self.assertEqual(
                probe.lib.b4n_export(probe.context, probe.output), -1
            )
            with self.assertRaises(ValueError):
                probe.once(ARMS[2], "multiply")
            probe.set(0, (0, 0, 2**200 - 1, 2**200 - 1))
            self.assertEqual(
                probe.once(ARMS[2], "multiply"), (0, (2**200 - 1) ** 2)
            )
        finally:
            probe.close()

    def test_supported_qos_request_records_success_or_failure_without_pinning(
        self,
    ):
        probe = NEONProbe(200, 2**200 - 1)
        try:
            report = probe.qos()
            self.assertEqual(report["requested_qos"], 0x21)
            self.assertGreaterEqual(report["qos_set"], 0)
            self.assertGreaterEqual(report["qos_get"], 0)
            if not report["thread_id_status"]:
                self.assertGreater(report["thread_id"], 0)
        finally:
            probe.close()

    def test_emitted_kernel_has_vector_products_and_scalar_control_does_not(
        self,
    ):
        text = (SOURCE / "build/neon/disassembly.txt").read_text()

        def body(name):
            match = re.search(
                r"_" + name + r":\n(.*?)(?=\n_[A-Za-z][^:\n]*:\n|\Z)",
                text,
                re.S,
            )
            self.assertIsNotNone(match, name)
            return match.group(1)

        vector_product = r"\b(?:umull|umlal)2?[^\n]*\bv[0-9]+"
        self.assertRegex(body("b4n_neon_product"), vector_product)
        self.assertNotRegex(body("b4n_scalar_product"), vector_product)


if __name__ == "__main__":
    unittest.main()
