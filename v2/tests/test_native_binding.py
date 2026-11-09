"""Build-time ABI descriptors preserve exact native identity admission."""

import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from v2.ecm import native_packet
from v2.ecm.native_packet import SOURCE, NativePacket
from v2.tests.test_native_correctness import deterministic_budget


@unittest.skipUnless(
    (SOURCE / "build/manifest.json").exists(), "scalar native build required"
)
class NativeBindingTests(unittest.TestCase):
    def test_first_load_does_not_construct_inline_cffi_parser(self):
        import cffi

        packet = NativePacket()
        self.addCleanup(packet.close)

        with patch.object(cffi, "FFI", side_effect=AssertionError("parser")):
            self.assertTrue(packet.available(1013, deterministic_budget()))

        self.assertIn("abi_binding", packet.manifest)

    def test_changed_or_missing_declared_descriptor_refuses_in_auto_mode(self):
        for missing in (False, True):
            with tempfile.TemporaryDirectory() as directory:
                target = Path(directory)
                for name in ("packet.c", "packet.h", "team.c", "curve.c"):
                    shutil.copy2(SOURCE / name, target / name)
                build = target / "build"
                build.mkdir()
                manifest = json.loads(
                    (SOURCE / "build/manifest.json").read_text()
                )
                for name in (
                    "manifest.json",
                    "libfactor_ecm.dylib",
                    manifest["abi_binding"]["file"],
                ):
                    shutil.copy2(SOURCE / "build" / name, build / name)
                descriptor = build / manifest["abi_binding"]["file"]
                if missing:
                    descriptor.unlink()
                else:
                    descriptor.write_bytes(
                        descriptor.read_bytes() + b"\n# changed"
                    )
                packet = NativePacket(mode="auto")
                try:
                    with patch.object(native_packet, "SOURCE", target):
                        with self.assertRaisesRegex(ValueError, "binding"):
                            packet.available(1013, deterministic_budget())
                    self.assertFalse(packet.unavailable)
                finally:
                    packet.close()
