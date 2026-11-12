"""Native prime coverage, certificate proofs and checked cache fallback."""

import hashlib
import json
import os
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from struct import pack

from v2.ecm import chain_data
from v2.ecm.native_curve import NativeCurvePacket
from v2.ecm.native_metadata import MAGIC, NativeMetadata
from v2.ecm.native_packet import SOURCE
from v2.execution import stage_jobs
from v2.execution.budget import Budget
from v2.execution.portfolio_config import (
    _execution_context,
    close_execution_context,
)
from v2.execution.schedules import SieveContext
from v2.tests.test_native_correctness import (
    FIXTURES,
    deterministic_budget,
    selected,
)


@unittest.skipUnless(
    (SOURCE / "build/manifest.json").exists(), "native build required"
)
class NativeMetadataTests(unittest.TestCase):
    def setUp(self):
        self.budget = deterministic_budget()
        self.packet = NativeCurvePacket()
        self.assertTrue(self.packet.available(1013, self.budget))
        self.metadata = NativeMetadata(self.packet.ffi, self.packet.lib)

    def tearDown(self):
        self.metadata.close()
        self.packet.close()
        self.budget.close()

    def values(self, metadata, lo, hi, bound=None):
        block, values = metadata.block(
            lo,
            hi,
            bound,
            None,
            self.budget,
            b1=2 if bound is None else bound,
            b2=hi - 1,
        )
        try:
            return list(values), metadata.encode(block)
        finally:
            self.assertEqual(self.packet.lib.b4_block_destroy(block), 0)

    def test_bulk_half_open_coverage_powers_and_word_edge(self):
        for lo, hi in (
            (2, 3),
            (3, 4),
            (4, 5),
            (49, 53),
            (100, 1001),
            (2**32 - 512, 2**32),
        ):
            with self.subTest(lo=lo, hi=hi):
                context = SieveContext(hi, segment_size=1024)
                primes, data = self.values(self.metadata, lo, hi, hi - 1)
                self.assertEqual(primes, context.prime_segment(lo, hi))
                self.assertEqual(len(data), 16 * len(primes))
                for index, prime in enumerate(primes):
                    power = prime
                    while power <= (hi - 1) // prime:
                        power *= prime
                    self.assertEqual(
                        int.from_bytes(
                            data[16 * index + 8 : 16 * index + 16], "little"
                        ),
                        power,
                    )

    def test_independent_catalogue_frontier_rejects_lost_guards(self):
        ffi, lib = self.packet.ffi, self.packet.lib
        scalar, code, masks, output, slots = next(
            row for row in chain_data.PRAC if row[2] != bytes(len(row[2]))
        )

        def create(certificate):
            body = (
                pack("<I", 1)
                + pack(
                    "<QHHBB",
                    scalar,
                    len(code),
                    len(certificate),
                    slots,
                    output,
                )
                + code
                + certificate
            )
            buffer = ffi.from_buffer("const unsigned char[]", body)
            try:
                return lib.b4_catalogue_create(buffer, len(body))
            finally:
                ffi.release(buffer)

        catalogue = create(masks)
        self.assertNotEqual(catalogue, ffi.NULL)
        lib.b4_catalogue_destroy(catalogue)
        # All-zero masks lose discarded intermediate coordinates even though
        # the unchanged bytecode still has a valid scalar certificate.
        self.assertEqual(create(bytes(len(masks))), ffi.NULL)

    def test_real_root_bulk_matches_python_canonical_jobs(self):
        cases = json.loads(FIXTURES.read_text())["cases"][:3]
        corpus = FIXTURES.parents[1] / "corpora/b4a_production_3p_v1.json"
        case = next(
            row
            for row in json.loads(corpus.read_text())["fixtures"]
            if row["population"] == "failed" and row["digits"] == 60
        )
        cases.append(dict(n=case["n"], seed=73, b1=2000, b2=81000))

        def run(case, kernel, *, checkpoint=None, pause=False):
            budget = Budget(seconds=10, cpu_seconds=20, rss_bytes=8 * 2**30)
            config = selected(kernel, case["b1"], case["b2"])
            if case["b1"] == 2000:
                config = replace(
                    config,
                    ecm_chain_mode="auto",
                    ecm_tiers=((2000, 81000, 9),),
                )
            context, programs = _execution_context({}, config, budget)
            source = programs or context
            packet = getattr(source, "native_packet", None)
            job = checkpoint or stage_jobs.new_job(
                "ecm", case["n"], case["seed"], case["b1"], case["b2"]
            )
            try:
                for _ in range(100000):
                    stage_jobs.advance_job(job, budget, source, config)
                    if job["done"]:
                        break
                    if pause and packet._read().next > 100:
                        break
                else:
                    self.fail("finite bulk job failed to finish")
                if packet is not None:
                    packet.sync(allow_result=True)
                    if case["b1"] == 2000:
                        self.assertIsNotNone(packet.blocks.metadata)
                return job
            finally:
                close_execution_context(source)
                budget.close()

        for case in cases:
            with self.subTest(n=case["n"], b1=case["b1"]):
                self.assertEqual(run(case, "c-gmp-full"), run(case, "python"))
                if case["b1"] == 2000:
                    paused = run(case, "c-gmp-full", pause=True)
                    self.assertFalse(paused["done"])
                    self.assertEqual(
                        run(case, "c-gmp-full", checkpoint=paused),
                        run(case, "python"),
                    )

    def inventory(self, directory, disposition):
        target = os.environ.get("FACTOR_NATIVE_CACHE_INVENTORY")
        if target:
            with Path(target).open("a") as stream:
                for path in Path(directory).rglob("*"):
                    if path.is_file():
                        stream.write(
                            json.dumps(
                                dict(
                                    path=str(path),
                                    size=path.stat().st_size,
                                    sha256=hashlib.sha256(
                                        path.read_bytes()
                                    ).hexdigest(),
                                    disposition=disposition,
                                )
                            )
                            + "\n"
                        )

    def test_cache_hit_truncation_incompatibility_and_forged_hash_fallback(
        self,
    ):
        parent = Path(
            os.environ.get(
                "FACTOR_NATIVE_TEST_CACHE", SOURCE / "build/cache-tests"
            )
        )
        parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=parent) as directory:
            metadata = NativeMetadata(
                self.packet.ffi, self.packet.lib, cache_directory=directory
            )
            try:
                expected = self.values(metadata, 3, 101)[0]
                self.assertEqual(self.values(metadata, 3, 101)[0], expected)
                self.assertGreater(metadata.stats["hits"], 0)
                identity = metadata._identity(
                    "prime-program-block",
                    lo=3,
                    hi=101,
                    bound=None,
                    b1=2,
                    b2=100,
                    program=None,
                )
                path = metadata._path(identity)
                original = path.read_bytes()
                for corruption in (
                    "truncated",
                    "incompatible",
                    "forged-coverage",
                ):
                    self.inventory(directory, "overwritten-for-" + corruption)
                    if corruption == "truncated":
                        raw = original[:9]
                    else:
                        offset = len(MAGIC)
                        length = int.from_bytes(
                            original[offset : offset + 4], "little"
                        )
                        header = json.loads(
                            original[offset + 4 : offset + 4 + length]
                        )
                        body = original[offset + 4 + length :]
                        if corruption == "incompatible":
                            header["identity"]["integer_width"] = 32
                        else:
                            body = body[16:]
                            header["bytes"] = len(body)
                            header["sha256"] = hashlib.sha256(body).hexdigest()
                        encoded = json.dumps(header).encode()
                        raw = MAGIC + pack("<I", len(encoded)) + encoded + body
                    path.write_bytes(raw)
                    self.assertEqual(
                        self.values(metadata, 3, 101)[0], expected
                    )
                self.assertGreater(metadata.stats["corrupt"], 0)
                self.assertEqual(metadata.stats["incompatible"], 1)
                self.assertEqual(metadata.stats["proof_failures"], 1)
                self.assertEqual(metadata.stats["mapped_bytes"], 0)
            finally:
                metadata.close()
                self.inventory(directory, "removed-by-test-cleanup")
