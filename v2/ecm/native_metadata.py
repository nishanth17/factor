"""Exact bulk preparation and optional data-only native metadata artifacts."""

import hashlib
import json
import os
import stat
import sys
import tempfile
from pathlib import Path
from struct import iter_unpack, pack

from . import chains
from .native_packet import source_identity

MAGIC = b"B4META1\n"
MAX_ARTIFACT = 2 * 2**20


def supported_context(context):
    """Keep arbitrary schedule providers on their original exact path."""
    from ..execution.schedules import ScheduleCache, SieveContext
    from .programs import ECMPrograms
    from .worker_state import BoundContext, DecodedPrograms

    while type(context) in (
        ECMPrograms,
        ScheduleCache,
        BoundContext,
        DecodedPrograms,
    ):
        if type(context) is BoundContext:
            context = context.source
        elif type(context) is DecodedPrograms:
            context = context.programs
        else:
            context = context.context
    return type(context) is SieveContext


class NativePlan:
    """The existing accepted plan identity, with proofs retained in C."""

    def __init__(self, identity, catalogue, unit_bytes):
        self.identity, self.catalogue = identity, catalogue
        self.unit_bytes = unit_bytes


class NativeValues:
    """Decode the immutable native index only for a canonical checkpoint."""

    def __init__(self, metadata, block):
        self.metadata, self.block = metadata, block
        self.count = int(metadata.lib.b4_block_count(block))

    def __len__(self):
        return self.count

    def __iter__(self):
        data = self.metadata.encode(self.block)
        return (prime for prime, _ in iter_unpack("<QQ", data))


class NativeMetadata:
    """One owner's verified catalogue and explicitly requested disk cache.

    Cache bytes contain no executable objects or native pointers. Every load
    binds format/source/program/bounds, checks its digest, then proves coverage
    and program semantics independently. Native blocks own copied data; the
    read buffer is released immediately. There is no mapping to outlive.
    """

    def __init__(self, ffi, lib, *, cache_directory=None):
        self.ffi, self.lib = ffi, lib
        self.catalogues = {}
        self.plans = {}
        self.max_record_bytes = 0
        self.used = 0
        self.cache_directory = (
            None
            if cache_directory is None
            else Path(cache_directory).resolve()
        )
        self.stats = dict(
            hits=0,
            misses=0,
            corrupt=0,
            incompatible=0,
            io_errors=0,
            proof_failures=0,
            written_bytes=0,
            read_bytes=0,
            mapped_bytes=0,
        )
        self.provenance = source_identity()

    def close(self):
        for catalogue in self.catalogues.values():
            self.lib.b4_catalogue_destroy(catalogue)
        self.catalogues.clear()
        self.plans.clear()
        self.used = 0

    def _identity(self, kind, **fields):
        return dict(
            format="factor-native-metadata-v1",
            kind=kind,
            integer_width=64,
            byte_order="little",
            host_byte_order=sys.byteorder,
            pointer_width=8 * self.ffi.sizeof("void *"),
            coverage="segmented-eratosthenes-half-open-v1",
            certificate="scalar-forward-ancestor-frontier-v1",
            provenance=self.provenance,
            **fields,
        )

    def _path(self, identity):
        key = hashlib.sha256(
            json.dumps(identity, sort_keys=True).encode()
        ).hexdigest()
        return self.cache_directory / (key + ".b4m")

    def _read(self, identity, *, digest=None):
        if self.cache_directory is None:
            return None
        try:
            fd = os.open(self._path(identity), os.O_RDONLY | os.O_NOFOLLOW)
            with os.fdopen(fd, "rb") as stream:
                info = os.fstat(stream.fileno())
                if (
                    not stat.S_ISREG(info.st_mode)
                    or info.st_size > MAX_ARTIFACT
                ):
                    raise ValueError("unsupported cache artifact geometry")
                raw = stream.read(MAX_ARTIFACT + 1)
            self.stats["read_bytes"] += len(raw)
        except FileNotFoundError:
            self.stats["misses"] += 1
            return None
        except OSError:
            self.stats["io_errors"] += 1
            return None
        try:
            if len(raw) < len(MAGIC) + 4 or raw[: len(MAGIC)] != MAGIC:
                raise ValueError("invalid cache format")
            offset = len(MAGIC)
            size = int.from_bytes(raw[offset : offset + 4], "little")
            offset += 4
            if size > 16384 or size > len(raw) - offset:
                raise ValueError("invalid cache header geometry")
            header = json.loads(raw[offset : offset + size])
            body = raw[offset + size :]
            if set(header) != {"identity", "sha256", "bytes"}:
                raise ValueError("invalid cache header fields")
            if header["identity"] != identity:
                self.stats["incompatible"] += 1
                return None
            actual = hashlib.sha256(body).hexdigest()
            if header["bytes"] != len(body) or header["sha256"] != actual:
                raise ValueError("corrupt cache payload")
            if digest is not None and actual != digest:
                raise ValueError("catalogue disagrees with immutable source")
            self.stats["hits"] += 1
            return body
        except (ValueError, TypeError, UnicodeError):
            self.stats["corrupt"] += 1
            return None

    def _write(self, identity, body):
        if self.cache_directory is None:
            return
        header = json.dumps(
            dict(
                identity=identity,
                bytes=len(body),
                sha256=hashlib.sha256(body).hexdigest(),
            ),
            sort_keys=True,
        ).encode()
        raw = MAGIC + pack("<I", len(header)) + header + body
        if len(raw) > MAX_ARTIFACT:
            raise ValueError("cache artifact exceeds finite codec geometry")
        temporary = None
        try:
            self.cache_directory.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(
                dir=self.cache_directory, prefix=".pending-", delete=False
            ) as stream:
                temporary = Path(stream.name)
                stream.write(raw)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, self._path(identity))
            self.stats["written_bytes"] += len(raw)
        except OSError:
            self.stats["io_errors"] += 1
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    def _catalogue(self, chain, budget):
        if chain in self.catalogues:
            return self.catalogues[chain]
        if chain:
            from . import chain_data

            rows = sorted(chain_data.PRAC)
            if len(rows) != 333:
                raise ValueError("incomplete accepted native catalogue")
            self.max_record_bytes = max(
                len(code) + len(masks) for _, code, masks, _, _ in rows
            )
            budget.consume(
                chains.CATALOG_PREPARATION_WORK
                + sum(len(code) + slots + 1 for _, code, _, _, slots in rows)
            )
            parts = [pack("<I", len(rows))]
            for scalar, code, masks, output, slots in rows:
                parts.extend(
                    (
                        pack(
                            "<QHHBB",
                            scalar,
                            len(code),
                            len(masks),
                            slots,
                            output,
                        ),
                        code,
                        masks,
                    )
                )
            body = b"".join(parts)
        else:
            body = pack("<I", 0)
        budget.consume(4096 + 4 * len(body))
        identity = self._identity(
            "catalogue", program=chains.CHAIN_VERSION if chain else None
        )
        cached = self._read(identity, digest=hashlib.sha256(body).hexdigest())
        if cached is not None:
            body = cached
        data = self.ffi.from_buffer("const unsigned char[]", body)
        try:
            catalogue = self.lib.b4_catalogue_create(data, len(body))
        finally:
            self.ffi.release(data)
        if catalogue == self.ffi.NULL:
            raise ValueError("native catalogue proof refused")
        self.catalogues[chain] = catalogue
        self.used += 4096 + 4 * len(body)
        if cached is None:
            self._write(identity, body)
        return catalogue

    def plan(self, job, context, budget):
        source = getattr(context, "chains", None)
        if (
            type(source) is not chains.ChainPlans
            or source.backend != "python-int"
            or job["b1"] not in source.bounds
            or not chains.supports_modulus(job["n"])
        ):
            return None
        identity = chains.identity(job["b1"], source.backend)
        budget.consume(1)
        if identity not in self.plans:
            catalogue = self._catalogue(True, budget)
            budget.consume(job["b1"])
            self.plans[identity] = NativePlan(
                identity, catalogue, 128 + 2 * self.max_record_bytes
            )
        return self.plans[identity]

    def encode(self, block):
        size = 16 * int(self.lib.b4_block_count(block))
        data = self.ffi.new("unsigned char[]", max(1, size))
        try:
            if self.lib.b4_block_encode(block, data, size):
                raise ValueError("native metadata encoding refused")
            return bytes(self.ffi.buffer(data, size))
        finally:
            self.ffi.release(data)

    def block(self, lo, hi, bound, plan, budget, *, b1, b2):
        catalogue = (
            plan.catalogue
            if plan is not None
            else self._catalogue(False, budget)
        )
        identity = self._identity(
            "prime-program-block",
            lo=lo,
            hi=hi,
            bound=bound,
            b1=b1,
            b2=b2,
            program=None if plan is None else plan.identity,
        )
        cached = self._read(identity)
        block = self.ffi.NULL
        if cached is not None:
            data = self.ffi.from_buffer("const unsigned char[]", cached)
            try:
                block = self.lib.b4_block_decode(
                    catalogue,
                    lo,
                    hi,
                    bound or 0,
                    int(plan is not None),
                    data,
                    len(cached),
                )
            finally:
                self.ffi.release(data)
            if block == self.ffi.NULL:
                self.stats["proof_failures"] += 1
                cached = None
        if block == self.ffi.NULL:
            block = self.lib.b4_block_generate(
                catalogue, lo, hi, bound or 0, int(plan is not None)
            )
        if block == self.ffi.NULL:
            raise ValueError("native exact bulk metadata refused")
        try:
            if not bound and self.lib.b4_block_prepare(block, b1, b2):
                raise ValueError("native bulk stage-two indices refused")
            if cached is None and self.cache_directory is not None:
                self._write(identity, self.encode(block))
            return block, NativeValues(self, block)
        except BaseException:
            self.lib.b4_block_destroy(block)
            raise
