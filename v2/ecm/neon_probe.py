"""Explicit ECM-local feasibility adapter; never used by the ECM executor."""

import hashlib
import json
import sys
from pathlib import Path

from ..common.arithmetic import BackendUnavailableError

SOURCE = Path(__file__).with_name("native")
ARMS = (
    "arithmetic-gmp-mpn",
    "arithmetic-scalar26",
    "arithmetic-neon26",
    "complete-production-mpz",
    "complete-hybrid-scalar26",
    "complete-hybrid-neon26",
    "layout-only-diagnostic",
)
WIDTHS = (200, 266, 333)
ENTRIES = 64


def file_sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(65536):
            digest.update(block)
    return digest.hexdigest()


class NEONProbe:
    """One private context with preallocated buffers and no worker threads."""

    def __init__(self, width, modulus):
        if width not in WIDTHS or type(width) is not int:
            raise ValueError("prototype widths are exactly 200/266/333")
        if (
            type(modulus) is not int
            or modulus.bit_length() != width
            or not modulus % 2
        ):
            raise ValueError("prototype requires an odd full-width modulus")
        if sys.implementation.name != "pypy" or sys.version_info[:2] != (
            3,
            11,
        ):
            raise RuntimeError("NEON prototype requires PyPy Python 3.11")
        self.ffi = self.lib = self.context = None
        self.buffers = []
        self.closed = False
        self.width, self.modulus = width, modulus
        directory = SOURCE / "build/neon"
        try:
            if (directory / "manifest.json").stat().st_size > 32768:
                raise ValueError("oversized NEON build manifest")
            self.manifest = json.loads(
                (directory / "manifest.json").read_text()
            )
            required = {
                "neon/" + name
                for name in (
                    "kernels.c",
                    "kernels.h",
                    "screen.c",
                    "screen.h",
                    "sink.c",
                )
            }
            required.update(("neon_build.py", "packet.c"))
            if set(self.manifest["source"]) != required:
                raise ValueError("incomplete NEON source identity")
            for name, pin in self.manifest["source"].items():
                path = (SOURCE / name).resolve()
                if (
                    not path.is_relative_to(SOURCE.resolve())
                    or file_sha(path) != pin
                ):
                    raise ValueError("NEON build source mismatch")
            library = directory / "libfactor_ecm_neon.dylib"
            dependency = Path(self.manifest["gmp_library"])
            if (
                library.stat().st_size > 512 * 1024
                or dependency.stat().st_size > 4 * 2**20
            ):
                raise ValueError("NEON dependency exceeds load bound")
            if (
                file_sha(library) != self.manifest["library_sha256"]
                or file_sha(dependency) != self.manifest["gmp_library_sha256"]
            ):
                raise ValueError("NEON dependency/build identity mismatch")
            import cffi

            self.ffi = cffi.FFI()
            self.ffi.cdef(
                "\n".join(
                    line
                    for line in (SOURCE / "neon/screen.h")
                    .read_text()
                    .splitlines()
                    if not line.startswith("#")
                )
            )
            self.lib = self.ffi.dlopen(str(library))
            if (
                self.lib.b4n_abi() != 1
                or self.manifest["abi"] != 1
                or self.ffi.string(self.lib.b4n_gmp_version()).decode()
                != self.manifest["gmp_runtime"]
            ):
                raise ValueError("NEON ABI/GMP identity mismatch")
            for declaration in (
                "unsigned char[168]",
                "b4n_value[4]",
                "b4n_output *",
                "b4n_observation *",
            ):
                self.buffers.append(self.ffi.new(declaration))
            self.transport, self.values, self.output, self.observation = (
                self.buffers
            )
            for i in range(4):
                self.values[i].data = self.transport + 42 * i
            self._values([modulus])
            self.context = self.lib.b4n_create(width, self.values[0])
            if self.context == self.ffi.NULL:
                self.context = None
                raise MemoryError("native NEON context allocation refused")
        except BaseException as error:
            self.close()
            if isinstance(error, (OSError, ImportError)):
                raise BackendUnavailableError(
                    "explicit NEON prototype build required"
                ) from error
            raise

    def _values(self, values):
        for i, value in enumerate(values):
            if type(value) is not int or not 0 <= value < 2**self.width:
                raise ValueError("unsigned prototype operand exceeds width")
            raw = value.to_bytes(max(1, (value.bit_length() + 7) // 8), "big")
            self.ffi.memmove(self.values[i].data, raw, len(raw))
            self.values[i].size = len(raw)

    def set(self, index, values, *, active=2):
        if self.closed:
            raise RuntimeError("NEON probe closed")
        if type(index) is not int or not 0 <= index < ENTRIES:
            raise ValueError("prototype entry exceeds fixed panel")
        if (
            type(active) is not int
            or active not in (1, 2)
            or len(values) != 2 * active
        ):
            raise ValueError(
                "prototype requires one or two independent operations"
            )
        self._values(values)
        if self.lib.b4n_set(self.context, index, self.values, active):
            raise ValueError("prototype input refused")

    def once(self, arm, operation, index=0):
        if (
            self.closed
            or arm not in ARMS
            or operation not in ("multiply", "square")
        ):
            raise ValueError("invalid prototype operation")
        if type(index) is not int or not 0 <= index < ENTRIES:
            raise ValueError("prototype entry exceeds fixed panel")
        if self.lib.b4n_once(
            self.context, ARMS.index(arm), operation == "square", index
        ):
            raise ValueError(
                "prototype result disagrees with GMP or was refused"
            )
        if self.lib.b4n_export(self.context, self.output):
            raise ValueError("prototype canonical export refused")
        return tuple(
            int.from_bytes(
                self.ffi.buffer(self.output.bytes[i], self.output.size[i])[:],
                "big",
            )
            for i in range(2)
        )

    def qos(self):
        self.lib.b4n_qos(self.observation)
        return {
            key: int(getattr(self.observation, key))
            for key in (
                "qos_set",
                "qos_get",
                "requested_qos",
                "observed_qos",
                "relative_priority",
                "thread_id",
                "thread_id_status",
            )
        }

    def batch(self, arm, operation, repeats, offset):
        """Stage A only: time C work, then validate results outside timing."""
        if (
            self.closed
            or arm not in ARMS
            or operation not in ("multiply", "square")
        ):
            raise ValueError("invalid prototype arm/operation")
        if (
            type(repeats) is not int
            or not 1 <= repeats <= 65536
            or type(offset) is not int
            or not 0 <= offset < ENTRIES
        ):
            raise ValueError("prototype batch exceeds finite bounds")
        status = self.lib.b4n_batch(
            self.context,
            ARMS.index(arm),
            operation == "square",
            repeats,
            offset,
            self.observation,
        )
        if status < 0:
            raise ValueError("prototype timing preflight refused")
        validated = self.lib.b4n_validate(self.context, self.observation) == 0
        row = {
            key: int(getattr(self.observation, key))
            for key in (
                "checksum",
                "pairs",
                "qos_set",
                "qos_get",
                "requested_qos",
                "observed_qos",
                "relative_priority",
                "thread_id",
                "thread_id_status",
            )
        }
        row.update(
            wall=float(self.observation.wall),
            cpu=float(self.observation.cpu),
            validated=validated,
            status=status,
            repeats=repeats,
            offset=offset,
            independent_operations=(
                0 if arm == ARMS[6] else 2 * int(self.observation.pairs)
            ),
        )
        return row

    def close(self):
        if self.context is not None and self.ffi is not None:
            self.lib.b4n_destroy(self.context)
        self.context = None
        if self.ffi is not None:
            for owner in reversed(self.buffers):
                self.ffi.release(owner)
        self.buffers.clear()
        if self.lib is not None:
            self.ffi.dlclose(self.lib)
        self.lib = self.ffi = None
        self.closed = True
