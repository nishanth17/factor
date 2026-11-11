"""Lazy, explicitly owned scalar ECM packets using canonical CFFI buffers."""

import hashlib
import json
import sys
import threading
from functools import lru_cache
from pathlib import Path

from ..common.arithmetic import BackendUnavailableError
from ..execution.resources import removed_arguments
from .chain_records import verify_frontier

NATIVE_MAX_WORKERS = 16
ABI = 1
IDENTITY = "b4a-c-gmp-packet-v1"
RESERVE = 2 * 2**20
SOURCE = Path(__file__).with_name("native")
# OpenMP may keep parked runtime threads until process exit. Retain its DSO;
# unloading the packet while a runtime worker can reference it is unsafe.
_OPENMP_BINDINGS = []


def poll_budget(budget):
    """Force a real poll on B12 budgets and preserve policy-wrapper checks."""
    method = getattr(budget, "publication", None)
    if method is None:
        method = getattr(budget, "poll", None)
    if getattr(method, "__self__", None) is budget:
        method()
    else:
        budget.consume(0)


class ProgramStore:
    """One bounded immutable CFFI program store per useful context/team."""

    def __init__(self):
        self.programs = {}
        self.used = 0

    def close(self, ffi):
        for owners in self.programs.values():
            for owner in owners:
                ffi.release(owner)
        self.programs.clear()
        self.used = 0


def _sha_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while raw := stream.read(65536):
            digest.update(raw)
    return digest.hexdigest()


@lru_cache(maxsize=1)
def source_identity():
    """Bind snapshots to the ABI, adapter and compiled source bytes."""
    digest = hashlib.sha256(Path(__file__).read_bytes())
    for name in ("packet.c", "packet.h", "team.c", "curve.c"):
        digest.update((SOURCE / name).read_bytes())
    if (SOURCE / "openmp.c").exists():
        digest.update((SOURCE / "openmp.c").read_bytes())
    for path in (
        Path(__file__).with_name("native_curve.py"),
        Path(__file__).with_name("native_curve_team.py"),
        Path(__file__).with_name("native_resources.py"),
        Path(__file__).parents[1] / "execution/ecm_native_full.py",
    ):
        digest.update(path.read_bytes())
    return IDENTITY + "/" + digest.hexdigest()


def _bytes(value):
    if type(value) is not int or value < 0:
        raise ValueError("native packet requires unsigned Python integers")
    return value.to_bytes(max(1, (value.bit_length() + 7) // 8), "big")


@removed_arguments(("memory_bytes",))
class NativePacket:
    """One non-reentrant adapter per context; no threads or callbacks.

    Native scratch survives slice returns. All program/input owners remain
    live until the packet ends. Publication occurs only after a final poll;
    interrupted native state is disposable and never enters a checkpoint.
    """

    identity = IDENTITY

    def __init__(
        self,
        *,
        mode="c-gmp",
        _variant="scalar",
        _binding=None,
        _programs=None,
    ):
        if mode not in ("c-gmp", "auto"):
            raise ValueError("invalid native packet mode")
        if _variant not in ("scalar", "pthread", "openmp"):
            raise ValueError("invalid native packet variant")
        self.mode = mode
        self.ffi = self.lib = self.context = None
        self.variant = _variant
        self._binding = _binding
        self.store = _programs or ProgramStore()
        self.owns_programs = _programs is None
        self.programs = self.store.programs
        self.buffers = []
        self.input_cache = [None] * 5
        self.output_views = {}
        self.unavailable = False
        self.closed = False
        self._lock = threading.Lock()

    def _load(self, budget):
        if self.closed:
            raise RuntimeError("native packet is closed")
        if self.unavailable:
            return False
        if self.lib is not None:
            return True
        if self._binding is not None:
            budget.consume(4096)
            ffi, lib, manifest = self._binding
            maximum = 0 if self.variant == "scalar" else NATIVE_MAX_WORKERS
            if (
                type(manifest.get("native_team_max_workers")) is not int
                or manifest["native_team_max_workers"] != maximum
                or int(lib.b4_team_max_workers()) != maximum
            ):
                raise ValueError("native team-count capability mismatch")
            self._allocate(ffi, lib, manifest)
            return True
        if sys.implementation.name != "pypy" or sys.version_info[:2] != (
            3,
            11,
        ):
            raise BackendUnavailableError(
                "native ECM requires PyPy Python 3.11"
            )

        # Fund hashing, context creation and unused setup before loading.
        budget.consume(1)
        directory = SOURCE / "build"
        if self.variant != "scalar":
            directory /= self.variant
        paths = [
            directory / "manifest.json",
            directory / "libfactor_ecm.dylib",
        ]
        paths.extend(
            SOURCE / name
            for name in ("packet.c", "packet.h", "team.c", "curve.c")
        )
        if (SOURCE / "openmp.c").exists():
            paths.append(SOURCE / "openmp.c")
        ffi = lib = None
        try:
            sizes = [path.stat().st_size for path in paths]
            if (
                sizes[0] > 16384
                or sizes[1] > 256 * 1024
                or max(sizes[2:]) > 65536
            ):
                raise ValueError("native load exceeds funded buffer sizes")
            budget.consume(4096 + sum(sizes))
            manifest = json.loads(paths[0].read_text())
            library = paths[1]
            dependency = Path(manifest["gmp_library"])
            if dependency.stat().st_size > 4 * 2**20:
                raise ValueError("native GMP dependency exceeds load bound")
            budget.consume(dependency.stat().st_size)
            if _sha_file(dependency) != manifest["gmp_library_sha256"]:
                raise ValueError("native GMP dependency identity mismatch")
            if self.variant == "openmp":
                runtime = (
                    Path(manifest["openmp"]["prefix"]) / "lib/libomp.dylib"
                )
                if runtime.stat().st_size > 4 * 2**20:
                    raise ValueError(
                        "native OpenMP dependency exceeds load bound"
                    )
                budget.consume(runtime.stat().st_size)
                if _sha_file(runtime) != manifest["openmp"]["runtime_sha256"]:
                    raise ValueError(
                        "native OpenMP dependency identity mismatch"
                    )
            if manifest["backend"] != self.variant:
                raise ValueError("native build backend identity mismatch")
            expected_sources = {
                "packet.c",
                "packet.h",
                "team.c",
                "curve.c",
            }
            if "whole_curve_openmp" in manifest:
                expected_sources.add("openmp.c")
            if set(manifest["source"]) != expected_sources:
                raise ValueError("incomplete native build source identity")
            if library.stat().st_size > 256 * 1024:
                raise ValueError("native library exceeds funded load size")
            if manifest["abi"] != ABI or any(
                hashlib.sha256((SOURCE / name).read_bytes()).hexdigest() != pin
                for name, pin in manifest["source"].items()
            ):
                raise ValueError("native build source/ABI mismatch")
            if (
                hashlib.sha256(library.read_bytes()).hexdigest()
                != manifest["library_sha256"]
            ):
                raise ValueError("native library identity mismatch")
            binding = manifest.get("abi_binding")
            if binding is None:
                import cffi

                ffi = cffi.FFI()
                ffi.cdef(
                    "\n".join(
                        line
                        for line in (SOURCE / "packet.h")
                        .read_text()
                        .splitlines()
                        if not line.startswith("#")
                    )
                )
            else:
                if binding.get("file") != "_factor_ecm_abi.abi":
                    raise ValueError("native binding name mismatch")
                path = directory / binding["file"]
                if not path.is_file() or path.stat().st_size > 128 * 1024:
                    raise ValueError("native binding missing or oversized")
                budget.consume(path.stat().st_size)
                raw = path.read_bytes()
                if hashlib.sha256(raw).hexdigest() != binding["sha256"]:
                    raise ValueError("native binding identity mismatch")
                # Execute the validated bytes, not a second path lookup.
                namespace = {}
                exec(compile(raw, str(path), "exec"), namespace)
                ffi = namespace["ffi"]
            if self.variant == "openmp" and _OPENMP_BINDINGS:
                ffi, lib, bound = _OPENMP_BINDINGS[0]
                if bound != manifest:
                    raise ValueError(
                        "cannot reload a different OpenMP build in one process"
                    )
            else:
                lib = ffi.dlopen(str(library))
                if self.variant == "openmp":
                    _OPENMP_BINDINGS.append((ffi, lib, manifest))
            if (
                ffi.string(lib.b4_gmp_version()).decode("ascii")
                != manifest["gmp_runtime"]
            ):
                raise ValueError("loaded native GMP version mismatch")
            if bool(lib.b4_openmp()) != (self.variant == "openmp"):
                raise ValueError("native team capability mismatch")
            if bool(lib.b4_pthread()) != (self.variant == "pthread"):
                raise ValueError("native pthread capability mismatch")
            if lib.b4_abi() != ABI:
                raise ValueError("native library ABI mismatch")
            if bool(lib.b4_limb_division()) != manifest.get(
                "limb_division", False
            ):
                raise ValueError("native limb build identity mismatch")
            if (
                "wide_gmp_remainder" in manifest
                and bool(lib.b4_wide_remainder())
                != manifest["wide_gmp_remainder"]
            ):
                raise ValueError("native wide-remainder identity mismatch")
            if (
                "coalesce_unit_yields" in manifest
                and bool(lib.b4_coalesced_unit_yields())
                != manifest["coalesce_unit_yields"]
            ):
                raise ValueError("native unit-yield build identity mismatch")
            if (
                "whole_curve_openmp" in manifest
                and bool(lib.b4_omp_whole_kernel())
                != manifest["whole_curve_openmp"]
            ):
                raise ValueError("native whole-curve OpenMP identity mismatch")
            if (
                "slice_actions" in manifest
                and int(lib.b4_slice_actions()) != manifest["slice_actions"]
            ):
                raise ValueError("native action quantum identity mismatch")
            maximum = 0 if self.variant == "scalar" else NATIVE_MAX_WORKERS
            if (
                type(manifest.get("native_team_max_workers")) is not int
                or manifest["native_team_max_workers"] != maximum
                or int(lib.b4_team_max_workers()) != maximum
            ):
                raise ValueError("native team-count capability mismatch")
            self._allocate(ffi, lib, manifest)
        except BaseException as error:
            if lib is not None and self.variant != "openmp":
                ffi.dlclose(lib)
            if isinstance(error, (ImportError, OSError)):
                if self.mode == "auto":
                    self.unavailable = True
                    return False
                backend = {
                    "openmp": "OpenMP team",
                    "pthread": "pthread team",
                    "scalar": "C/GMP",
                }[self.variant]
                raise BackendUnavailableError(
                    f"native ECM {backend} backend unavailable; "
                    "build v2.ecm.native.build explicitly"
                ) from error
            raise

        return True

    def _allocate(self, ffi, lib, manifest):
        context = lib.b4_create()
        if context == ffi.NULL:
            raise MemoryError("native context allocation failed")
        buffers = []
        try:
            for declaration in (
                "unsigned char[1024]",
                "b4_value[5]",
                "b4_record[16]",
                "b4_output *",
            ):
                buffers.append(ffi.new(declaration))
            transport, inputs = buffers[:2]
            output = buffers[3]
            for index, offset in enumerate((0, 128, 256, 384, 512)):
                inputs[index].data = transport + offset
            views = {
                name: ffi.buffer(getattr(output, name), 128)
                for name in ("x", "z", "guard")
            }
        except BaseException:
            lib.b4_destroy(context)
            for owner in reversed(buffers):
                ffi.release(owner)
            raise
        self.ffi, self.lib, self.context = ffi, lib, context
        self.buffers = buffers
        self.transport, self.inputs, self.records_buffer, self.output = buffers
        self.output_views = views
        self.manifest = manifest
        self.build = dict(
            manifest, gmp=ffi.string(lib.b4_gmp_version()).decode()
        )

    def available(self, n, budget):
        if not self._lock.acquire(blocking=False):
            raise RuntimeError("native packet already active")
        try:
            return self._available(n, budget)
        finally:
            self._lock.release()

    def _available(self, n, budget):
        if self.closed:
            raise RuntimeError("native packet is closed")
        eligible = (
            type(n) is int and 3 <= n and n % 2 and n.bit_length() <= 1024
        )
        if not eligible:
            if self.mode == "auto":
                return False
            raise ValueError("native ECM requires an odd 3..1024-bit modulus")
        return self._load(budget)

    def _values(self, values):
        """Make raw diagnostic ABI owners; the caller must release them."""
        ffi = self.ffi
        owners = [
            ffi.new("unsigned char[]", _bytes(value)) for value in values
        ]
        records = [ffi.new("b4_value *") for _ in values]
        for record, owner, value in zip(records, owners, values):
            record.data = owner
            record.size = max(1, (value.bit_length() + 7) // 8)
        return owners, records

    def _inputs(self, values):
        for index, value in enumerate(values):
            if type(value) is not int or value < 0:
                raise ValueError(
                    "native packet requires unsigned Python integers"
                )
            if self.input_cache[index] == value:
                continue
            raw = _bytes(value)
            capacity = 512 if index == 4 else 128
            if len(raw) > capacity:
                raise ValueError("native input exceeds its canonical buffer")
            # Cache equality is only a transport optimization. C still binds
            # every begin to the explicitly supplied modulus/a24 and point.
            self.ffi.memmove(self.inputs[index].data, raw, len(raw))
            self.inputs[index].size = len(raw)
            self.input_cache[index] = value
        return [self.inputs[index] for index in range(len(values))]

    def _export(self):
        output = self.output
        if self.lib.b4_result(self.context, output):
            raise ValueError("native output is not committed")
        values = []
        for name in ("x", "z", "guard"):
            size = getattr(output, name + "_size")
            if not 1 <= size <= 128:
                raise ValueError("invalid native output length")
            raw = self.output_views[name][:size]
            if len(raw) > 1 and raw[0] == 0:
                raise ValueError("noncanonical native output")
            values.append(int.from_bytes(raw, "big"))
        return values, int(output.completed)

    def _finish(self, budget):
        while True:
            poll_budget(budget)
            result = self.lib.b4_step(self.context, 5, 128)
            if result < 0:
                raise ValueError("native packet step refused")
            if result:
                break
        poll_budget(budget)
        return self._export()

    def _begin(self, n, budget):
        if not self._lock.acquire(blocking=False):
            raise RuntimeError("native packet already active")
        try:
            if not self._available(n, budget):
                raise BackendUnavailableError("native packet was not admitted")
        except BaseException:
            self._lock.release()
            raise

    def end(self):
        """Discard private state after a scalar call or drained team region."""
        if self.context is not None:
            self.lib.b4_abort(self.context)
        self._lock.release()

    def prepare_ladder(self, scalar, point, n, a24, budget):
        """Fund transport and hold private input buffers until end()."""
        if type(scalar) is not int or not 0 <= scalar < 2**4096:
            raise ValueError("native scalar exceeds 4096 bits")
        if not isinstance(point, (tuple, list)) or len(point) != 2:
            raise ValueError("native point requires two canonical coordinates")
        self._begin(n, budget)
        try:
            budget.consume(7 * ((n.bit_length() + 7) // 8) + 512)
            values = self._inputs([n, *point, a24, scalar])
            if self.lib.b4_ladder(self.context, ABI, *values):
                raise ValueError("native ladder input refused")
        except BaseException:
            self.end()
            raise

    def prepare_records(self, entries, point, n, a24, budget):
        """Bind verified immutable records; mutable GMP state stays private."""
        if not 0 < len(entries) <= 16:
            raise ValueError("native packet exceeds batch16")
        if not isinstance(point, (tuple, list)) or len(point) != 2:
            raise ValueError("native point requires two canonical coordinates")
        self._begin(n, budget)
        try:
            selected = []
            for _, _, action, _, _ in entries:
                key = (action.record, action.masks)
                if key not in self.programs:
                    code, masks = action.record.code, action.masks
                    size = 512 + len(code) + len(masks)
                    if (
                        len(self.programs) >= 333
                        or self.store.used + size > RESERVE // 2
                    ):
                        raise MemoryError("native program retention refused")
                    budget.consume(size)
                    verify_frontier(*key)
                    code_owner = self.ffi.new("unsigned char[]", code)
                    try:
                        mask_owner = self.ffi.new("unsigned char[]", masks)
                    except BaseException:
                        self.ffi.release(code_owner)
                        raise
                    self.programs[key] = (code_owner, mask_owner)
                    self.store.used += size
                selected.append((action, self.programs[key]))
            budget.consume(
                7 * ((n.bit_length() + 7) // 8) + 128 * len(entries)
            )
            for target, (action, (code, masks)) in zip(
                self.records_buffer, selected
            ):
                target.code, target.masks = code, masks
                target.code_size = len(action.record.code)
                target.masks_size = len(action.masks)
                target.slots, target.output = (
                    action.record.slots,
                    action.record.output,
                )
            values = self._inputs([n, *point, a24])
            if self.lib.b4_records(
                self.context, ABI, *values, self.records_buffer, len(entries)
            ):
                raise ValueError("native program input refused")
        except BaseException:
            self.end()
            raise

    @staticmethod
    def checked(result, n, count=None):
        values, completed = result
        if any(value >= n for value in values):
            raise ValueError("native output exceeds modulus")
        if count is None:
            return tuple(values[:2])
        if not 1 <= completed <= count:
            raise ValueError("invalid native program output")
        return tuple(values[:2]), values[2], completed

    def ladder(self, scalar, point, n, a24, budget):
        self.prepare_ladder(scalar, point, n, a24, budget)
        try:
            return self.checked(self._finish(budget), n)
        finally:
            self.end()

    def records(self, entries, point, n, a24, budget):
        self.prepare_records(entries, point, n, a24, budget)
        try:
            return self.checked(self._finish(budget), n, len(entries))
        finally:
            self.end()

    def close(self):
        if not self._lock.acquire(blocking=False):
            raise RuntimeError("cannot close an active native packet")
        try:
            if self.context is not None:
                self.lib.b4_destroy(self.context)
            self.context = None
            if self.owns_programs and self.ffi is not None:
                self.store.close(self.ffi)
            self.output_views.clear()
            for owner in reversed(self.buffers):
                self.ffi.release(owner)
            self.buffers.clear()
            self.input_cache = [None] * 5
            self.transport = self.inputs = self.records_buffer = (
                self.output
            ) = None
            if (
                self.lib is not None
                and self.variant != "openmp"
                and self._binding is None
            ):
                self.ffi.dlclose(self.lib)
            self.lib = self.ffi = None
            self.closed = True
        finally:
            self._lock.release()
