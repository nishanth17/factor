"""Explicit project-local build; never invoked by the factoring library."""

import argparse
import hashlib
import json
import platform
import re
import subprocess
import sys
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cc", default="/usr/bin/clang")
    parser.add_argument("--gmp-prefix", type=Path, required=True)
    parser.add_argument("--openmp-prefix", type=Path)
    parser.add_argument("--pthread", action="store_true")
    parser.add_argument("--vector-report", action="store_true")
    parser.add_argument("--limb-division", action="store_true")
    parser.add_argument("--coalesce-unit-yields", action="store_true")
    parser.add_argument("--extended-quantum", action="store_true")
    parser.add_argument("--split-team-signals", action="store_true")
    parser.add_argument("--optimization", choices=("O2", "O3"), default="O2")
    parser.add_argument("--cpu-target")
    args = parser.parse_args()
    if args.pthread and args.openmp_prefix is not None:
        parser.error("select pthread or OpenMP, not both")
    if args.split_team_signals and not args.pthread:
        parser.error("split team signals require pthreads")
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("build requires PyPy implementing Python 3.11")

    required = [
        args.gmp_prefix / "include/gmp.h",
        args.gmp_prefix / "lib/libgmp.10.dylib",
    ]
    if args.openmp_prefix is not None:
        required += [
            args.openmp_prefix / "include/omp.h",
            args.openmp_prefix / "lib/libomp.dylib",
        ]
    for path in required:
        if not path.is_file():
            backend = "OpenMP" if "omp" in path.name else "C/GMP"
            parser.error(f"requested {backend} backend unavailable: {path}")

    import cffi

    source = Path(__file__).resolve().parent
    output = source / "build"
    if args.openmp_prefix is not None:
        output /= "openmp"
    elif args.pthread:
        output /= "pthread"
    output.mkdir(parents=True, exist_ok=True)
    library = output / "libfactor_ecm.dylib"
    flags = [
        "-std=c11",
        "-" + args.optimization,
        "-Wall",
        "-Wextra",
        "-Werror",
    ]
    if args.cpu_target is not None:
        flags.append("-mcpu=" + args.cpu_target)
    if args.limb_division:
        flags.append("-DB4_MPN_DIVISION")
    if args.coalesce_unit_yields:
        flags.append("-DB4_COALESCE_UNIT_YIELDS")
    if args.extended_quantum:
        flags.append("-DB4_EXTENDED_QUANTUM")
    if args.split_team_signals:
        flags.append("-DB4_SPLIT_SIGNALS")
    if args.vector_report:
        flags += ["-Rpass=loop-vectorize", "-Rpass-missed=loop-vectorize"]
    omp_flags = []
    if args.pthread:
        omp_flags = ["-DB4_PTHREAD", "-pthread"]
    if args.openmp_prefix is not None:
        omp_flags = [
            "-DB4_OPENMP",
            "-Xpreprocessor",
            "-fopenmp",
            "-I" + str(args.openmp_prefix / "include"),
            "-L" + str(args.openmp_prefix / "lib"),
            "-Wl,-rpath," + str(args.openmp_prefix / "lib"),
            "-lomp",
        ]
    command = [
        args.cc,
        *flags,
        *omp_flags,
        "-dynamiclib",
        str(source / "packet.c"),
        str(source / "team.c"),
        "-I" + str(args.gmp_prefix / "include"),
        "-L" + str(args.gmp_prefix / "lib"),
        "-Wl,-rpath," + str(args.gmp_prefix / "lib"),
        "-lgmp",
        "-o",
        str(library),
    ]
    subprocess.run(command, check=True)
    ffi = cffi.FFI()
    ffi.cdef(
        "unsigned b4_abi(void); unsigned b4_openmp(void); "
        "unsigned b4_pthread(void); "
        "const char *b4_gmp_version(void);"
    )
    lib = ffi.dlopen(str(library))
    try:
        if lib.b4_abi() != 1:
            raise RuntimeError("built ECM packet has an incompatible ABI")
        gmp_runtime = ffi.string(lib.b4_gmp_version()).decode("ascii")
        if bool(lib.b4_openmp()) != (args.openmp_prefix is not None):
            raise RuntimeError(
                "built ECM team capability differs from request"
            )
        if bool(lib.b4_pthread()) != args.pthread:
            raise RuntimeError("built pthread capability differs from request")
    finally:
        ffi.dlclose(lib)
    declarations = "\n".join(
        line
        for line in (source / "packet.h").read_text().splitlines()
        if not line.startswith("#")
    )
    builder = cffi.FFI()
    builder.cdef(declarations)
    builder.set_source("_factor_ecm_abi", None)
    binding = output / "_factor_ecm_abi.abi"
    builder.emit_python_code(str(binding))
    manifest = dict(
        abi=1,
        backend="pthread"
        if args.pthread
        else ("openmp" if args.openmp_prefix is not None else "scalar"),
        compiler=subprocess.check_output([args.cc, "--version"], text=True),
        command=command,
        target=platform.machine(),
        platform=platform.platform(),
        python=sys.version,
        cffi=cffi.__version__,
        gmp_prefix=str(args.gmp_prefix.resolve()),
        gmp_runtime=gmp_runtime,
        gmp_header=str((args.gmp_prefix / "include/gmp.h").resolve()),
        gmp_header_sha256=sha(args.gmp_prefix / "include/gmp.h"),
        gmp_header_version=".".join(
            re.search(
                r"#define\s+" + name + r"\s+(\d+)",
                (args.gmp_prefix / "include/gmp.h").read_text(),
            )[1]
            for name in (
                "__GNU_MP_VERSION",
                "__GNU_MP_VERSION_MINOR",
                "__GNU_MP_VERSION_PATCHLEVEL",
            )
        ),
        linked_libraries=subprocess.check_output(
            ["/usr/bin/otool", "-L", str(library)], text=True
        ),
        gmp_library=str((args.gmp_prefix / "lib/libgmp.10.dylib").resolve()),
        gmp_library_sha256=sha(args.gmp_prefix / "lib/libgmp.10.dylib"),
        source={
            name: sha(source / name)
            for name in ("packet.c", "packet.h", "team.c", "curve.c")
        },
        library_sha256=sha(library),
        abi_binding=dict(file=binding.name, sha256=sha(binding)),
        limb_division=args.limb_division,
        coalesce_unit_yields=args.coalesce_unit_yields,
        slice_actions=256 if args.extended_quantum else 128,
        split_team_signals=args.split_team_signals,
        optimization=args.optimization,
        cpu_target=args.cpu_target,
        openmp=None
        if args.openmp_prefix is None
        else dict(
            prefix=str(args.openmp_prefix.resolve()),
            header_sha256=sha(args.openmp_prefix / "include/omp.h"),
            runtime_sha256=sha(args.openmp_prefix / "lib/libomp.dylib"),
        ),
    )
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(library)


if __name__ == "__main__":
    main()
