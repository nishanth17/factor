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
    parser.add_argument("--vector-report", action="store_true")
    args = parser.parse_args()
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
    output.mkdir(parents=True, exist_ok=True)
    library = output / "libfactor_ecm.dylib"
    flags = ["-std=c11", "-O2", "-Wall", "-Wextra", "-Werror"]
    if args.vector_report:
        flags += ["-Rpass=loop-vectorize", "-Rpass-missed=loop-vectorize"]
    omp_flags = []
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
    finally:
        ffi.dlclose(lib)
    manifest = dict(
        abi=1,
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
        source={name: sha(source / name) for name in ("packet.c", "packet.h")},
        library_sha256=sha(library),
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
