"""Build only the optional, ECM-local two-lane NEON feasibility screen."""

import argparse
import json
import platform
import subprocess
import sys
from pathlib import Path

from .build import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmp-prefix", type=Path, required=True)
    parser.add_argument("--cc", default="/usr/bin/clang")
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("NEON build requires PyPy Python 3.11")
    if sys.platform != "darwin" or platform.machine() != "arm64":
        raise RuntimeError("this optional prototype targets ARM64 macOS")
    import cffi

    root = Path(__file__).resolve().parent
    source, output = root / "neon", root / "build/neon"
    output.mkdir(parents=True, exist_ok=True)
    flags = [
        "-std=c11",
        "-O2",
        "-Wall",
        "-Wextra",
        "-Werror",
        "-fno-lto",
        "-fno-vectorize",
        "-fno-slp-vectorize",
        "-I" + str(args.gmp_prefix / "include"),
    ]
    commands, objects = [], []
    for name in ("kernels", "screen", "sink"):
        obj = output / (name + ".o")
        command = [
            args.cc,
            *flags,
            "-c",
            str(source / (name + ".c")),
            "-o",
            str(obj),
        ]
        subprocess.run(command, check=True)
        commands.append(command)
        objects.append(str(obj))
    library = output / "libfactor_ecm_neon.dylib"
    command = [
        args.cc,
        "-dynamiclib",
        "-fno-lto",
        *objects,
        "-L" + str(args.gmp_prefix / "lib"),
        "-Wl,-rpath," + str(args.gmp_prefix / "lib"),
        "-lgmp",
        "-o",
        str(library),
    ]
    subprocess.run(command, check=True)
    commands.append(command)
    disassembly = subprocess.check_output(
        ["/usr/bin/otool", "-tvV", str(library)]
    )
    (output / "disassembly.txt").write_bytes(disassembly)

    ffi = cffi.FFI()
    ffi.cdef("unsigned b4n_abi(void); const char *b4n_gmp_version(void);")
    lib = ffi.dlopen(str(library))
    try:
        if lib.b4n_abi() != 1:
            raise RuntimeError("incompatible NEON prototype ABI")
        gmp = ffi.string(lib.b4n_gmp_version()).decode("ascii")
    finally:
        ffi.dlclose(lib)
    pins = {
        "neon/" + name: sha(source / name)
        for name in (
            "kernels.c",
            "kernels.h",
            "screen.c",
            "screen.h",
            "sink.c",
        )
    }
    pins["neon_build.py"] = sha(Path(__file__))
    pins["packet.c"] = sha(root / "packet.c")
    manifest = dict(
        identity="b4a-ecm-neon-v3",
        abi=1,
        compiler=subprocess.check_output([args.cc, "--version"], text=True),
        commands=commands,
        target=platform.machine(),
        platform=platform.platform(),
        python=sys.version,
        cffi=cffi.__version__,
        gmp_runtime=gmp,
        gmp_library=str((args.gmp_prefix / "lib/libgmp.10.dylib").resolve()),
        gmp_library_sha256=sha(args.gmp_prefix / "lib/libgmp.10.dylib"),
        source=pins,
        library_sha256=sha(library),
        disassembly_sha256=sha(output / "disassembly.txt"),
        codegen_review="pending; inspect scalar and explicit NEON symbols",
        auto_vectorization=False,
        lto=False,
        native_threads=1,
    )
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(library)


if __name__ == "__main__":
    main()
