"""Build bounded GMP-ECM reference inputs with pinned upstream C sources.

Only scratch allocation capacity changes. All decoding uses the upstream
implementation with assertions, isolated under process/time/storage limits.
Generated output is then checked by two independent integer interpreters.
"""

import argparse
import hashlib
import json
import platform
import resource
import shutil
import struct
import subprocess
import sys
import time
from pathlib import Path

from ....common import prime_sieve
from ...support.paths import (
    BENCHMARK_ROOT,
    REPOSITORY_ROOT,
    source_path,
)
from . import c6_chains

ROOT = REPOSITORY_ROOT
UPSTREAM = BENCHMARK_ROOT / "inputs/upstream/c6_gmp_ecm"
PIN = "8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e"


def require_runtime():
    if platform.python_implementation() != "PyPy" or sys.version_info[:2] != (
        3,
        11,
    ):
        raise RuntimeError("requires PyPy implementing Python 3.11")


def limits():
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    resource.setrlimit(resource.RLIMIT_FSIZE, (16 * 1024**2,) * 2)


def run(command, directory, name, *, input_text=None):
    started = time.perf_counter()
    with (directory / name).open("w") as log:
        process = subprocess.Popen(
            command,
            cwd=directory,
            text=True,
            stdin=subprocess.PIPE if input_text is not None else None,
            stdout=log,
            stderr=subprocess.STDOUT,
            preexec_fn=limits,
        )
        pending = input_text
        try:
            while True:
                try:
                    process.communicate(input=pending, timeout=0.1)
                    break
                except subprocess.TimeoutExpired:
                    pending = None
                    rss = subprocess.check_output(
                        ["ps", "-o", "rss=", "-p", str(process.pid)],
                        text=True,
                    ).strip()
                    if (
                        rss and int(rss) > 512 * 1024
                    ) or time.perf_counter() - started > 60:
                        raise RuntimeError(
                            "upstream process exceeded finite allowance"
                        )
            if process.returncode:
                raise RuntimeError("upstream command failed; see " + name)
        except BaseException:
            process.kill()
            process.communicate()
            raise
    return time.perf_counter() - started


def prepare(directory):
    directory.mkdir(parents=True, exist_ok=False)
    for name in ("LucasChainGen.c", "LucasChainGen.h", "LCG_macros.h"):
        shutil.copyfile(UPSTREAM / name, directory / name)
    header = directory / "LucasChainGen.h"
    original = source_path(header).read_text()
    assert original.count("MAX_CODE_OR_PRIME_COUNT 6000000") == 1
    header.write_text(
        original.replace(
            "MAX_CODE_OR_PRIME_COUNT 6000000", "MAX_CODE_OR_PRIME_COUNT 4096"
        )
    )

    source = (source_path(UPSTREAM / "ecm.c")).read_text()
    start = source.index("/* PBMcL additions")
    stop = source.index("/* end PBMcL additions */")
    declarations = source[start:stop]
    start = source.index("/* functions for using optimal")
    stop = source.index("/* Input: x is initial point")
    decoder = source[start:stop]
    wrapper = r"""
int main(void) {
    uint64_t p, code;
    while (scanf("%" SCNu64 " %" SCNu64, &p, &code) == 2) {
        assert(p >= 2 && p <= 2000);
        chain_element chain[64] = {{0}};
        chain[0].value = 1;
        chain[1] = (chain_element){2, 0, 0, 0};
        chain[2] = (chain_element){3, 0, 1, 1};
        uint8_t length = generate_Lucas_chain(p, code, chain);
        assert(length < 64);
        printf("%" PRIu64 " %u\n", p, length);
        for (unsigned i=1; i<=length; i++)
            printf("%" PRIu64 " %u %u %u\n", chain[i].value,
                   chain[i].comp_offset_1, chain[i].comp_offset_2,
                   chain[i].dif_offset);
    }
    return 0;
}
"""
    (directory / "decode.c").write_text(
        "#include <stdio.h>\n#include <stdint.h>\n#include <inttypes.h>\n"
        "#include <assert.h>\n#define ASSERT assert\n"
        + declarations
        + decoder
        + wrapper
    )
    return {
        "generator_compile_seconds": run(
            ["cc", "-O2", "-pthread", "LucasChainGen.c", "-o", "generate"],
            directory,
            "compile-generator.txt",
        ),
        "decoder_compile_seconds": run(
            ["cc", "-O2", "decode.c", "-o", "decode"],
            directory,
            "compile-decoder.txt",
        ),
    }


def generate(directory):
    generation = run(
        ["./generate", "-B1", "2000", "-nT", "1"], directory, "generator.txt"
    )
    raw = (source_path(directory / "Lchain_codes.dat")).read_bytes()
    primes = list(prime_sieve.prime_sieve(2001))
    assert len(raw) == 8 * (len(primes) - 4)
    codes = [0] * 4 + list(struct.unpack("=" + "Q" * (len(primes) - 4), raw))
    assert all(codes[4:])
    decoding = run(
        ["./decode"],
        directory,
        "decoded.txt",
        input_text="".join(f"{p} {c}\n" for p, c in zip(primes, codes)),
    )
    rows, lines = (
        [],
        iter(
            (source_path(directory / "decoded.txt")).read_text().splitlines()
        ),
    )
    for prime, code in zip(primes, codes):
        scalar, count = map(int, next(lines).split())
        assert scalar == prime and count < 64
        rows.append(
            dict(
                prime=prime,
                code=code,
                elements=[
                    list(map(int, next(lines).split())) for _ in range(count)
                ],
            )
        )
    assert next(lines, None) is None
    started = time.perf_counter()
    c6_chains.decode_rows(rows)
    verification = time.perf_counter() - started
    return rows, dict(
        generation_seconds=generation,
        decoding_seconds=decoding,
        verification_seconds=verification,
        code_bytes=len(raw),
        code_sha256=hashlib.sha256(raw).hexdigest(),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    require_runtime()
    from .c6_study import performance_window

    with performance_window("upstream generation"):
        timings = prepare(args.scratch.resolve())
        rows, costs = generate(args.scratch.resolve())
        timings.update(costs)
        data = dict(
            upstream_commit=PIN,
            bound=2000,
            capacity=4096,
            threads=1,
            source_sha256={
                p.name: hashlib.sha256(source_path(p).read_bytes()).hexdigest()
                for p in UPSTREAM.iterdir()
                if p.is_file()
            },
            code_sha256=costs["code_sha256"],
            records=rows,
        )
        c6_chains.DATA.write_text(
            json.dumps(data, separators=(",", ":")) + "\n"
        )
        (args.scratch / "costs.json").write_text(json.dumps(timings, indent=2))
        print(json.dumps(timings, indent=2))


if __name__ == "__main__":
    main()
