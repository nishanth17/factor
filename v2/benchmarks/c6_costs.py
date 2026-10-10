"""C6 construction, storage, generation and interpreter diagnostics."""

import argparse
import json
import resource
import shutil
import statistics
import sys
import time
from pathlib import Path

from . import build_c6_inputs as builder
from . import c6_chains as chains
from . import c6_study as study
from .p41_campaign import PYTHON_BACKEND, gmp_backend
from .prac_oracle import affine_multiply, matches


def measured(function):
    captures = []
    for warmup, count in ((3, 9), (5, 18), (8, 27)):
        started, runs = time.perf_counter(), 0
        while time.perf_counter() - started < warmup:
            function()
            runs += 1
        samples = []
        for _ in range(count):
            start = time.perf_counter()
            function()
            samples.append(time.perf_counter() - start)
        row = dict(
            warmup_seconds=time.perf_counter() - started - sum(samples),
            validated_warmups=runs,
            samples=samples,
            median=statistics.median(samples),
            relative_iqr=study.spread(samples),
        )
        captures.append(row)
        if row["relative_iqr"] <= 0.15:
            break
    return dict(attempts=captures, accepted=captures[-1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--upstream-build", type=Path, required=True)
    args = parser.parse_args()
    builder.require_runtime()
    settings = study.protocol()
    with study.performance_window("construction and interpreter diagnostics"):
        catalog = chains.load_lucas()
        report = dict(
            protocol=settings,
            runtime=sys.version,
            generation=[],
            construction={},
            interpreters={},
            storage={},
        )
        for i in range(9):
            directory = args.output.parent / f"generator-repeat-{i}"
            directory.mkdir()
            for name in ("generate", "decode"):
                shutil.copyfile(args.upstream_build / name, directory / name)
                (directory / name).chmod(0o700)
            records, costs = builder.generate(directory.resolve())
            expected = json.loads(chains.DATA.read_text())["records"]
            if records != expected:
                raise AssertionError("upstream regeneration mismatch")
            report["generation"].append(costs)
        for method in ("checked", "compact", "lucas", "rolling"):

            def construct():
                program = chains.build_program(2000, method, PYTHON_BACKEND)
                assert len(program) == 303
                return program

            report["construction"][method] = measured(construct)
            program = construct()
            actions = {id(a): a for row in program for a in row[2:]}.values()
            if method == "checked":
                report["storage"][method] = dict(
                    unique_records=len(list(actions)),
                    instructions=sum(len(a.instructions) for a in actions),
                )
            else:
                report["storage"][method] = dict(
                    unique_records=len(list(actions)),
                    instructions=sum(len(a.record.code) // 4 for a in actions),
                    code_bytes=sum(len(a.record.code) for a in actions),
                    max_registers=max(a.record.slots for a in actions),
                )
        report["construction"]["load_decode_verify"] = measured(
            chains.load_lucas
        )

        # Compare exactly the same eligible chains: no search or cherry-picked
        # arithmetic subset within the recognized CF family.
        subset = [
            c
            for c in catalog.values()
            if chains.continued_fraction_bits(c) is not None
        ]
        report["cf_primes"] = [c.scalar for c in subset]
        for backend_name, backend in (
            ("int", PYTHON_BACKEND),
            ("gmp", gmp_backend()),
        ):
            modulus, curve_a, point = 1009, 6, (3, 293)
            expected = [
                affine_multiply(c.scalar, point, modulus, curve_a)
                for c in subset
            ]
            for layout in ("compact", "rolling", "three"):
                if layout == "three":
                    actions = [
                        chains.ThreePointExecutor(c, backend) for c in subset
                    ]
                else:
                    actions = [
                        chains.Executor(
                            chains.compact(c, rolling=layout == "rolling"),
                            backend,
                        )
                        for c in subset
                    ]
                n, a24 = backend.integer(modulus), backend.integer(2)
                start_point = backend.integer(3), backend.integer(1)

                def execute():
                    for action, target in zip(actions, expected):
                        actual = action(start_point, n, a24)
                        if not matches(actual, target, modulus):
                            raise AssertionError("CF interpreter mismatch")

                report["interpreters"][backend_name + "/" + layout] = measured(
                    execute
                )
        report["maxrss_bytes_macos"] = resource.getrusage(
            resource.RUSAGE_SELF
        ).ru_maxrss
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report["storage"], indent=2))


if __name__ == "__main__":
    main()
