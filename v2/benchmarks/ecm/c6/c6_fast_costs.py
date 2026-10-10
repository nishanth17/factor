"""Separate optimized C6 construction and identical-chain layout costs."""

import argparse
import hashlib
import json
import resource
import sys
from pathlib import Path

from ...support.paths import (
    source_path,
)
from ...support.prac_oracle import affine_multiply, matches
from ..p41.p41_campaign import PYTHON_BACKEND, gmp_backend
from . import build_c6_cf_inputs as cf_builder
from . import c6_cf, c6_chains, c6_fast, c6_fast_study, c6_study
from .build_c6_inputs import require_runtime
from .c6_costs import measured


class ThreePoint(c6_fast.FastRecord):
    """Algorithm 1 with the same arithmetic and factor-coverage masks."""

    def __init__(self, chain, backend):
        bits = c6_chains.continued_fraction_bits(chain)
        if bits is None:
            raise ValueError("not an identical-arithmetic CF record")
        super().__init__(c6_chains.compact(chain), backend)
        self.bits = bits
        self.run = self.three

    def three(self, point, n, a24):
        add, double = self.backend.ecm.point_add, self.backend.ecm.point_double
        factor, masks = c6_fast._factor, self.masks
        a = point
        product = factor(1, a, masks[0])
        b = double(*a, n, a24)
        product = factor(product, b, masks[1])
        c = add(*a, *b, *a, n)
        product = factor(product, c, masks[2])
        for bit, mask in zip(self.bits, masks[3:]):
            if bit:
                a, b, c = a, c, add(*a, *c, *b, n)
            else:
                a, b, c = b, c, add(*b, *c, *a, n)
            product = factor(product, c, mask)
        return c, product


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require_runtime()
    settings = c6_fast_study.protocol()
    selection = json.loads(source_path(c6_fast_study.SELECTION).read_text())
    cf_path = cf_builder.PROTOCOL.parent / "c6_cf_selection.json"
    cf_selection = json.loads(source_path(cf_path).read_text())
    cf_settings = json.loads(source_path(cf_builder.PROTOCOL).read_text())
    for name, digest in cf_settings["sha256"].items():
        if (
            hashlib.sha256(
                (source_path(c6_study.ROOT / name)).read_bytes()
            ).hexdigest()
            != digest
        ):
            raise ValueError("CF frozen source changed: " + name)
    for backend, arms in cf_selection["candidates"].items():
        selection["candidates"][backend] += arms
    report = dict(
        protocol=settings,
        source_sha256=hashlib.sha256(
            source_path(Path(__file__)).read_bytes()
        ).hexdigest(),
        selection=selection,
        runtime=sys.version,
        cf_protocol=cf_settings,
        cf_generation=[],
        construction={},
        layouts={},
        storage={},
    )
    with c6_study.performance_window("C6 optimized costs and layouts"):
        expected_cf = json.loads(source_path(c6_cf.DATA).read_text())
        for index in range(9):
            data, costs = cf_builder.generate(
                (
                    c6_study.ROOT
                    / "v2/audit/results/c6-fast"
                    / ("costs-" + args.output.stem)
                    / f"cf-regenerate-{index}"
                ).resolve()
            )
            if data != expected_cf:
                raise AssertionError("CF regeneration changed")
            report["cf_generation"].append(costs)
        report["load_cf_verify"] = measured(c6_cf.load_catalog)
        expected = json.loads(source_path(c6_fast.DATA).read_text())

        def generate():
            actual = c6_fast.generate_catalog()
            if actual != expected:
                raise AssertionError("precomputation changed")

        report["catalog_precomputation"] = measured(generate)
        catalog = c6_fast.load_catalog()
        report["load_all_verify"] = measured(c6_fast.load_catalog)
        for family, records in catalog.items():
            report["storage"][family] = dict(
                records=len(records),
                bytecode_bytes=sum(
                    len(c.record.code) for c in records.values()
                ),
                mask_bytes=sum(len(c.masks) for c in records.values()),
                point_registers=max(c.record.slots for c in records.values()),
            )
        for name, arms in selection["candidates"].items():
            backend = PYTHON_BACKEND if name == "int" else gmp_backend()
            for arm in arms:

                def construct():
                    if arm.startswith("fast/cf/"):
                        _, _, mode, batch = arm.split("/")
                        program = c6_cf.build_program(
                            2000, mode, backend, int(batch)
                        )
                    else:
                        program = c6_fast_study.construct(arm, backend)
                    if len(program.entries) != 303:
                        raise AssertionError("wrong full-stage program")
                    return program

                report["construction"][name + "/" + arm] = measured(construct)
                program = construct()
                actions = [row[2] for row in program.entries]
                report["storage"][name + "/" + arm] = dict(
                    source_bytes=sum(a.source_bytes for a in actions),
                    code_bytes=sum(len(a.record.code) for a in actions),
                    instructions=sum(len(a.record.code) // 4 for a in actions),
                    additions=sum(
                        op[3] != c6_chains.DOUBLE
                        for a in actions
                        for op in c6_fast.instructions(a.record)
                    ),
                    doublings=sum(
                        op[3] == c6_chains.DOUBLE
                        for a in actions
                        for op in c6_fast.instructions(a.record)
                    ),
                    point_registers=max(a.record.slots for a in actions),
                    coverage_bytes=sum(len(a.masks) for a in actions),
                    guarded_coordinates=sum(
                        sum(bool(m & 1) + bool(m & 2) for m in a.masks)
                        for a in actions
                    ),
                )
        upstream = c6_chains.load_lucas()
        subset = [
            c
            for c in upstream.values()
            if c6_chains.continued_fraction_bits(c) is not None
        ]
        report["identical_arithmetic_primes"] = [c.scalar for c in subset]
        expected = [
            affine_multiply(c.scalar, (3, 293), 1009, 6) for c in subset
        ]
        for name, backend in (("int", PYTHON_BACKEND), ("gmp", gmp_backend())):
            for layout in ("compact", "rolling", "three"):
                actions = [
                    ThreePoint(c, backend)
                    if layout == "three"
                    else c6_fast.FastRecord(
                        c6_chains.compact(c, rolling=layout == "rolling"),
                        backend,
                    )
                    for c in subset
                ]
                n, a24 = backend.integer(1009), backend.integer(2)
                point = backend.integer(3), backend.integer(1)

                def execute():
                    for action, target in zip(actions, expected):
                        if not matches(action(point, n, a24), target, 1009):
                            raise AssertionError("independent layout mismatch")

                report["layouts"][name + "/" + layout] = measured(execute)
        report["maxrss_bytes_macos"] = resource.getrusage(
            resource.RUSAGE_SELF
        ).ru_maxrss
        with args.output.open("x") as output:
            json.dump(report, output, indent=2)
            output.write("\n")
        c6_fast_study.protocol()


if __name__ == "__main__":
    main()
