"""Isolated B4 control and bounded C6 kernels; production is untouched."""

import hashlib
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

from . import c6_cf, c6_fast
from .c6_chains import DOUBLE

PIN = Path(__file__).parent / "inputs/controls/c6_b4_ecm_a521573.py.txt"
PIN_HASH = "326f04ff19c07d85225603aa9247ee830e5f9c4b2dc7a72ee231da461ec7fa50"


def load_backend(backend):
    """Load the exact B4 ECM module under an isolated package name."""
    raw = PIN.read_bytes()
    if hashlib.sha256(raw).hexdigest() != PIN_HASH:
        raise ValueError("B4 source control changed")
    module = ModuleType("v2._c6_b4_control")
    module.__package__ = "v2"
    sys.modules[module.__name__] = module
    exec(compile(raw, str(PIN), "exec"), module.__dict__)
    return SimpleNamespace(
        ecm=module, prac=backend.prac, gcd=backend.gcd, integer=backend.integer
    )


def point_add(px, pz, qx, qz, rx, rz, n):
    """B4's two early addition reductions, with identical canonical exits."""
    u = (px - pz) * (qx + qz) % n
    v = (px + pz) * (qx - qz) % n
    total, difference = u + v, u - v
    return rz * total * total % n, rx * difference * difference % n


def point_double(px, pz, n, a24):
    """B4's two early square reductions, with identical canonical exits."""
    total, difference = px + pz, px - pz
    aa = total * total % n
    bb = difference * difference % n
    delta = aa - bb
    return aa * bb % n, delta * (bb + a24 * delta) % n


def double_add(point, other, difference, n, a24):
    """Share P's sum/difference for independent 2P and P+Q operations."""
    px, pz = point
    qx, qz = other
    rx, rz = difference
    total, delta = px + pz, px - pz
    aa, bb = total * total % n, delta * delta % n
    u, v = delta * (qx + qz) % n, total * (qx - qz) % n
    added_total, added_difference = u + v, u - v
    delta = aa - bb
    return (
        (aa * bb % n, delta * (bb + a24 * delta) % n),
        (
            rz * added_total * added_total % n,
            rx * added_difference * added_difference % n,
        ),
    )


def fusion_plan(operations):
    """Fuse only adjacent independent D/A pairs; retain both coverage masks."""
    plan, index = [], 0
    while index < len(operations):
        first = operations[index]
        if index + 1 < len(operations):
            second = operations[index + 1]
            double, add = (
                (first, second) if first[3] == DOUBLE else (second, first)
            )
            reads = (second[1],) if second[3] == DOUBLE else second[1:4]
            if (
                double[3] == DOUBLE
                and add[3] != DOUBLE
                and double[1] in add[1:3]
                and first[0] not in reads
                and first[0] != second[0]
            ):
                other = add[2] if double[1] == add[1] else add[1]
                plan.append(
                    (
                        2,
                        double[0],
                        double[1],
                        other,
                        add[3],
                        double[4],
                        add[0],
                        add[4],
                    )
                )
                index += 2
                continue
        dest, left, right, difference, mask = first
        plan.append(
            (
                int(difference != DOUBLE),
                dest,
                left,
                right,
                difference,
                mask,
                0,
                0,
            )
        )
        index += 1
    return tuple(plan)


class FusedRecord(c6_fast.FastRecord):
    """Independent adjacent-pair fusion over already certified operations."""

    def __init__(self, certified, backend):
        super().__init__(certified, backend)
        self.plan = fusion_plan(self.operations)
        self.run = self.fused

    def fused(self, point, n, a24):
        points = [None] * self.record.slots
        points[0] = point
        product = c6_fast._factor(1, point, self.masks[0])
        for kind, dest, left, right, diff, mask, adest, amask in self.plan:
            if kind == 0:
                result = point_double(*points[left], n, a24)
            elif kind == 1:
                result = point_add(
                    *points[left], *points[right], *points[diff], n
                )
            else:
                result, added = double_add(
                    points[left], points[right], points[diff], n, a24
                )
                product = c6_fast._factor(product, added, amask)
                points[adest] = added
            product = c6_fast._factor(product, result, mask)
            points[dest] = result
        return points[self.record.output], product


def build_program(family, mode, backend, batch):
    """Preserve strict recovery kernels; vary only the certified fast path."""
    if mode not in ("late", "reduced", "fused"):
        raise ValueError("unknown combined kernel")
    program = (
        c6_cf.build_program(2000, "tuple", backend, batch)
        if family == "cf"
        else c6_fast.build_program(2000, family, "tuple", backend, batch)
    )
    if mode == "late":
        return program
    reduced = SimpleNamespace(
        ecm=SimpleNamespace(point_add=point_add, point_double=point_double),
        gcd=backend.gcd,
    )
    entries = []
    for prime, power, action, unit in program.entries:
        certified = c6_fast.CertifiedRecord(action.record, action.masks)
        if mode == "fused":
            replacement = FusedRecord(certified, backend)
        else:
            replacement = c6_fast.FastRecord(certified, backend)
            # The strict property retains the original checked backend. The
            # bound interpreter receives only residue-equivalent fast kernels.
            runner = c6_fast.FastRecord(certified, reduced)
            replacement.run = runner.interpret
        entries.append((prime, power, replacement, unit))
    return c6_fast.Program(entries, backend, batch)
