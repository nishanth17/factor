"""C6 optimized, certificate-guarded chain execution; no production routing.

Factor coverage follows polynomial divisibility, independently of curve order.
Doubling carries both input-coordinate factors into output Z; differential
addition carries difference-Z factors into output X and difference-X factors
into output Z. A verified frontier therefore covers every old X/Z check.
"""

import json
from dataclasses import dataclass
from pathlib import Path

from ....common import prime_sieve, utils
from ....ecm import prac
from ...support.paths import (
    BENCHMARK_ROOT,
)
from . import c6_chains as reference

DATA = BENCHMARK_ROOT / "inputs/controls/c6_fast_records.json"
MODES = ("tuple", "calls", "inline")
MAX_SOURCE_BYTES = 128 * 1024
MAX_PROGRAM_SOURCE_BYTES = 8 * 1024 * 1024
MAX_BATCH = 64


def instructions(record):
    code = record.code
    result = []
    for start in range(0, len(code), 4):
        stop = start + 4
        result.append(tuple(code[start:stop]))
    return tuple(result)


def frontier(record):
    """Find coordinate leaves in the polynomial factor-propagation DAG."""
    reference.verify(record)
    versions = [None] * record.slots
    versions[0] = (0, 1)
    propagated = set()
    for index, (dest, left, right, difference) in enumerate(
        instructions(record), 1
    ):
        parents = (
            versions[left]
            if difference == reference.DOUBLE
            else versions[difference]
        )
        propagated.update(parents)
        versions[dest] = 2 * index, 2 * index + 1
    # The result coordinates are certified at the block boundary, or by
    # the next record's certificate. Avoid multiplying them at every prime.
    propagated.update(versions[record.output])
    masks = bytes(
        (
            int(2 * i not in propagated)
            | (int(2 * i + 1 not in propagated) << 1)
        )
        for i in range(len(record.code) // 4 + 1)
    )
    verify_frontier(record, masks)
    return masks


def verify_frontier(record, masks):
    """Independently prove coverage using forward ancestor bitsets.

    Guards together with the two output coordinates must cover every
    intermediate coordinate, including discarded points and the input.
    Checking the next record therefore also certifies this record's output.
    """
    reference.verify(record)
    count = len(record.code) // 4
    if type(masks) is not bytes or len(masks) != count + 1:
        raise ValueError("wrong coverage certificate size")
    if any(mask > 3 for mask in masks):
        raise ValueError("invalid coordinate mask")
    states = [None] * record.slots
    states[0] = (1, 2)
    covered = (1 if masks[0] & 1 else 0) | (2 if masks[0] & 2 else 0)
    for index, (dest, left, right, difference) in enumerate(
        instructions(record), 1
    ):
        x, z = 1 << (2 * index), 1 << (2 * index + 1)
        if difference == reference.DOUBLE:
            previous_x, previous_z = states[left]
            z |= previous_x | previous_z
        else:
            known_x, known_z = states[difference]
            x |= known_z
            z |= known_x
        if masks[index] & 1:
            covered |= x
        if masks[index] & 2:
            covered |= z
        states[dest] = x, z
    output_x, output_z = states[record.output]
    covered |= output_x | output_z
    if covered != (1 << (2 * (count + 1))) - 1:
        raise ValueError("certificate loses an intermediate factor")
    return True


def _factor(value, point, mask):
    if mask & 1:
        value *= point[0]
    if mask & 2:
        value *= point[1]
    return value


def generated(record, masks, backend, mode):
    """Compile bounded straight-line code from verified numeric operands.

    The inline mode preserves the pinned kernels' expressions and operation
    order exactly. No source fragment is accepted from a record or caller.
    """
    lines = ["def run(point, n, a24):", "    x0, z0 = point", "    guard = 1"]

    def guard(slot, mask):
        if mask & 1:
            lines.append(f"    guard *= x{slot}")
        if mask & 2:
            lines.append(f"    guard *= z{slot}")

    guard(0, masks[0])
    for index, (dest, left, right, difference) in enumerate(
        instructions(record), 1
    ):
        if mode == "calls":
            if difference == reference.DOUBLE:
                expression = f"double(x{left}, z{left}, n, a24)"
            else:
                expression = (
                    f"add(x{left}, z{left}, x{right}, z{right}, "
                    f"x{difference}, z{difference}, n)"
                )
            lines.append(f"    x{dest}, z{dest} = {expression}")
        elif difference == reference.DOUBLE:
            lines += [
                f"    total, difference = x{left} + z{left}, "
                f"x{left} - z{left}",
                "    sum_squared = total * total",
                "    difference_squared = difference * difference",
                "    delta = sum_squared - difference_squared",
                f"    x{dest}, z{dest} = "
                "(sum_squared * difference_squared % n,",
                "        delta * (difference_squared + a24 * delta) % n)",
            ]
        else:
            lines += [
                f"    u = (x{left} - z{left}) * (x{right} + z{right})",
                f"    v = (x{left} + z{left}) * (x{right} - z{right})",
                "    total, difference = u + v, u - v",
                f"    x{dest}, z{dest} = (z{difference} * total * total % n,",
                f"        x{difference} * difference * difference % n)",
            ]
        guard(dest, masks[index])
    lines.append(f"    return (x{record.output}, z{record.output}), guard")
    source = "\n".join(lines) + "\n"
    if len(source.encode()) > MAX_SOURCE_BYTES:
        raise ValueError("generated function exceeds source cap")
    namespace = dict(
        add=backend.ecm.point_add,
        double=backend.ecm.point_double,
        __builtins__={},
    )
    exec(compile(source, "<verified-c6-chain>", "exec"), namespace)
    return namespace["run"], len(source.encode())


@dataclass(frozen=True)
class CertifiedRecord:
    """An immutable scalar record whose factor coverage is checked once."""

    record: reference.Record
    masks: bytes

    def __post_init__(self):
        verify_frontier(self.record, self.masks)


class FastRecord:
    """A verified record, coverage proof and finite checked recovery."""

    def __init__(self, record, backend, mode="tuple", masks=None):
        if mode not in MODES:
            raise ValueError("unknown fast executor")
        if type(record) is CertifiedRecord:
            if masks is not None:
                raise ValueError("cannot replace a certified mask")
            certified = record
        else:
            certified = CertifiedRecord(
                record, frontier(record) if masks is None else masks
            )
        self.record, self.backend = certified.record, backend
        self.masks = certified.masks
        self._strict = None
        record = self.record
        self.operations = tuple(
            (*operation, mask)
            for operation, mask in zip(instructions(record), self.masks[1:])
        )
        self.source_bytes = 0
        self.run = self.interpret
        if mode != "tuple":
            self.run, self.source_bytes = generated(
                record, self.masks, backend, mode
            )

    @property
    def strict(self):
        if self._strict is None:
            self._strict = reference.Executor(self.record, self.backend)
        return self._strict

    def interpret(self, point, n, a24):
        points = [None] * self.record.slots
        points[0] = point
        product = _factor(1, point, self.masks[0])
        add = self.backend.ecm.point_add
        double = self.backend.ecm.point_double
        for dest, left, right, difference, mask in self.operations:
            if difference == reference.DOUBLE:
                result = double(*points[left], n, a24)
            else:
                result = add(
                    *points[left], *points[right], *points[difference], n
                )
            product = _factor(product, result, mask)
            points[dest] = result
        return points[self.record.output], product

    def __call__(self, point, n, a24):
        point = point[0] % n, point[1] % n
        result, product = self.run(point, n, a24)
        if self.backend.gcd(product * result[0] * result[1], n) == 1:
            return result
        return self.strict(point, n, a24)


class Program:
    """Batch coverage GCDs, with one finite strict replay per bad block."""

    def __init__(self, entries, backend, batch):
        utils.require_integer(batch, "batch", 1)
        if batch > MAX_BATCH or len(entries) > reference.MAX_RECORDS:
            raise ValueError("fast program exceeds bounds")
        unique = {id(row[2]): row[2] for row in entries}.values()
        if (
            sum(action.source_bytes for action in unique)
            > MAX_PROGRAM_SOURCE_BYTES
        ):
            raise ValueError("generated program exceeds source cap")
        self.entries = tuple(entries)
        self.backend, self.batch = backend, batch

    def __call__(self, point, n, a24, extra=None):
        if extra is None:
            extra = {}
        for name in (
            "guard_batches",
            "block_replays",
            "record_replay_allowance",
            "prime_power_replays",
            "prime_units_replayed",
        ):
            extra.setdefault(name, 0)
        point = point[0] % n, point[1] % n
        index, entries, gcd = 0, self.entries, self.backend.gcd
        while index < len(entries):
            start, original, product = index, point, 1
            stop = min(index + self.batch, len(entries))
            while index < stop:
                point, guard = entries[index][2].run(point, n, a24)
                product = product * guard % n
                index += 1
                if not product:
                    break
            extra["guard_batches"] += 1
            if gcd(product * point[0] * point[1], n) == 1:
                continue

            # Reconstruct the original block using the unchanged checked
            # implementation. It extracts any proper factor before a point
            # is discarded and performs its bounded ladder/prime-unit retry.
            # Even a saturated product cannot hide factors split across n.
            from .c6_study import stage_one

            extra["block_replays"] += 1
            extra["record_replay_allowance"] += index - start
            strict = tuple(
                (
                    prime,
                    power,
                    action.strict,
                    reference.Executor(unit.record, self.backend),
                )
                for prime, power, action, unit in entries[start:index]
            )
            point, factor = stage_one(
                original, n, a24, strict, self.backend, extra
            )
            if point is None or factor is not None:
                return point, factor
        return point, None


def ladder_chain(scalar):
    """An independent binary chain for same-executor arithmetic controls."""
    utils.require_integer(scalar, "scalar", 1)
    if scalar >= 2**32:
        raise ValueError("scalar exceeds record cap")
    if scalar == 1:
        return prac.Chain(1, (), 0, "binary")
    operations = [(0, 0, -1)]
    q, r = 0, 1
    for bit in bin(scalar)[3:]:
        operations.append((q, r, 0))
        total = len(operations)
        doubled = r if bit == "1" else q
        operations.append((doubled, doubled, -1))
        if bit == "1":
            q, r = total, len(operations)
        else:
            q, r = len(operations), total
    chain = prac.Chain(scalar, tuple(operations), q, "binary")
    prac.verify_chain(chain)
    return chain


def generate_catalog():
    """Precompute all required prime powers <=2000, never full-lcm search."""
    lucas = reference.load_lucas()
    data = dict(schema=1, bound=2000, families={})
    for family in ("ladder", "prac", "lucas"):
        records = {}
        for prime in prime_sieve.prime_sieve(2001):
            power = prime
            while power <= 2000:
                if family == "lucas":
                    chain = reference.compose(lucas[prime], power)
                elif family == "prac":
                    chain = prac.get_chain(power)
                else:
                    chain = ladder_chain(power)
                record = reference.compact(chain)
                masks = frontier(record)
                records[str(power)] = dict(
                    code=record.code.hex(),
                    output=record.output,
                    slots=record.slots,
                    masks=masks.hex(),
                )
                power *= prime
        if len(records) > reference.MAX_RECORDS:
            raise ValueError("too many precomputed records")
        data["families"][family] = records
    prac.clear_cache()
    return data


def load_catalog(path=DATA, families=None):
    with Path(path).open("rb") as stream:
        raw = stream.read(reference.MAX_FILE_BYTES + 1)
    if len(raw) > reference.MAX_FILE_BYTES:
        raise ValueError("precomputed catalog exceeds file cap")
    data = json.loads(raw)
    if data["schema"] != 1 or data["bound"] != 2000:
        raise ValueError("unknown precomputed catalog")
    if set(data["families"]) != {"ladder", "prac", "lucas"}:
        raise ValueError("wrong chain families")
    result = {}
    for family, rows in data["families"].items():
        if families is not None and family not in families:
            continue
        if len(rows) > reference.MAX_RECORDS:
            raise ValueError("too many records")
        result[family] = {}
        for scalar, row in rows.items():
            record = reference.Record(
                int(scalar),
                bytes.fromhex(row["code"]),
                row["output"],
                row["slots"],
            )
            masks = bytes.fromhex(row["masks"])
            result[family][int(scalar)] = CertifiedRecord(record, masks)
    return result


def build_program(bound, family, mode, backend, batch=16, catalog=None):
    utils.require_integer(bound, "bound", 2)
    if bound > reference.MAX_BOUND:
        raise ValueError("fast C6 supports B1 <= 2000")
    if catalog is None:
        catalog = load_catalog(families=(family,))
    if family not in catalog:
        raise ValueError("unknown precomputed chain family")
    records, actions, units, entries = catalog[family], {}, {}, []
    for prime in prime_sieve.prime_sieve(bound + 1):
        power = utils.prime_power(prime, bound)
        if power not in actions:
            actions[power] = FastRecord(records[power], backend, mode)
        if prime not in units:
            units[prime] = records[prime]
        entries.append((prime, power, actions[power], units[prime]))
    return Program(entries, backend, batch)
