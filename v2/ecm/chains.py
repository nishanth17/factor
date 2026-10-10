"""Finite run-owned C6 PRAC/Lucas plans for atomic production chunks."""

import hashlib
import json
from collections import OrderedDict
from pathlib import Path
from types import SimpleNamespace

from ..common import arithmetic, prime_sieve, utils
from . import core as ecm
from . import prac
from .chain_records import (
    DOUBLE,
    Executor,
    Record,
    _factor,
    instructions,
    verify_frontier,
)

CATALOG = (
    Path(__file__).parents[1]
    / "benchmarks/inputs/controls/c6_fast_records.json"
)
CATALOG_SHA256 = (
    "b316003a7238c1240f52dd28f48f507931ef7f4061362333dc436f016d531ed9"
)
CATALOG_BYTES = 202461
CHAIN_VERSION = "ecm-c6-chunks-v1/" + CATALOG_SHA256
BATCH = 16
SCRATCH_BYTES = 4 * 1024**2
PLAN_BYTES = 4 * 1024**2
MIN_MEMORY_BYTES = SCRATCH_BYTES + PLAN_BYTES
MIN_MODULUS = 10**39
MAX_MODULUS = 10**80


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


def supports_modulus(n):
    """Keep inputs outside the C6 measured size band on the B4 ladder."""
    return MIN_MODULUS <= n < MAX_MODULUS


def identity(bound, backend):
    """Include schedule, catalog, kernel, guard and recovery policy."""
    engine = arithmetic.get_backend(backend)
    family = "prac-reduced" if backend == "python-int" else "lucas-tuple"
    return f"{CHAIN_VERSION}/{bound}/{family}/{engine.identity}/batch16"


class Action:
    """Own certified operations and an eagerly verified strict interpreter."""

    def __init__(self, record, masks, backend):
        verify_frontier(record, masks)
        self.record, self.masks = record, masks
        self.operations = tuple(
            (*operation, mask)
            for operation, mask in zip(instructions(record), masks[1:])
        )
        self.strict = Executor(record, backend)
        native = backend.name == "python-int"
        self.add = point_add if native else ecm.point_add
        self.double = point_double if native else ecm.point_double

    def run(self, point, n, a24):
        points = [None] * self.record.slots
        points[0] = point
        product = _factor(1, point, self.masks[0])
        add, double = self.add, self.double
        for dest, left, right, difference, mask in self.operations:
            if difference == DOUBLE:
                result = double(*points[left], n, a24)
            else:
                result = add(
                    *points[left], *points[right], *points[difference], n
                )
            product = _factor(product, result, mask)
            points[dest] = result
        return points[self.record.output], product

    @property
    def strict_work(self):
        # Covers point operations, both coordinate GCDs and one checked
        # ladder retry; these conservative work units are not elapsed time.
        return (
            5 * len(self.operations) + 10 * self.record.scalar.bit_length() + 8
        )


class ChainPlan:
    """Own a bound plan without retaining curve points or pending products."""

    def __init__(self, bound, backend, budget):
        utils.require_integer(bound, "chain bound", 2)
        if bound > 2000:
            raise ValueError("ECM chains support B1 <= 2000")
        self.identity = identity(bound, backend)
        self.bound, self.backend = bound, backend
        family = "prac" if backend == "python-int" else "lucas"
        engine = SimpleNamespace(
            name=backend,
            ecm=ecm,
            prac=prac,
            gcd=arithmetic.gcd,
        )
        # Reading, hashing and parsing are paid even if verification fails.
        budget.consume(CATALOG_BYTES)
        with CATALOG.open("rb") as stream:
            raw = stream.read(CATALOG_BYTES + 1)
        if (
            len(raw) != CATALOG_BYTES
            or hashlib.sha256(raw).hexdigest() != CATALOG_SHA256
        ):
            raise ValueError("ECM chain catalog identity mismatch")
        data = json.loads(raw)
        if data["schema"] != 1 or data["bound"] != 2000:
            raise ValueError("incompatible ECM chain catalog")
        rows = data["families"][family]
        if len(rows) != 333:
            raise ValueError("incomplete ECM chain catalog")
        actions = {}
        owned = 4096
        for scalar, row in rows.items():
            scalar = int(scalar)
            code, masks = (
                bytes.fromhex(row["code"]),
                bytes.fromhex(row["masks"]),
            )
            steps = len(code) // 4
            # Scalar proof, independent coverage, strict verification and
            # decoded tuples all consume the shared cumulative allowance.
            budget.consume(4 * steps + row["slots"] + 1)
            record = Record(scalar, code, row["output"], row["slots"])
            actions[scalar] = Action(record, masks, engine)
            owned += 4096 + 512 * steps
        if owned > PLAN_BYTES:
            raise MemoryError("verified ECM chain plan exceeds its reserve")
        budget.consume(bound)
        primes = prime_sieve.prime_sieve(bound + 1)
        self.entries = {}
        for prime in primes:
            power = utils.prime_power(prime, bound)
            action, unit = actions[power], actions[prime]
            remaining, exponent = power, 0
            while remaining > 1:
                remaining //= prime
                exponent += 1
            fast_work = 3 * len(action.operations) + 3
            reserve = (
                fast_work + action.strict_work + exponent * unit.strict_work
            )
            self.entries[prime] = (power, action, unit, reserve)
        self.owned_bytes = owned

    def execute(self, powers, point, n, a24, budget):
        """Certify at most 16 records or replay their saved input once."""
        if not 0 < len(powers) <= BATCH:
            raise ValueError("chain chunk exceeds batch16")
        entries = []
        previous = 1
        for prime, power in powers:
            if prime <= previous or prime not in self.entries:
                raise ValueError("invalid chain chunk schedule")
            entry = self.entries[prime]
            if power != entry[0]:
                raise ValueError("chain chunk power disagrees with bound")
            entries.append((prime, *entry))
            previous = prime
        # No arithmetic or mutation precedes this reservation. Recovery may
        # inspect every record and every repeated-prime unit, including a
        # checked ladder retry. Unused recovery credit is not refunded.
        budget.consume(1 + sum(entry[4] for entry in entries))
        original, product = tuple(point), 1
        completed = 0
        for _, _, action, _, _ in entries:
            point, guard = action.run(point, n, a24)
            product = product * guard % n
            completed += 1
            if not product:
                break
        if arithmetic.gcd(product * point[0] * point[1], n) == 1:
            return point, None, False

        # The aggregate can equal n with different coordinate factors. Do
        # not treat that as failure before strict replay has inspected them.
        point = original
        for prime, power, action, unit, _ in entries[:completed]:
            saved = point
            try:
                point = action.strict(point, n, a24)
                if arithmetic.gcd(point[1], n) != n:
                    continue
            except prac.NonunitPointError as result:
                if result.factor is not None:
                    return None, result.factor, True
            point, remaining = saved, power
            while remaining > 1:
                try:
                    point = unit.strict(point, n, a24)
                except prac.NonunitPointError as result:
                    return None, result.factor, True
                if arithmetic.gcd(point[1], n) == n:
                    return None, None, True
                remaining //= prime
        # A unit-certified prefix that ended early cannot certify the tail.
        # Returning None requests the existing saved-chunk replay path.
        return (point if completed == len(entries) else None), None, True


class ChainPlans:
    """Own an LRU with simultaneous construction/recovery scratch reserved.

    Evict before allocating a replacement. A refusal publishes no partial
    plan; already charged verification stays consumed. Resume starts empty.
    """

    def __init__(self, memory_bytes, backend, bounds):
        utils.require_integer(memory_bytes, "chain memory", MIN_MEMORY_BYTES)
        self.memory_bytes, self.backend = memory_bytes, backend
        self.bounds = frozenset(bounds)
        self.plans = OrderedDict()
        self.used_bytes = SCRATCH_BYTES
        self.hits = self.misses = self.evictions = 0

    def get(self, bound, budget):
        if bound not in self.bounds or not 2 <= bound <= 2000:
            return None
        key = identity(bound, self.backend)
        budget.consume(1)
        if key in self.plans:
            self.hits += 1
            self.plans.move_to_end(key)
            return self.plans[key]
        while self.plans and self.used_bytes + PLAN_BYTES > self.memory_bytes:
            _, old = self.plans.popitem(last=False)
            self.used_bytes -= old.owned_bytes
            self.evictions += 1
            # Dropping the cache key alone leaves this local owning the old
            # plan while its replacement is allocated. Release both owners.
            del old
        self.misses += 1
        plan = ChainPlan(bound, self.backend, budget)
        self.plans[key] = plan
        self.used_bytes += plan.owned_bytes
        return plan


def verify_progress(job, backend, verifier):
    """Check the certified chunk prefix before rehydrating a chain snapshot."""
    phase = job["phase"]
    expected_chunks = 19
    if phase in ("setup", "stage_one", "replay"):
        cursor, powers = job["cursor"], job["powers"]
        if len(powers) > BATCH:
            raise ValueError("checkpoint chain chunk exceeds batch16")
        consumed = list(verifier.primes(2, cursor["left"]))
        consumed += cursor["values"][: cursor["index"]]
        if len(powers) > len(consumed):
            raise ValueError("chain chunk exceeds consumed prime prefix")
        first_pending = len(consumed) - len(powers)
        expected = [
            [prime, utils.prime_power(prime, job["b1"])]
            for prime in consumed[first_pending:]
        ]
        if powers != expected:
            raise ValueError("checkpoint chain powers disagree with schedule")
        committed = len(consumed) - len(powers)
        expected_chunks = (committed + BATCH - 1) // BATCH + int(
            phase == "replay"
        )
    if expected_chunks:
        if job.get("chain_identity") != identity(job["b1"], backend):
            raise ValueError("missing certified chain prefix identity")
        if job.get("chain_chunks") != expected_chunks:
            raise ValueError("chain cursor disagrees with committed chunks")
    elif "chain_identity" in job:
        raise ValueError("unexpected certified chain prefix")
    points = [] if phase == "setup" else [job["value"]]
    if phase == "replay":
        index, power = job["replay_index"], job["replay_power"]
        utils.require_integer(index, "chain replay index", 0)
        utils.require_integer(power, "chain replay power", 1)
        if index > len(job["powers"]):
            raise ValueError("chain replay index exceeds chunk")
        if index == len(job["powers"]):
            if power != 1:
                raise ValueError("completed chain replay has pending power")
        else:
            prime, target = job["powers"][index]
            if power >= target or target % power:
                raise ValueError("chain replay power exceeds pending scalar")
            remaining = power
            while remaining > 1 and remaining % prime == 0:
                remaining //= prime
            if remaining != 1:
                raise ValueError("chain replay power has an unrelated prime")
        points.append(job["replay_value"])
    for point in points:
        if not isinstance(point, list) or len(point) != 2:
            raise ValueError("invalid chain checkpoint point")
        for coordinate in point:
            utils.require_integer(coordinate, "chain coordinate", 0)
            if coordinate >= job["n"]:
                raise ValueError("noncanonical chain checkpoint coordinate")
        if arithmetic.gcd(point[1], job["n"]) != 1:
            raise ValueError("uncertified chain checkpoint continuation")
