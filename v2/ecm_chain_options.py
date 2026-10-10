"""Explicit Lucas/CF alternatives to the accepted batch-16 PRAC plan.

The default executor and its legacy identity remain in ecm_chains. Optional
plans reuse its certified actions, atomic reservations and strict recovery.
"""

import hashlib
import json
from types import SimpleNamespace

from . import arithmetic, ecm, prac, prime_sieve, utils
from . import ecm_chains as core

CF_CATALOG = core.CATALOG.with_name("b3_cf_records.json")
CF_SHA256 = "fd1e650180c1c576fe9bd95e64221f3d19bb0a730a772d356f0876d1e91284ae"
CF_BYTES = 62430


def identity(bound, backend, family):
    """Pin the optional family without changing legacy PRAC/Lucas IDs."""
    if family not in ("lucas", "cf"):
        raise ValueError("optional ECM family must be lucas or cf")
    digest = CF_SHA256 if family == "cf" else core.CATALOG_SHA256
    return core.identity(bound, backend) + f"/optional-{family}/{digest}/v1"


class ChainPlan(core.ChainPlan):
    """Use accepted scalar/frontier records with unchanged batch16 recovery."""

    def __init__(self, bound, backend, budget, family):
        utils.require_integer(bound, "chain bound", 2)
        if bound > 2000:
            raise ValueError("ECM chains support B1 <= 2000")
        self.identity = identity(bound, backend, family)
        self.bound, self.backend = bound, backend
        path, digest, size = (
            (CF_CATALOG, CF_SHA256, CF_BYTES)
            if family == "cf"
            else (core.CATALOG, core.CATALOG_SHA256, core.CATALOG_BYTES)
        )
        budget.consume(size)
        with path.open("rb") as stream:
            raw = stream.read(size + 1)
        if len(raw) != size or hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError("optional ECM catalog identity mismatch")
        data = json.loads(raw)
        if data["schema"] != 1 or data["bound"] != 2000:
            raise ValueError("incompatible optional ECM catalog")
        rows = data["families"][family]
        if len(rows) != 333:
            raise ValueError("incomplete optional ECM catalog")
        engine = SimpleNamespace(
            name=backend, ecm=ecm, prac=prac, gcd=arithmetic.gcd
        )
        actions, owned = {}, 4096
        for scalar, row in rows.items():
            scalar = int(scalar)
            code = bytes.fromhex(row["code"])
            masks = bytes.fromhex(row["masks"])
            steps = len(code) // 4
            # Keep the production verification/reservation contract even
            # though these records were independently certified in C6.
            budget.consume(4 * steps + row["slots"] + 1)
            record = core.Record(scalar, code, row["output"], row["slots"])
            actions[scalar] = core.Action(record, masks, engine)
            owned += 4096 + 512 * steps
        if owned > core.PLAN_BYTES:
            raise MemoryError("optional ECM plan exceeds its reserve")
        budget.consume(bound)
        self.entries = {}
        for prime in prime_sieve.prime_sieve(bound + 1):
            power = utils.prime_power(prime, bound)
            action, unit = actions[power], actions[prime]
            remaining, exponent = power, 0
            while remaining > 1:
                remaining //= prime
                exponent += 1
            reserve = (
                3 * len(action.operations)
                + 3
                + action.strict_work
                + exponent * unit.strict_work
            )
            self.entries[prime] = (power, action, unit, reserve)
        self.owned_bytes = owned


class ChainPlans(core.ChainPlans):
    """Keep optional records in the same finite, invocation-owned LRU."""

    def __init__(self, memory_bytes, backend, bounds, family):
        super().__init__(memory_bytes, backend, bounds)
        identity(2000, backend, family)
        self.family = family

    def get(self, bound, budget):
        if bound not in self.bounds or not 2 <= bound <= 2000:
            return None
        key = identity(bound, self.backend, self.family)
        budget.consume(1)
        if key in self.plans:
            self.hits += 1
            self.plans.move_to_end(key)
            return self.plans[key]
        while (
            self.plans
            and self.used_bytes + core.PLAN_BYTES > self.memory_bytes
        ):
            _, old = self.plans.popitem(last=False)
            self.used_bytes -= old.owned_bytes
            self.evictions += 1
            del old
        self.misses += 1
        plan = ChainPlan(bound, self.backend, budget, self.family)
        self.plans[key] = plan
        self.used_bytes += plan.owned_bytes
        return plan


def verify_progress(job, backend, verifier, family):
    """Verify optional identity and the unchanged batch16 prefix proof."""
    if "chain_identity" in job:
        if job["chain_identity"] != identity(job["b1"], backend, family):
            raise ValueError("incompatible optional ECM chain progress")
        job = dict(job, chain_identity=core.identity(job["b1"], backend))
    core.verify_progress(job, backend, verifier)
