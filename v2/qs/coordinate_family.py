"""Exact coordinate view using the production root recurrence unchanged."""

import json

from ..common import arithmetic, utils
from .families import (
    PolynomialFamily,
    _checked_resources,
    _checksum,
    _identity,
)
from .sign_schedule import SignSchedule

FAMILY_SCHEMA = "c10a07-coordinate-family-v1"
CHECKPOINT_BYTES = 8192
HELPER_BYTES = 65536


class CoordinateFamily:
    """Permute/negate CRT coordinates once, then use ordinary adjacent steps.

    The private kernel keeps matching term/prime order. Public A primes remain
    sorted for assignment identity and checkpoint validation. A translated face
    changes initial signs, not the initial ordinal or the root recurrence.
    """

    def __init__(
        self,
        base,
        a_primes,
        *,
        schedule,
        budget,
        memory_bytes,
        parity=False,
        root_update_policy="standard",
        root_update_limit=None,
    ):
        if not isinstance(
            schedule, SignSchedule
        ) or schedule.factor_count != len(a_primes):
            raise ValueError("schedule does not match this A family")
        extra = 1024 + 512 * len(a_primes)
        if extra > HELPER_BYTES or memory_bytes < extra:
            raise MemoryError("coordinate view exceeds its memory allowance")
        budget.consume(len(a_primes) + 1)
        kernel = PolynomialFamily(
            base,
            a_primes,
            budget=budget,
            memory_bytes=memory_bytes - extra,
            parity=parity,
        )
        order = (0,) + schedule.variable_order
        terms = kernel.terms
        kernel.terms = tuple(
            -terms[i] if schedule.negative_mask & (1 << i) else terms[i]
            for i in order
        )
        kernel.a_primes = tuple(a_primes[i] for i in order)
        kernel.identity = _checksum([kernel.identity, schedule.mapping()])
        kernel._configure_root_updates(root_update_policy, root_update_limit)
        kernel.workspace_bytes += extra
        self._kernel = kernel
        self.a_primes = a_primes
        self.schedule = schedule

    def __getattr__(self, name):
        if name == "_kernel":
            raise AttributeError(name)
        return getattr(self._kernel, name)

    @property
    def next_index(self):
        return self._kernel.next_index

    @next_index.setter
    def next_index(self, value):
        self._kernel.next_index = value

    @property
    def current(self):
        return self._kernel.current

    @current.setter
    def current(self, value):
        self._kernel.current = value

    @property
    def budget(self):
        return self._kernel.budget

    @budget.setter
    def budget(self, value):
        self._kernel.budget = value

    def checkpoint(self):
        payload = dict(
            schema="c10a01-coordinate-family-v2"
            if self.parity
            else FAMILY_SCHEMA,
            base_identity=self.base_identity,
            a_primes=list(self.a_primes),
            schedule=self.schedule.mapping(),
            next_index=self.next_index,
            resources=dict(
                work_used=self.budget.used,
                wall_used=self.budget.wall_used,
                cpu_used=self.budget.cpu_used,
            ),
        )
        if self.parity:
            payload["parity"] = True
        if self.root_update_policy != "standard":
            payload["schema"] = "c10a01-coordinate-family-v3"
            payload["backend"] = arithmetic.backend_for(self.base.n).identity
            payload["root_update_policy"] = self.root_update_policy
            payload["root_update_limit"] = self.root_update_limit
        _checked_resources(payload["resources"])
        if len(json.dumps(payload).encode()) > CHECKPOINT_BYTES:
            raise MemoryError("coordinate family checkpoint exceeds its cap")
        return dict(payload=payload, sha256=_checksum(payload))

    @classmethod
    def from_checkpoint(
        cls, base, checkpoint, *, budget, memory_bytes, expected_schedule=None
    ):
        if not isinstance(checkpoint, dict) or set(checkpoint) != {
            "payload",
            "sha256",
        }:
            raise ValueError("invalid coordinate family checkpoint")
        payload = checkpoint["payload"]
        if not isinstance(payload, dict) or set(payload) - {
            "parity",
            "root_update_policy",
            "root_update_limit",
            "backend",
        } != {
            "schema",
            "base_identity",
            "a_primes",
            "schedule",
            "next_index",
            "resources",
        }:
            raise ValueError("invalid coordinate family payload")
        schedule = SignSchedule.from_mapping(payload["schedule"])
        if expected_schedule is not None and schedule != expected_schedule:
            raise ValueError("saved sign schedule must be preserved")
        primes = payload["a_primes"]
        if (
            not isinstance(primes, list)
            or len(primes) != schedule.factor_count
        ):
            raise ValueError("invalid coordinate family A primes")
        for prime in primes:
            utils.require_integer(prime, "saved A prime", 3)
            if prime not in base.primes:
                raise ValueError("saved A prime is outside the base")
        index = utils.require_integer(payload["next_index"], "next_index", 0)
        if index > 1 << (len(primes) - 1):
            raise ValueError("coordinate cursor exceeds its family")
        resources = _checked_resources(payload["resources"])
        parity = payload.get("parity", False)
        schema = "c10a01-coordinate-family-v2" if parity else FAMILY_SCHEMA
        if type(parity) is not bool or (not parity and "parity" in payload):
            raise ValueError("invalid coordinate parity policy")
        policy = payload.get("root_update_policy", "standard")
        limit = payload.get("root_update_limit")
        fields = {"root_update_policy", "root_update_limit"}
        if fields.intersection(payload):
            if not fields.issubset(payload) or policy == "standard":
                raise ValueError("invalid coordinate root policy")
            schema = "c10a01-coordinate-family-v3"
            if (
                payload.get("backend")
                != arithmetic.backend_for(base.n).identity
            ):
                raise ValueError("coordinate root-table backend mismatch")
        elif "backend" in payload:
            raise ValueError("legacy coordinate schema cannot add a backend")
        if payload["schema"] != schema or payload[
            "base_identity"
        ] != _identity(base):
            raise ValueError("coordinate family target/convention mismatch")
        if len(json.dumps(payload).encode()) > CHECKPOINT_BYTES or checkpoint[
            "sha256"
        ] != _checksum(payload):
            raise ValueError("coordinate family integrity/size mismatch")
        if (
            budget.used < resources["work_used"]
            or budget.wall_used < resources["wall_used"]
            or budget.cpu_used < resources["cpu_used"]
        ):
            raise ValueError("resume must retain charged resources")
        family = cls(
            base,
            tuple(primes),
            schedule=schedule,
            parity=parity,
            root_update_policy=policy,
            root_update_limit=limit,
            budget=budget,
            memory_bytes=memory_bytes,
        )
        if index > family.root_update_limit:
            raise ValueError("coordinate cursor exceeds root-update prefix")
        if index:
            family.current = family._make_step(index - 1, None)
        family.next_index = index
        return family
