"""Exact coordinate view using the production root recurrence unchanged."""

import json

from ..common import utils
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

    def __init__(self, base, a_primes, *, schedule, budget, memory_bytes):
        if not isinstance(
            schedule, SignSchedule
        ) or schedule.factor_count != len(a_primes):
            raise ValueError("schedule does not match this A family")
        extra = 1024 + 512 * len(a_primes)
        if extra > HELPER_BYTES or memory_bytes < extra:
            raise MemoryError("coordinate view exceeds its memory allowance")
        budget.consume(len(a_primes) + 1)
        kernel = PolynomialFamily(
            base, a_primes, budget=budget, memory_bytes=memory_bytes - extra
        )
        order = (0,) + schedule.variable_order
        terms = kernel.terms
        kernel.terms = tuple(
            -terms[i] if schedule.negative_mask & (1 << i) else terms[i]
            for i in order
        )
        kernel.a_primes = tuple(a_primes[i] for i in order)
        kernel.identity = _checksum([kernel.identity, schedule.mapping()])
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
            schema=FAMILY_SCHEMA,
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
        if not isinstance(payload, dict) or set(payload) != {
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
        if payload["schema"] != FAMILY_SCHEMA or payload[
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
            budget=budget,
            memory_bytes=memory_bytes,
        )
        if index:
            family.current = family._make_step(index - 1, None)
        family.next_index = index
        return family
