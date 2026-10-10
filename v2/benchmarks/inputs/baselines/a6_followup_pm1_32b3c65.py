"""Opt-in finite p-1 campaigns with exact, verified bound continuation."""

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from math import gcd, isqrt, prod

from . import constants, utils
from .budget import Budget, BudgetExhaustedError
from .factor import FactorizationResult
from .schedules import RATIO_VERSION, SieveContext, prime_power_ratio
from .stage_jobs import prime_cursor

CHECKPOINT_VERSION = 1
EXECUTION_VERSION = "pm1-campaign-v1"


@dataclass(frozen=True)
class PM1Config:
    """Predeclare inclusive (B1, B2) rungs for one explicit base.

    Raising B1 invalidates stage-two arithmetic. Equal-B1 rungs append only
    the new B2 interval. Saturation stops this base, including later rungs.
    Work units belong to this API and are distinct from portfolio job units.
    """

    bounds: tuple = ((constants.PM1_B1, constants.PM1_B2),)
    chunk_size: int = 16
    gcd_batch: int = constants.GCD_BATCH_SIZE
    recovery_limit: int = constants.RHO_RECOVERY_LIMIT
    segment_size: int = 1024
    memory_bytes: int = 8_388_608
    max_input_bits: int = 4096

    def __post_init__(self):
        for name, minimum in (
            ("chunk_size", 1),
            ("gcd_batch", 1),
            ("recovery_limit", 0),
            ("segment_size", 1),
            ("memory_bytes", 65_536),
            ("max_input_bits", 2),
        ):
            utils.require_integer(getattr(self, name), name, minimum)
        if max(self.chunk_size, self.gcd_batch) > 256:
            raise ValueError("arithmetic chunks must not exceed 256")
        bounds = tuple(tuple(pair) for pair in self.bounds)
        if not 1 <= len(bounds) <= 64:
            raise ValueError("declare between one and 64 bound rungs")
        previous = (1, 1)
        for pair in bounds:
            if len(pair) != 2:
                raise ValueError("rungs need B1 and B2")
            b1, b2 = pair
            utils.require_integer(b1, "B1", 2)
            utils.require_integer(b2, "B2", b1)
            if b1 < previous[0] or b2 < previous[1] or pair == previous:
                raise ValueError("campaign bounds must increase monotonically")
            previous = pair
        object.__setattr__(self, "bounds", bounds)
        if self.memory_bytes - self.workspace_reserve < 8192:
            raise MemoryError("p-1 state exceeds configured workspace cap")

    @property
    def workspace_reserve(self):
        """Conservatively include replay, cache, bigint and JSON copies."""
        coordinate = 128 + self.max_input_bits // 8
        return (
            16_384
            + 256 * self.segment_size
            + 16 * coordinate * (self.chunk_size + self.gcd_batch + 64)
            + 256 * self.chunk_size * self.bounds[-1][0].bit_length()
        )


@dataclass(frozen=True)
class PM1Run:
    """A validated split or unresolved input, plus cumulative resources."""

    result: FactorizationResult
    divisor: int | None
    reason: str
    work_used: int
    wall_seconds: float
    cpu_seconds: float
    verification_work: int
    checkpoint: dict


def _initial(n, base, config):
    b1, b2 = config.bounds[0]
    return {
        "n": n,
        "base": base,
        "rung": 0,
        "b1": b1,
        "b2": b2,
        "old_b1": 1,
        "phase": "setup",
        "value": base % n,
        "cursor": prime_cursor(2, b1 + 1),
        "powers": [],
        "terms": [],
        "product": 1,
        "recovery": 0,
        "steps": 0,
        "execution_work": 0,
        "done": False,
        "factor": None,
        "reason": None,
    }


def _finish(state, reason, divisor=None):
    if divisor is not None and not utils.valid_divisor(divisor, state["n"]):
        raise ValueError("p-1 produced an invalid divisor")
    state.update(done=True, factor=divisor, reason=reason)


def _load_segment(state, context, budget, config):
    cursor = state["cursor"]
    left = cursor["next"]
    right = min(cursor["hi"], left + 2 * config.segment_size)
    budget.consume(config.segment_size + len(context.base_primes))
    values = context.prime_segment(left, right)
    cursor.update(left=left, next=right, values=values, index=0)


def _next_rung(state, budget, config):
    budget.consume()
    rung = state["rung"] + 1
    if rung == len(config.bounds):
        _finish(state, "exhausted")
        return
    old_b1, old_b2 = state["b1"], state["b2"]
    b1, b2 = config.bounds[rung]
    state.update(rung=rung, b1=b1, b2=b2, recovery=0)
    if b1 == old_b1:
        # All prior batches were checked. Preserve a**last_prime and its
        # residue-specific gap cache only while the stage-one value is fixed.
        state.update(
            phase="stage_two", cursor=prime_cursor(old_b2 + 1, b2 + 1)
        )
    else:
        state.update(
            old_b1=old_b1,
            phase="fill",
            cursor=prime_cursor(2, b1 + 1),
            powers=[],
            previous_prime=0,
            stage_two_value=1,
            gap_powers={},
        )


def _step(state, budget, context, config):
    """Reserve one finite action before changing any serializable state."""
    n, phase = state["n"], state["phase"]
    if phase == "setup":
        budget.consume(n.bit_length())
        divisor = gcd(state["base"], n)
        if utils.valid_divisor(divisor, n):
            _finish(state, "factor_found", divisor)
        elif divisor != 1 or n == 2:
            _finish(state, "nonunit" if divisor != 1 else "exhausted")
        elif n % 2 == 0:
            _finish(state, "factor_found", 2)
        else:
            state["phase"] = "fill"
        return

    if phase == "fill":
        cursor = state["cursor"]
        if len(state["powers"]) == config.chunk_size:
            budget.consume()
            state["phase"] = "power"
        elif cursor["index"] < len(cursor["values"]):
            count = min(
                config.chunk_size - len(state["powers"]),
                len(cursor["values"]) - cursor["index"],
            )
            budget.consume(count)
            start = cursor["index"]
            stop = start + count
            for prime in cursor["values"][start:stop]:
                ratio = prime_power_ratio(prime, state["old_b1"], state["b1"])
                if ratio != 1:
                    state["powers"].append([prime, ratio])
            cursor["index"] += count
        elif cursor["next"] < cursor["hi"]:
            _load_segment(state, context, budget, config)
        else:
            budget.consume()
            if state["powers"]:
                state["phase"] = "power"
            else:
                state.update(
                    phase="stage_two",
                    cursor=prime_cursor(state["b1"] + 1, state["b2"] + 1),
                    previous_prime=0,
                    stage_two_value=1,
                    gap_powers={},
                )
        return

    if phase == "power":
        # Reserve a conservative exponent-bit allowance before even building
        # the chunk product. A refused action retains the original residue.
        budget.consume(
            1 + sum(power.bit_length() for _, power in state["powers"])
        )
        scalar = prod(power for _, power in state["powers"])
        value = pow(state["value"], scalar, n)
        divisor = gcd(value - 1, n)
        if utils.valid_divisor(divisor, n):
            _finish(state, "factor_found", divisor)
        elif divisor == n:
            state.update(
                phase="replay",
                replay_value=state["value"],
                replay_index=0,
                replay_power=1,
            )
        else:
            state.update(value=value, powers=[], phase="fill")
        return

    if phase == "replay":
        position = state["replay_index"]
        if state["recovery"] >= config.recovery_limit or position == len(
            state["powers"]
        ):
            budget.consume()
            _finish(state, "saturated")
            return
        prime, power = state["powers"][position]
        budget.consume(prime.bit_length() + 1)
        state["replay_value"] = pow(state["replay_value"], prime, n)
        state["replay_power"] *= prime
        state["recovery"] += 1
        divisor = gcd(state["replay_value"] - 1, n)
        if utils.valid_divisor(divisor, n):
            _finish(state, "factor_found", divisor)
        elif divisor == n:
            _finish(state, "saturated")
        elif state["replay_power"] == power:
            state.update(replay_index=position + 1, replay_power=1)
        return

    if phase == "gcd":
        budget.consume()
        divisor = gcd(state["product"], n)
        if utils.valid_divisor(divisor, n):
            _finish(state, "factor_found", divisor)
        elif divisor == n:
            state.update(phase="term_replay", recovery_index=0)
        else:
            state.update(terms=[], product=1, phase="stage_two")
        return

    if phase == "term_replay":
        position = state["recovery_index"]
        if state["recovery"] >= config.recovery_limit or position == len(
            state["terms"]
        ):
            budget.consume()
            _finish(state, "saturated")
            return
        budget.consume()
        divisor = gcd(state["terms"][position], n)
        state["recovery_index"] += 1
        state["recovery"] += 1
        if utils.valid_divisor(divisor, n):
            _finish(state, "factor_found", divisor)
        return

    if phase != "stage_two":
        raise ValueError("unknown p-1 execution phase")
    cursor = state["cursor"]
    if len(state["terms"]) == config.gcd_batch:
        budget.consume()
        state["phase"] = "gcd"
    elif cursor["index"] < len(cursor["values"]):
        count = min(
            config.gcd_batch - len(state["terms"]),
            len(cursor["values"]) - cursor["index"],
        )
        start = cursor["index"]
        stop = start + count
        primes = cursor["values"][start:stop]
        previous = state["previous_prime"]
        gaps = []
        for prime in primes:
            gaps.append(prime - previous)
            previous = prime
        # Cache hits change CPU cost, never the allowance or resume identity.
        budget.consume(sum(gap.bit_length() + 1 for gap in gaps))
        cache = state["gap_powers"]
        residue, value = state["stage_two_value"], state["value"]
        for gap in gaps:
            key = str(gap)
            if key not in cache:
                if len(cache) == 64:
                    cache.clear()
                cache[key] = pow(value, gap, n)
            residue = residue * cache[key] % n
            term = (residue - 1) % n
            state["terms"].append(term)
            state["product"] = state["product"] * term % n
        cursor["index"] += count
        state.update(previous_prime=previous, stage_two_value=residue)
    elif cursor["next"] < cursor["hi"]:
        _load_segment(state, context, budget, config)
    elif state["terms"]:
        budget.consume()
        state["phase"] = "gcd"
    else:
        _next_rung(state, budget, config)


def _advance(state, budget, context, config):
    before = budget.used
    _step(state, budget, context, config)
    state["execution_work"] += budget.used - before
    state["steps"] += 1


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _restore(checkpoint, n, base, config):
    """Check finite metadata before charged deterministic reconstruction."""
    try:
        encoded = _canonical(checkpoint["payload"])
        if len(encoded.encode()) > config.memory_bytes // 2:
            raise ValueError("checkpoint exceeds storage cap")
        if (
            hashlib.sha256(encoded.encode()).hexdigest()
            != checkpoint["sha256"]
        ):
            raise ValueError("checkpoint checksum mismatch")
        payload = json.loads(encoded)
        if (
            type(payload["version"]) is not int
            or payload["version"] != CHECKPOINT_VERSION
            or payload["execution"] != EXECUTION_VERSION
            or payload["backend"] != "python-int"
            or payload["schedule"] != RATIO_VERSION
            or _canonical(payload["config"]) != _canonical(asdict(config))
            or type(payload["n"]) is not int
            or payload["n"] != n
            or type(payload["base"]) is not int
            or payload["base"] != base
        ):
            raise ValueError("incompatible p-1 checkpoint")
        utils.require_integer(payload["work_used"], "work_used", 0)
        for name in ("wall_used", "cpu_used"):
            value = payload[name]
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0
            ):
                raise ValueError("invalid cumulative p-1 time")
        overhead = checkpoint["serialization_overhead"]
        if len(overhead) != 2 or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0
            for value in overhead
        ):
            raise ValueError("invalid serialization time")
        if (
            hashlib.sha256(_canonical(overhead).encode()).hexdigest()
            != checkpoint["overhead_sha256"]
        ):
            raise ValueError("serialization time checksum mismatch")
        payload["wall_used"] += overhead[0]
        payload["cpu_used"] += overhead[1]
        state = payload["state"]
        for name in ("steps", "execution_work"):
            utils.require_integer(state[name], name, 0)
        if (
            not state["steps"]
            <= state["execution_work"]
            <= payload["work_used"]
        ):
            raise ValueError("invalid p-1 progress allowance")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("malformed p-1 checkpoint") from error
    return payload


def factorize_pm1_bounded(
    n, *, base=2, config=None, budget=None, checkpoint=None, max_actions=None
):
    """Run a finite one-base campaign; resume only an identical assignment.

    Resume reconstructs all committed actions before reusing arithmetic and
    charges that verification to the cumulative budget. Small repeated grants
    can be spent entirely on verification; allow replay plus progress.
    Changing bounds/base/config requires a fresh explicitly declared campaign.
    No primality claim is made: even a successful split remains unclassified.
    """
    utils.require_integer(n, "n", 2)
    utils.require_integer(base, "base", 2)
    if type(n) is not int or type(base) is not int:
        raise TypeError("bounded p-1 uses Python integers")
    config = PM1Config() if config is None else config
    budget = Budget() if budget is None else budget
    if max(n.bit_length(), base.bit_length()) > config.max_input_bits:
        raise ValueError("input or base exceeds configured bit cap")
    if budget.used or budget.prior_wall or budget.prior_cpu:
        raise ValueError("supply an unused Budget with total allowances")
    if max_actions is not None:
        utils.require_integer(max_actions, "max_actions", 0)
    state = _initial(n, base, config)
    saved = None
    if checkpoint is not None:
        payload = _restore(checkpoint, n, base, config)
        if budget.work_limit < payload["work_used"]:
            raise ValueError("budget is below consumed work")
        budget.used = payload["work_used"]
        budget.prior_wall = payload["wall_used"]
        budget.prior_cpu = payload["cpu_used"]
        saved = payload["state"]

    reason, verification_work = "paused", 0
    try:
        budget.consume(0)
        budget.consume((isqrt(config.bounds[-1][1]) + 1) // 2)
        context = SieveContext(
            config.bounds[-1][1] + 1,
            segment_size=config.segment_size,
            memory_bytes=config.memory_bytes - config.workspace_reserve,
        )
        if saved is not None:
            started = budget.used
            try:
                for _ in range(saved["steps"]):
                    if state["done"]:
                        raise ValueError(
                            "checkpoint continues past a terminal outcome"
                        )
                    _advance(state, budget, context, config)
                # JSON numeric types are part of canonical arithmetic state;
                # Python equality would accept 0.0 == 0 and False == 0.
                if _canonical(state) != _canonical(saved):
                    raise ValueError("corrupt p-1 arithmetic progress")
            finally:
                verification_work = budget.used - started
            state = saved
        actions = 0
        while not state["done"] and (
            max_actions is None or actions < max_actions
        ):
            _advance(state, budget, context, config)
            actions += 1
        budget.consume(0)
        reason = state["reason"] if state["done"] else "paused"
    except BudgetExhaustedError as error:
        reason = str(error)
        if saved is not None:
            # Until verification finishes, never publish its partial replay
            # as execution progress or trust the retained arithmetic result.
            state = saved

    divisor = state["factor"] if saved is None or state is saved else None
    if (
        reason in ("work_limit", "wall_limit", "cpu_limit", "cancelled")
        and saved is not None
    ):
        divisor = None
    remaining = (divisor, n // divisor) if divisor is not None else (n,)
    result = FactorizationResult(n, 1, (), tuple(sorted(remaining)))
    payload = {
        "version": CHECKPOINT_VERSION,
        "execution": EXECUTION_VERSION,
        "backend": "python-int",
        "schedule": RATIO_VERSION,
        "n": n,
        "base": base,
        "config": asdict(config),
        "state": state,
        "work_used": budget.used,
        "wall_used": budget.wall_used,
        "cpu_used": budget.cpu_used,
    }
    encoded = _canonical(payload)
    if len(encoded.encode()) > config.memory_bytes // 2:
        raise MemoryError("serialized p-1 checkpoint exceeds storage cap")
    snapshot = {
        "payload": json.loads(encoded),
        "sha256": hashlib.sha256(encoded.encode()).hexdigest(),
    }
    overhead = [
        budget.wall_used - payload["wall_used"],
        budget.cpu_used - payload["cpu_used"],
    ]
    snapshot["serialization_overhead"] = overhead
    snapshot["overhead_sha256"] = hashlib.sha256(
        _canonical(overhead).encode()
    ).hexdigest()
    return PM1Run(
        result,
        divisor,
        reason,
        budget.used,
        budget.wall_used,
        budget.cpu_used,
        verification_work,
        snapshot,
    )
