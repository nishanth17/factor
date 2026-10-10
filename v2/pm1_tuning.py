"""Opt-in bounded p-1 challengers with a separate checkpoint identity."""

from dataclasses import dataclass
from math import gcd

from . import pm1_bounded as control
from . import utils


@dataclass(frozen=True)
class PM1TuningConfig(control.PM1Config):
    """Select finite exponent/gap/paired variants without changing defaults.

    Bit caps count the sum of factor bit lengths, an upper bound on the
    exponent length. Paired records cover only primes in the declared range;
    singletons use their exact ordinary relation, without relocation.
    """

    execution_version = "pm1-tuning-v1"

    chunk_bits: int = 0
    gap_mode: str = "cached"
    gap_entries: int = 64
    wheel: int = 0

    def __post_init__(self):
        utils.require_integer(self.chunk_bits, "chunk_bits", 0)
        utils.require_integer(self.gap_entries, "gap_entries", 1)
        utils.require_integer(self.wheel, "wheel", 0)
        if self.chunk_bits and not 32 <= self.chunk_bits <= 4096:
            raise ValueError("chunk_bits must be zero or between 32 and 4096")
        if self.gap_entries > 256 or self.gap_mode not in (
            "cached",
            "recurrence",
        ):
            raise ValueError("invalid bounded gap mode or capacity")
        if self.wheel not in (0, 30, 210):
            raise ValueError("wheel must be zero, 30 or 210")
        super().__post_init__()
        if (
            self.chunk_bits
            and self.bounds[-1][0].bit_length() > self.chunk_bits
        ):
            raise ValueError("a single prime power exceeds the bit cap")

    @property
    def workspace_reserve(self):
        coordinate = 128 + self.max_input_bits // 8
        return (
            super().workspace_reserve
            + 8 * coordinate * (self.wheel + 2 + self.gap_entries)
            + 512 * (self.wheel + 2 * self.gcd_batch)
        )

    def advance(self, state, budget, context):
        old_b1 = state["b1"]
        phase = state["phase"]
        if phase == "fill" and self.chunk_bits:
            _fill_bits(state, budget, context, self)
        elif self.wheel and phase in (
            "stage_two",
            "wheel_terms",
            "gcd",
            "term_replay",
        ):
            _paired_step(state, budget, context, self)
        elif self.gap_mode == "recurrence" and phase == "stage_two":
            _even_gaps(state, budget, context, self)
        else:
            control._step(state, budget, context, self)

        if state["b1"] != old_b1:
            # Only a fixed completed stage-one residue owns these tables.
            # Increased B1 invalidates them even if a residue happens to agree.
            for key in tuple(state):
                if key.startswith(("wheel_", "even_")):
                    del state[key]


def _fill_bits(state, budget, context, config):
    cursor = state["cursor"]
    if len(state["powers"]) >= config.chunk_size or cursor["index"] >= len(
        cursor["values"]
    ):
        control._step(state, budget, context, config)
        return

    bit_count = sum(power.bit_length() for _, power in state["powers"])
    additions, scanned, full = [], 0, False
    for prime in cursor["values"][cursor["index"] :]:
        ratio = control.prime_power_ratio(prime, state["old_b1"], state["b1"])
        if ratio != 1:
            bits = ratio.bit_length()
            if bit_count + bits > config.chunk_bits:
                full = True
                break
            additions.append([prime, ratio])
            bit_count += bits
        scanned += 1
        if len(state["powers"]) + len(additions) == config.chunk_size:
            break

    # Scratch ratios have not changed progress. Reserve scanned candidates
    # before publishing either the chunk or its cursor.
    budget.consume(max(1, scanned))
    state["powers"].extend(additions)
    cursor["index"] += scanned
    if full or not scanned or bit_count == config.chunk_bits:
        state["phase"] = "power"


def _even_gaps(state, budget, context, config):
    cursor = state["cursor"]
    if len(state["terms"]) == config.gcd_batch or cursor["index"] >= len(
        cursor["values"]
    ):
        control._step(state, budget, context, config)
        return

    count = min(
        config.gcd_batch - len(state["terms"]),
        len(cursor["values"]) - cursor["index"],
    )
    primes = cursor["values"][cursor["index"] : cursor["index"] + count]
    previous, gaps = state["previous_prime"], []
    for prime in primes:
        gaps.append(prime - previous)
        previous = prime
    needed = max(
        (
            gap // 2
            for gap in gaps
            if gap % 2 == 0 and gap // 2 <= config.gap_entries
        ),
        default=0,
    )
    powers = state.get("even_powers", [])
    growth = max(0, needed - len(powers))
    budget.consume(sum(gap.bit_length() + 1 for gap in gaps) + growth)

    n, value = state["n"], state["value"]
    if growth:
        square = value * value % n
        while len(powers) < needed:
            powers.append((powers[-1] if powers else 1) * square % n)
        state["even_powers"] = powers
    residue, cache = state["stage_two_value"], state["gap_powers"]
    for gap in gaps:
        if gap % 2 == 0 and gap // 2 <= len(powers):
            power = powers[gap // 2 - 1]
        else:
            key = str(gap)
            if key not in cache:
                if len(cache) == 64:
                    cache.clear()
                cache[key] = pow(value, gap, n)
            power = cache[key]
        residue = residue * power % n
        term = (residue - 1) % n
        state["terms"].append(term)
        state["product"] = state["product"] * term % n
    cursor["index"] += count
    state.update(previous_prime=previous, stage_two_value=residue)


def _center(prime, wheel):
    return (prime + wheel // 2) // wheel * wheel


def _pair_setup(state, budget, config):
    n, value, wheel = state["n"], state["value"], config.wheel
    budget.consume(2 * n.bit_length() + 2 * wheel + 4)
    divisor = gcd(value, n)
    if divisor != 1:
        control._finish(
            state,
            "factor_found" if utils.valid_divisor(divisor, n) else "nonunit",
            divisor if utils.valid_divisor(divisor, n) else None,
        )
        return

    inverse = pow(value, -1, n)
    forward, backward = [1], [1]
    for _ in range(wheel // 2):
        forward.append(forward[-1] * value % n)
        backward.append(backward[-1] * inverse % n)
    state.update(
        wheel_forward=forward,
        wheel_backward=backward,
        wheel_step=forward[-1] * forward[-1] % n,
        wheel_inverse_step=backward[-1] * backward[-1] % n,
        wheel_inverse=inverse,
        wheel_pending=[],
        wheel_records=[],
        wheel_record_index=0,
        wheel_term_primes=[],
    )


def _compile_center(state, budget, config):
    pending = state["wheel_pending"]
    center = _center(pending[0], config.wheel)
    previous = state.get("wheel_center")
    distance = 0 if previous is None else (center - previous) // config.wheel
    cost = (
        2 * center.bit_length()
        if previous is None
        else 2 * max(1, distance.bit_length())
    )
    budget.consume(cost + len(pending) + 4)

    n = state["n"]
    if previous is None:
        giant = pow(state["value"], center, n)
        inverse = pow(state["wheel_inverse"], center, n)
    else:
        step, inverse_step = state["wheel_step"], state["wheel_inverse_step"]
        if distance != 1:
            step = pow(step, distance, n)
            inverse_step = pow(inverse_step, distance, n)
        giant = state["wheel_giant"] * step % n
        inverse = state["wheel_inverse_giant"] * inverse_step % n
    records = {}
    for prime in pending:
        records.setdefault(abs(prime - center), []).append(prime)
    state.update(
        wheel_center=center,
        wheel_giant=giant,
        wheel_inverse_giant=inverse,
        wheel_records=list(records.values()),
        wheel_record_index=0,
        wheel_pending=[],
        phase="wheel_terms",
    )


def _paired_terms(state, budget, config):
    records, start = state["wheel_records"], state["wheel_record_index"]
    count = min(config.gcd_batch - len(state["terms"]), len(records) - start)
    if not count:
        budget.consume()
        state["phase"] = "gcd" if state["terms"] else "stage_two"
        return
    budget.consume(2 * count + 1)

    n, center, giant = state["n"], state["wheel_center"], state["wheel_giant"]
    trace = giant + state["wheel_inverse_giant"]
    forward, backward = state["wheel_forward"], state["wheel_backward"]
    for primes in records[start : start + count]:
        offset = primes[0] - center
        if len(primes) == 2:
            distance = abs(offset)
            # This is A^-center times the two direct prime relations.
            # Multiplication by this unit preserves their product's GCD.
            term = (trace - forward[distance] - backward[distance]) % n
        else:
            baby = forward[offset] if offset >= 0 else backward[-offset]
            term = (giant * baby - 1) % n
        state["terms"].append(term)
        state["wheel_term_primes"].append(primes)
        state["product"] = state["product"] * term % n
    state["wheel_record_index"] += count
    if state["wheel_record_index"] == len(records):
        state.update(wheel_records=[], wheel_record_index=0, phase="stage_two")


def _paired_replay(state, budget, config):
    index = state["recovery_index"]
    records = state["wheel_term_primes"]
    if state["recovery"] >= config.recovery_limit or index == len(records):
        budget.consume()
        control._finish(state, "saturated")
        return

    position = state.get("wheel_recovery_prime", 0)
    prime = records[index][position]
    budget.consume(prime.bit_length() + 1)
    divisor = gcd(pow(state["value"], prime, state["n"]) - 1, state["n"])
    state["recovery"] += 1
    position += 1
    if position == len(records[index]):
        state.update(recovery_index=index + 1, wheel_recovery_prime=0)
    else:
        state["wheel_recovery_prime"] = position
    if utils.valid_divisor(divisor, state["n"]):
        control._finish(state, "factor_found", divisor)


def _paired_step(state, budget, context, config):
    phase = state["phase"]
    if phase == "gcd":
        control._step(state, budget, context, config)
        if state["phase"] == "stage_two":
            state["wheel_term_primes"] = []
        elif state["phase"] == "term_replay":
            state["wheel_recovery_prime"] = 0
        return
    if phase == "term_replay":
        _paired_replay(state, budget, config)
        return
    if phase == "wheel_terms":
        _paired_terms(state, budget, config)
        return
    if "wheel_forward" not in state:
        _pair_setup(state, budget, config)
        return
    if len(state["terms"]) == config.gcd_batch:
        budget.consume()
        state["phase"] = "gcd"
        return
    if state["wheel_records"]:
        budget.consume()
        state["phase"] = "wheel_terms"
        return

    cursor, pending = state["cursor"], state["wheel_pending"]
    if cursor["index"] < len(cursor["values"]):
        prime = cursor["values"][cursor["index"]]
        center = _center(prime, config.wheel)
        if pending and center != _center(pending[0], config.wheel):
            _compile_center(state, budget, config)
            return
        count = 0
        for prime in cursor["values"][cursor["index"] :]:
            if _center(prime, config.wheel) != center:
                break
            count += 1
        budget.consume(count + 1)
        pending.extend(
            cursor["values"][cursor["index"] : cursor["index"] + count]
        )
        cursor["index"] += count
    elif cursor["next"] < cursor["hi"]:
        control._load_segment(state, context, budget, config)
    elif pending:
        _compile_center(state, budget, config)
    elif state["terms"]:
        budget.consume()
        state["phase"] = "gcd"
    else:
        control._next_rung(state, budget, config)
