"""Serializable atomic steps for rho, p-1, and Montgomery ECM candidates."""

import random
from bisect import bisect_right
from math import gcd, isqrt, prod

from . import ecm, utils


def prime_cursor(lo, hi):
    """Create a half-open cursor whose buffered segment survives a pause."""
    return {"left": lo, "next": lo, "hi": hi, "values": [], "index": 0}


def peek_prime(cursor, context, budget):
    """Return the next prime without committing consumption.

    Segment generation is charged before replacing the buffer. At most one
    segment's values is retained; an exhausted cursor returns None.
    """
    while cursor["index"] >= len(cursor["values"]):
        left = cursor["next"]
        if left >= cursor["hi"]:
            return None
        right = min(cursor["hi"], left + 2 * context.segment_size)
        budget.consume(context.segment_size + len(context.base_primes))
        segment = getattr(context, "prime_segment", None)
        values = (
            list(context.primes(left, right))
            if segment is None
            else segment(left, right)
        )
        cursor.update(left=left, next=right, values=values, index=0)
    return cursor["values"][cursor["index"]]


def take_prime(cursor):
    """Commit a prime only after its arithmetic action was reserved."""
    cursor["index"] += 1


def new_job(kind, n, seed, b1=0, b2=0):
    """Create one reproducible walk, base, or curve assignment."""
    return {
        "kind": kind,
        "n": n,
        "seed": seed,
        "b1": b1,
        "b2": b2,
        "phase": "setup",
        "factor": None,
        "done": False,
        "cursor": prime_cursor(2, b1 + 1) if b1 >= 2 else None,
        "powers": [],
        "terms": [],
        "product": 1,
        "recovery": 0,
    }


def _finish(job, divisor=None):
    """Store a validated outcome, never accepting saturation as success."""
    if divisor is not None and not utils.valid_divisor(divisor, job["n"]):
        raise ValueError("candidate produced an invalid divisor")
    job.update(done=True, factor=divisor)


def _rho_step(job, budget, config):
    """Advance one Brent batch or one bounded saturation recovery step."""
    n = job["n"]
    if job["phase"] == "setup":
        budget.consume()
        generator = random.Random(job["seed"])
        job.update(
            y=generator.randrange(1, n),
            offset=generator.randrange(1, n),
            length=1,
            used=0,
            advance=0,
            position=0,
            phase="advance",
        )
    allowance = config.rho_evaluations - job["used"]
    if allowance <= 0:
        _finish(job)
        return
    if job["phase"] == "recover":
        if job["recovery"] >= min(config.recovery_limit, job["count"]):
            _finish(job)
            return
        budget.consume(2)
        job["saved"] = (job["saved"] ** 2 + job["offset"]) % n
        job["used"] += 1
        job["recovery"] += 1
        divisor = gcd(abs(job["x"] - job["saved"]), n)
        if utils.valid_divisor(divisor, n):
            _finish(job, divisor)
        return
    if job["phase"] == "advance":
        count = min(
            config.rho_batch, job["length"] - job["advance"], allowance
        )
        budget.consume(count)
        if job["advance"] == 0:
            job["x"] = job["y"]
        y, offset = job["y"], job["offset"]
        # Keep arithmetic in locals within the already-reserved atomic batch.
        # The committed boundary state is identical to per-iteration writes.
        for _ in range(count):
            y = (y * y + offset) % n
        job["y"] = y
        job["used"] += count
        job["advance"] += count
        if job["advance"] == job["length"]:
            job["phase"] = "batch"
        return
    count = min(config.rho_batch, job["length"] - job["position"], allowance)
    budget.consume(count + 1)
    job["saved"] = job["y"]
    y, offset, x = job["y"], job["offset"], job["x"]
    product = 1
    for _ in range(count):
        y = (y * y + offset) % n
        product = product * abs(x - y) % n
    job["y"] = y
    job["used"] += count
    job["position"] += count
    divisor = gcd(product, n)
    if utils.valid_divisor(divisor, n):
        _finish(job, divisor)
    elif divisor == n:
        job.update(phase="recover", count=count)
    elif job["position"] == job["length"]:
        job.update(
            length=2 * job["length"], position=0, advance=0, phase="advance"
        )


def _apply(job, value, scalar):
    """Apply a stage-one exponent to a residue or projective point."""
    if job["kind"] == "pm1":
        return pow(value, scalar, job["n"])
    return list(ecm.scalar_multiply(scalar, *value, job["n"], job["a24"]))


def _state_gcd(job, value):
    """Expose a nonunit without affine conversion."""
    term = value - 1 if job["kind"] == "pm1" else value[1]
    return gcd(term, job["n"])


def _stage_one(job, budget, context, config):
    """Accumulate a moderate chunk, then replay saturation at prime units."""
    if job["phase"] == "replay":
        position = job["replay_index"]
        if position >= len(job["powers"]):
            _finish(job)
            return
        prime, power = job["powers"][position]
        budget.consume(prime.bit_length() + 1)
        job["replay_value"] = _apply(job, job["replay_value"], prime)
        job["replay_power"] *= prime
        divisor = _state_gcd(job, job["replay_value"])
        if utils.valid_divisor(divisor, job["n"]):
            _finish(job, divisor)
        elif divisor == job["n"]:
            _finish(job)
        elif job["replay_power"] == power:
            job.update(replay_index=position + 1, replay_power=1)
        return
    prime = (
        None
        if len(job["powers"]) >= config.chunk_size
        else peek_prime(job["cursor"], context, budget)
    )
    if prime is not None and len(job["powers"]) < config.chunk_size:
        cursor = job["cursor"]
        count = min(
            config.chunk_size - len(job["powers"]),
            len(cursor["values"]) - cursor["index"],
        )
        budget.consume(count)
        for _ in range(count):
            prime = cursor["values"][cursor["index"]]
            job["powers"].append([prime, utils.prime_power(prime, job["b1"])])
            take_prime(cursor)
        return
    if not job["powers"]:
        job.update(
            phase="stage_two_setup",
            cursor=prime_cursor(job["b1"] + 1, job["b2"] + 1),
        )
        return
    scalar = prod(power for _, power in job["powers"])
    budget.consume(scalar.bit_length() + 1)
    start = job["value"]
    value = _apply(job, start, scalar)
    divisor = _state_gcd(job, value)
    if utils.valid_divisor(divisor, job["n"]):
        _finish(job, divisor)
    elif divisor == job["n"]:
        job.update(
            phase="replay",
            replay_index=0,
            replay_power=1,
            replay_value=start,
        )
    else:
        job.update(value=value, powers=[])


def _batch_check(job, budget):
    """Replay each retained term on saturation, charging every GCD."""
    if not job["terms"]:
        return
    if job["phase"] == "term_replay":
        if job["recovery"] >= len(job["terms"]):
            _finish(job)
            return
        budget.consume()
        divisor = gcd(job["terms"][job["recovery"]], job["n"])
        job["recovery"] += 1
        if utils.valid_divisor(divisor, job["n"]):
            _finish(job, divisor)
        return
    budget.consume()
    divisor = gcd(job["product"], job["n"])
    if utils.valid_divisor(divisor, job["n"]):
        _finish(job, divisor)
    elif divisor == job["n"]:
        job.update(phase="term_replay", recovery=0)
    else:
        job.update(terms=[], product=1)


def _stage_two(job, budget, context, config):
    """Stream p-1 gaps or ECM even baby steps; retain one GCD batch."""
    n = job["n"]
    if job["phase"] == "term_replay":
        _batch_check(job, budget)
        return
    if len(job["terms"]) == config.gcd_batch:
        _batch_check(job, budget)
        return
    prime = peek_prime(job["cursor"], context, budget)
    if prime is None:
        _batch_check(job, budget)
        if prime is None and not job["terms"] and not job["done"]:
            _finish(job)
        return
    if job["kind"] == "pm1":
        cursor = job["cursor"]
        count = min(
            config.gcd_batch - len(job["terms"]),
            len(cursor["values"]) - cursor["index"],
        )
        start = cursor["index"]
        stop = start + count
        primes = cursor["values"][start:stop]
        previous = job["previous_prime"]
        gaps = []
        for candidate in primes:
            gaps.append(candidate - previous)
            previous = candidate
        budget.consume(sum(gap.bit_length() + 1 for gap in gaps))
        cache = job.setdefault("gap_powers", {})
        value, residue = job["value"], job["stage_two_value"]
        product, terms = job["product"], job["terms"]
        for candidate, gap in zip(primes, gaps):
            key = str(gap)  # String keys survive a JSON roundtrip unchanged.
            if key not in cache:
                if len(cache) >= 64:
                    cache.clear()
                cache[key] = pow(value, gap, n)
            residue = residue * cache[key] % n
            term = (residue - 1) % n
            terms.append(term)
            product = product * term % n
        cursor["index"] = stop
        job.update(
            stage_two_value=residue, product=product, previous_prime=previous
        )
        return
    elif job["distance"] < 2:
        budget.consume(prime.bit_length())
        term = ecm.scalar_multiply(prime, *job["value"], n, job["a24"])[1]
    else:
        step = 2 * job["distance"]
        if prime > job["center"] + step:
            budget.consume(2)
            giant = list(
                ecm.point_add(
                    *job["giant"], *job["baby"][-1], *job["previous"], n
                )
            )
            job.update(
                previous=job["giant"],
                giant=giant,
                center=job["center"] + step,
            )
            return
        cursor = job["cursor"]
        end = min(
            cursor["index"] + config.gcd_batch - len(job["terms"]),
            bisect_right(cursor["values"], job["center"] + step),
        )
        budget.consume(2 * (end - cursor["index"]))
        giant_x, giant_z = job["giant"]
        baby, terms, product = job["baby"], job["terms"], job["product"]
        for position in range(cursor["index"], end):
            candidate = cursor["values"][position]
            index = (candidate - job["center"]) // 2
            bx, bz = baby[index]
            term = (giant_x * bz - bx * giant_z) % n
            terms.append(term)
            product = product * term % n
        cursor["index"] = end
        job["product"] = product
        return
    job["terms"].append(term)
    job["product"] = job["product"] * term % n
    take_prime(job["cursor"])


def advance_job(job, budget, context, config):
    """Commit one bounded candidate action; exhaustion preserves its state."""
    if job["done"]:
        return
    if job["kind"] == "rho":
        _rho_step(job, budget, config)
        return
    if job["phase"] == "setup":
        budget.consume(job["n"].bit_length())
        if job["kind"] == "pm1":
            base = 2 + job["seed"] % max(1, job["n"] - 3)
            divisor = gcd(base, job["n"])
            job["base"] = base
            if divisor != 1:
                _finish(job, divisor if divisor < job["n"] else None)
                return
            job["value"] = base
        else:
            sigma = random.Random(job["seed"]).randrange(6, 2**63)
            setup = ecm.setup_curve(job["n"], sigma)
            job["sigma"] = sigma
            if setup.factor is not None or setup.retry:
                _finish(job, setup.factor)
                return
            job.update(value=list(setup.point), a24=setup.a24)
        job["phase"] = "stage_one"
        return
    if job["phase"] in ("stage_one", "replay"):
        _stage_one(job, budget, context, config)
        return
    if job["phase"] == "stage_two_setup":
        if job["kind"] == "pm1":
            budget.consume()
            job.update(previous_prime=0, stage_two_value=1, phase="stage_two")
            return
        center = job["b1"] if job["b1"] % 2 else job["b1"] - 1
        distance = min(isqrt(job["b2"]), (center - 1) // 2)
        if distance < 2:
            budget.consume()
            job.update(distance=distance, phase="stage_two")
            return
        budget.consume(2)
        job.update(
            center=center,
            distance=distance,
            baby=[
                None,
                list(ecm.point_double(*job["value"], job["n"], job["a24"])),
            ],
            phase="baby_steps",
        )
        return
    if job["phase"] == "baby_steps":
        index = len(job["baby"])
        budget.consume(2)
        if index == 2:
            point = ecm.point_double(*job["baby"][1], job["n"], job["a24"])
        else:
            point = ecm.point_add(
                *job["baby"][-1],
                *job["baby"][1],
                *job["baby"][-2],
                job["n"],
            )
        job["baby"].append(list(point))
        if index == job["distance"]:
            job["phase"] = "giant_setup"
        return
    if job["phase"] == "giant_setup":
        budget.consume(2 * job["center"].bit_length())
        job.update(
            giant=list(
                ecm.scalar_multiply(
                    job["center"], *job["value"], job["n"], job["a24"]
                )
            ),
            previous=list(
                ecm.scalar_multiply(
                    job["center"] - 2 * job["distance"],
                    *job["value"],
                    job["n"],
                    job["a24"],
                )
            ),
            phase="stage_two",
        )
        return
    _stage_two(job, budget, context, config)
