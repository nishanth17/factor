"""Bounded +/- continuation with curve-private points and scalar recovery."""

from bisect import bisect_left
from math import gcd as integer_gcd

from . import ecm, utils
from .arithmetic import gcd


def verify_progress(job, config, verifier):
    """Check paired cursor, certificate, table and replay storage on resume.

    As with stage one, the checkpoint checksum protects arithmetic state;
    these additional checks bind executable records to the exact schedule.
    No table normalization assumes a unit over a composite modulus.
    """
    from struct import pack

    from .budget import Budget
    from .ecm_programs import ProgramBlock, pair_coverage

    phase = job["phase"]
    if phase in ("setup", "stage_one", "replay", "stage_two_setup"):
        return
    if phase not in (
        "pair_baby",
        "pair_terms",
        "pair_term_replay",
        "pair_scalar_replay",
    ):
        raise ValueError("invalid paired phase")
    wheel = config.ecm_pair_wheel
    distance, n = (
        wheel // 2 if wheel else config.ecm_pair_distance,
        job["n"],
    )
    origin = wheel or job["b1"] - (job["b1"] % 2 == 0)
    if job["distance"] != distance or not (
        origin <= job["center"] <= max(origin, job["b2"] + distance)
        and (
            job["center"] == origin
            if not distance
            else (job["center"] - origin) % (2 * distance) == 0
        )
    ):
        raise ValueError("invalid paired recurrence metadata")

    def point(value):
        if (
            not isinstance(value, list)
            or len(value) != 2
            or any(
                type(coordinate) is not int or not 0 <= coordinate < n
                for coordinate in value
            )
        ):
            raise ValueError("invalid paired point")

    point(job["value"])
    if wheel:
        from .ecm_wheel import verify_table

        verify_table(job, wheel, point)
    elif distance:
        baby = job["baby"]
        if not (2 <= len(baby) <= distance // 2 + 1) or baby[0] is not None:
            raise ValueError("invalid paired baby table")
        for value in baby[1:]:
            point(value)
        if phase != "pair_baby":
            if len(baby) != distance // 2 + 1:
                raise ValueError("incomplete paired baby table")
            for key in ("pair_step", "giant", "previous"):
                point(job[key])

    records, index, cursor = (
        job["pair_records"],
        job["pair_index"],
        job["cursor"],
    )
    if wheel and (
        cursor["left"] != job["b1"] + 1
        and cursor["left"] % wheel != wheel // 2
        or cursor["next"] not in (job["b1"] + 1, job["b2"] + 1)
        and cursor["next"] % wheel != wheel // 2
    ):
        raise ValueError("paired window splits a wheel cell")
    utils.require_integer(index, "paired record index", 0)
    if (not records and index) or (records and not index < len(records)):
        raise ValueError("invalid paired program position")
    if records:
        if any(
            type(value) is not int for record in records for value in record
        ):
            raise ValueError("noncanonical paired coverage integer")
        if cursor["index"] != len(cursor["values"]):
            raise ValueError("paired program disagrees with prime cursor")
        block = ProgramBlock(
            cursor["left"],
            cursor["next"],
            None,
            b"".join(pack("<Q", prime) for prime in cursor["values"]),
            b"",
        )
        coverage = pair_coverage(
            block,
            b1=job["b1"],
            b2=job["b2"],
            distance=distance,
            memory_bytes=4096 + 512 * len(cursor["values"]),
            budget=Budget(
                work_limit=config.segment_size + 1,
                seconds=None,
                cpu_seconds=None,
            ),
            wheel=wheel,
        )
        if records != [list(record) for record in coverage.records()]:
            raise ValueError("corrupt paired coverage records")

    terms, replay = job["terms"], job["term_records"]
    if len(terms) != len(replay) or len(terms) > config.gcd_batch:
        raise ValueError("invalid paired replay storage")
    product = 1
    for term, record in zip(terms, replay):
        utils.require_integer(term, "paired term", 0)
        if term >= n or len(record) != 4:
            raise ValueError("invalid paired replay record")
        center, offset, minus, plus = record
        for value in record:
            utils.require_integer(value, "paired certificate", 0)
        if offset:
            if not (
                distance
                and (
                    integer_gcd(offset, wheel) == 1
                    if wheel
                    else offset % 2 == 0
                )
                and offset <= distance
                and center >= origin
                and (center - origin) % (2 * distance) == 0
                and (not minus or minus == center - offset)
                and (not plus or plus == center + offset)
            ):
                raise ValueError("invalid paired replay certificate")
        elif (minus, plus) != (center, 0):
            raise ValueError("invalid direct replay certificate")
        if not (minus or plus):
            raise ValueError("empty paired replay certificate")
        for prime in (minus, plus):
            if prime and not (
                job["b1"] < prime <= job["b2"]
                and list(verifier.primes(prime, prime + 1)) == [prime]
            ):
                raise ValueError("invalid paired replay prime")
        product = product * term % n
    if job["product"] != product:
        raise ValueError("invalid paired product")
    if phase in ("pair_term_replay", "pair_scalar_replay"):
        utils.require_integer(job["recovery"], "paired recovery", 0)
        if not terms or job["recovery"] > len(terms):
            raise ValueError("invalid paired recovery position")
        if phase == "pair_scalar_replay" and (
            job["recovery"] == len(terms)
            or job["recovery_side"] not in (0, 1)
            or type(job["recovery_side"]) is not int
        ):
            raise ValueError("invalid paired scalar recovery position")


def _finish(job, divisor=None):
    if divisor is not None and not utils.valid_divisor(divisor, job["n"]):
        raise ValueError("paired ECM produced an invalid divisor")
    job.update(done=True, factor=divisor)


def _check(job, budget):
    """Recover a saturated product, including two factors in a single pair."""
    n = job["n"]
    if job["phase"] == "pair_scalar_replay":
        record = job["term_records"][job["recovery"]]
        prime = record[2 + job["recovery_side"]]
        budget.consume(prime.bit_length() + 1 if prime else 1)
        divisor = (
            gcd(ecm.scalar_multiply(prime, *job["value"], n, job["a24"])[1], n)
            if prime
            else 1
        )
        if utils.valid_divisor(divisor, n):
            _finish(job, divisor)
        elif job["recovery_side"] == 0:
            job["recovery_side"] = 1
        else:
            job.update(phase="pair_term_replay", recovery=job["recovery"] + 1)
        return

    if job["phase"] == "pair_term_replay":
        if job["recovery"] == len(job["terms"]):
            _finish(job)
            return
        budget.consume()
        divisor = gcd(job["terms"][job["recovery"]], n)
        if utils.valid_divisor(divisor, n):
            _finish(job, divisor)
        elif divisor == n:
            # One cross-product can cover opposite prime-order components.
            # That term's GCD cannot separate them; replay both certificates.
            job.update(phase="pair_scalar_replay", recovery_side=0)
        else:
            job["recovery"] += 1
        return

    budget.consume()
    divisor = gcd(job["product"], n)
    if utils.valid_divisor(divisor, n):
        _finish(job, divisor)
    elif divisor == n:
        job.update(phase="pair_term_replay", recovery=0)
    else:
        job.update(terms=[], term_records=[], product=1)


def advance(job, budget, context, config, peek_prime):
    """Commit one charged table, program, recurrence, product or replay action.

    D is even; odd centers are separated by 2D. The baby table contains
    [2]Q through [D]Q, while a separate [2D]Q advances giants. No points or
    products enter the reusable store. Segment boundaries may split a pair.
    """
    n, phase, distance = job["n"], job["phase"], config.ecm_pair_distance
    if phase == "stage_two_setup":
        budget.consume(2 if distance else 1)
        job.update(
            distance=distance,
            center=job["b1"] - (job["b1"] % 2 == 0),
            pair_records=[],
            pair_index=0,
            term_records=[],
            phase="pair_baby" if distance else "pair_terms",
        )
        if distance:
            job["baby"] = [
                None,
                list(ecm.point_double(*job["value"], n, job["a24"])),
            ]
        return

    if phase == "pair_baby":
        index = len(job["baby"])
        if index > distance // 2:
            budget.consume(2 + 2 * job["center"].bit_length())
            job.update(
                pair_step=list(
                    ecm.point_double(*job["baby"][-1], n, job["a24"])
                ),
                giant=list(
                    ecm.scalar_multiply(
                        job["center"], *job["value"], n, job["a24"]
                    )
                ),
                previous=list(
                    ecm.scalar_multiply(
                        job["center"] - 2 * distance,
                        *job["value"],
                        n,
                        job["a24"],
                    )
                ),
                phase="pair_terms",
            )
            return
        budget.consume(2)
        point = (
            ecm.point_double(*job["baby"][1], n, job["a24"])
            if index == 2
            else ecm.point_add(
                *job["baby"][-1], *job["baby"][1], *job["baby"][-2], n
            )
        )
        job["baby"].append(list(point))
        return

    if phase in ("pair_term_replay", "pair_scalar_replay"):
        _check(job, budget)
        return
    if len(job["terms"]) == config.gcd_batch:
        _check(job, budget)
        return

    if job["pair_index"] == len(job["pair_records"]):
        if peek_prime(job["cursor"], context, budget) is None:
            if job["terms"]:
                _check(job, budget)
            else:
                _finish(job)
            return
        records = context.coverage(
            job["cursor"],
            b1=job["b1"],
            b2=job["b2"],
            distance=distance,
            budget=budget,
        )
        job.update(pair_records=records, pair_index=0)
        job["cursor"]["index"] = len(job["cursor"]["values"])
        return

    _products(job, budget, config, distance)


def _products(job, budget, config, distance, *, wheel=None):
    """Execute certified terms with shared recurrence and recovery."""
    n = job["n"]
    record = job["pair_records"][job["pair_index"]]
    center, offset, minus, plus = record
    if offset and job["center"] < center:
        budget.consume(2)
        # With a wheel, the first giant is [W]Q and its predecessor is O.
        # Differential addition with O as the difference degenerates;
        # doubling is the exact first transition to [2W]Q.
        giant = (
            list(ecm.point_double(*job["giant"], n, job["a24"]))
            if (wheel and job["center"] == wheel)
            else list(
                ecm.point_add(
                    *job["giant"], *job["pair_step"], *job["previous"], n
                )
            )
        )
        job.update(
            previous=job["giant"],
            giant=giant,
            center=job["center"] + 2 * distance,
        )
        return

    if offset:
        if job["center"] != center:
            raise ValueError("paired centers must advance monotonically")
        # Reserve a whole same-center slice, bounded by the GCD batch and
        # this decoded segment. Its points and recovery certificates stay
        # private, while local arithmetic avoids per-term dispatcher calls.
        start = job["pair_index"]
        end = start
        records = job["pair_records"]
        limit = min(len(records), start + config.gcd_batch - len(job["terms"]))
        while end < limit and records[end][0] == center and records[end][1]:
            end += 1
        budget.consume(2 * (end - start))
        giant_x, giant_z = job["giant"]
        product = job["product"]
        for position in range(start, end):
            entry = records[position]
            baby_index = (
                bisect_left(job["baby_offsets"], entry[1])
                if wheel
                else entry[1] // 2
            )
            bx, bz = job["baby"][baby_index]
            term = (giant_x * bz - bx * giant_z) % n
            job["terms"].append(term)
            job["term_records"].append(entry)
            product = product * term % n
        job.update(product=product, pair_index=end)
    else:
        budget.consume(center.bit_length() + 1)
        term = ecm.scalar_multiply(center, *job["value"], n, job["a24"])[1]
        job["terms"].append(term)
        job["term_records"].append(record)
        job["product"] = job["product"] * term % n
        job["pair_index"] += 1
    if job["pair_index"] == len(job["pair_records"]):
        job.update(pair_records=[], pair_index=0)
