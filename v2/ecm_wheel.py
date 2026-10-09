"""Complete bounded wheel cells and sparse curve-private baby tables.

Each prime maps to its nearest multiple of W. Program boundaries lie halfway
between centers, so sieve boundaries cannot split a pair. W fits inside the
accepted prime-segment bound; no B2-sized plan or pending-prime queue exists.
This is one distance set, without relocation or overlapping-window matching.
"""

from math import gcd

from . import ecm, ecm_paired, utils


def peek_prime(cursor, context, budget, wheel):
    """Fill a bounded block ending between cells, including partial tails."""
    while cursor["index"] >= len(cursor["values"]):
        left = cursor["next"]
        if left >= cursor["hi"]:
            return None
        right = min(cursor["hi"], left + 2 * context.segment_size)
        if right < cursor["hi"]:
            right -= (right - wheel // 2) % wheel
        if right <= left:
            raise ValueError("wheel cell exceeds the prime segment bound")
        values = context.program_segment(left, right, budget)
        cursor.update(left=left, next=right, values=values, index=0)
    return cursor["values"][cursor["index"]]


def verify_table(job, wheel, point):
    """Verify sparse ordering, partial construction and canonical points."""
    index = utils.require_integer(job["baby_next"], "baby cursor", 3)
    if index % 2 == 0 or index > wheel // 2 + 2:
        raise ValueError("invalid sparse baby cursor")
    expected = [
        offset for offset in range(1, index, 2) if gcd(offset, wheel) == 1
    ]
    offsets = job["baby_offsets"]
    if (
        any(type(offset) is not int for offset in offsets)
        or offsets != expected
        or len(job["baby"]) != len(offsets)
    ):
        raise ValueError("invalid sparse baby offsets")
    for value in job["baby"]:
        point(value)
    for key in ("wheel_previous", "wheel_last", "pair_step"):
        point(job[key])
    if job["phase"] != "pair_baby":
        if index <= wheel // 2:
            raise ValueError("incomplete sparse baby table")
        for key in ("giant", "previous"):
            point(job[key])


def advance(job, budget, context, config):
    """Execute one charged action using the accepted arithmetic contracts."""
    n, phase, wheel = job["n"], job["phase"], config.ecm_pair_wheel
    distance = wheel // 2
    if phase == "stage_two_setup":
        budget.consume(2)
        job.update(
            distance=distance,
            center=max(wheel, ((job["b1"] + 1 + distance) // wheel) * wheel),
            pair_records=[],
            pair_index=0,
            term_records=[],
            baby=[job["value"][:]],
            baby_offsets=[1],
            baby_next=3,
            wheel_previous=job["value"][:],
            wheel_last=job["value"][:],
            pair_step=list(ecm.point_double(*job["value"], n, job["a24"])),
            phase="pair_baby",
        )
        return

    if phase == "pair_baby":
        offset = job["baby_next"]
        if offset > distance:
            center = job["center"]
            budget.consume(wheel.bit_length() + 2 * center.bit_length() + 1)
            step = list(
                ecm.scalar_multiply(wheel, *job["value"], n, job["a24"])
            )
            job.update(
                pair_step=step,
                giant=list(
                    ecm.scalar_multiply(center, *job["value"], n, job["a24"])
                ),
                previous=list(
                    ecm.scalar_multiply(
                        center - wheel, *job["value"], n, job["a24"]
                    )
                ),
                phase="pair_terms",
            )
            return
        # Generate odd multiples with two scratch points. Discarding a
        # noncoprime distance saves storage, but its arithmetic is still paid.
        budget.consume(2)
        value = list(
            ecm.point_add(
                *job["wheel_last"],
                *job["pair_step"],
                *job["wheel_previous"],
                n,
            )
        )
        if gcd(offset, wheel) == 1:
            job["baby"].append(value)
            job["baby_offsets"].append(offset)
        job.update(
            wheel_previous=job["wheel_last"],
            wheel_last=value,
            baby_next=offset + 2,
        )
        return

    if phase in ("pair_term_replay", "pair_scalar_replay") or (
        len(job["terms"]) == config.gcd_batch
    ):
        ecm_paired._check(job, budget)
        return

    if job["pair_index"] == len(job["pair_records"]):
        if peek_prime(job["cursor"], context, budget, wheel) is None:
            if job["terms"]:
                ecm_paired._check(job, budget)
            else:
                ecm_paired._finish(job)
            return
        records = context.coverage(
            job["cursor"],
            b1=job["b1"],
            b2=job["b2"],
            distance=distance,
            budget=budget,
            wheel=wheel,
        )
        job.update(pair_records=records, pair_index=0)
        job["cursor"]["index"] = len(job["cursor"]["values"])
        return

    ecm_paired._products(job, budget, config, distance, wheel=wheel)
