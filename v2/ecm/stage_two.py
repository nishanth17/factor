"""Worker-private immutable indices for the existing ECM continuation."""

from bisect import bisect_right
from dataclasses import dataclass

from ..execution import stage_jobs
from ..execution.worker_contract import integer
from . import core as ecm


@dataclass(frozen=True)
class IndexBlock:
    offsets: tuple
    ends: tuple
    centers: tuple
    first: int
    last: int


class StageTwoSchedule:
    """Retain a funded prefix of exact prime-to-baby-index metadata.

    Keys belong to one compatible internal prime provider. Curve points and
    products never enter this cache. Refused blocks use the original path.
    """

    def __init__(self, segment_size, memory_bytes):
        self.segment_size = integer(segment_size, "segment size", 1)
        self.memory_bytes = integer(
            memory_bytes, "index memory", 0, 64 * 2**20
        )
        self.used_bytes = 4096 + 256 * segment_size
        self.blocks = {}
        self.hits = self.misses = self.refusals = self.build_bytes = 0

    def take(self, cursor, b1, b2, distance, center, limit, budget):
        """Return a group end and the immutable full block, or decline."""
        values, start = cursor["values"], cursor["index"]
        lo, hi = cursor["left"], cursor["next"]
        if distance < 2 or not 0 <= start < len(values) or limit < 1:
            raise ValueError("invalid stage-two index cursor")
        if len(values) > self.segment_size or hi - lo > 2 * self.segment_size:
            raise ValueError("stage-two index block exceeds scratch")
        budget.consume(1)
        key = (lo, hi, b1, b2, distance)
        block = self.blocks.get(key)
        if block is None:
            self.misses += 1
            step = 2 * distance
            groups = min(len(values), (hi - lo) // step + 2)
            reserve = 512 + 64 * len(values) + 256 * groups
            if (
                len(self.blocks) >= 4096
                or self.used_bytes + reserve > self.memory_bytes
            ):
                self.refusals += 1
                return None

            # Reserve the entire build before allocating; interruption keeps
            # that conservative charge and never publishes a partial block.
            budget.consume(len(values) + groups + 1)
            poll = getattr(budget, "poll", None)
            origin = b1 if b1 % 2 else b1 - 1
            offsets, ends, centers = [], [], []
            previous = lo - 1
            for position, prime in enumerate(values):
                if position % 256 == 0:
                    poll() if poll is not None else budget.consume(0)
                if not (
                    type(prime) is int
                    and previous < prime < hi
                    and prime >= lo
                    and b1 < prime <= b2
                    and prime % 2
                ):
                    raise ValueError("incompatible internal prime block")
                group = (prime - origin - 1) // step
                selected_center = origin + group * step
                offset = (prime - selected_center) // 2
                if not 1 <= offset <= distance:
                    raise ValueError("invalid even-baby index")
                if not centers or selected_center != centers[-1]:
                    if centers:
                        ends.append(position)
                    centers.append(selected_center)
                offsets.append(offset)
                previous = prime
            ends.append(len(values))
            block = IndexBlock(
                tuple(offsets),
                tuple(ends),
                tuple(centers),
                values[0],
                values[-1],
            )
            self.blocks[key] = block
            self.used_bytes += reserve
            self.build_bytes += reserve
        else:
            self.hits += 1

        if (
            len(values) != len(block.offsets)
            or values[0] != block.first
            or values[-1] != block.last
        ):
            raise ValueError("stage-two prime provider changed")
        group = bisect_right(block.ends, start)
        if block.centers[group] != center:
            raise ValueError("stage-two giant cursor changed")
        end = min(start + limit, block.ends[group])
        return end, block.offsets


def advance_prepared_job(job, budget, context, config):
    """Use prepared indices only in the unchanged classical ECM branch."""
    if getattr(getattr(context, "native_packet", None), "full_ecm", False):
        return stage_jobs.advance_job(job, budget, context, config)
    if (
        job["kind"] != "ecm"
        or job["phase"] != "stage_two"
        or job["distance"] < 2
        or config.ecm_pair_distance is not None
        or config.ecm_pair_wheel is not None
    ):
        return stage_jobs.advance_job(job, budget, context, config)
    if len(job["terms"]) == config.gcd_batch:
        stage_jobs._batch_check(job, budget)
        return
    prime = stage_jobs.peek_prime(job["cursor"], context, budget)
    if prime is None:
        stage_jobs._batch_check(job, budget)
        if not job["terms"] and not job["done"]:
            stage_jobs._finish(job)
        return

    n, step = job["n"], 2 * job["distance"]
    if prime > job["center"] + step:
        budget.consume(2)
        giant = list(
            getattr(context, "point_add", ecm.point_add)(
                *job["giant"], *job["baby"][-1], *job["previous"], n
            )
        )
        job.update(
            previous=job["giant"], giant=giant, center=job["center"] + step
        )
        return

    cursor = job["cursor"]
    selected = context.stage_two_schedule.take(
        cursor,
        job["b1"],
        job["b2"],
        job["distance"],
        job["center"],
        config.gcd_batch - len(job["terms"]),
        budget,
    )
    if selected is None:
        # peek_prime already owns any newly generated buffer; the original
        # step sees that same cursor and does not charge generation twice.
        return stage_jobs._stage_two(job, budget, context, config)
    end, offsets = selected
    budget.consume(2 * (end - cursor["index"]))
    giant_x, giant_z = job["giant"]
    baby, terms, product = job["baby"], job["terms"], job["product"]
    for position in range(cursor["index"], end):
        bx, bz = baby[offsets[position]]
        term = (giant_x * bz - bx * giant_z) % n
        terms.append(term)
        product = product * term % n
    cursor["index"] = end
    job["product"] = product
