"""Exact density, bounded periodic scores and fused scanning.

Selectors are explicit; automatic adoption requires separate acceptance.
The packed_queue component changes sparse event storage only; scores, event
order, paid work and conservative reservations remain identical.
"""

from ..common.arithmetic import gcd
from .power_sieve import prime_power_roots

KERNELS = (
    "density",
    "presieve",
    "presieve_powers",
    "unroll4",
    "scan_fused",
    "parity_two",
    "packed_queue",
)
SCAN_CHUNK = 64


def components(marking):
    """Decode a canonical, duplicate-free experimental kernel combination."""
    if marking in ("dense", "sparse", "bucket", None):
        return ()
    if not isinstance(marking, str):
        raise ValueError("unknown marking kernel")
    parts = tuple(marking.split("+"))
    if not parts or parts != tuple(k for k in KERNELS if k in parts):
        raise ValueError("unknown, repeated or unordered marking kernels")
    if "presieve" in parts and "presieve_powers" in parts:
        raise ValueError("choose one periodic score initializer")
    return parts


def combine(parts):
    """Use the registry order for stable configuration and receipt names."""
    parts = tuple(parts)
    if len(set(parts)) != len(parts) or any(p not in KERNELS for p in parts):
        raise ValueError("invalid kernel combination")
    if "presieve" in parts and "presieve_powers" in parts:
        raise ValueError("choose one periodic score initializer")
    return "+".join(k for k in KERNELS if k in parts) or "sparse"


def workspace(marking, width, count):
    """Reserve periods and simultaneous build, rotation and scan scratch."""
    parts = components(marking)
    reserve = 4096 if "parity_two" in parts else 0
    if "presieve" in parts or "presieve_powers" in parts:
        reserve += 32768 + 8 * width * (160 + (count + 7) // 8)
    if "presieve_powers" in parts:
        # Distinct prime powers <= width number at most width - 1. Include
        # their tuples, sort scratch, selected keys and old/new cache sets.
        # The lifter itself keeps only two <=64-root lists at a time.
        reserve += 65536 + 1024 * width
    if "scan_fused" in parts:
        reserve += 32768 + 8 * width
    return reserve


def parity_two_initial(collector, maximum):
    """Return certified universal 2-power levels, score and recovery bit.

    Periodic initializers already supply these levels and keep their own
    exact path. A sparse-queue slice may have removed 2 from this root view.
    No higher roots or conservative overmarks are replaced by this helper.
    """
    cfg = collector.config
    roots = collector._roots
    if (
        cfg.score_backend != "list"
        or cfg.score_policy not in ("powers", "fixed")
        or cfg.small_prime_cutoff > 2
        or collector._presieve is not None
        or not roots
        or roots[0].prime != 2
        or not roots[0].all_positions
    ):
        return 0, 0, 0
    polynomial = collector.polynomial
    collector.budget.consume(polynomial.n_prime.bit_length() + 1)
    if polynomial.a % 4 != 2 or polynomial.b % 2 != 1:
        return 0, 0, 0
    residue = polynomial.n_prime % 8
    if residue not in (1, 5):
        return 0, 0, 0
    levels = min(2 if residue == 1 else 1, max(0, maximum.bit_length() - 1))
    bit = int(levels > 0 and cfg.division == "bucket")
    return levels, levels * collector._logs[0], bit


class Presieve:
    """Periodic first-level scores and bits, with exact power tails."""

    def __init__(self):
        self.key, self.scores, self.hits = None, (), ()
        self.indices = frozenset()
        self.builds = self.replays = 0

    def apply(self, collector, lo, hi, maximum, buckets):
        if collector.config.score_backend != "list":
            return frozenset()
        period, indices = 1, []
        for index, roots in enumerate(collector._roots):
            prime = roots.prime
            if prime < collector.config.small_prime_cutoff:
                continue
            # The lifter emits no first level above the norm bound. Exclude
            # its capped universal branch, whose weight is not one log(p).
            if prime > maximum or (roots.all_positions and prime > 64):
                continue
            if period * prime > collector.config.block_width:
                break
            period *= prime
            indices.append(index)
        key = collector.polynomial, tuple(indices), buckets
        if self.key != key:
            collector.budget.consume(period + len(indices) + 1)
            scores, hits = [0] * period, [0] * period
            for index in indices:
                roots = collector._roots[index]
                residues = (0,) if roots.all_positions else roots.roots
                step = 1 if roots.all_positions else roots.prime
                for root in residues:
                    positions = range(root, period, step)
                    collector.budget.consume(len(positions) + 1)
                    for position in positions:
                        scores[position] += collector._logs[index]
                        if buckets:
                            hits[position] |= 1 << index
            self.key, self.scores, self.hits = key, scores, hits
            self.indices = frozenset(indices)
            self.builds += 1
        width, offset = hi - lo, lo % period
        collector.budget.consume(width + 1)
        repeats = (offset + width + period - 1) // period
        end = offset + width
        collector._scores[:width] = (self.scores * repeats)[offset:end]
        initial = (self.hits * repeats)[offset:end]
        words = getattr(collector, "_word_hits", None)
        if words is not None and collector._hits is words:
            words.assign_prefix(initial, collector.budget)
        else:
            collector._hits[:width] = initial
        self.replays += 1
        return self.indices


class PowerPresieve:
    """Bounded LCM period for a prefix of exact prime-power mark levels."""

    def __init__(self):
        self.key, self.scores, self.hits = None, (), ()
        self.levels = frozenset()
        self.builds = self.replays = 0

    def apply(self, collector, lo, hi, maximum, buckets):
        if collector.config.score_backend != "list":
            return frozenset()
        cap = min(maximum, collector.config.block_width)
        key = collector.polynomial, cap, buckets
        if self.key != key:
            steps = []
            for index, roots in enumerate(collector._roots):
                prime = roots.prime
                if prime > cap:
                    break
                collector.budget.consume(1)
                if prime < collector.config.small_prime_cutoff:
                    continue
                if not roots.roots and not roots.all_positions:
                    continue
                # The large universal branch emits a norm-dependent upper
                # bound at modulus one, not ordinary log(p) contributions.
                if roots.all_positions and prime > 64:
                    continue
                modulus = prime
                while modulus <= cap:
                    collector.budget.consume(1)
                    steps.append((modulus, index))
                    modulus *= prime

            collector.budget.consume(len(steps) + 1)
            steps.sort()
            period, limits = 1, {}
            for modulus, index in steps:
                collector.budget.consume(1)
                candidate = period // gcd(period, modulus) * modulus
                if candidate > collector.config.block_width:
                    break
                period = candidate
                limits[index] = modulus

            collector.budget.consume(period + len(limits) + 1)
            scores, hits, levels = [0] * period, [0] * period, set()
            for index, limit in limits.items():
                weight = collector._logs[index]
                marks = prime_power_roots(
                    collector.polynomial,
                    collector._roots[index],
                    limit,
                    weight,
                    collector.budget,
                    inverses=collector._lift_inverses,
                    shared_inverses=collector._shared_lift_inverses,
                    root_transfers=collector._power_transfers,
                )
                for modulus, residues, power_weight in marks:
                    # Full-period lifting can reach its singular-root cap
                    # before window-pruned lifting. Retain only its normal
                    # levels; the original stream supplies any fallback.
                    if power_weight != weight:
                        break
                    collector.budget.consume(len(residues) + 1)
                    levels.add((index, modulus))
                    first_level = modulus == collector._roots[index].prime
                    for root in residues:
                        positions = range(root, period, modulus)
                        collector.budget.consume(len(positions) + 1)
                        for position in positions:
                            scores[position] += weight
                            if first_level and buckets:
                                hits[position] |= 1 << index

            # Publish only a complete cache. Cancellation during construction
            # leaves the old period intact and every attempted action charged.
            collector.budget.consume(len(levels) + 1)
            frozen_levels = frozenset(levels)
            self.key, self.scores, self.hits, self.levels = (
                key,
                scores,
                hits,
                frozen_levels,
            )
            self.builds += 1

        width, offset = hi - lo, lo % len(self.scores)
        collector.budget.consume(width + 1)
        repeats = (offset + width + len(self.scores) - 1) // len(self.scores)
        end = offset + width
        collector._scores[:width] = (self.scores * repeats)[offset:end]
        initial = (self.hits * repeats)[offset:end]
        words = getattr(collector, "_word_hits", None)
        if words is not None and collector._hits is words:
            words.assign_prefix(initial, collector.budget)
        else:
            collector._hits[:width] = initial
        self.replays += 1
        return self.levels


def add_unrolled(scores, hits, offset, width, step, weight, bit):
    """Write one already-reserved root in four-stride groups and exact tails.

    Keep every write in the scalar root's order. There are no extra budget
    polls or merged reservations: refusal/resume retains the original atomic
    boundary. Sparse roots are handled by the caller before this function.
    """
    if bit:
        while offset + 3 * step < width:
            second, third, fourth = (
                offset + step,
                offset + 2 * step,
                offset + 3 * step,
            )
            hits[offset] |= bit
            scores[offset] += weight
            hits[second] |= bit
            scores[second] += weight
            hits[third] |= bit
            scores[third] += weight
            hits[fourth] |= bit
            scores[fourth] += weight
            offset += 4 * step
        if offset < width:
            hits[offset] |= bit
            scores[offset] += weight
            offset += step
            if offset < width:
                hits[offset] |= bit
                scores[offset] += weight
                offset += step
                if offset < width:
                    hits[offset] |= bit
                    scores[offset] += weight
    else:
        while offset + 3 * step < width:
            scores[offset] += weight
            scores[offset + step] += weight
            scores[offset + 2 * step] += weight
            scores[offset + 3 * step] += weight
            offset += 4 * step
        if offset < width:
            scores[offset] += weight
            offset += step
            if offset < width:
                scores[offset] += weight
                offset += step
                if offset < width:
                    scores[offset] += weight


def word_candidates(scores, threshold, width, budget=None):
    """Yield ordered hits using guarded packed-lane threshold comparisons.

    A lane's high bit prevents a borrow reaching its neighbor. Compare the
    low bits against threshold's low bits, then combine the original high
    bits. This also handles scores above 127; a wrapping signed-byte offset
    would lose such hits. Wider thresholds use wider lanes, without changing
    the original scores needed by exact candidate refinement.
    """
    if threshold <= 0:
        yield from range(width)
        return
    lane_bytes = max(1, (threshold.bit_length() + 7) // 8)
    lane_bits = 8 * lane_bytes
    high = 1 << (lane_bits - 1)
    lanes = max(1, SCAN_CHUNK // lane_bytes)
    if width <= 0:
        return
    full_count = min(lanes, width)
    ones = ((1 << (lane_bits * full_count)) - 1) // ((1 << lane_bits) - 1)
    full_high_mask = high * ones
    full_low_mask = (high - 1) * ones
    full_subtrahend = (threshold & (high - 1)) * ones
    for start in range(0, width, lanes):
        if budget is not None:
            budget.consume(0)
        stop = min(width, start + lanes)
        count = stop - start
        if lane_bytes == 1:
            data = (
                scores[start:stop]
                if isinstance(scores, bytearray)
                else bytes(
                    min(scores[i], threshold) for i in range(start, stop)
                )
            )
            packed = int.from_bytes(data, "little")
        else:
            # Clipping only this temporary scan view is exact: score >= T
            # iff min(score,T) >= T. Stored scores and thresholds are intact.
            packed = int.from_bytes(
                b"".join(
                    min(scores[i], threshold).to_bytes(lane_bytes, "little")
                    for i in range(start, stop)
                ),
                "little",
            )
        if count == full_count:
            high_mask, low_mask = full_high_mask, full_low_mask
            subtrahend = full_subtrahend
        else:
            ones = ((1 << (lane_bits * count)) - 1) // ((1 << lane_bits) - 1)
            high_mask, low_mask = high * ones, (high - 1) * ones
            subtrahend = (threshold & (high - 1)) * ones
        low_pass = (((packed & low_mask) | high_mask) - subtrahend) & high_mask
        selected = (
            low_pass & packed
            if threshold & high
            else low_pass | (packed & high_mask)
        )
        while selected:
            bit = selected & -selected
            yield start + (bit.bit_length() - 1) // lane_bits
            selected ^= bit


def counters(collector):
    """Expose period construction/replay apart from semantic hit counts."""
    table = collector._presieve
    if table is None:
        return {}
    result = {
        "kernel_presieve_builds": table.builds,
        "kernel_presieve_replays": table.replays,
    }
    if isinstance(table, PowerPresieve):
        result["kernel_presieve_power_levels"] = len(table.levels)
    return result
