"""Exact density, periodic first-level scores and bounded fused scanning.

Selectors are explicit; automatic adoption requires separate acceptance.
"""

KERNELS = ("density", "presieve", "scan_fused")
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
    return parts


def combine(parts):
    """Use the registry order for stable configuration and receipt names."""
    parts = tuple(parts)
    if len(set(parts)) != len(parts) or any(p not in KERNELS for p in parts):
        raise ValueError("invalid kernel combination")
    return "+".join(k for k in KERNELS if k in parts) or "sparse"


def workspace(marking, width, count):
    """Reserve periods and simultaneous build, rotation and scan scratch."""
    parts = components(marking)
    reserve = 0
    if "presieve" in parts:
        reserve += 32768 + 8 * width * (160 + (count + 7) // 8)
    if "scan_fused" in parts:
        reserve += 32768 + 8 * width
    return reserve


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
        collector._hits[:width] = (self.hits * repeats)[offset:end]
        self.replays += 1
        return self.indices


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
    for start in range(0, width, lanes):
        if budget is not None:
            budget.consume(0)
        stop = min(width, start + lanes)
        count = stop - start
        if lane_bytes == 1 and isinstance(scores, bytearray):
            packed = int.from_bytes(scores[start:stop], "little")
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
        ones = ((1 << (lane_bits * count)) - 1) // ((1 << lane_bits) - 1)
        high_mask = high * ones
        low_mask = (high - 1) * ones
        low_pass = (
            ((packed & low_mask) | high_mask) - (threshold & (high - 1)) * ones
        ) & high_mask
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
    return {
        "kernel_presieve_builds": table.builds,
        "kernel_presieve_replays": table.replays,
    }
