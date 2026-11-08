"""Static scalar specializations of the generic collector implementations.

These checked-in functions remove disabled grouping, power-period, fused-scan
and packed-queue branches. Selection occurs at collector construction. The
commented generic implementations remain the arithmetic/accounting oracle;
changes here require matching buffer, paid-prefix and resume checks. There is
no runtime AST transformation or source compilation.
"""

from .sieve_collector import (
    MAX_WINDOW_WIDTH,
    BudgetExhaustedError,
    PollingBudget,
    SieveResult,
    _CollectionProgress,
    checked_position,
    kernels,
    require_published,
)
from .sparse_queue import ENTRY_BYTES, array, clear


def direct_scalar(self, lo, hi, stats):
    """Mark complete root hits in bounded prime chunks, without F(x)."""
    width = hi - lo
    self._omitted_allowance = 0
    self.budget.consume(width)
    two_levels, initial_score, initial_hit = (0, 0, 0)
    powers = self.config.score_policy in ("powers", "fixed")
    parity_two = "parity_two" in self.config._kernels
    if parity_two:
        lower, maximum = self._bounds(lo, hi)
        two_levels, initial_score, initial_hit = kernels.parity_two_initial(
            self, maximum
        )
    period_initializes = (
        powers
        and self._presieve is not None
        and (self.config.score_backend == "list")
    )
    words = self._word_hits if self._hits is self._word_hits else None
    if words is not None:
        words.begin(initial_hit, self.budget)
    if not period_initializes:
        for index in range(width):
            self._scores[index] = initial_score
            if words is None:
                self._hits[index] = initial_hit
    if not parity_two:
        lower, maximum = self._bounds(lo, hi)
    self._block_norm_bits = (lower.bit_length(), maximum.bit_length())
    if two_levels:
        stats["parity_two_blocks"] = stats.get("parity_two_blocks", 0) + 1
        stats["parity_two_seeded_hits"] = (
            stats.get("parity_two_seeded_hits", 0) + two_levels * width
        )
    lower_threshold = max(0, lower.bit_length() - 1)
    lower_threshold -= (self._candidate_bound - 1).bit_length()
    if (
        self.config.score_policy == "adaptive"
        and lower_threshold <= 0
        and (self.config.threshold_extra == 0)
        and (self.config.division != "bucket")
    ):
        stats["blocks"] += 1
        stats["threshold_min"] = 0
        stats["skipped_score_blocks"] = (
            stats.get("skipped_score_blocks", 0) + 1
        )
        return 0
    weights = self._logs if powers else self._weights(maximum)
    omitted = 0
    buckets = (
        self.config.marking == "bucket" or self.config.division == "bucket"
    )
    chunk = self.config.metadata_chunk
    unroll = "unroll4" in self.config._kernels
    density = (
        "density" in self.config._kernels or unroll
    ) and self.config.score_backend == "list"
    presieved = (
        self._presieve.apply(self, lo, hi, maximum, buckets)
        if powers and self._presieve is not None
        else frozenset()
    )
    for start in range(0, len(self._roots), chunk):
        for index in range(start, min(start + chunk, len(self._roots))):
            roots, weight = (self._roots[index], weights[index])
            prime = roots.prime
            recovery_mask = self._recovery_mask(index) if buckets else 0
            skipped = prime < self.config.small_prime_cutoff
            if skipped:
                omitted += weight * maximum.bit_length() if powers else weight
            if powers and (not skipped):
                first_level = True
                marks = (
                    self._power_marks(index, maximum, lo, hi, stats)
                    if self.config.power_plan_bytes
                    else self._power_roots(
                        self.polynomial,
                        roots,
                        maximum,
                        weight,
                        self.budget,
                        lo=lo,
                        hi=hi,
                        inverses=self._lift_inverses,
                        shared_inverses=self._shared_lift_inverses,
                        root_transfers=self._power_transfers,
                    )
                )
                for modulus, lifted, power_weight in marks:
                    if (
                        index == 0
                        and two_levels
                        and (power_weight == weight)
                        and (modulus in (2, 4))
                        and (modulus <= 1 << two_levels)
                    ):
                        self.budget.consume(1)
                        if first_level:
                            stats["root_hits"] += width
                        stats["power_hits"] += width
                        first_level = False
                        continue
                    for root in lifted:
                        if first_level and index in presieved or False:
                            hits = range((root - lo) % modulus, width, modulus)
                            self.budget.consume(len(hits) + 1)
                            if first_level:
                                stats["root_hits"] += len(hits)
                            stats["power_hits"] += len(hits)
                            continue
                        if first_level and buckets and (words is not None):
                            hits = range((root - lo) % modulus, width, modulus)
                            self.budget.consume(len(hits) + 1)
                            stats["root_hits"] += len(hits)
                            stats["power_hits"] += len(hits)
                            words.add(index, hits, self._scores, power_weight)
                            continue
                        if density:
                            offset = (root - lo) % modulus
                            if modulus >= width:
                                count = int(offset < width)
                                self.budget.consume(count + 1)
                                if first_level:
                                    stats["root_hits"] += count
                                stats["power_hits"] += count
                                if count:
                                    if first_level and buckets:
                                        self._hits[offset] |= recovery_mask
                                    self._scores[offset] += power_weight
                                continue
                            hits = range(offset, width, modulus)
                            self.budget.consume(len(hits) + 1)
                            stats["power_hits"] += len(hits)
                            if first_level:
                                stats["root_hits"] += len(hits)
                            if unroll:
                                kernels.add_unrolled(
                                    self._scores,
                                    self._hits,
                                    offset,
                                    width,
                                    modulus,
                                    power_weight,
                                    recovery_mask
                                    if first_level and buckets
                                    else 0,
                                )
                                continue
                            if first_level and buckets:
                                bit = recovery_mask
                                for hit in hits:
                                    self._hits[hit] |= bit
                                    self._scores[hit] += power_weight
                            else:
                                for hit in hits:
                                    self._scores[hit] += power_weight
                            continue
                        hits = range((root - lo) % modulus, width, modulus)
                        self.budget.consume(len(hits) + 1)
                        if first_level:
                            stats["root_hits"] += len(hits)
                        stats["power_hits"] += len(hits)
                        if (
                            first_level
                            and buckets
                            and (self.config.score_backend == "list")
                        ):
                            scores, hit_masks = (self._scores, self._hits)
                            for hit in hits:
                                hit_masks[hit] |= recovery_mask
                                scores[hit] += power_weight
                        else:
                            if first_level and buckets:
                                for hit in hits:
                                    self._hits[hit] |= recovery_mask
                            self._add(hits, power_weight, stats)
                    first_level = False
                continue
            residues = (0,) if roots.all_positions else roots.roots
            step = 1 if roots.all_positions else prime
            for root in residues:
                offset = (root - lo) % step
                hits = range(offset, width, step)
                self.budget.consume(len(hits) + 1)
                stats["root_hits"] += len(hits)
                if buckets:
                    if words is not None:
                        words.add(index, hits)
                    else:
                        for hit in hits:
                            self._hits[hit] |= recovery_mask
                if skipped or powers or self.config.marking == "bucket":
                    continue
                if (
                    self.config.marking == "sparse"
                    or bool(self.config._kernels)
                ) and step >= width:
                    hits = (offset,) if offset < width else ()
                self._add(hits, weight, stats)
    if not powers and self.config.marking == "bucket":
        for offset in range(width):
            bits = self._hits[offset]
            while bits:
                bit = bits & -bits
                index = bit.bit_length() - 1
                if self._roots[index].prime >= self.config.small_prime_cutoff:
                    self._add((offset,), weights[index], stats)
                bits ^= bit
    threshold = (
        self._lower_score(lower)
        if self.config.score_policy == "fixed"
        else max(0, lower.bit_length() - 1)
    )
    self._omitted_allowance = omitted
    threshold -= (
        self._residual_score()
        if self.config.score_policy == "fixed"
        else (self._candidate_bound - 1).bit_length()
    ) + omitted
    threshold = max(0, threshold + self.config.threshold_extra)
    if self.config.score_backend != "list":
        threshold = self._clip_threshold(threshold)
    stats["blocks"] += 1
    stats["threshold_min"] = min(stats["threshold_min"], threshold)
    stats["threshold_max"] = max(stats["threshold_max"], threshold)
    return threshold


def collect_scalar(self, lo, hi, *, snapshot):
    """Share validation, traversal and publication for both results."""
    require_published(self)
    checked_position(lo)
    checked_position(hi)
    if not 0 <= hi - lo <= MAX_WINDOW_WIDTH:
        raise ValueError("window must be ordered and at most 1000000 wide")
    stats = dict.fromkeys(
        (
            "blocks",
            "scanned",
            "candidates",
            "root_hits",
            "power_hits",
            "division_primes",
            "division_steps",
            "zeros",
            "composite_residuals",
            "admitted_atoms",
            "matches",
            "duplicates",
            "evictions",
            "dropped_partials",
            "slice_bytes",
            "slice_allocations",
            "threshold_max",
            "plan_builds",
            "plan_replays",
            "plan_refusals",
        ),
        0,
    )
    if self._graph is not None:
        stats.update(dict.fromkeys(self._double_counters, 0))
    stats["threshold_min"] = 2**63
    position, reason, divisor = (lo, "complete", None)
    span_budget = (
        self.budget if isinstance(self.budget, PollingBudget) else None
    )
    span_ledger = None if span_budget is None else span_budget.budget
    reporting = (
        self._reporting
        if getattr(self.config, "extension_policy", None)
        else None
    )
    kernel_before = kernels.counters(self)
    batcher = None
    if self._checked is not None and self.config.batch_policy == "shared-v2":
        from .cofactor_policy import CofactorBatch

        batcher = CofactorBatch(self.config.batch_size)
    try:
        if (
            (self.config.power_plan_bytes or self.config.sparse_queue_bytes)
            and lo < hi
            and (
                self._plan_window is None
                or lo < self._plan_window[0]
                or hi > self._plan_window[1]
            )
        ):
            self._set_plan_interval(lo, hi)
        while position < hi:
            block_lo = position
            block_hi = min(hi, block_lo + self.config.block_width)
            if getattr(self.config, "extension_policy", None):
                block_hi = min(
                    hi,
                    (position // self.config.block_width + 1)
                    * self.config.block_width,
                )
                if position % self.config.block_width == 0:
                    self._candidate_bound = self._reporting.decide()
            if batcher is not None:
                batcher.records.clear()
                batcher.hints.clear()
            threshold = self._sieve(block_lo, block_hi, stats)
            if self.config.division == "resieve":
                if not self._resieve(block_lo, block_hi, threshold, stats):
                    reason = "memory_limit"
                    break
            while position < block_hi:
                offset = position - block_lo
                selected = self._scores[offset] >= threshold
                if not selected and span_budget:
                    affordable = span_ledger.work_limit - span_ledger.used
                    amount = span_budget._remaining
                    if amount > affordable:
                        amount = affordable
                    if amount < 1:
                        amount = 1
                    stop = min(block_hi - block_lo, offset + amount)
                    end = offset + 1
                    while end < stop and self._scores[end] < threshold:
                        end += 1
                    amount = end - offset
                    span_budget.consume_span(amount)
                    position += amount
                    stats["scanned"] += amount
                    if reporting is not None:
                        reporting.record(amount, 0, 0)
                    continue
                if selected:
                    if not self._piecewise_passes(offset):
                        self.budget.consume()
                        atom, divisor = (None, None)
                    elif batcher is None:
                        atom, divisor = self._divide(position, offset, stats)
                    else:
                        atom, divisor = batcher.recover(
                            self, position, offset, block_hi, threshold, stats
                        )
                    if divisor is not None:
                        position += 1
                        stats["scanned"] += 1
                        if getattr(self.config, "extension_policy", None):
                            self._reporting.record(1, 1, 0)
                        reason = "factor_found"
                        break
                    if atom is not None:
                        refusal = self._admit(atom, stats)
                        if refusal:
                            reason = refusal
                            break
                else:
                    self.budget.consume()
                position += 1
                stats["scanned"] += 1
                if getattr(self.config, "extension_policy", None):
                    self._reporting.record(1, int(selected), 0)
            if reason != "complete":
                break
            if (
                getattr(self.config, "extension_policy", None)
                and position % self.config.block_width == 0
            ):
                self._candidate_bound = self._reporting.decide()
    except BudgetExhaustedError:
        reason = self.budget.reason
    if stats["blocks"] == 0:
        stats["threshold_min"] = 0
    self._resieved.clear()
    self._scratch_bytes = 0
    stats["plan_bytes"] = self._plan_bytes
    stats.update(
        {
            key: value - kernel_before.get(key, 0)
            for key, value in kernels.counters(self).items()
        }
    )
    if not snapshot:
        return _CollectionProgress(
            divisor,
            position,
            reason,
            stats,
            self._workspace + self._scratch_peak_bytes,
        )
    return SieveResult(
        tuple(self._atoms.values()),
        tuple(self._full),
        tuple(self._combined),
        self.partial_ids,
        divisor,
        position,
        reason,
        stats,
        self._workspace + self._scratch_peak_bytes,
    )


def build_native(collector, window, byte_cap):
    """Publish a whole schedule or use exact native fallback after refusal.

    Unlike a family-wide CRT scheduler, eligibility is relative to a working
    block. Every sparse root is enumerated once across this polynomial's
    interval. Higher powers retain the production Hensel and singular logic.
    """
    lo, hi, maximum = window
    width = collector.config.block_width
    roots = collector._roots
    split = collector._sparse_split
    count = (hi - lo + width - 1) // width
    reserve = 65536 + 1024 * count + 256 + 16 * split
    collector.budget.consume(len(roots) + count + 1)
    if reserve > byte_cap or count > 256:
        return None
    if collector._sparse_singular:
        return None
    word_max = (1 << 8 * array("I").itemsize) - 1
    weight_max = max(1, maximum.bit_length()) * max(collector._logs)
    if max(width - 1, len(roots) - 1, weight_max << 1 | 1) > word_max:
        return None
    queues = [array("I") for _ in range(count)]
    budget = collector.budget
    grouped = isinstance(budget, PollingBudget)
    for index in range(split, len(roots)):
        root = roots[index]
        marks = collector._power_roots(
            collector.polynomial,
            root,
            maximum,
            collector._logs[index],
            collector.budget,
            lo=lo,
            hi=hi,
            inverses=collector._lift_inverses,
            shared_inverses=collector._shared_lift_inverses,
            root_transfers=collector._power_transfers,
        )
        for modulus, residues, weight in marks:
            encoded_weight = weight << 1 | (modulus == root.prime)
            for residue in residues:
                hits = range(lo + (residue - lo) % modulus, hi, modulus)
                hit_count = len(hits)
                if not hit_count:
                    continue
                if reserve + ENTRY_BYTES * hit_count > byte_cap:
                    return None
                reserve += ENTRY_BYTES * hit_count
                block, offset = divmod(hits.start - lo, width)
                block_step, offset_step = (
                    divmod(modulus, width) if hit_count > 1 else (0, 0)
                )
                emitted = 0
                while emitted < hit_count:
                    if grouped:
                        remaining = budget._remaining or 1
                        affordable = max(
                            1,
                            (budget.budget.work_limit - budget.budget.used)
                            // 4,
                        )
                        chunk = min(hit_count - emitted, remaining, affordable)
                        if chunk == 1:
                            budget.consume(4)
                        else:
                            budget.consume_actions(chunk, 4)
                    else:
                        chunk = 1
                        budget.consume(4)
                    for _ in range(chunk):
                        queues[block].extend((offset, index, encoded_weight))
                        block += block_step
                        offset += offset_step
                        if offset >= width:
                            block += 1
                            offset -= width
                    emitted += chunk
    collector.budget.consume(count + 1)
    return (split, tuple(queues), reserve)


def sieve_native(collector, lo, hi, stats):
    """Keep dense marking and exact recovery, avoiding empty sparse visits."""
    cfg = collector.config
    window = collector._plan_window
    if (
        cfg.score_policy != "powers"
        or cfg.score_backend != "list"
        or cfg.division != "bucket"
        or cfg.small_prime_cutoff
        or cfg.threshold_extra
        or (not window)
        or (lo - window[0]) % cfg.block_width
        or (hi != min(window[1], lo + cfg.block_width))
        or collector._sparse_square
    ):
        return collector._sieve_direct(lo, hi, stats)
    key = (
        collector.polynomial.identity,
        window,
        cfg.block_width,
        cfg.sparse_queue_bytes,
    )
    if getattr(collector, "_sparse_queue_key", None) != key:
        clear(collector)
        cap = min(
            cfg.sparse_queue_bytes,
            max(0, cfg.memory_limit - collector._workspace),
        )
        collector._workspace += cap
        collector._queue_workspace = cap
        try:
            schedule = build_native(collector, window, cap)
            collector._sparse_dense_roots = (
                collector._roots[: schedule[0]] if schedule is not None else ()
            )
        except BaseException:
            clear(collector)
            raise
        retained = schedule[2] if schedule is not None else 0
        collector._workspace -= cap - retained
        collector._queue_workspace = retained
        collector._sparse_queue_key = key
        collector._sparse_queue = schedule
        if schedule is not None:
            collector._sparse_queue_peak = max(
                getattr(collector, "_sparse_queue_peak", 0), retained
            )
        stats["queue_builds" if schedule else "queue_refusals"] = 1
    schedule = collector._sparse_queue
    if schedule is None:
        return collector._sieve_direct(lo, hi, stats)
    split, queues, reserve = schedule
    entries = queues[(lo - window[0]) // cfg.block_width]
    roots = collector._roots
    try:
        collector._roots = collector._sparse_dense_roots
        threshold = collector._sieve_direct(lo, hi, stats)
    finally:
        collector._roots = roots
    hits = len(entries) // 3
    collector.budget.consume(hits + 1)
    first_count = 0
    masks = collector._hit_masks
    mask_count = len(masks)
    words = getattr(collector, "_word_hits", None)
    if collector._hits is not words:
        words = None
    for entry in range(0, len(entries), 3):
        offset, index, weight = (
            entries[entry],
            entries[entry + 1],
            entries[entry + 2],
        )
        collector._scores[offset] += weight >> 1
        if weight & 1:
            if words is not None:
                words.add_one(index, offset)
            else:
                mask = masks[index] if index < mask_count else 1 << index
                collector._hits[offset] |= mask
            first_count += 1
    stats["power_hits"] += hits
    stats["root_hits"] += first_count
    stats["queue_hits"] = stats.get("queue_hits", 0) + hits
    stats["queue_reserved_peak"] = max(
        stats.get("queue_reserved_peak", 0), reserve
    )
    stats["sparse_prime_visits_avoided"] = (
        stats.get("sparse_prime_visits_avoided", 0) + len(roots) - split
    )
    return threshold
