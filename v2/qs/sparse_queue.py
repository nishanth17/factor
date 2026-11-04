"""Bounded exact sparse hit schedules; derived state is rebuilt on resume."""

from array import array

from ..execution.work_budget import PollingBudget
from .power_sieve import prime_power_roots

HELPER_BYTES = 32 * 2**20
# Four copies of a three-word entry cover buffer growth and coexistence.
ENTRY_BYTES = 12 * array("I").itemsize


def build(collector, window, byte_cap):
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
    # Sparse singular roots may trigger conservative universal overmarks
    # which depend on the observation window. Keep the native route for them.
    if collector._sparse_singular:
        return None
    word_max = (1 << (8 * array("I").itemsize)) - 1
    weight_max = max(1, maximum.bit_length()) * max(collector._logs)
    if max(width - 1, len(roots) - 1, (weight_max << 1) | 1) > word_max:
        return None
    queues = [array("I") for _ in range(count)]
    budget = collector.budget
    grouped = isinstance(budget, PollingBudget)
    for index in range(split, len(roots)):
        root = roots[index]
        marks = prime_power_roots(
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
            encoded_weight = (weight << 1) | (modulus == root.prime)
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
                        # A polling callback may change the shared work
                        # limit. Keep its first event separate, then group
                        # only the remaining already-polled actions.
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
                        # Advance quotient/remainder by addition. Sparse
                        # moduli visit each block at most once; emitted words
                        # and the root traversal order remain unchanged.
                        queues[block].extend((offset, index, encoded_weight))
                        block += block_step
                        offset += offset_step
                        if offset >= width:
                            block += 1
                            offset -= width
                    emitted += chunk
    collector.budget.consume(count + 1)
    return split, tuple(queues), reserve


def clear(collector):
    """Release derived scratch only; retained relations remain untouched."""
    collector._sparse_queue_key = None
    collector._sparse_queue = None
    collector._sparse_dense_roots = ()
    collector._workspace -= getattr(collector, "_queue_workspace", 0)
    collector._queue_workspace = 0


def sieve(collector, lo, hi, stats):
    """Keep dense marking and exact recovery, avoiding empty sparse visits."""
    cfg = collector.config
    window = collector._plan_window
    if (
        cfg.score_policy != "powers"
        or cfg.score_backend != "list"
        or cfg.division != "bucket"
        or cfg.small_prime_cutoff
        or cfg.threshold_extra
        or not window
        or (lo - window[0]) % cfg.block_width
        or hi != min(window[1], lo + cfg.block_width)
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
        # Reserve construction coexistence before allocating any queue. The
        # published queue and every subsequent atom compete for this same cap.
        cap = min(
            cfg.sparse_queue_bytes,
            max(0, cfg.memory_limit - collector._workspace),
        )
        collector._workspace += cap
        collector._queue_workspace = cap
        try:
            schedule = build(collector, window, cap)
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
