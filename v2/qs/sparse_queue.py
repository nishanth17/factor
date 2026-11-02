"""Bounded exact sparse hit schedules; derived state is rebuilt on resume."""

from math import isqrt

from .power_sieve import prime_power_roots

HELPER_BYTES = 32 * 2**20
ENTRY_BYTES = 384


def build(collector, window, byte_cap):
    """Publish a whole schedule or use exact native fallback after refusal.

    Unlike a family-wide CRT scheduler, eligibility is relative to a working
    block. Every sparse root is enumerated once across this polynomial's
    interval. Higher powers retain the production Hensel and singular logic.
    """
    lo, hi, maximum = window
    width = collector.config.block_width
    roots = collector._roots
    split = next(
        (i for i, r in enumerate(roots) if r.prime >= width), len(roots)
    )
    count = (hi - lo + width - 1) // width
    reserve = 65536 + 1024 * count
    collector.budget.consume(len(roots) + count + 1)
    if reserve > byte_cap or count > 256:
        return None
    # Sparse singular roots may trigger conservative universal overmarks
    # which depend on the observation window. Keep the native route for them.
    if any(collector.polynomial.n_prime % r.prime == 0 for r in roots[split:]):
        return None
    queues = [[] for _ in range(count)]
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
        )
        for modulus, residues, weight in marks:
            for residue in residues:
                hits = range(lo + (residue - lo) % modulus, hi, modulus)
                first = modulus == root.prime
                for position in hits:
                    collector.budget.consume(4)
                    if reserve + ENTRY_BYTES > byte_cap:
                        return None
                    reserve += ENTRY_BYTES
                    block, offset = divmod(position - lo, width)
                    queues[block].append((offset, index, weight, first))
    collector.budget.consume(count + 1)
    return split, tuple(tuple(q) for q in queues), reserve


def clear(collector):
    """Release derived scratch only; retained relations remain untouched."""
    collector._sparse_queue_key = None
    collector._sparse_queue = None
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
        or isqrt(collector.polynomial.n_prime) ** 2
        == collector.polynomial.n_prime
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
        collector._roots = roots[:split]
        threshold = collector._sieve_direct(lo, hi, stats)
    finally:
        collector._roots = roots
    collector.budget.consume(len(entries) + 1)
    first_count = 0
    for offset, index, weight, first in entries:
        collector._scores[offset] += weight
        if first:
            collector._hits[offset] |= 1 << index
            first_count += 1
    stats["power_hits"] += len(entries)
    stats["root_hits"] += first_count
    stats["queue_hits"] = stats.get("queue_hits", 0) + len(entries)
    stats["queue_reserved_peak"] = max(
        stats.get("queue_reserved_peak", 0), reserve
    )
    stats["sparse_prime_visits_avoided"] = (
        stats.get("sparse_prime_visits_avoided", 0) + len(roots) - split
    )
    return threshold
