"""Shape bounds for bounded-yield ownership and shared recovery scratch."""


def graph_reserve(plan, choices):
    # Token tuples/integers, hash-table slots, lifting handles and per-handle
    # cycle-index growth. Original token objects are shared through XORs.
    return (
        1024
        + 256 * (len(plan.tokens) + len(plan.vector) + len(plan.provenance))
        + 64 * choices
    )


def base_batch_reserve(base, size, budget):
    budget.consume(len(base.primes) + 1)
    product_bits = sum(p.bit_length() for p in base.primes)
    return 2 * (65536 + 16 * product_bits) + size * 8192


def reserve_batch(collector, lo, hi):
    """Grow a scratch high-water mark before a physical block is processed."""
    _, maximum = collector._bounds(lo, hi)
    bits = (
        collector.polynomial.supported_a.bit_length()
        + max(1, maximum).bit_length()
    )
    product, support = 1, 0
    for prime in collector.factor_base.primes:
        collector.budget.consume(product.bit_length() + prime.bit_length())
        product *= prime
        if product.bit_length() > bits:
            break
        support += 1
    needed = collector._batch_base + collector.config.batch_size * (
        512 * support + 16 * bits
    )
    delta = max(0, needed - collector._batch_reserve)
    if (
        collector._workspace + delta + collector._scratch_bytes
        > collector.config.memory_limit
    ):
        raise MemoryError("bounded batch scratch exceeds memory cap")
    collector._batch_reserve += delta
    collector._workspace += delta
