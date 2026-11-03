"""Shape bounds for bounded-yield ownership and shared recovery scratch."""


def integer_reserve(value):
    """Bound a private integer header, digit array and allocator alignment.

    PyPy stores large integers in machine-word digits, using at least 31 bits
    per digit. Sixteen bytes per digit and 256 bytes of fixed overhead cover
    both digit storage and its object/array owners on the supported runtime.
    This is a conservative logical bound, not a sys.getsizeof/RSS estimate.
    """
    return 256 + 16 * ((abs(value).bit_length() + 30) // 31)


def polynomial_reserve(polynomial):
    # The arena owns one immutable context and one index entry per identity.
    # Shared n/multiplier integers are conservatively counted per context.
    return 1024 + sum(
        integer_reserve(getattr(polynomial, name))
        for name in ("n", "multiplier", "a", "b", "c", "square_coefficient")
    )


def row_reserve(leaves, base):
    # Packed factors have a separate slab owner. These charges cover the row
    # object, leaf/identity tuples, residual metadata and materialized scalars.
    return (
        4096
        + 512 * len(leaves)
        + 128 * sum(len(a.residual_primes) for a in leaves)
        + 2 * integer_reserve(base.n_prime)
    )


def pack_owned(collector, operation):
    """Keep allocated slabs charged when a later packing step refuses."""
    arena = collector._arena
    before = arena.allocated_bytes
    try:
        return operation()
    finally:
        # A failed plan may leave the current slab pinned by the allocator.
        # Keep that allocation charged; do not hide it in an uncommitted atom.
        delta = arena.allocated_bytes - before
        collector._workspace += delta
        collector._slab_bytes += delta


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
    if bits <= getattr(collector, "_batch_bits", 0):
        return
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
    # The required support and scalar storage are monotone in this bound.
    # A smaller norm cannot exceed the already owned high-water reservation.
    collector._batch_bits = bits
