"""Bounded lossless relation rows with legacy admission/work accounting."""

MAGIC = b"QSP1"
MAX_BITS = 32768


def norm_bits(config, n):
    """Bound a packet norm and its distinct-prime support conservatively."""
    return 2 + max(
        n.bit_length() + 20,
        2
        * (
            config.factor_count * config.base_bound.bit_length()
            + config.half_width.bit_length()
            + 2
        ),
    )


def legacy_length(row):
    """Exact compact-JSON length without constructing its string."""
    gray, position, sign, factors, residual, primes, opaque, extended = row
    factors_size = 2 + max(0, len(factors) - 1)
    factors_size += sum(3 + len(str(p)) + len(str(e)) for p, e in factors)
    primes_size = (
        2 + max(0, len(primes) - 1) + sum(len(str(p)) for p in primes)
    )
    return (
        9
        + len(str(gray))
        + len(str(position))
        + len(str(sign))
        + factors_size
        + len(str(residual))
        + primes_size
        + (4 if opaque is None else len(str(opaque)))
        + (4 if extended else 5)
    )


def _put(output, value):
    if type(value) is not int or value < 0 or value.bit_length() > MAX_BITS:
        raise ValueError("packed integer exceeds its finite bounds")
    while value >= 128:
        output.append((value & 127) | 128)
        value >>= 7
    output.append(value)


def _row_bytes(row):
    if not isinstance(row, list) or len(row) != 8:
        raise ValueError("invalid packed relation shape")
    gray, position, sign, factors, residual, primes, opaque, extended = row
    if (
        type(position) is not int
        or type(sign) is not int
        or sign not in (-1, 1)
    ):
        raise ValueError("invalid packed position/sign")
    if type(extended) is not bool:
        raise ValueError("invalid packed relation kind")
    output = bytearray()
    _put(output, gray)
    _put(output, 2 * position if position >= 0 else -2 * position - 1)
    output.append(int(sign < 0))
    _put(output, len(factors))
    for prime, exponent in factors:
        _put(output, prime)
        _put(output, exponent)
    _put(output, residual)
    _put(output, len(primes))
    for prime in primes:
        _put(output, prime)
    if opaque is not None and (type(opaque) is not int or opaque < 0):
        raise ValueError("invalid packed opaque cofactor")
    _put(output, 0 if opaque is None else opaque + 1)
    output.append(int(extended))
    return output


class PackedRows:
    """Reuse one finite buffer; published bytes never alias mutable scratch."""

    def __init__(self, capacity):
        if type(capacity) is not int or not 4096 <= capacity <= 16 * 2**20:
            raise ValueError("invalid packed buffer capacity")
        self.capacity = capacity
        self.buffer = bytearray(capacity)
        self.reset()

    def reset(self):
        self.position, self.count = 8, 0
        self.buffer[:4] = MAGIC

    def __len__(self):
        return self.count

    def append(self, row):
        encoded = _row_bytes(row)
        if self.position + len(encoded) > self.capacity:
            raise MemoryError("packed publication exceeds its prepaid buffer")
        self.buffer[self.position : self.position + len(encoded)] = encoded
        self.position += len(encoded)
        self.count += 1

    def finish(self):
        self.buffer[4:8] = self.count.to_bytes(4, "little")
        return bytes(memoryview(self.buffer)[: self.position])


def decode(raw, *, max_rows, max_exponents, max_bits, legacy_limit):
    """Reject malformed dimensions before allocating unbounded row objects."""
    if (
        type(raw) is not bytes
        or len(raw) < 8
        or raw[:4] != MAGIC
        or len(raw) > legacy_limit
        or not 1 <= max_bits <= MAX_BITS
    ):
        raise ValueError("invalid packed packet header or allowance")
    count = int.from_bytes(raw[4:8], "little")
    if count > max_rows:
        raise ValueError("packed atom cap exceeded")
    offset = 8

    def take():
        nonlocal offset
        value, shift = 0, 0
        while True:
            if offset >= len(raw) or shift >= max_bits:
                raise ValueError("truncated or oversized packed integer")
            member = raw[offset]
            offset += 1
            value |= (member & 127) << shift
            if not member & 128:
                if (shift and member == 0) or value.bit_length() > max_bits:
                    raise ValueError(
                        "noncanonical or oversized packed integer"
                    )
                return value
            shift += 7

    rows, size = [], 2
    for _ in range(count):
        gray, position, negative, exponents = take(), take(), take(), take()
        if negative > 1 or exponents > max_exponents:
            raise ValueError("invalid packed sign or exponent dimension")
        factors = [[take(), take()] for _ in range(exponents)]
        residual, endpoints = take(), take()
        if endpoints > 3:
            raise ValueError("packed residual arity exceeds its bound")
        primes = [take() for _ in range(endpoints)]
        opaque, extended = take(), take()
        if extended > 1:
            raise ValueError("invalid packed relation kind")
        row = [
            gray,
            position // 2 if not position & 1 else -(position // 2) - 1,
            -1 if negative else 1,
            factors,
            residual,
            primes,
            None if opaque == 0 else opaque - 1,
            bool(extended),
        ]
        size += legacy_length(row) + 1
        if size > legacy_limit:
            raise ValueError("packed rows exceed legacy admission allowance")
        rows.append(row)
    if offset != len(raw):
        raise ValueError("trailing packed relation bytes")
    return rows


def writer_for(state, capacity):
    """Scratch lasts only until the owned pool's quiescent finish."""
    writer = state.get("_qs_writer")
    if writer is None or writer.capacity != capacity:
        state.pop("_qs_writer", None)
        writer = state["_qs_writer"] = PackedRows(capacity)
    writer.reset()
    return writer
