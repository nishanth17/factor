"""Integer bounds for fixed-point base-two sieve scores."""

SCORE_SCALE = 32
MANTISSA_BITS = 10
MANTISSA_UNIT = 1 << MANTISSA_BITS


def _exact_bounds(value):
    """Bound 32*log2(value) using a finite integer power."""
    lower = (value**SCORE_SCALE).bit_length() - 1
    return lower, lower + bool(value & (value - 1))


_MANTISSA = tuple(
    _exact_bounds(value)
    for value in range(MANTISSA_UNIT, 2 * MANTISSA_UNIT + 1)
)


def log_bounds(value):
    """Return conservative integer lower/upper scores for a positive integer.

    Large values use their top eleven bits. Adjacent mantissa bins bound the
    exact value, so neither floating-point arithmetic nor a large power is
    needed per candidate. The table includes its upper endpoint.
    """
    if value <= 0:
        raise ValueError("log score requires a positive integer")
    shift = max(0, value.bit_length() - 1 - MANTISSA_BITS)
    if not shift:
        return _exact_bounds(value)
    mantissa = value >> shift
    lower, upper = _MANTISSA[mantissa - MANTISSA_UNIT]
    if value != mantissa << shift:
        upper = _MANTISSA[mantissa + 1 - MANTISSA_UNIT][1]
    return SCORE_SCALE * shift + lower, SCORE_SCALE * shift + upper
