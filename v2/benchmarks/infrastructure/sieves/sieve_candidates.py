"""Diagnostic wheel/pre-sieve/bitset candidates; production stays readable."""

from bisect import bisect_left
from math import isqrt

from ....common import prime_sieve

RESIDUES = (1, 7, 11, 13, 17, 19, 23, 29)
RESIDUE_INDEX = {value: index for index, value in enumerate(RESIDUES)}
PATTERN = bytes(
    int(all((2 * i + 1) % p for p in (3, 5, 7))) for i in range(105)
)


def wheel_thirty(lo, hi):
    """Exact wheel-30 control with explicit Python per-strike indexing."""
    lo = max(lo, 2)
    if hi <= lo:
        return []
    first = lo // 30 * 8 + bisect_left(RESIDUES, lo % 30)
    stop = hi // 30 * 8 + bisect_left(RESIDUES, hi % 30)
    flags = bytearray(b"\x01") * (stop - first)

    for prime in prime_sieve.small_sieve(isqrt(hi - 1) + 1):
        if prime <= 5:
            continue
        start = max(prime * prime, ((lo + prime - 1) // prime) * prime)
        for value in range(start, hi, prime):
            residue = RESIDUE_INDEX.get(value % 30)
            if residue is not None:
                flags[value // 30 * 8 + residue - first] = 0

    result = [prime for prime in (2, 3, 5) if lo <= prime < hi]
    for index, survives in enumerate(flags):
        slot = first + index
        value = slot // 8 * 30 + RESIDUES[slot % 8]
        if survives and value >= 7:
            result.append(value)

    return result


def presieved(lo, hi, *, segment_size=4096):
    """Copy a phase-correct odd 3*5*7 pattern, then slice-mark other primes."""
    lo = max(lo, 2)
    if hi <= lo:
        return []
    result = [prime for prime in (2, 3, 5, 7) if lo <= prime < hi]
    bases = prime_sieve.small_sieve(isqrt(hi - 1) + 1)

    for left in range(max(lo, 3) | 1, hi, 2 * segment_size):
        right = min(hi, left + 2 * segment_size)
        size = (right - left + 1) // 2
        phase = (left // 2) % len(PATTERN)
        stop = phase + size
        flags = bytearray((PATTERN * ((size + phase) // 105 + 1))[phase:stop])

        for prime in bases:
            if prime <= 7:
                continue
            if prime * prime >= right:
                break
            start = max(prime * prime, ((left + prime - 1) // prime) * prime)
            if start % 2 == 0:
                start += prime
            index = (start - left) // 2
            if index < size:
                flags[index::prime] = b"\x00" * (
                    (size - 1 - index) // prime + 1
                )

        result.extend(
            left + 2 * i for i in range(size) if flags[i] and left + 2 * i > 7
        )

    return sorted(result)


def integer_bitset(lo, hi):
    """Geometric bit masks mark odd multiples; include full extraction cost."""
    lo = max(lo, 2)
    if hi <= lo:
        return []
    left = max(lo, 3) | 1
    size = max(0, (hi - left + 1) // 2)
    flags = (1 << size) - 1

    for prime in prime_sieve.small_sieve(isqrt(hi - 1) + 1):
        if prime == 2:
            continue
        start = max(prime * prime, ((left + prime - 1) // prime) * prime)
        if start % 2 == 0:
            start += prime
        index = (start - left) // 2
        if index < size:
            count = (size - 1 - index) // prime + 1
            strikes = ((1 << (prime * count)) - 1) // ((1 << prime) - 1)
            flags &= ~(strikes << index)

    result = [2] if lo <= 2 < hi else []
    while flags:
        bit = flags & -flags
        result.append(left + 2 * (bit.bit_length() - 1))
        flags ^= bit
    return result
