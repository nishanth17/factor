"""Prime sieves with a uniform half-open endpoint contract.

prime_sieve(hi), small_sieve(hi), and sieve_of_atkin(hi) emit p < hi.
segmented_sieve(lo, hi) emits lo <= p < hi. All marking state is local.
The original mod-60 Atkin recurrences/tables are retained with exact integer
arithmetic, corrected emission, and shortened final segments.
"""

from math import isqrt

from .. import constants
from . import utils

UNDER_60 = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59)
ATKIN_RESIDUES = (1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59)

# Solutions of 4*f²+g²=d (mod 60), f<=15, g<=30.
DFG_1 = [
    [1, 0, 1],
    [1, 0, 11],
    [1, 0, 19],
    [1, 0, 29],
    [1, 2, 15],
    [1, 3, 5],
    [1, 3, 25],
    [1, 5, 9],
    [1, 5, 21],
    [1, 7, 15],
    [1, 8, 15],
    [1, 10, 9],
    [1, 10, 21],
    [1, 12, 5],
    [1, 12, 25],
    [1, 13, 15],
    [13, 1, 3],
    [13, 1, 27],
    [13, 4, 3],
    [13, 4, 27],
    [13, 6, 7],
    [13, 6, 13],
    [13, 6, 17],
    [13, 6, 23],
    [13, 9, 7],
    [13, 9, 13],
    [13, 9, 17],
    [13, 9, 23],
    [13, 11, 3],
    [13, 11, 27],
    [13, 14, 3],
    [13, 14, 27],
    [17, 2, 1],
    [17, 2, 11],
    [17, 2, 19],
    [17, 2, 29],
    [17, 7, 1],
    [17, 7, 11],
    [17, 7, 19],
    [17, 7, 29],
    [17, 8, 1],
    [17, 8, 11],
    [17, 8, 19],
    [17, 8, 29],
    [17, 13, 1],
    [17, 13, 11],
    [17, 13, 19],
    [17, 13, 29],
    [29, 1, 5],
    [29, 1, 25],
    [29, 4, 5],
    [29, 4, 25],
    [29, 5, 7],
    [29, 5, 13],
    [29, 5, 17],
    [29, 5, 23],
    [29, 10, 7],
    [29, 10, 13],
    [29, 10, 17],
    [29, 10, 23],
    [29, 11, 5],
    [29, 11, 25],
    [29, 14, 5],
    [29, 14, 25],
    [37, 2, 9],
    [37, 2, 21],
    [37, 3, 1],
    [37, 3, 11],
    [37, 3, 19],
    [37, 3, 29],
    [37, 7, 9],
    [37, 7, 21],
    [37, 8, 9],
    [37, 8, 21],
    [37, 12, 1],
    [37, 12, 11],
    [37, 12, 19],
    [37, 12, 29],
    [37, 13, 9],
    [37, 13, 21],
    [41, 2, 5],
    [41, 2, 25],
    [41, 5, 1],
    [41, 5, 11],
    [41, 5, 19],
    [41, 5, 29],
    [41, 7, 5],
    [41, 7, 25],
    [41, 8, 5],
    [41, 8, 25],
    [41, 10, 1],
    [41, 10, 11],
    [41, 10, 19],
    [41, 10, 29],
    [41, 13, 5],
    [41, 13, 25],
    [49, 0, 7],
    [49, 0, 13],
    [49, 0, 17],
    [49, 0, 23],
    [49, 1, 15],
    [49, 4, 15],
    [49, 5, 3],
    [49, 5, 27],
    [49, 6, 5],
    [49, 6, 25],
    [49, 9, 5],
    [49, 9, 25],
    [49, 10, 3],
    [49, 10, 27],
    [49, 11, 15],
    [49, 14, 15],
    [53, 1, 7],
    [53, 1, 13],
    [53, 1, 17],
    [53, 1, 23],
    [53, 4, 7],
    [53, 4, 13],
    [53, 4, 17],
    [53, 4, 23],
    [53, 11, 7],
    [53, 11, 13],
    [53, 11, 17],
    [53, 11, 23],
    [53, 14, 7],
    [53, 14, 13],
    [53, 14, 17],
    [53, 14, 23],
]

# Solutions of 3*f²+g²=d (mod 60), f<=10, g<=30.
DFG_2 = [
    [7, 1, 2],
    [7, 1, 8],
    [7, 1, 22],
    [7, 1, 28],
    [7, 3, 10],
    [7, 3, 20],
    [7, 7, 10],
    [7, 7, 20],
    [7, 9, 2],
    [7, 9, 8],
    [7, 9, 22],
    [7, 9, 28],
    [19, 1, 4],
    [19, 1, 14],
    [19, 1, 16],
    [19, 1, 26],
    [19, 5, 2],
    [19, 5, 8],
    [19, 5, 22],
    [19, 5, 28],
    [19, 9, 4],
    [19, 9, 14],
    [19, 9, 16],
    [19, 9, 26],
    [31, 3, 2],
    [31, 3, 8],
    [31, 3, 22],
    [31, 3, 28],
    [31, 5, 4],
    [31, 5, 14],
    [31, 5, 16],
    [31, 5, 26],
    [31, 7, 2],
    [31, 7, 8],
    [31, 7, 22],
    [31, 7, 28],
    [43, 1, 10],
    [43, 1, 20],
    [43, 3, 4],
    [43, 3, 14],
    [43, 3, 16],
    [43, 3, 26],
    [43, 7, 4],
    [43, 7, 14],
    [43, 7, 16],
    [43, 7, 26],
    [43, 9, 10],
    [43, 9, 20],
]

# Solutions of 3*f²-g²=d (mod 60), f<=10, g<=30. Seeds may be negative.
DFG_3 = [
    [11, 0, 7],
    [11, 0, 13],
    [11, 0, 17],
    [11, 0, 23],
    [11, 2, 1],
    [11, 2, 11],
    [11, 2, 19],
    [11, 2, 29],
    [11, 3, 4],
    [11, 3, 14],
    [11, 3, 16],
    [11, 3, 26],
    [11, 5, 2],
    [11, 5, 8],
    [11, 5, 22],
    [11, 5, 28],
    [11, 7, 4],
    [11, 7, 14],
    [11, 7, 16],
    [11, 7, 26],
    [11, 8, 1],
    [11, 8, 11],
    [11, 8, 19],
    [11, 8, 29],
    [23, 1, 10],
    [23, 1, 20],
    [23, 2, 7],
    [23, 2, 13],
    [23, 2, 17],
    [23, 2, 23],
    [23, 3, 2],
    [23, 3, 8],
    [23, 3, 22],
    [23, 3, 28],
    [23, 4, 5],
    [23, 4, 25],
    [23, 6, 5],
    [23, 6, 25],
    [23, 7, 2],
    [23, 7, 8],
    [23, 7, 22],
    [23, 7, 28],
    [23, 8, 7],
    [23, 8, 13],
    [23, 8, 17],
    [23, 8, 23],
    [23, 9, 10],
    [23, 9, 20],
    [47, 1, 4],
    [47, 1, 14],
    [47, 1, 16],
    [47, 1, 26],
    [47, 2, 5],
    [47, 2, 25],
    [47, 3, 10],
    [47, 3, 20],
    [47, 4, 1],
    [47, 4, 11],
    [47, 4, 19],
    [47, 4, 29],
    [47, 6, 1],
    [47, 6, 11],
    [47, 6, 19],
    [47, 6, 29],
    [47, 7, 10],
    [47, 7, 20],
    [47, 8, 5],
    [47, 8, 25],
    [47, 9, 4],
    [47, 9, 14],
    [47, 9, 16],
    [47, 9, 26],
    [59, 0, 1],
    [59, 0, 11],
    [59, 0, 19],
    [59, 0, 29],
    [59, 1, 2],
    [59, 1, 8],
    [59, 1, 22],
    [59, 1, 28],
    [59, 4, 7],
    [59, 4, 13],
    [59, 4, 17],
    [59, 4, 23],
    [59, 5, 4],
    [59, 5, 14],
    [59, 5, 16],
    [59, 5, 26],
    [59, 6, 7],
    [59, 6, 13],
    [59, 6, 17],
    [59, 6, 23],
    [59, 9, 2],
    [59, 9, 8],
    [59, 9, 22],
    [59, 9, 28],
]


def small_sieve(hi):
    """Exact wheel-six Eratosthenes: list p < hi with sliced byte marking."""
    utils.require_integer(hi, "hi")
    if hi <= 3:
        return [2] if hi == 3 else []
    # Slot i represents (3*i+1)|1, the alternating residues 6*k-1/6*k+1.
    # The final slot count excludes hi exactly, without filtering each prime.
    size = 2 * (hi // 6) + int(hi % 6 > 1)
    flags = bytearray(b"\x01") * size
    flags[0] = 0  # The first wheel slot is 1, not a prime.

    for index in range(1, isqrt(hi - 1) // 3 + 1):
        if flags[index]:
            prime = (3 * index + 1) | 1
            step = 2 * prime
            # Two progressions mark multiples in both surviving residues.
            starts = (
                prime * prime // 3,
                (prime * prime + 4 * prime - 2 * prime * (index & 1)) // 3,
            )
            for start in starts:
                if start < size:
                    count = (size - 1 - start) // step + 1
                    flags[start::step] = b"\x00" * count

    return [2, 3] + [
        (3 * index + 1) | 1 for index in range(1, size) if flags[index]
    ]


def segmented_sieve(lo, hi, *, segment_size=constants.LOWER_SEGMENT_SIZE):
    """Materialize primes in [lo, hi); only the marking buffer is segmented.

    This API does not cap output or cache base primes across calls.
    Such scheduling/context work belongs to Phase 2.
    """
    utils.require_integer(lo, "lo")
    utils.require_integer(hi, "hi")
    utils.require_integer(segment_size, "segment_size", 1)
    if hi <= max(lo, 2):
        return []
    lo = max(lo, 2)
    base_primes = small_sieve(isqrt(hi - 1) + 1)
    result = [2] if lo <= 2 < hi else []
    # Index i represents left+2*i. A stride of p skips even multiples.
    start = max(lo, 3) | 1

    for left in range(start, hi, 2 * segment_size):
        right = min(hi, left + 2 * segment_size)
        size = (right - left + 1) // 2
        flags = bytearray(b"\x01") * size

        for prime in base_primes[1:]:
            if prime * prime >= right:
                break
            first = max(prime * prime, ((left + prime - 1) // prime) * prime)
            if first % 2 == 0:
                first += prime
            index = (first - left) // 2
            if index < size:
                count = (size - 1 - index) // prime + 1
                flags[index::prime] = b"\x00" * count

        result.extend(left + 2 * i for i in range(size) if flags[i])

    return result


def _enumerate_quadratic_1(residue, f, g, start, width, rows):
    """Toggle 4*x²+y²=60*k+d using x steps of 15 and y steps of 30."""
    x, y0, end = f, g, start + width
    k0 = (4 * f * f + g * g - residue) // 60
    while k0 < end:
        k0 += 2 * x + 15
        x += 15
    while True:
        x -= 15
        k0 -= 2 * x + 15
        if x <= 0:
            return
        while k0 < start:
            k0 += y0 + 15
            y0 += 30
        k, y = k0, y0
        while k < end:
            rows[residue][(k - start) >> 5] ^= 1 << ((k - start) & 31)
            k += y + 15
            y += 30


def _enumerate_quadratic_2(residue, f, g, start, width, rows):
    """Toggle 3*x²+y²=60*k+d using x steps of 10 and y steps of 30."""
    x, y0, end = f, g, start + width
    k0 = (3 * f * f + g * g - residue) // 60
    while k0 < end:
        k0 += x + 5
        x += 10
    while True:
        x -= 10
        k0 -= x + 5
        if x <= 0:
            return
        while k0 < start:
            k0 += y0 + 15
            y0 += 30
        k, y = k0, y0
        while k < end:
            rows[residue][(k - start) >> 5] ^= 1 << ((k - start) & 31)
            k += y + 15
            y += 30


def _enumerate_quadratic_3(residue, f, g, start, width, rows):
    """Toggle 3*x²-y²=60*k+d; preserve signed k seeds and require y<x."""
    x, y0, end = f, g, start + width
    k0 = (3 * f * f - g * g - residue) // 60

    while True:
        while k0 >= end:
            if x <= y0:
                return
            k0 -= y0 + 15
            y0 += 30

        k, y = k0, y0
        while k >= start and y < x:
            rows[residue][(k - start) >> 5] ^= 1 << ((k - start) & 31)
            k -= y + 15
            y += 30
        k0 += x + 5
        x += 10


def sieve_of_atkin(hi):
    """Corrected segmented mod-60 Atkin with private residue rows."""
    utils.require_integer(hi, "hi")
    result = [prime for prime in UNDER_60 if prime < hi]
    if hi <= 61:
        return result
    root = isqrt(hi - 1)
    base_primes = small_sieve(root + 1)
    block_size = 60 * root
    last_k = (hi - 1) // 60

    for start in range(1, last_k + 1, block_size):
        width = min(block_size, last_k - start + 1)
        rows = {
            residue: [0] * ((width + 31) // 32) for residue in ATKIN_RESIDUES
        }

        for table, enumerate_points in (
            (DFG_1, _enumerate_quadratic_1),
            (DFG_2, _enumerate_quadratic_2),
            (DFG_3, _enumerate_quadratic_3),
        ):
            for residue, f, g in table:
                enumerate_points(residue, f, g, start, width, rows)

        for prime in base_primes:
            if prime < 7:
                continue
            square = prime * prime
            inverse = utils.modular_inverse(60, square)
            for residue in ATKIN_RESIDUES:
                # Solve 60*(start+offset)+residue == 0 modulo prime².
                offset = -(60 * start + residue) * inverse % square
                for index in range(offset, width, square):
                    rows[residue][index >> 5] &= ~(1 << (index & 31))

        for index in range(width):
            scaled = 60 * (start + index)

            for residue in ATKIN_RESIDUES:
                candidate = scaled + residue
                if candidate >= hi:
                    return result
                if rows[residue][index >> 5] & (1 << (index & 31)):
                    result.append(candidate)

    return result


def prime_sieve(hi):
    """List primes below hi; retain Eratosthenes as the default backend."""
    utils.require_integer(hi, "hi")
    if hi <= constants.SMALL_THRESHOLD:
        return [prime for prime in UNDER_60 if prime < hi]
    if hi <= constants.ERAT_THRESHOLD:
        return small_sieve(hi)
    return segmented_sieve(2, hi)
