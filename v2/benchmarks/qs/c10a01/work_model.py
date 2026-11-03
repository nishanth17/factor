"""Conservative analytic postprocessing charges for the frozen capacities."""

from .configuration import COUNTS, ROWS

RATE = 8_000_000_000
SETUP = 10_000_000_000


def charges(digits):
    """Bound all passes without assuming sparse fill or cache hits."""
    r, c = ROWS[digits], COUNTS[digits] + 1
    # This integer bound exceeds log2(10**digits), without floating arithmetic.
    bits = (digits * 3322 + 999) // 1000
    words = 1 + (c + 63) // 64 + (r + 63) // 64
    atom = 513 * 512 + 2 * 64**2
    preparation = r * (257 + 256 * (atom + c))
    compact = ((c + 63) // 64 + 1) * (c + r * c)
    compact += ((c + 63) // 64) * r * c + r
    filtering = 5 * r * c * words + 2 * r + compact
    solving = r * (c + 1) * (r + c + 1)
    # One dependency can select all rows. Include cache-key and extraction
    # charges; every returned dependency is bounded by the original row count.
    extraction = r * (r * c * (bits + 2) + 1)
    per_pass = preparation + filtering + solving + extraction
    # Shared full/combined rows only grow; unchanged counts do not re-prepare.
    # Include initial empty preparation and one final forced attempt.
    passes = r + 2
    return dict(
        bits=bits,
        rows=r,
        columns=c,
        atom=atom,
        preparation=preparation,
        filtering=filtering,
        solving=solving,
        extraction=extraction,
        passes=passes,
        postprocessing=passes * per_pass,
        largest_atomic=max(r * c * (bits + 1), r * c * words),
    )


def work_limit(digits, seconds):
    """Rate headroom for collection plus separately bounded postprocessing."""
    return 4 * RATE * seconds + SETUP + charges(digits)["postprocessing"]
