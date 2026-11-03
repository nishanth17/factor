"""Bounded centered-Gray root tables; construction policies are exact peers."""

from ..common import utils

ROOT_UPDATE_POLICIES = ("standard", "cached", "batch")


def root_update_workspace(base_count, factor_count, limit):
    """Reserve table/list coexistence and per-prime batch-inversion scratch."""
    coordinates = (limit - 1).bit_length()
    if not coordinates:
        return 0
    return 8192 + base_count * (64 + 128 * coordinates) + 1024 * factor_count


def build_root_updates(family, policy, limit):
    """Build only coordinates consumed by the declared adjacent Gray prefix.

    For t_j=(A/q_j)*gamma_j, 2*t_j/A = 2*gamma_j/q_j modulo p.
    The batch path uses the already cached A inverse, including A=2*prod(q)
    for parity families. No additional inversions or unbounded caches occur.
    Coordinate permutation and signed terms must already have been applied.
    """
    if policy not in ROOT_UPDATE_POLICIES:
        raise ValueError("unknown root-update policy")
    utils.require_integer(limit, "root-update limit", 1)
    if limit > family.count:
        raise ValueError("root-update limit exceeds the family")
    if policy == "standard" or limit == 1:
        return (), 0
    extra = root_update_workspace(
        len(family.base.entries), len(family.terms), limit
    )
    if family.workspace_bytes + extra > family.memory_bytes:
        raise MemoryError("root-update table exceeds memory_bytes")
    count = (limit - 1).bit_length()
    family.budget.consume((count + 1) * (family.a.bit_length() + 1))
    rows = [[] for _ in range(count)]
    gammas = tuple(
        term // (family.a // prime)
        for term, prime in zip(family.terms, family.a_primes)
    )
    for entry, inverse in zip(family.base.entries, family.inverses):
        prime = entry.prime
        family.budget.consume(
            (len(family.terms) + count + 1)
            * (family.a.bit_length() + prime.bit_length() ** 2)
        )
        if inverse is None:
            for row in rows:
                row.append(None)
        elif policy == "cached":
            for index, row in enumerate(rows, 1):
                row.append(2 * family.terms[index] * inverse % prime)
        else:
            prefixes = [1]
            for factor in family.a_primes:
                prefixes.append(prefixes[-1] * factor % prime)
            # Multiplying A^{-1} by 2 gives (A/2)^{-1}, not its half.
            reciprocal = inverse * (2 if family.parity else 1) % prime
            for index in range(len(family.terms) - 1, -1, -1):
                factor_inverse = reciprocal * prefixes[index] % prime
                if 1 <= index <= count:
                    rows[index - 1].append(
                        2 * gammas[index] * factor_inverse % prime
                    )
                reciprocal = reciprocal * family.a_primes[index] % prime
    return tuple(tuple(row) for row in rows), extra
