"""Unpruned Hensel traversal selected once at collector construction.

The generic power_sieve implementation remains the opt-in pruning oracle.
Keep arithmetic, paid prefixes and singular/window behavior identical.
"""

from ..common.arithmetic import pow
from .power_sieve import MAX_LIFT_ROOTS


def prime_power_roots(
    polynomial,
    roots,
    maximum,
    log,
    budget,
    *,
    lo=None,
    hi=None,
    inverses=None,
    shared_inverses=None,
    root_transfers=None,
):
    """Yield modulus/root/ceil-log marks covering every p-adic valuation.

    Lift nonsingular roots uniquely. Singular roots branch only within a
    fixed workspace cap; if that cap is reached, an upper-bound mark at the
    original modulus covers every remaining valuation. False candidates are
    safe; omitted smooth candidates are forbidden. Inputs are certified by
    the collector. A root list has at most 64 entries and work is charged
    before computing any next level.
    """
    if (lo is None) != (hi is None):
        raise ValueError("both window endpoints are required")
    prime = roots.prime
    initial = roots.roots
    if roots.all_positions:
        if prime > MAX_LIFT_ROOTS:
            yield 1, (0,), maximum.bit_length() * log
            return
        initial = tuple(range(prime))

    modulus, current = prime, initial
    target_bits = None
    transferable = (
        root_transfers is not None
        and prime > 2
        and not roots.all_positions
        and len(initial) == 2
    )

    while current and modulus <= maximum:
        modulus_bits = modulus.bit_length()
        if lo is not None:
            budget.consume(len(current) * (modulus_bits + 1))
            # Every higher lift is in this residue class. A class without
            # a hit in this window cannot acquire one at a higher power.
            # A complete period contains every residue; retain the original
            # tuple instead of filtering and copying it in that case.
            if hi - lo < modulus:
                current = tuple(
                    root for root in current if lo + (root - lo) % modulus < hi
                )
                if not current:
                    break

        yield modulus, current, log
        # Python integers are exact. The temporary is bounded by maximum's
        # bit length plus one prime and replaces a division at every level.
        next_modulus = modulus * prime
        if next_modulus > maximum:
            break
        if target_bits is None:
            target_bits = polynomial.n_prime.bit_length()
        cost = target_bits + 4 * modulus_bits + 1
        if transferable and len(current) == 2:
            shifted = root_transfers.shifted(
                polynomial, prime, next_modulus, initial, cost, budget
            )
            if shifted is not None:
                current, modulus = shifted, next_modulus
                continue
        lifted = []

        for root in current:
            budget.consume(cost)

            # The derivative modulo p is constant within an original root
            # class. Reuse its inverse (or singular zero) across all lifts.
            key = prime, root % prime
            inverse = None if inverses is None else inverses.get(key)
            if inverse is None:
                derivative = 2 * (polynomial.a * root + polynomial.b) % prime
                inverse = (
                    shared_inverses.inverse(prime, derivative)
                    if shared_inverses is not None
                    else pow(derivative, -1, prime)
                    if derivative
                    else 0
                )
                if inverses is not None:
                    inverses[key] = inverse
            # Derived lift residues are bounded by maximum, rather than the
            # public position cap. Their invariant is F(root) == 0 mod modulus.
            quotient = (
                (
                    (polynomial.a * root + 2 * polynomial.b) * root
                    + polynomial.c
                )
                // modulus
                % prime
            )
            if inverse:
                correction = -quotient * inverse % prime
                lifted.append(root + modulus * correction)
            elif quotient == 0:
                if len(lifted) + prime > MAX_LIFT_ROOTS:
                    # Earlier levels have already added positive scores.
                    # This universal remaining allowance can overestimate.
                    yield prime, initial, maximum.bit_length() * log
                    return

                lifted.extend(
                    root + modulus * offset for offset in range(prime)
                )

        complete_pair = transferable and len(current) == len(lifted) == 2
        current = tuple(lifted)
        if complete_pair:
            root_transfers.remember(
                polynomial, modulus, next_modulus, current, budget
            )
        modulus = next_modulus
