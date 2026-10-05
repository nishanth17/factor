"""Bounded exact Hensel lifting for an optional prime-power sieve."""

MAX_LIFT_ROOTS = 64


def prime_power_roots(
    polynomial, roots, maximum, log, budget, *, lo=None, hi=None, inverses=None
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

    while current and modulus <= maximum:
        if lo is not None:
            budget.consume(len(current) * (modulus.bit_length() + 1))
            # Every higher lift is in this residue class. A class without
            # a hit in this window cannot acquire one at a higher power.
            current = tuple(
                root for root in current if lo + (root - lo) % modulus < hi
            )
            if not current:
                break

        yield modulus, current, log
        if modulus > maximum // prime:
            break
        lifted = []

        for root in current:
            budget.consume(
                polynomial.n_prime.bit_length() + 4 * modulus.bit_length() + 1
            )

            derivative = 2 * (polynomial.a * root + polynomial.b) % prime
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
            if derivative:
                key = prime, root % prime
                inverse = None if inverses is None else inverses.get(key)
                if inverse is None:
                    inverse = pow(derivative, -1, prime)
                    if inverses is not None:
                        inverses[key] = inverse
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

        current = tuple(lifted)
        modulus *= prime
