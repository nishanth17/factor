"""One bounded parent CRT template; derived state is never checkpointed."""

from math import prod

from .families import MAX_A_FACTORS
from .polynomial import Polynomial
from .sign_schedule import schedule_for


class CoefficientCache:
    """Reuse selected roots, schedule and CRT terms within one packet merge.

    Terms become reusable only after the same per-coordinate reservation as
    uncached reconstruction. Every polynomial retains those reservations, so
    a cache hit cannot change refusal/cancellation boundaries or paid work.
    """

    def __init__(self, base, primes, config):
        if (
            not isinstance(primes, tuple)
            or not 1 <= len(primes) <= MAX_A_FACTORS
        ):
            raise ValueError("finite immutable A factors required")
        self.base, self.primes, self.config = base, primes, config
        self.parity = (
            config.polynomial_policy == "parity2" and base.n_prime % 4 == 1
        )
        self.a = prod(primes) * (2 if self.parity else 1)
        # The factor base already owns this checked immutable column map.
        # Avoid a full-base membership scan for each of up to64 Gray steps.
        self.roots = tuple(
            base.entries[base._columns[prime] - 1].square_roots[0]
            for prime in primes
        )
        self.schedule = schedule_for(
            len(primes),
            "swap_high"
            if config.sign_policy == "swap_high" and len(primes) >= 3
            else "control",
        )
        self.terms = [None] * len(primes)

    def polynomial_at(self, base, primes, ordinal, config, budget):
        if (
            base is not self.base
            or config is not self.config
            or primes != self.primes
        ):
            raise ValueError("coefficient template belongs to another context")
        a = self.a
        mask = self.schedule.mask(ordinal)
        b = a // 2 if self.parity else 0
        for index, prime in enumerate(primes):
            budget.consume(a.bit_length() + prime.bit_length() ** 2)
            term = self.terms[index]
            if term is None:
                quotient = a // prime
                term = (
                    quotient * pow(quotient, -1, prime) * self.roots[index] % a
                )
                self.terms[index] = term
            b += -term if mask & (1 << index) else term
        b = (b + a // 2) % a - a // 2
        return Polynomial(base.n, base.multiplier, a, b)
