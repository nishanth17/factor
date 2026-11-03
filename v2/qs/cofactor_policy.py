"""Finite exact cofactor strategies and synchronous bounded batch hints.

SQUFOF follows the integer recurrence in Gower--Wagstaff, section 3.3.
Method order and limits are experimental controls, not portable optima.
Every hint is checked as a proper divisor before endpoint certification.
"""

from collections import deque
from math import gcd, isqrt

from ..common import utils
from ..rho.brent import RhoStats, factorize_rho
from .large_primes import split_two_primes


def squfof(number, budget, *, steps=1024, queue_limit=64):
    """Bound both directions and the small-form queue; return a proper split."""
    utils.require_integer(number, "cofactor", 4)
    utils.require_integer(steps, "SQUFOF steps", 1)
    utils.require_integer(queue_limit, "SQUFOF queue", 1)
    if number >= 2**80 or steps > 65536 or queue_limit > 256:
        raise ValueError("SQUFOF domain or effort exceeds its cap")
    if number % 2 == 0:
        return 2
    root = isqrt(number)
    if root * root == number:
        return root
    discriminant = 2 * number if number % 4 == 1 else number
    if discriminant != number:
        root = isqrt(discriminant)
    p, old_q, q = root, 1, discriminant - root * root
    if q == 0:
        divisor = gcd(root, number)
        return divisor if utils.valid_divisor(divisor, number) else None
    small_limit = isqrt(isqrt(64 * discriminant))
    queue = deque()
    cost = 8 * discriminant.bit_length() ** 2
    for iteration in range(steps):
        budget.consume(cost)
        quotient = (root + p) // q
        next_p = quotient * q - p
        if q <= small_limit:
            small = q // 2 if q % 2 == 0 else q
            if q % 2 == 0 or q <= small_limit // 2:
                if len(queue) == queue_limit:
                    return None
                queue.append((small, p % small))
        old_q, q, p = q, old_q + quotient * (p - next_p), next_p
        if iteration % 2:
            continue
        square = isqrt(q)
        if square * square != q:
            continue
        budget.consume(len(queue))
        match = next(
            (
                i
                for i, (a, b) in enumerate(queue)
                if a == square and (p - b) % square == 0
            ),
            None,
        )
        if match is not None:
            if square == 1:
                return None
            for _ in range(match + 1):
                queue.popleft()
            continue

        old_q = square
        p += square * ((root - p) // square)
        q, remainder = divmod(discriminant - p * p, old_q)
        if remainder or q <= 0:
            raise ArithmeticError("SQUFOF inverse-square invariant failed")
        for _ in range(steps):
            budget.consume(cost)
            quotient = (root + p) // q
            next_p = quotient * q - p
            if p == next_p:
                divisor = gcd(q // 2 if q % 2 == 0 else q, number)
                return (
                    divisor if utils.valid_divisor(divisor, number) else None
                )
            old_q, q, p = q, old_q + quotient * (p - next_p), next_p
            if q <= 0:
                raise ArithmeticError("SQUFOF reduction invariant failed")
        return None
    return None


def batch_divisors(residuals, budget, *, maximum=256):
    """Shared-factor hints using P mod R², with exact whole-GCD recovery."""
    if not isinstance(residuals, tuple) or len(residuals) > maximum:
        raise ValueError("cofactor batch exceeds its count cap")
    if any(type(r) is not int or not 2 <= r < 2**80 for r in residuals):
        raise ValueError("cofactor batch exceeds its integer domain")
    product = 1
    for value in residuals:
        budget.consume((product.bit_length() + 1) * value.bit_length())
        product *= value
    result = []
    for value in residuals:
        budget.consume(
            product.bit_length() * value.bit_length() + value.bit_length() ** 2
        )
        quotient = (product % (value * value)) // value
        divisor = gcd(value, quotient)
        if divisor == value:
            # Repeated complete residuals and endpoints shared on both sides
            # can saturate the batch GCD. The bounded descent must not label
            # the whole residual a successful split or a prime.
            divisor = None
            for other in residuals:
                budget.consume(value.bit_length() * other.bit_length())
                trial = gcd(value, other)
                if utils.valid_divisor(trial, value):
                    divisor = trial
                    break
        result.append(divisor if utils.valid_divisor(divisor, value) else None)
    return tuple(result)


class KnownEndpoints:
    """A bounded product of certified endpoints, independent of residual repeats."""

    def __init__(self, maximum=256):
        self.maximum = maximum
        self.primes, self.seen, self.product = [], set(), 1

    def add(self, primes, budget):
        for prime in primes:
            if prime not in self.seen and len(self.primes) < self.maximum:
                budget.consume(
                    (self.product.bit_length() + 1) * prime.bit_length()
                )
                self.product *= prime
                self.primes.append(prime)
                self.seen.add(prime)

    def divisor(self, residual, budget):
        if not self.primes or residual in self.seen:
            # A previously certified prime cannot yield a proper split. In
            # particular, avoid a saturated product GCD followed by a scan.
            budget.consume(1)
            return None
        budget.consume(self.product.bit_length() * residual.bit_length())
        divisor = gcd(self.product, residual)
        if utils.valid_divisor(divisor, residual):
            return divisor
        if divisor == residual:
            for prime in self.primes:
                budget.consume(residual.bit_length())
                if residual % prime == 0 and prime < residual:
                    return prime
        return None


def split_cofactor(
    residual,
    bound,
    budget,
    checked,
    policy,
    *,
    hint=None,
    square_checked=False,
):
    """One resolved finite schedule, with endpoint proofs shared privately."""
    if policy == "reference":
        return split_two_primes(residual, bound, budget, known_composite=True)
    bits, evaluations = residual.bit_length(), 0
    budget.consume(bits**2)
    divisor = hint
    if divisor is None and not square_checked:
        root = isqrt(residual)
        if root * root == residual:
            divisor = root
    if divisor is not None and not utils.valid_divisor(divisor, residual):
        raise ValueError("invalid cofactor hint")
    if divisor is None and policy in ("portfolio-v2", "ecm-v2"):
        divisor = squfof(residual, budget, steps=1024)
    if divisor is None:
        # Two short walks preserve a deterministic inexpensive first try;
        # the second, larger attempt is paid only after the first fails.
        for seed, limit, batch in ((7, 512, 16), (11, 2048, 32)):
            budget.consume(limit * bits**2)
            stats = RhoStats()
            divisor = factorize_rho(
                residual,
                seed=seed,
                max_attempts=1,
                max_evaluations=limit,
                batch_size=batch,
                recovery_limit=min(batch, 64),
                stats=stats,
                _known_composite=True,
            )
            evaluations += stats.evaluations
            budget.consume(0)
            if divisor is not None:
                break
    if divisor is None and policy == "ecm-v2":
        from ..ecm.core import factorize_ecm

        # Existing arithmetic, one finite curve, and an explicit conservative
        # reservation. This fallback is selectable; it is not a new default
        # stage before every small residual.
        budget.consume(100000 * bits**2)
        divisor = factorize_ecm(
            residual,
            b1=50,
            b2=200,
            seed=17,
            max_curves=1,
            _known_composite=True,
        )
        budget.consume(0)
    if divisor is None:
        return (), "split_failed", evaluations
    if not utils.valid_divisor(divisor, residual):
        raise ArithmeticError("cofactor method returned an invalid divisor")
    pair = tuple(sorted((int(divisor), int(residual // divisor))))
    if pair[-1] > bound:
        return (), "endpoint_bound", evaluations
    if any(not checked.prime(p, budget) for p in pair):
        return (), "not_two_primes", evaluations
    if pair[0] * pair[1] != residual:
        raise ArithmeticError("cofactor reconstruction failed")
    return pair, "square" if pair[0] == pair[1] else "double", evaluations


class CofactorBatch:
    """Disposable lookahead within one sieve block; no persistent pending work.

    Exact division happens once for at most batch_size score candidates.
    Publication remains in physical position order. An interruption discards
    only derived scratch; already charged recovery is not refunded on resume.
    """

    def __init__(self, size):
        self.size, self.records, self.hints = size, {}, {}

    def recover(self, collector, position, offset, block_hi, threshold, stats):
        if position not in self.records:
            self.records.clear()
            self.hints.clear()
            original = collector._divide_double

            def capture(worker, x, value, remaining, exponents, counters):
                return (x, value, remaining, tuple(exponents)), None

            try:
                collector._divide_double = capture
                block_lo = position - offset
                for x in range(position, block_hi):
                    collector.budget.consume(1)
                    if collector._scores[x - block_lo] < threshold:
                        continue
                    self.records[x] = collector._divide(x, x - block_lo, stats)
                    if len(self.records) == self.size:
                        break
            finally:
                collector._divide_double = original
            eligible = [
                (x, result[0][2])
                for x, result in self.records.items()
                if result[0] is not None
                and result[0][2] > 1
                and result[0][2] <= collector.config.large_product_bound
            ]
            if len(eligible) > 1:
                divisors = batch_divisors(
                    tuple(int(r) for _, r in eligible), collector.budget
                )
                self.hints = dict(zip((x for x, _ in eligible), divisors))
                stats["cofactor_batches"] = (
                    stats.get("cofactor_batches", 0) + 1
                )
        record, divisor = self.records.pop(position)
        if record is None:
            return None, divisor
        return collector._divide_double(
            collector, *record, stats, hint=self.hints.pop(position, None)
        )
