"""Finite exact cofactor strategies and synchronous bounded batch hints.

SQUFOF follows the integer recurrence in Gower--Wagstaff, section 3.3.
Method order and limits are experimental controls, not portable optima.
Every hint is checked as a proper divisor before endpoint certification.
"""

from collections import OrderedDict, deque
from math import gcd, isqrt

from ..common import utils
from ..rho.brent import RhoStats, factorize_rho
from .large_primes import split_two_primes


def squfof(number, budget, *, steps=1024, queue_limit=64):
    """Bound both directions and the queue; return a proper split."""
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
    if any(type(r) is not int or not 2 <= r < 2**128 for r in residuals):
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
    """A bounded product of certified endpoints, independent of repeats."""

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
    dispatch=False,
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
    if divisor is None and dispatch and residual * 4 >= bound * bound:
        # The endpoint envelope forces balanced factors here. A tiny Fermat
        # window is exact; exhausting it says nothing about compositeness.
        middle = isqrt(residual) + 1
        for _ in range(32):
            budget.consume(4 * bits**2)
            delta = middle * middle - residual
            root = isqrt(delta)
            if root * root == delta and utils.valid_divisor(
                middle - root, residual
            ):
                divisor = middle - root
                break
            middle += 1
    if (
        divisor is None
        and policy in ("portfolio-v2", "ecm-v2")
        and (not dispatch or 32 < bits < 80)
    ):
        divisor = squfof(residual, budget, steps=512 if dispatch else 1024)
    if divisor is None:
        # Two short walks preserve a deterministic inexpensive first try;
        # the second, larger attempt is paid only after the first fails.
        schedule = (
            ((7, 128, 8), (11, 512, 16))
            if dispatch and bits <= 32
            else ((7, 512, 16), (11, 2048, 32))
        )
        for seed, limit, batch in schedule:
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


class SplitMemo:
    """Finite exact residual/policy cache, including paid failed attempts."""

    def __init__(self, capacity):
        if type(capacity) is not int or not 0 <= capacity <= 4096:
            raise ValueError("invalid split memo capacity")
        self.capacity, self.values = capacity, OrderedDict()

    def get(self, key, budget):
        budget.consume(1)
        result = self.values.get(key)
        if result is not None:
            self.values.move_to_end(key)
        return result

    def remember(self, key, value):
        if not self.capacity:
            return
        self.values[key] = value
        self.values.move_to_end(key)
        if len(self.values) > self.capacity:
            self.values.popitem(last=False)


def split_residual(
    residual,
    bound,
    budget,
    checked,
    policy,
    *,
    arity=3,
    hint=None,
    dispatch=True,
    memo=None,
):
    """Finite at-most-three-prime decomposition with exact endpoint proof.

    A failed two-prime classification can still contain one composite factor.
    A separate bounded walk decomposes it; no PRP result certifies an endpoint.
    """
    if arity not in (2, 3) or not 2 <= residual <= bound**arity:
        return (), "endpoint_bound", 0
    if hint is not None and not utils.valid_divisor(hint, residual):
        raise ValueError("invalid residual split hint")
    key = residual, bound, policy, arity, dispatch
    previous = memo.get(key, budget) if memo is not None else None
    if previous is not None and (previous[0] or hint is None):
        factors, reason = previous
        if factors and (
            _product(factors) != residual
            or any(p > bound or not checked.prime(p, budget) for p in factors)
        ):
            raise ValueError("cached cofactor proof changed")
        return factors, reason, 0

    def finished(factors, reason, evaluations):
        if memo is not None:
            memo.remember(key, (factors, reason))
        return factors, reason, evaluations

    pair, reason, evaluations = (), "not_two_primes", 0
    if residual <= bound**2:
        pair, reason, evaluations = split_cofactor(
            residual,
            bound,
            budget,
            checked,
            policy,
            hint=hint,
            dispatch=dispatch,
        )
    if pair or arity == 2:
        return finished(pair, reason, evaluations)
    pending, factors = [residual], []
    for attempt in range(2 * arity - 1):
        if not pending:
            break
        value = pending.pop()
        if checked.prime(value, budget):
            if value > bound:
                return finished((), "endpoint_bound", evaluations)
            factors.append(value)
            continue
        if len(factors) + len(pending) + 2 > arity:
            return finished((), "not_two_primes", evaluations)
        root = isqrt(value)
        divisor = hint if attempt == 0 else None
        if divisor is None and root * root == value:
            divisor = root
        if divisor is None:
            budget.consume(2048 * value.bit_length() ** 2)
            stats = RhoStats()
            divisor = factorize_rho(
                value,
                seed=23 + attempt,
                max_attempts=1,
                max_evaluations=2048,
                batch_size=32,
                recovery_limit=32,
                stats=stats,
                _known_composite=True,
            )
            evaluations += stats.evaluations
            budget.consume(0)
        if not utils.valid_divisor(divisor, value):
            return finished((), "split_failed", evaluations)
        pending.extend((int(divisor), int(value // divisor)))
    if pending or _product(factors) != residual:
        return finished((), "split_failed", evaluations)
    result = tuple(sorted(factors))
    return finished(result, "double", evaluations)


def _product(values):
    product = 1
    for value in values:
        product *= value
    return product


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
                return (x, value, remaining, exponents), None

            try:
                collector._divide_double = capture
                block_lo = position - offset
                positions = []
                for x in range(position, block_hi):
                    collector.budget.consume(1)
                    if collector._scores[x - block_lo] < threshold:
                        continue
                    if not collector._piecewise_passes(x - block_lo):
                        continue
                    positions.append(x)
                    if len(positions) == self.size:
                        break
                dense = getattr(collector.config, "smooth_batch_size", 0)
                collector._smooth_hints = {}
                if (
                    dense
                    and len(positions) >= dense
                    and collector.config.division != "resieve"
                ):
                    from .recovery_batch import BaseRecovery

                    values = tuple(
                        abs(collector.polynomial.value(x)) for x in positions
                    )
                    if all(values):
                        helper = getattr(collector, "_base_recovery", None)
                        if helper is None:
                            try:
                                helper = BaseRecovery(
                                    collector.factor_base,
                                    collector._batch_reserve // 2,
                                    collector.budget,
                                )
                            except MemoryError:
                                helper = None
                            collector._base_recovery = helper
                        if helper is not None:
                            residuals = helper.residuals(
                                values, collector.budget
                            )
                            collector._smooth_hints = dict(
                                zip(positions, residuals)
                            )
                for x in positions:
                    self.records[x] = collector._divide(x, x - block_lo, stats)
            finally:
                collector._divide_double = original
                collector._smooth_hints = {}
            eligible = [
                (x, result[0][2])
                for x, result in self.records.items()
                if result[0] is not None
                and result[0][2] > 1
                and result[0][2] <= collector.config.large_product_bound
            ]
            if len(eligible) > 1:
                residuals = tuple(int(r) for _, r in eligible)
                if getattr(collector.config, "recovery_workers", 0):
                    from .recovery_batch import PendingHints

                    pending = PendingHints(residuals, collector.budget)
                    # Safe independent parent work can overlap the hint task.
                    try:
                        for value in residuals[:16]:
                            collector._checked.prime(value, collector.budget)
                    finally:
                        divisors = pending.result()
                else:
                    divisors = batch_divisors(residuals, collector.budget)
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
