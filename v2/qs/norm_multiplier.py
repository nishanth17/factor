"""Finite actual-polynomial previews for an opt-in coupled multiplier score."""

from dataclasses import dataclass
from fractions import Fraction
from math import prod

from ..common import arithmetic, utils
from .assignment_stream import AssignmentStream
from .count_policy import CountPolicyRefusedError, endpoint_for_count
from .factor_base import build_factor_base, checked_target
from .multiplier import (
    WIDE_PRIME_BOUND,
    _log_weight,
    select_coupled_multiplier,
)
from .polynomial import Polynomial
from .sign_schedule import schedule_for

SHORTLIST_SIZE = 4
PREVIEW_FAMILIES = 4
PREVIEW_POLYNOMIALS = 4
PREVIEW_POSITIONS = 17
PREVIEW_RECEIPT_BYTES = 65536


class NormPolicyRefusedError(CountPolicyRefusedError):
    """No shortlisted candidate has a valid finite preview population."""


@dataclass(frozen=True)
class NormMultiplierChoice:
    """Exact selection plus small receipts; sampled bases never escape."""

    multiplier: int
    divisor: int | None
    scores: tuple
    density_scores: tuple
    previews: tuple
    refused: tuple


def preview_positions(half_width):
    """Fixed integer positions, including both ends of the actual window."""
    return tuple(
        sorted(
            {
                -half_width + 2 * half_width * j // (PREVIEW_POSITIONS - 1)
                for j in range(PREVIEW_POSITIONS)
            }
        )
    )


def preview_polynomials(base, primes, config, budget):
    """Construct exact selected CRT polynomials without any root table.

    This computes the same centered B as the declared sign/Gray schedule.
    It avoids charging/building unused full-factor-base recurrence state.
    """
    parity = config.polynomial_policy == "parity2" and base.n_prime % 4 == 1
    a = arithmetic.backend_for(base.n).integer(prod(primes))
    a *= 2 if parity else 1
    entries = {
        entry.prime: entry for entry in base.entries if entry.prime in primes
    }
    terms = []
    for prime in primes:
        budget.consume(a.bit_length() + prime.bit_length() ** 2)
        quotient = a // prime
        gamma = (
            utils.modular_inverse(quotient, prime)
            * entries[prime].square_roots[0]
        ) % prime
        terms.append(quotient * gamma)
    schedule = schedule_for(
        len(primes),
        "swap_high"
        if config.sign_policy == "swap_high" and len(primes) >= 3
        else "control",
    )
    ordinals = sorted(
        {
            (config.gray_limit - 1) * j // (PREVIEW_POLYNOMIALS - 1)
            for j in range(PREVIEW_POLYNOMIALS)
        }
    )
    result = []
    for ordinal in ordinals:
        budget.consume((len(primes) + 1) * (a.bit_length() + 1))
        mask = schedule.mask(ordinal)
        raw = a // 2 if parity else 0
        raw += sum(
            -term if mask & (1 << j) else term for j, term in enumerate(terms)
        )
        b = (raw + a // 2) % a - a // 2
        result.append(Polynomial(base.n, base.multiplier, a, b))
    return tuple(result)


def _preview(n, multiplier, config, seed, budget, memory_bytes):
    """Own one base and population at a time, then return scalar receipts."""
    bound = config.base_bound
    if hasattr(config, "base_count"):
        bound, divisor = endpoint_for_count(
            n,
            multiplier,
            config.base_count,
            budget,
            ceiling=bound,
            memory_bytes=memory_bytes,
        )
        if divisor is not None:
            return divisor, None
    setup = build_factor_base(
        n,
        multiplier=multiplier,
        bound=bound,
        budget=budget,
        memory_bytes=memory_bytes,
    )
    if setup.divisor is not None:
        return setup.divisor, None
    base = setup.factor_base
    parity = config.polynomial_policy == "parity2" and base.n_prime % 4 == 1
    try:
        stream = AssignmentStream(
            base,
            config.half_width * (2 if parity else 1),
            factor_count=config.factor_count,
            family_count=PREVIEW_FAMILIES,
            pool_size=config.pool_size,
            seed=seed,
            policy=config.assignment_policy,
            budget=budget,
            memory_bytes=memory_bytes,
        )
    except ValueError as error:
        raise NormPolicyRefusedError(str(error)) from error
    positions = preview_positions(config.half_width)
    logs, samples, correction = 0, 0, Fraction(0)
    coefficients = []
    for index in range(len(stream)):
        primes = stream[index]
        budget.consume(len(primes) * 8192)
        correction += sum(
            (
                Fraction(_log_weight(p), p - 1)
                for p in primes
                if p < WIDE_PRIME_BOUND
            ),
            Fraction(0),
        )
        polynomials = preview_polynomials(base, primes, config, budget)
        coefficients.append(int(polynomials[0].a))
        for polynomial in polynomials:
            for position in positions:
                # Bound exact powering for log_weight before its allocation.
                bits = (
                    base.n_prime.bit_length()
                    + polynomial.a.bit_length()
                    + 2 * config.half_width.bit_length()
                    + 4
                )
                budget.consume(256 * bits + bits * bits)
                value = abs(polynomial.value(position))
                if value == 0:
                    divisor = arithmetic.gcd(polynomial.u_value(position), n)
                    if utils.valid_divisor(divisor, n):
                        return divisor, None
                logs += _log_weight(max(1, value))
                samples += 1
    return None, dict(
        multiplier=multiplier,
        bound=bound,
        base_count=len(base.entries),
        coefficients=tuple(coefficients),
        log_sum=logs,
        samples=samples,
        correction_numerator=correction.numerator,
        correction_denominator=correction.denominator * len(stream),
        assignment_identity=stream.identity,
    )


def select_norm_multiplier(n, *, config, seed, budget, memory_bytes):
    """Rerank four density nominees using fixed actual-A/norm samples.

    Scores are an exact rational heuristic in units of 1/256 binary log.
    Returned floors are display values; exact numerator/denominator receipts
    determine ranking and ties. This predicts neither yield nor setup time.
    """
    checked_target(n, 1)
    if config.mode != "siqs" or config.assignment_policy == "reference":
        raise ValueError("norm multiplier policy requires streaming SIQS")
    utils.require_integer(seed, "seed", 0)
    if seed >= 2**64:
        raise ValueError("norm selector seed exceeds 64 bits")
    scratch = 262144 + 1024 * (
        2 * n.bit_length() + 4 * config.half_width.bit_length() + 128
    )
    if scratch > memory_bytes:
        raise MemoryError("norm preview scratch exceeds memory_bytes")
    available = memory_bytes - scratch
    density = select_coupled_multiplier(
        n,
        budget=budget,
        memory_bytes=available,
    )
    if density.divisor is not None:
        return NormMultiplierChoice(
            density.multiplier, density.divisor, (), density.scores, (), ()
        )
    shortlist = sorted(density.scores, key=lambda item: (-item[1], item[0]))[
        :SHORTLIST_SIZE
    ]
    previews, refused, scored = [], [], []
    for multiplier, old_score in shortlist:
        try:
            divisor, receipt = _preview(
                n,
                multiplier,
                config,
                seed,
                budget,
                available,
            )
        except (CountPolicyRefusedError, NormPolicyRefusedError) as error:
            refused.append((multiplier, str(error)))
            continue
        if divisor is not None:
            return NormMultiplierChoice(
                multiplier,
                divisor,
                (),
                density.scores,
                tuple(previews),
                tuple(refused),
            )
        score = (
            old_score
            + Fraction(_log_weight(multiplier), 2)
            - Fraction(receipt["log_sum"], receipt["samples"])
            - Fraction(
                receipt["correction_numerator"],
                receipt["correction_denominator"],
            )
        )
        receipt["score_numerator"] = score.numerator
        receipt["score_denominator"] = score.denominator
        previews.append(receipt)
        scored.append((multiplier, score))
    if not scored:
        raise NormPolicyRefusedError(
            "all four norm multiplier previews refused: " + repr(refused)
        )
    selected = max(scored, key=lambda item: (item[1], -item[0]))[0]
    scores = tuple(
        (k, score.numerator // score.denominator) for k, score in scored
    )
    return NormMultiplierChoice(
        selected, None, scores, density.scores, tuple(previews), tuple(refused)
    )
