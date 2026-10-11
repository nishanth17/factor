"""Prepare disjoint confirmation inputs after a separate recipe freeze."""

import math
from collections import Counter

from ...suites.build_phase_two_corpus import certified_prime
from .round2_prime_proofs import MAX_NODES, UniformPrimeSource, verify_proofs


def strata(balanced_count, uneven_count):
    """Specify exact decimal bands independently of generated outcomes."""
    if any(
        type(c) is not int or not 1 <= c <= 8
        for c in (balanced_count, uneven_count)
    ):
        raise ValueError("confirmation stratum counts must be 1..8")
    for digits in (30, 40):
        yield digits, digits // 2, "balanced", balanced_count
        for small in (8, 10, 12, 14 if digits == 30 else 16):
            yield digits, small, f"uneven_{small}", uneven_count


def validate(corpus, excluded):
    """Verify primes, strata, products and top-level factor disjointness."""
    verify_proofs(corpus["certificates"])
    expected = {
        (digits, kind): (small, count)
        for digits, small, kind, count in strata(**corpus["counts"])
    }
    old_numbers = {f["n"] for f in excluded}
    factors_used = {p for f in excluded for p, _e in f["factors"]}
    identifiers, numbers, counts = set(), set(), Counter()
    for fixture in corpus["fixtures"]:
        identity = (fixture["digits"], fixture["kind"])
        if identity not in expected:
            raise ValueError("unexpected confirmation stratum")
        small, _count = expected[identity]
        factors = fixture["factors"]
        if (
            len(factors) != 2
            or any(
                type(p) is not int
                or type(e) is not int
                or e != 1
                or str(p) not in corpus["certificates"]
                for p, e in factors
            )
            or factors[0][0] == factors[1][0]
            or math.prod(p**e for p, e in factors) != fixture["n"]
            or len(str(fixture["n"])) != fixture["digits"]
            or len(str(min(p for p, _e in factors))) != small
            or fixture["small_digits"] != small
        ):
            raise ValueError("invalid certified confirmation factorization")
        if (
            fixture["id"] in identifiers
            or fixture["n"] in numbers | old_numbers
            or any(p in factors_used for p, _e in factors)
        ):
            raise ValueError(
                "confirmation overlaps revealed inputs or factors"
            )
        identifiers.add(fixture["id"])
        numbers.add(fixture["n"])
        factors_used.update(p for p, _e in factors)
        counts[identity] += 1
    if counts != Counter(
        {key: count for key, (_s, count) in expected.items()}
    ):
        raise ValueError("incomplete confirmation strata")


def generate(seed, excluded, *, balanced_count=2, uneven_count=2):
    """Generate only after committing a candidate and a new finite recipe.

    Uniform candidate sampling covers every factor up to twenty digits.
    Larger partners retain the disclosed Pocklington bias. Failed proof work
    aborts the entire attempt; only declared decimal/disjointness rejection
    may draw a new pair. This helper does not select or freeze a policy.
    """
    recipe = tuple(strata(balanced_count, uneven_count))
    source = UniformPrimeSource(seed)
    certificates = source.certificates
    old_numbers = {f["n"] for f in excluded}
    factors_used = {p for f in excluded for p, _e in f["factors"]}
    fixtures, numbers = [], set()

    def prime(digits):
        if digits <= 20:
            return source.prime(digits)
        lower, upper = 10 ** (digits - 1), 10**digits - 1
        for _ in range(1000):
            source.budget.consume(0)
            # At most 32 decimal digits are requested. Sixteen reserved nodes
            # cover that recursive Pocklington chain without overshooting the
            # certificate allowance inside the older generator.
            if len(certificates) + 16 > MAX_NODES:
                raise RuntimeError("confirmation proof storage exhausted")
            bits = source.generator.randint(
                lower.bit_length(), upper.bit_length()
            )
            value = certified_prime(bits, source.generator, certificates)
            source.budget.consume(0)
            if lower <= value <= upper:
                return value
        raise RuntimeError("larger-prime decimal rejection exhausted")

    for digits, small, kind, count in recipe:
        for index in range(count):
            for _ in range(1000):
                p, q = prime(small), prime(digits - small)
                n = p * q
                if (
                    p != q
                    and len(str(n)) == digits
                    and n not in old_numbers | numbers
                    and p not in factors_used
                    and q not in factors_used
                ):
                    break
            else:
                raise RuntimeError("confirmation pair rejection exhausted")
            fixtures.append(
                dict(
                    id=f"c3r2fresh_{digits}_{kind}_{index}",
                    kind=kind,
                    n=n,
                    factors=sorted([(p, 1), (q, 1)]),
                    digits=digits,
                    small_digits=small,
                )
            )
            numbers.add(n)
            factors_used.update((p, q))

    counts = dict(balanced_count=balanced_count, uneven_count=uneven_count)
    corpus = dict(
        seed=seed,
        fixtures=fixtures,
        certificates=certificates,
        counts=counts,
        bias=(
            "Uniform odd candidates for factors <=20 digits; larger partners "
            "use nonuniform Pocklington p=kq+1; joint decimal and factor-"
            "disjointness rejection conditions the population"
        ),
        generation_work=source.budget.used,
        generation_draws=1_000_000 - source.generator.remaining,
    )
    source.budget.consume(0)
    validate(corpus, excluded)
    source.budget.consume(0)
    corpus["generation_wall"] = source.budget.wall_used
    corpus["generation_cpu"] = source.budget.cpu_used
    return corpus
