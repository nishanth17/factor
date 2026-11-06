"""Bounded disjoint confirmation inputs, generated only after policy freeze."""

import argparse
import fcntl
import hashlib
import json
import math
import random
import resource
import signal
import time
from pathlib import Path

DIGITS = (30, 40, 50, 60, 70, 80, 100)


def smaller_sizes(digits):
    """Prespecified feasible factor strata; labels never enter dispatch."""
    plans = {
        30: (10, 12, 14, 15),
        40: (10, 15, 19, 20),
        50: (10, 15, 20, 24),
        60: (10, 18, 22, 27),
        70: (10, 18, 24, 30),
        80: (10, 18, 24, 32),
        100: (10, 18, 24, 33),
    }
    if digits not in plans:
        raise ValueError("unfrozen size band")
    return plans[digits]


def confirmation_plan(include_ecm=False):
    """Declare twelve additional ECM moduli before any generation."""
    if type(include_ecm) is not bool:
        raise ValueError("ECM confirmation selection must be Boolean")
    plans = [
        (f"fresh-{digits}-small{small}", digits, small, "certified_enriched")
        for digits in DIGITS
        for small in smaller_sizes(digits)
    ]
    if include_ecm:
        plans.extend(
            (
                f"fresh-ecm-{digits}-small{small}",
                digits,
                small,
                "selected_ecm_enriched",
            )
            for digits in (60, 80, 100)
            for small in (15, 17, 18, 20)
        )
    return plans


def uniform_panel(seed, primary_count=36, auxiliary_count=12):
    """Draw without conditioning on factors, parity or primality."""
    if not 12 <= primary_count <= 100 or not 4 <= auxiliary_count <= 100:
        raise ValueError("finite independent-input counts required")
    generator = random.Random(seed)
    return [
        dict(
            id=f"uniform-{digits}-{index:03}",
            n=generator.randrange(10 ** (digits - 1), 10**digits),
            digits=digits,
            population="uniform",
        )
        for digits in DIGITS
        for index in range(
            primary_count if digits in (60, 80, 100) else auxiliary_count
        )
    ]


def build(seed, *, seconds=120, draws=200_000, include_ecm=False):
    """Generate independent certified pairs within a finite draw/CPU lane."""
    from ...suites.build_phase_two_corpus import (
        certified_prime,
        verify_certificates,
    )

    if not 0 < seconds <= 120 or not 1000 <= draws <= 200_000:
        raise ValueError("corpus allowance exceeds the campaign allocation")
    plans = confirmation_plan(include_ecm)
    started = time.monotonic()

    class LimitedRandom(random.Random):
        used = 0

        def randrange(self, *args):
            self.used += 1
            if self.used > draws or time.monotonic() - started > seconds:
                raise TimeoutError("finite confirmation generation exhausted")
            return super().randrange(*args)

    generator = LimitedRandom(seed)
    certificates, fixtures, used_primes, moduli = {}, [], set(), set()

    def prime(digits):
        bits = (10**digits - 1).bit_length()
        while True:
            value = certified_prime(bits, generator, certificates)
            if len(str(value)) == digits and value not in used_primes:
                used_primes.add(value)
                return value

    for identity, digits, small, population in plans:
        while True:
            p = prime(small)
            q = prime(digits - small)
            n = p * q
            if len(str(n)) == digits and n not in moduli:
                break
        moduli.add(n)
        fixtures.append(
            dict(
                id=identity,
                n=n,
                digits=digits,
                factors=[[p, 1], [q, 1]],
                smaller_digits=small,
                population=population,
            )
        )
    verify_certificates(certificates)
    assert all(
        math.prod(p**e for p, e in row["factors"]) == row["n"]
        for row in fixtures
    )
    return dict(
        schema=1,
        generation_seed=seed,
        fixtures=fixtures,
        certificates=certificates,
        random_draws=generator.used,
        generation_wall=time.monotonic() - started,
        distribution=(
            "Pocklington-conditioned enriched primes; not uniform integers"
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy-freeze", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=2026110901)
    arguments = parser.parse_args()
    policy = json.loads(arguments.policy_freeze.read_text())
    if policy.get("status") != "frozen-for-confirmation":
        raise ValueError("fresh generation requires a frozen selected policy")
    if arguments.output.exists():
        raise ValueError("confirmation inputs are immutable")
    with open("/private/tmp/factor-performance.lock", "r+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        resource.setrlimit(resource.RLIMIT_CPU, (120, 121))
        signal.alarm(120)
        corpus = build(
            arguments.seed,
            include_ecm=policy.get("include_ecm_confirmation", False),
        )
        corpus["uniform"] = uniform_panel(
            arguments.seed + 1,
            policy.get("primary_uniform_count", 36),
            policy.get("auxiliary_uniform_count", 12),
        )
        excluded = set(policy["exclusion_moduli"])
        values = [row["n"] for row in corpus["fixtures"] + corpus["uniform"]]
        if len(values) != len(set(values)) or excluded.intersection(values):
            raise ValueError("fresh panel overlaps prior exposure or itself")
        corpus["policy_freeze_sha256"] = hashlib.sha256(
            arguments.policy_freeze.read_bytes()
        ).hexdigest()
        with arguments.output.open("x") as stream:
            json.dump(corpus, stream, indent=2, allow_nan=False)
        print(
            json.dumps(
                dict(
                    enriched=len(corpus["fixtures"]),
                    uniform=len(corpus["uniform"]),
                    draws=corpus["random_draws"],
                )
            )
        )


if __name__ == "__main__":
    main()
