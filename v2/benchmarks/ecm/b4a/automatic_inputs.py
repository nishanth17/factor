"""Bounded certified, factor-disjoint inputs for automatic ECM confirmation."""

import argparse
import hashlib
import json
import math
import random
import time
from pathlib import Path

from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)

BANDS = (30, 40, 50, 60, 70, 80, 100)


def exclusions(directory):
    """Inventory earlier inputs and top-level factors without factoring."""
    numbers, factors, files = set(), set(), {}
    for path in sorted(Path(directory).glob("*.json")):
        raw = path.read_bytes()
        data = json.loads(raw)
        files[path.name] = hashlib.sha256(raw).hexdigest()
        pending = [data]
        while pending:
            value = pending.pop()
            if isinstance(value, dict):
                if type(value.get("n")) is int:
                    numbers.add(value["n"])
                for entry in value.get("factors") or ():
                    if (
                        isinstance(entry, (tuple, list))
                        and len(entry) == 2
                        and type(entry[0]) is int
                    ):
                        factors.add(entry[0])
                pending.extend(value.values())
            elif isinstance(value, list):
                pending.extend(value)
    return numbers, factors, files


def build(directory, *, seed=202611130412, seconds=45, draws=200000):
    """Generate every planned look before exposing any factoring outcome."""
    if not 0 < seconds <= 90 or not 1000 <= draws <= 200000:
        raise ValueError("finite generation bounds required")
    started = time.monotonic()
    old_numbers, old_factors, files = exclusions(directory)

    class LimitedRandom(random.Random):
        used = 0

        def randrange(self, *args):
            self.used += 1
            if self.used > draws or time.monotonic() - started > seconds:
                raise TimeoutError("automatic input generation exhausted")
            return super().randrange(*args)

    generator = LimitedRandom(seed)
    proofs, seen, primes = {}, set(old_numbers), set(old_factors)

    def prime(digits):
        bits = (10**digits - 1).bit_length()
        while True:
            if len(proofs) > 20000:
                raise ValueError("certificate node allowance exhausted")
            value = certified_prime(bits, generator, proofs)
            if len(str(value)) == digits and value not in primes:
                primes.add(value)
                return value

    def semiprime(digits, small, identity, population):
        while True:
            p, q = prime(small), prime(digits - small)
            n = p * q
            if len(str(n)) == digits and n not in seen:
                seen.add(n)
                return dict(
                    id=identity,
                    n=n,
                    digits=digits,
                    smaller_factor_digits=small,
                    population=population,
                    factors=sorted([[p, 1], [q, 1]]),
                )

    fixtures = []
    for digits in BANDS:
        smalls = (13, 14, 15) if digits == 30 else (15, 18, digits // 2)
        for index in range(9):
            small = smalls[index % 3]
            fixtures.append(
                semiprime(
                    digits,
                    small,
                    f"auto-{digits}-{index}",
                    "balanced-quota" if index % 3 == 2 else "ecm-relevant",
                )
            )
    integration = [
        semiprime(d, d // 2, f"auto-integration-{d}", "balanced-integration")
        for d in (60, 70)
    ]
    controls = []
    for digits in (30, 60, 80, 100):
        value = prime(digits)
        controls.append(
            dict(
                id=f"auto-prime-{digits}",
                n=value,
                population="prime-control",
                factors=[[value, 1]],
            )
        )
    uniform_generator = random.Random(seed + 1)
    uniform = [
        dict(
            id=f"auto-uniform-{d}-{i}",
            digits=d,
            n=uniform_generator.randrange(10 ** (d - 1), 10**d),
            population="uniform-unfiltered",
        )
        for d in BANDS
        for i in range(6)
    ]
    values = [r["n"] for r in fixtures + integration + controls + uniform]
    if len(set(values)) != len(values) or set(values) & old_numbers:
        raise ValueError("new input overlap; no outcome-dependent resampling")
    verify_certificates(proofs)
    for row in fixtures + integration + controls:
        if math.prod(p**e for p, e in row["factors"]) != row["n"]:
            raise ArithmeticError(
                "generated factorization did not reconstruct"
            )
    result = dict(
        schema="b4a-automatic-confirmation-v1",
        generation_seed=seed,
        fixtures=fixtures,
        integration=integration,
        controls=controls,
        uniform=uniform,
        certificates=proofs,
        generation_wall=time.monotonic() - started,
        random_draws=generator.used,
        excluded_inputs=len(old_numbers),
        excluded_top_factors=len(old_factors),
        exclusion_files=files,
        distribution=(
            "Constructive nonuniform Pocklington primes, conditioned "
            "on exact decimal factor/product sizes and disjointness. "
            "Known factors and strata never enter execution selection. "
            "Nine inputs per band form three predeclared looks; "
            "large balanced inputs target quota exhaustion without "
            "assuming future outcomes. Uniform controls are unsifted. "
            "This is not a random-integer population/p90 corpus."
        ),
    )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--exclusions", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=202611130412)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("refusing to overwrite generated evidence")
    result = build(args.exclusions, seed=args.seed)
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(
        json.dumps(
            dict(
                fixtures=len(result["fixtures"]),
                integration=len(result["integration"]),
                draws=result["random_draws"],
                wall=result["generation_wall"],
            )
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
