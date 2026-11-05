"""Finite fresh inputs for cost calibration and policy confirmation."""

import argparse
import fcntl
import hashlib
import json
import math
import os
import random
import resource
import signal
import sys
import time
from pathlib import Path

CALIBRATION_DIGITS = (55, 75, 81, 90)


def enriched_plan():
    """Separate direct semiprimes from preprocessed cofactor coverage cases."""
    rows = []
    for digits, child, exponent in ((60, 55, 16), (80, 75, 16), (100, 90, 33)):
        for index, small in enumerate((15, 18, 17, 20)):
            direct = index < 2
            rows.append(
                dict(
                    id=f"ecm-{digits}-{index}",
                    digits=digits,
                    cofactor_digits=digits if direct else child,
                    smaller_cofactor_digits=small,
                    two_exponent=0 if direct else exponent,
                    population="ecm_direct" if direct else "ecm_recursive",
                    algorithm_seed=7919 if index in (0, 3) else 104729,
                )
            )
    return rows


def build(phase, seed, *, excluded=(), seconds=90, draws=200000):
    """Prove candidate primality with independent recursive certificates."""
    from ...suites.build_phase_two_corpus import (
        certified_prime,
        verify_certificates,
    )

    if phase not in ("calibration", "confirmation"):
        raise ValueError("unknown input phase")
    if type(seed) is not int or seed < 0:
        raise ValueError("nonnegative generation seed required")
    if not 0 < seconds <= 90 or not 1000 <= draws <= 200000:
        raise ValueError("finite generation allowance required")
    started = time.monotonic()

    class LimitedRandom(random.Random):
        used = 0

        def randrange(self, *args):
            self.used += 1
            if self.used > draws or time.monotonic() - started > seconds:
                raise TimeoutError("targeted corpus generation exhausted")
            return super().randrange(*args)

    generator = LimitedRandom(seed)
    proofs, seen, fixtures = {}, set(), []

    def prime(digits):
        bits = (10**digits - 1).bit_length()
        while True:
            value = certified_prime(bits, generator, proofs)
            if len(str(value)) == digits and value not in seen:
                seen.add(value)
                return value

    if phase == "calibration":
        for digits in CALIBRATION_DIGITS:
            value = prime(digits)
            fixtures.append(
                dict(
                    id=f"cost-{digits}",
                    n=value,
                    digits=digits,
                    population="certified_prime_fixed_work",
                    factors=[[value, 1]],
                    algorithm_seed=7919,
                )
            )
        uniform, controls = [], []
    else:
        for plan in enriched_plan():
            while True:
                small = plan["smaller_cofactor_digits"]
                p, q = prime(small), prime(plan["cofactor_digits"] - small)
                child = p * q
                n = child * 2 ** plan["two_exponent"]
                if (
                    len(str(child)) == plan["cofactor_digits"]
                    and len(str(n)) == plan["digits"]
                ):
                    break
            factors = [[p, 1], [q, 1]]
            if plan["two_exponent"]:
                factors.append([2, plan["two_exponent"]])
                proofs["2"] = dict(kind="trial")
            fixtures.append(dict(plan, n=n, factors=sorted(factors)))
        # No factor/parity sifting or outcome-dependent resampling.
        uniform = []
        uniform_rng = random.Random(seed + 1)
        for digits in (30, 40, 50, 60, 70, 80, 100):
            count = 12 if digits in (60, 80, 100) else 1
            for index in range(count):
                uniform.append(
                    dict(
                        id=f"uniform-{digits}-{index:02}",
                        n=uniform_rng.randrange(
                            10 ** (digits - 1), 10**digits
                        ),
                        digits=digits,
                        population="uniform",
                        algorithm_seed=7919 if index % 2 == 0 else 104729,
                    )
                )
        controls = []
        for digits in (30, 40, 50, 70):
            value = prime(digits)
            controls.extend(
                (
                    dict(
                        id=f"prime-{digits}",
                        n=value,
                        digits=digits,
                        factors=[[value, 1]],
                        population="prime_control",
                        algorithm_seed=7919,
                    ),
                    dict(
                        id=f"power-{digits}",
                        n=10 ** (digits - 1),
                        digits=digits,
                        factors=[[2, digits - 1], [5, digits - 1]],
                        population="exposed_power_control",
                        algorithm_seed=7919,
                    ),
                )
            )
        proofs["5"] = dict(kind="trial")
    verify_certificates(proofs)
    fresh = (
        fixtures
        + uniform
        + [r for r in controls if r["population"] == "prime_control"]
    )
    values = [row["n"] for row in fresh]
    if len(set(values)) != len(values) or set(values).intersection(excluded):
        raise ValueError("fresh corpus overlaps earlier exposure or itself")
    for row in fixtures + controls:
        if math.prod(p**e for p, e in row["factors"]) != row["n"]:
            raise ValueError("input factorization does not reconstruct")
    return dict(
        schema="b12-f4-targeted-corpus-v1",
        phase=phase,
        generation_seed=seed,
        fixtures=fixtures,
        uniform=uniform,
        controls=controls,
        certificates=proofs,
        random_draws=generator.used,
        generation_wall=time.monotonic() - started,
        distribution=(
            "Pocklington-conditioned enriched inputs; "
            "separately unconditioned uniform integers"
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase", choices=("calibration", "confirmation"), required=True
    )
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--policy-freeze", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    policy = json.loads(args.policy_freeze.read_text())
    expected = (
        "frozen-for-calibration"
        if args.phase == "calibration"
        else "frozen-for-confirmation"
    )
    if policy["status"] != expected or args.output.exists():
        raise ValueError(
            "input generation requires a new output and frozen policy"
        )
    descriptor = int(os.environ["FACTOR_DIAGNOSTIC_LOCK_FD"])
    actual, expected_lock = (
        os.fstat(descriptor),
        os.stat("/private/tmp/factor-performance.lock"),
    )
    if (actual.st_dev, actual.st_ino) != (
        expected_lock.st_dev,
        expected_lock.st_ino,
    ):
        raise RuntimeError("shared lock identity changed")
    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_CPU, (95, 96))
    signal.alarm(98)
    data = build(args.phase, args.seed, excluded=policy["exclusion_moduli"])
    data["policy_sha256"] = hashlib.sha256(
        args.policy_freeze.read_bytes()
    ).hexdigest()
    with args.output.open("x") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
    print(
        json.dumps(
            dict(
                fixtures=len(data["fixtures"]),
                uniform=len(data["uniform"]),
                controls=len(data["controls"]),
                draws=data["random_draws"],
                sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),
            )
        )
    )


if __name__ == "__main__":
    main()
