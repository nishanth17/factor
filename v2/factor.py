#!/usr/bin/env pypy3
"""Complete or partial Python 3 factorization with explicit certainty."""

import argparse
import json
import sys
from dataclasses import dataclass, replace
from functools import lru_cache
from math import prod
from pathlib import Path
from typing import Tuple

if __package__ in (None, ""):
    # Support both `python -m v2.factor` and the existing direct script entry.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    __package__ = "v2"

from . import constants, ecm, pollard_rho, prime_sieve, utils


@dataclass(frozen=True)
class PrimeFactor:
    """A terminal factor and multiplicity with its primality evidence."""

    value: int
    exponent: int
    certainty: utils.Primality


@dataclass(frozen=True)
class FactorizationResult:
    """No cofactor may disappear, including on an algorithm's failure."""

    original: int
    sign: int
    factors: Tuple[PrimeFactor, ...]
    remaining: Tuple[int, ...]

    def __post_init__(self):
        if self.reconstruct() != self.original:
            raise ValueError("factorization does not reconstruct its input")

    @property
    def complete(self):
        """All cofactors were split or classified; this is not a proof flag."""
        return not self.remaining

    @property
    def proven(self):
        """True only for a complete result whose factors are all proven."""
        return self.complete and all(
            factor.certainty is utils.Primality.PROVEN
            for factor in self.factors
        )

    def reconstruct(self):
        return (
            self.sign
            * prod(f.value**f.exponent for f in self.factors)
            * (prod(self.remaining))
        )


@lru_cache(maxsize=8)
def _trial_primes(bound):
    """Reuse a bounded exact prime list between factorization calls."""
    return tuple(prime_sieve.prime_sieve(bound + 1))


def factorize_bf(n, *, bound=constants.TRIAL_BOUND):
    """Return (factor/exponent pairs, remainder) using exact trial division."""
    utils.require_integer(n, minimum=1)
    utils.require_integer(bound, "bound", 2)
    factors = []
    if n == 1:
        return factors, n
    # Refresh the exact bound only when division reduces n.
    root = utils.isqrt(n)

    for prime in _trial_primes(bound):
        # Every smaller prime was exhausted, so this remainder is prime.
        if prime > root:
            if n > 1:
                factors.append((n, 1))
                n = 1
            break

        exponent = 0
        while n % prime == 0:
            n //= prime
            exponent += 1
        if exponent:
            factors.append((prime, exponent))
            root = utils.isqrt(n)

    return factors, n


def factorize(
    n,
    verbose=False,
    level=3,
    *,
    seed=None,
    trial_bound=constants.TRIAL_BOUND,
    rho_attempts=constants.RHO_ATTEMPTS,
    rho_evaluations=constants.RHO_EVALUATIONS,
    ecm_curves=constants.MAX_CURVES_ECM,
    ecm_b1=constants.ECM_B1,
    ecm_b2=constants.ECM_B2,
    primality_rounds=constants.PRIMALITY_ROUNDS,
):
    """Factor n, retaining unresolved composites and probable-prime labels.

    Zero has no finite prime factorization and is rejected. Negative inputs
    retain their sign. Levels 3/2/1/0 enable trial+rho+ECM/rho+ECM/ECM/none.
    Local algorithm allowances are explicit; a global time/memory scheduler
    and default p-1 integration remain Phase 2 work.
    """
    utils.require_integer(n)
    if n == 0:
        raise ValueError("zero has no finite prime factorization")
    utils.require_integer(level, "level", 0)
    if level > 3:
        raise ValueError("level must be between 0 and 3")
    utils.require_integer(trial_bound, "trial_bound", 2)
    utils.require_integer(rho_attempts, "rho_attempts", 0)
    utils.require_integer(rho_evaluations, "rho_evaluations", 0)
    utils.require_integer(ecm_curves, "ecm_curves", 0)
    utils.require_integer(ecm_b1, "ecm_b1", 2)
    utils.require_integer(ecm_b2, "ecm_b2", ecm_b1)
    utils.require_integer(primality_rounds, "primality_rounds", 1)

    # Preserve Random's seed contract without allocating it on exact paths.
    if seed is not None and not isinstance(
        seed, (int, float, str, bytes, bytearray)
    ):
        raise TypeError(
            "seed must be None, int, float, str, bytes or bytearray"
        )

    original = n
    sign = -1 if n < 0 else 1
    n = abs(n)
    if n == 1:
        return FactorizationResult(original, sign, (), ())

    generator = (
        utils.resolve_rng(seed) if n >= utils.DETERMINISTIC_LIMIT else None
    )
    initial = utils.classify_prime(
        n, tolerance=primality_rounds, rng=generator
    )
    if initial is not utils.Primality.COMPOSITE:
        # A prime must not scan the whole trial table before its proof/test.
        return FactorizationResult(
            original, sign, (PrimeFactor(n, 1, initial),), ()
        )

    counts = {}
    classifications = {n: initial}
    remaining = []
    if level >= 3:
        trial_factors, n = factorize_bf(n, bound=trial_bound)
        for prime, exponent in trial_factors:
            counts[prime] = counts.get(prime, 0) + exponent
            classifications[prime] = utils.Primality.PROVEN

    pending = [n] if n > 1 else []

    while pending:
        cofactor = pending.pop()
        classification = classifications.get(cofactor)
        if classification is None:
            classification = utils.classify_prime(
                cofactor, tolerance=primality_rounds, rng=generator
            )
            classifications[cofactor] = classification

        if classification is not utils.Primality.COMPOSITE:
            counts[cofactor] = counts.get(cofactor, 0) + 1
            continue

        if level >= 1:
            # An exact square splits directly, avoiding random root searches.
            root = utils.isqrt(cofactor)
            if root * root == cofactor:
                pending.extend((root, root))
                continue

        divisor = None
        if level >= 1 and generator is None:
            # Exact preprocessing and proven primes need no random state.
            generator = utils.resolve_rng(seed)
        if level >= 2 and cofactor <= constants.SIZE_THRESHOLD_RHO:
            divisor = pollard_rho.factorize_rho(
                cofactor,
                verbose,
                rng=generator,
                max_attempts=rho_attempts,
                max_evaluations=rho_evaluations,
                _known_composite=True,
            )

        if level >= 1 and not utils.valid_divisor(divisor, cofactor):
            divisor = ecm.factorize_ecm(
                cofactor,
                verbose,
                rng=generator,
                max_curves=ecm_curves,
                b1=ecm_b1,
                b2=ecm_b2,
                _known_composite=True,
            )

        if utils.valid_divisor(divisor, cofactor):
            # Both children remain pending until independently classified.
            pending.extend((divisor, cofactor // divisor))
        else:
            # Failure preserves the cofactor for result reconstruction.
            remaining.append(cofactor)

    factors = tuple(
        PrimeFactor(prime, exponent, classifications[prime])
        for prime, exponent in sorted(counts.items())
    )
    return FactorizationResult(
        original, sign, factors, tuple(sorted(remaining))
    )


def print_factorization(n, result):
    """Format an honest complete/partial result, including probable labels."""
    if n != result.original:
        raise ValueError("result belongs to a different input")
    terms = ["-1"] if result.sign < 0 else []
    for factor in result.factors:
        term = f"{factor.value}^{factor.exponent}"
        if factor.certainty is utils.Primality.PROBABLE:
            term += " (probable prime)"
        terms.append(term)

    terms.extend(f"{cofactor} (unresolved)" for cofactor in result.remaining)
    expression = " * ".join(terms) if terms else "1"
    return f"{n} = {expression}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", nargs="?", type=int)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--bounded", action="store_true")
    parser.add_argument(
        "--method",
        choices=("auto", "qs", "mpqs", "siqs", "sss", "sssf"),
        default="auto",
        help="select a named engine after preprocessing, or auto dispatch",
    )
    qs_options = parser.add_argument_group(
        "QS / MPQS / SIQS search and fallback"
    )
    qs_options.add_argument(
        "--siqs",
        action="store_true",
        help="enable SIQS after rho/p-1/ECM in the bounded auto portfolio",
    )
    qs_options.add_argument(
        "--qs-base-bound",
        "--siqs-base-bound",
        type=int,
        help="factor-base prime bound (default: 1000)",
    )
    qs_options.add_argument(
        "--qs-half-width",
        "--siqs-half-width",
        type=int,
        help="half-width of each polynomial interval (default: 512)",
    )
    qs_options.add_argument(
        "--qs-max-half-width", "--siqs-max-half-width", type=int
    )
    qs_options.add_argument(
        "--qs-factor-count", "--siqs-factor-count", type=int
    )
    qs_options.add_argument(
        "--qs-family-count",
        "--siqs-family-count",
        type=int,
        help="finite MPQS/SIQS family allowance; QS uses one polynomial",
    )
    qs_options.add_argument("--qs-pool-size", "--siqs-pool-size", type=int)
    qs_options.add_argument(
        "--qs-assignment-policy",
        "--siqs-assignment-policy",
        choices=("reference", "nearest", "flyer"),
    )
    qs_options.add_argument(
        "--qs-polynomials-per-family",
        "--siqs-polynomials-per-family",
        type=int,
    )
    qs_options.add_argument(
        "--qs-residual-bound", "--siqs-residual-bound", type=int
    )
    parser.add_argument("--sss-base-bound", type=int, default=1000)
    parser.add_argument("--sss-rounds", type=int, default=256)
    parser.add_argument("--work-limit", type=int)
    parser.add_argument("--seconds", type=float, default=30)
    parser.add_argument("--cpu-seconds", type=float, default=30)
    parser.add_argument("--memory-mib", type=int)
    parser.add_argument("--fermat-steps", type=int, default=0)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--resume", type=Path)
    parser.add_argument(
        "--ecm-curves", type=int, default=constants.MAX_CURVES_ECM
    )
    args = parser.parse_args()
    use_sss = args.method in ("sss", "sssf")
    use_qs = args.siqs or args.method in ("qs", "mpqs", "siqs")
    if args.siqs and args.method != "auto":
        parser.error("--siqs is an auto fallback; use --method siqs alone")
    qs_parameters = {
        name: getattr(args, "qs_" + name)
        for name in (
            "base_bound",
            "half_width",
            "max_half_width",
            "factor_count",
            "family_count",
            "pool_size",
            "assignment_policy",
            "polynomials_per_family",
        )
        if getattr(args, "qs_" + name) is not None
    }
    if not use_qs and (qs_parameters or args.qs_residual_bound is not None):
        parser.error("QS parameters require --siqs or --method qs/mpqs/siqs")

    if args.memory_mib is None:
        args.memory_mib = 80 if use_sss or use_qs else 8
    if args.work_limit is None:
        args.work_limit = 200_000_000 if use_sss or use_qs else 2_000_000
    try:
        checkpoint = None
        if args.resume:
            if args.resume.stat().st_size > args.memory_mib * 1024 * 1024:
                raise ValueError("checkpoint file exceeds configured cap")
            checkpoint = json.loads(args.resume.read_text())
        number = args.n
        if number is None:
            number = (
                checkpoint["payload"]["state"]["original"]
                if checkpoint
                else int(input("Enter number: "))
            )

        if args.bounded or args.resume or args.checkpoint or use_sss or use_qs:
            from .budget import Budget
            from .portfolio import PortfolioConfig, factorize_bounded
            from .qs.sss import SSSConfig

            parameters = dict(
                ecm_tiers=(
                    (constants.ECM_B1, constants.ECM_B2, args.ecm_curves),
                ),
                memory_bytes=args.memory_mib * 1024 * 1024,
                fermat_steps=args.fermat_steps,
            )
            if use_qs:
                from .qs import SIQSConfig

                # Leave headroom for the parent, schedules and packed resume
                # state; the portfolio validates simultaneous ownership.
                qs_parameters["memory_bytes"] = max(
                    0, args.memory_mib * 1024 * 1024 - 16 * 1024 * 1024
                )
                qs_parameters["mode"] = (
                    "siqs" if args.method == "auto" else args.method
                )
                qs_config = SIQSConfig(**qs_parameters)
                if args.qs_residual_bound is not None:
                    qs_config = replace(
                        qs_config,
                        collector=replace(
                            qs_config.collector,
                            residual_bound=args.qs_residual_bound,
                        ),
                    )
                parameters["siqs"] = qs_config
                if args.method in ("qs", "mpqs", "siqs"):
                    parameters.update(
                        rho_attempts=0, pm1_attempts=0, ecm_tiers=()
                    )

            if use_sss:
                parameters.update(
                    rho_attempts=0,
                    pm1_attempts=0,
                    ecm_tiers=(),
                    sss=SSSConfig(
                        mode=args.method,
                        base_bound=args.sss_base_bound,
                        search_rounds=args.sss_rounds,
                        memory_bytes=max(
                            0, args.memory_mib * 1024 * 1024 - 16 * 1024 * 1024
                        ),
                    ),
                )

            config = PortfolioConfig(**parameters)
            run = factorize_bounded(
                number,
                seed=args.seed,
                config=config,
                checkpoint=checkpoint,
                budget=Budget(
                    work_limit=args.work_limit,
                    seconds=args.seconds,
                    cpu_seconds=args.cpu_seconds,
                ),
            )
            result = run.result
            if args.checkpoint:
                args.checkpoint.write_text(
                    json.dumps(run.checkpoint, indent=2) + "\n"
                )
            if args.verbose:
                print(f"Portfolio: {run.reason}, work={run.work_used}")
                for event in run.events:
                    # The checkpoint stage names the shared job implementation;
                    # CLI diagnostics identify the selected polynomial engine.
                    if event["stage"] == "siqs":
                        event = dict(event, stage=config.siqs.mode)
                    print(event)
        else:
            result = factorize(
                number,
                args.verbose,
                seed=args.seed,
                ecm_curves=args.ecm_curves,
            )
    except (TypeError, ValueError, OSError, KeyError, MemoryError) as error:
        parser.error(str(error))

    print(print_factorization(number, result))
    return 0 if result.complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
