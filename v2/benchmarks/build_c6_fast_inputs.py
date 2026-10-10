"""Freeze optimized C6 records and fresh independently certified inputs."""

import argparse
import hashlib
import json
import random
import signal
from pathlib import Path

from . import c6_fast, c6_study
from .build_c6_inputs import require_runtime
from .build_phase_two_corpus import certified_prime, verify_certificates

HELDOUT = Path(__file__).parent / "inputs/corpora/c6_fast_heldout.json"
PROTOCOL = Path(__file__).parent / "inputs/controls/c6_fast_protocol.json"
GENERATION_SEED = 6104602


class LimitedRandom(random.Random):
    """Bound all random draws used by certificate construction/filtering."""

    def __init__(self, seed):
        super().__init__(seed)
        self.draws = 0

    def getrandbits(self, count):
        self.draws += 1
        if self.draws > 500000:
            raise RuntimeError("certificate generation draw cap")
        return super().getrandbits(count)


def heldout():
    generator, certificates = LimitedRandom(GENERATION_SEED), {}

    def prime(digits):
        bits = (10**digits - 1).bit_length()
        for _ in range(64):
            value = certified_prime(bits, generator, certificates)
            if 10 ** (digits - 1) <= value < 10**digits:
                return value
        raise RuntimeError("decimal-prime candidate cap")

    fixtures = []
    for digits in (40, 50, 60, 70, 80):
        for shape in ("balanced", "small10"):
            first_digits = digits // 2 if shape == "balanced" else 10
            for _ in range(64):
                factors = [prime(first_digits), prime(digits - first_digits)]
                n = factors[0] * factors[1]
                if len(str(n)) == digits and factors[0] != factors[1]:
                    break
            else:
                raise RuntimeError("composite candidate cap")
            fixtures.append(
                dict(
                    id=f"{shape}_{digits}d",
                    n=n,
                    digits=digits,
                    factors=factors,
                    factor_digits=[len(str(p)) for p in factors],
                )
            )
    required = {}

    def retain(value):
        proof = certificates[str(value)]
        required[str(value)] = proof
        if proof["kind"] != "trial":
            retain(proof["q"])

    for fixture in fixtures:
        for value in fixture["factors"]:
            retain(value)
    verify_certificates(required)
    return dict(
        generation_seed=GENERATION_SEED,
        random_draws=generator.draws,
        certificates=required,
        fixtures=fixtures,
        seeds=[56839, 64758],
    )


def freeze():
    original = c6_study.protocol()
    paths = list(original["sha256"])
    paths += [
        "v2/benchmarks/" + name
        for name in (
            "c6_fast.py",
            "c6_fast_study.py",
            "build_c6_fast_inputs.py",
            "inputs/controls/c6_fast_records.json",
            "inputs/corpora/c6_fast_heldout.json",
        )
    ]
    settings = dict(
        baseline_commit=original["baseline_commit"],
        strict_control_commit="6d54a14",
        original_protocol="c6_protocol.json",
        sha256={
            path: hashlib.sha256(
                (c6_study.ROOT / path).read_bytes()
            ).hexdigest()
            for path in paths
        },
        b1=2000,
        b2=147396,
        curves=8,
        seconds_per_attempt=20,
        seeds=[41001, 48920],
        confirmation_seeds=[56839, 64758],
        cases=original["cases"],
        families=["ladder", "prac", "lucas"],
        modes=list(c6_fast.MODES),
        batches=[1, 16, 64],
        backends=["int", "gmp"],
        sampling=[[3, 9], [5, 18], [8, 27]],
        max_relative_iqr=0.15,
        bootstrap_seed=6104603,
        bootstrap_replicates=4000,
        screen="Paired reused full stages on the original training cohort. "
        "Select the lowest stable median paired runtime ratio "
        "per family/backend, "
        "tie within 1% favors tuple then calls then inline and smaller batch. "
        "Retain all rejected/unstable captures.",
        confirmation="Freeze winners before running fresh held-out input/seed "
        "stages and complete campaigns. Positive paired 95% interval, stable "
        "samples, no correctness failure, no >5-point class completion loss. "
        "No universal percentage floor; report cold and amortized costs.",
        search="Bounded offline CF search only if optimized common-executor "
        "results show credible remaining opportunity against the ladder. "
        "Any new family receives independent integer/projective/coverage "
        "verification and a separately frozen confirmation protocol.",
        limits=dict(
            records_per_family=512,
            steps=512,
            point_registers=16,
            batch_records=64,
            source_bytes_per_record=c6_fast.MAX_SOURCE_BYTES,
            source_bytes_per_program=c6_fast.MAX_PROGRAM_SOURCE_BYTES,
            catalog_bytes=1048576,
            generator_seconds=60,
            generator_draws=500000,
            worker_seconds=1200,
            phase_seconds=2700,
        ),
    )
    PROTOCOL.write_text(json.dumps(settings, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", action="store_true")
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    require_runtime()
    with c6_study.performance_window("optimized C6 input freeze"):
        if args.generate:

            def expired(signum, frame):
                raise RuntimeError("input generation exceeded 60 seconds")

            previous = [
                signal.signal(sig, expired)
                for sig in (signal.SIGALRM, signal.SIGPROF)
            ]
            for timer in (signal.ITIMER_REAL, signal.ITIMER_PROF):
                signal.setitimer(timer, 60)
            try:
                for path, data in (
                    (c6_fast.DATA, c6_fast.generate_catalog()),
                    (HELDOUT, heldout()),
                ):
                    raw = json.dumps(
                        data, sort_keys=True, separators=(",", ":")
                    )
                    if len(raw.encode()) > 1048576:
                        raise RuntimeError("generated input exceeds byte cap")
                    with path.open("x") as output:
                        output.write(raw + "\n")
                    print(path.name, len(raw) + 1)
            finally:
                for timer in (signal.ITIMER_REAL, signal.ITIMER_PROF):
                    signal.setitimer(timer, 0)
                for sig, handler in zip(
                    (signal.SIGALRM, signal.SIGPROF), previous
                ):
                    signal.signal(sig, handler)
        if args.freeze:
            freeze()
            print("Frozen", PROTOCOL.name)


if __name__ == "__main__":
    main()
