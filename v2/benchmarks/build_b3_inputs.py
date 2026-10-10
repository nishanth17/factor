"""Freeze B3 mainline controls and disjoint certified production cohorts."""

import argparse
import hashlib
import json
import signal
import subprocess
from math import prod
from pathlib import Path

from .a6_pm1 import require_runtime
from .b3_production import BASELINE, CORPUS, PROTOCOL, ROOT, SCOPES
from .build_c6_fast_inputs import LimitedRandom
from .build_phase_two_corpus import certified_prime, verify_certificates

BASE_COMMIT = "76f06b0110941017b926caa51d341db6221f470b"


def write_new(path, data):
    with path.open("x") as stream:
        json.dump(data, stream, sort_keys=True, indent=2)
        stream.write("\n")


def build_corpus():
    certificates, fixtures = {}, []
    for split, seed in (("training", 93001), ("confirmation", 100921)):
        generator = LimitedRandom(seed)

        def prime(digits):
            for _ in range(64):
                value = certified_prime(
                    (10**digits - 1).bit_length(), generator, certificates
                )
                if len(str(value)) == digits:
                    return value
            raise RuntimeError("prime decimal-size draw cap")

        for digits in (40, 50, 60, 70, 80):
            for shape in ("balanced", "small10"):
                target = digits // 2 if shape == "balanced" else 10
                for _ in range(64):
                    factors = [prime(target), prime(digits - target)]
                    n = prod(factors)
                    if len(str(n)) == digits and factors[0] != factors[1]:
                        break
                else:
                    raise RuntimeError("composite decimal-size draw cap")
                fixtures.append(
                    dict(
                        id=f"{shape}_{digits}d",
                        split=split,
                        n=n,
                        factors=factors,
                        factor_digits=list(
                            map(lambda p: len(str(p)), factors)
                        ),
                    )
                )
        for index, sizes in enumerate(((6, 6, 8), (6, 10, 20))):
            factors = [prime(size) for size in sizes]
            fixtures.append(
                dict(
                    id=f"recursive_{index}",
                    split=split,
                    n=prod(factors),
                    factors=factors,
                    factor_digits=list(sizes),
                )
            )
    verify_certificates(certificates)
    old_values = set()
    for path in CORPUS.parent.glob("*.json"):
        data = json.loads(path.read_text())
        if isinstance(data, dict):
            for fixture in data.get("fixtures", []):
                if isinstance(fixture, dict) and "n" in fixture:
                    old_values.add(fixture["n"])
    values = [fixture["n"] for fixture in fixtures]
    if len(set(values)) != len(values) or set(values) & old_values:
        raise ValueError("cohorts are not fresh and disjoint")
    return dict(
        generation_seeds=[93001, 100921],
        fixtures=fixtures,
        certificates=certificates,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    require_runtime()
    if any(path.exists() for path in (BASELINE, CORPUS, PROTOCOL)):
        raise ValueError("B3 freeze is immutable; choose new input names")
    signal.signal(
        signal.SIGALRM,
        lambda *unused: (_ for _ in ()).throw(
            TimeoutError("B3 generation cap")
        ),
    )
    signal.alarm(60)
    try:
        corpus = build_corpus()
    finally:
        signal.alarm(0)
    paths = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", BASE_COMMIT, "v2"], text=True
    ).splitlines()
    paths = [
        path
        for path in paths
        if path.endswith(".py")
        and not any(
            part in ("tests", "benchmarks", "audit")
            for part in Path(path).parts
        )
    ]
    sources = {
        path: subprocess.check_output(
            ["git", "show", BASE_COMMIT + ":" + path], text=True
        )
        for path in paths
    }
    write_new(
        BASELINE,
        dict(
            commit=BASE_COMMIT,
            source=sources,
            sha256={
                p: hashlib.sha256(s.encode()).hexdigest()
                for p, s in sources.items()
            },
        ),
    )
    write_new(CORPUS, corpus)
    frozen_paths = [ROOT / path for path in paths]
    frozen_paths += [
        ROOT / "v2/ecm_chains.py",
        ROOT / "v2/ecm_chain_records.py",
        Path(__file__),
        Path(__file__).with_name("b3_production.py"),
        BASELINE,
        CORPUS,
        ROOT / "v2/benchmarks/inputs/controls/c6_fast_records.json",
    ]
    write_new(
        PROTOCOL,
        dict(
            base_commit=BASE_COMMIT,
            sha256={
                str(p.relative_to(ROOT)): hashlib.sha256(
                    p.read_bytes()
                ).hexdigest()
                for p in frozen_paths
            },
            training_seeds=[31001, 38920],
            confirmation_seeds=[46839, 54758],
            scopes=list(SCOPES),
            backends=["python-int", "gmpy2-mpz"],
            b1=2000,
            b2=147396,
            curves=8,
            chunk=16,
            work=2_000_000,
            seconds=20,
            cpu_seconds=20,
            memory_bytes=33554432,
            chain_bytes=8388608,
            program_bytes=262144,
            segment_size=128,
            sampling=[[3, 9], [5, 18], [8, 27]],
            max_relative_iqr=0.15,
            bootstrap_seed=193001,
            bootstrap_replicates=4000,
            selection=(
                "No new chain/kernel/search choices. On training compare "
                "fresh1 and 2/4/8/16 curves with preparation charged per "
                "campaign. Production reuse routing remains explicit "
                "opt-in, B1=2000, chunk16, tier >=8 curves. Confirm all "
                "scopes on untouched inputs; retain ladder defaults on "
                "inconclusive/full-run losses."
            ),
            acceptance=(
                "Zero invalid result/coordinate failure. Positive paired "
                "median saving with 95% bootstrap interval, CPU and both "
                "chronological halves positive; stable samples and <=5 "
                "percentage-point completion regression in each declared "
                "fixture class. No universal percentage floor. Do not "
                "extend stable inconclusive comparisons. Claims limited to "
                "confirmed scope."
            ),
            accounting=(
                "Timers include actual plan construction, verification, "
                "conversion, recovery, failed curves, output validation and "
                "checkpoint JSON. Atomic strict/prime-unit replay reserves "
                "are charged on every chain chunk; cache misses rebuild "
                "under the same cumulative allowance. Certificates and "
                "independent affine target construction precede warmed "
                "timing and belong to separate cold totals. Stage reuse "
                "means a fresh plan per fixture/seed campaign, not an "
                "uncharged global cache. Eviction is a low-level "
                "2000/1999/2000 stress, outside enabled routing."
            ),
            limits=(
                "60-second/500000-draw corpus generation; 1200-second "
                "worker; 2700-second phase. No online chain search. Cold "
                "nine processes per arm/backend kept separate; no "
                "instrumented profile is performance evidence."
            ),
        ),
    )
    print(
        json.dumps(
            dict(fixtures=len(corpus["fixtures"]), source_modules=len(sources))
        )
    )


if __name__ == "__main__":
    main()
