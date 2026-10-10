"""Freeze B4 mainline sources and independently certified bounded inputs."""

import hashlib
import json
import random
import subprocess

from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from ...support.paths import (
    BENCHMARK_ROOT,
    REPOSITORY_ROOT,
    source_path,
)

ROOT = REPOSITORY_ROOT
INPUTS = BENCHMARK_ROOT / "inputs"


def digest(path):
    return hashlib.sha256(source_path(path).read_bytes()).hexdigest()


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as output:
        json.dump(value, output, indent=2)
        output.write("\n")


def build():
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    paths = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", revision, "v2"],
        cwd=ROOT,
        text=True,
    ).splitlines()
    paths = [
        path
        for path in paths
        if path.endswith(".py")
        and not path.startswith(("v2/tests/", "v2/benchmarks/"))
    ]
    sources = {
        path: subprocess.check_output(
            ["git", "show", f"{revision}:{path}"], cwd=ROOT, text=True
        )
        for path in paths
    }
    baseline = INPUTS / "baselines/b4_mainline.json"
    write_new(
        baseline,
        dict(
            revision=revision,
            source=sources,
            sha256={
                name: hashlib.sha256(source.encode()).hexdigest()
                for name, source in sources.items()
            },
        ),
    )

    generator = random.Random(2026100942)
    certificates, fixtures = {}, []
    for split in ("training", "held_out"):
        for case, small, bits in (
            ("64", 22, 64),
            ("128", 26, 128),
            ("256", 28, 256),
            ("329", 30, 329),
        ):
            for index in range(2):
                p = certified_prime(small, generator, certificates)
                q = certified_prime(bits - small + 1, generator, certificates)
                fixtures.append(
                    dict(
                        id=f"{split}_{case}_{index}",
                        split=split,
                        case=case,
                        n=p * q,
                        factors=sorted([[p, 1], [q, 1]]),
                    )
                )
        p = certified_prime(17, generator, certificates)
        fixtures.append(
            dict(
                id=f"{split}_power",
                split=split,
                case="power",
                n=p**3,
                factors=[[p, 3]],
            )
        )
    # Balanced controls measure finite stage work without requiring a split.
    for bits in (64, 128, 256, 329):
        p = certified_prime(bits // 2, generator, certificates)
        q = certified_prime(bits - bits // 2 + 1, generator, certificates)
        fixtures.append(
            dict(
                id=f"stage_{bits}",
                split="stages",
                case=str(bits),
                n=p * q,
                factors=sorted([[p, 1], [q, 1]]),
            )
        )
    verify_certificates(certificates)
    corpus = INPUTS / "corpora/b4_corpus.json"
    write_new(
        corpus,
        dict(
            generation_seed=2026100942,
            certificates=certificates,
            fixtures=fixtures,
        ),
    )
    write_new(
        INPUTS / "controls/b4_protocol.json",
        dict(
            schema=1,
            revision=revision,
            controls_sha256=digest(baseline),
            corpus_sha256=digest(corpus),
            seeds=[7, 19, 41],
            stage_seeds=[7, 19, 41],
            kernel_bits=[64, 96, 128, 166, 200, 256, 329, 512, 768, 1024],
            kernel_sigmas=[7, 19, 41],
            kernel_scalar=2**127 + 2026100942,
            cases={
                "64": [200, 20000, 16],
                "128": [1000, 50000, 16],
                "256": [1000, 50000, 16],
                "329": [1000, 50000, 16],
                "power": [50, 2000, 8],
            },
            stage_bounds=[1000, 50000],
            stage_curves=3,
            work_limit=5000000,
            seconds=30,
            cpu_seconds=30,
            memory_bytes=16 * 2**20,
            max_input_bits=331,
            sampling=[[3, 9], [5, 31], [8, 63]],
            relative_iqr_limit=0.15,
            backends=["python-int", "gmpy2-mpz"],
            candidate_limit=5,
            selection="Fastest stable whole-factoring arm per backend on the "
            "pooled training cohort, provided no declared class loses >5 "
            "completion points. Select even when gain is below promotion "
            "threshold, for a bounded fresh confirmation; no further tuning.",
            promotion="Zero correctness failures; same finite limits; fresh "
            "complete-run median improves >=10% with conditional 95% timing "
            "interval above zero, or completion improves >=10 percentage "
            "points with fixture-cluster interval above zero; no declared "
            "class loses >5 completion points. Unstable/censored/inconclusive "
            "evidence retains baseline. Stage/kernel wins cannot promote.",
            timing_scope="Setup, conversion, normalization, recovery, "
            "result reconstruction and checkpoint serialization included. "
            "Kernel diagnostics, stage timings, full factoring, cold starts "
            "and instrumented profiles are reported separately.",
            machine_lock="/private/tmp/factor-performance.lock",
            profile_status="Baseline breakdown must precede candidate freeze.",
        ),
    )


if __name__ == "__main__":
    build()
