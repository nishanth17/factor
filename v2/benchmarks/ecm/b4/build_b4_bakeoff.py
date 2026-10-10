"""Generate disjoint certified controls for the requested bakeoff."""

import random
import signal
import subprocess
from pathlib import Path

from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from . import b4_bakeoff as bakeoff
from . import b4_common as common
from . import b4_kernels as kernels
from . import b4_study as study
from .build_b4_inputs import write_new

GENERATION_SEED = 2026100957


def build():
    """Create new inputs once; a finite alarm bounds prime generation."""
    common.require_runtime()
    study.verify_freeze()
    original, old = common.inputs()
    if any(
        path.exists()
        for path in (bakeoff.CORPUS, bakeoff.PROTOCOL, bakeoff.FREEZE)
    ):
        raise ValueError("bakeoff inputs already exist")
    generator = random.Random(GENERATION_SEED)
    certificates, fixtures = {}, []
    used = {fixture["n"] for fixture in old["fixtures"]}
    signal.alarm(120)
    try:
        for split in ("screen", "confirmation"):
            for case, small, bits in (
                ("64", 22, 64),
                ("128", 26, 128),
                ("256", 28, 256),
                ("329", 30, 329),
                ("power", 17, 51),
            ):
                for index in range(1 if case == "power" else 2):
                    for _ in range(100):
                        p = certified_prime(small, generator, certificates)
                        if case == "power":
                            number, factors = p**3, [[p, 3]]
                        else:
                            q = certified_prime(
                                bits - small + 1, generator, certificates
                            )
                            number, factors = p * q, sorted([[p, 1], [q, 1]])
                        if number not in used:
                            break
                    else:
                        raise ValueError(
                            "bounded distinct input search failed"
                        )
                    used.add(number)
                    fixtures.append(
                        dict(
                            id=f"bakeoff_{split}_{case}_{index}",
                            split=split,
                            case=case,
                            n=number,
                            factors=factors,
                        )
                    )
        verify_certificates(certificates)
    finally:
        signal.alarm(0)
    write_new(
        bakeoff.CORPUS,
        dict(
            generation_seed=GENERATION_SEED,
            certificates=certificates,
            fixtures=fixtures,
        ),
    )
    protocol = dict(original)
    protocol.update(
        schema=2,
        parent_protocol_sha256=common.digest(common.PROTOCOL),
        corpus_sha256=common.digest(bakeoff.CORPUS),
        seeds=[11, 29, 53],
        order_seed=GENERATION_SEED,
        sample_cpu_seconds=0.5,
        sampling=[[3, 9], [5, 27], [8, 63]],
        window_seconds=2700,
        selection="Fastest positive, stable paired complete-cohort gain per "
        "backend with no >5-point class completion regression. "
        "The independent "
        "confirmation tests only that frozen winner against its baseline.",
        promotion="No universal percentage floor. Require zero correctness "
        "failures, matched finite budgets, stable baseline/candidate/paired "
        "ratio IQR <=15%, fresh paired 95% timing interval above zero, "
        "positive "
        "CPU gain, positive effects in both chronological halves and positive "
        "screen direction. Scope claims to the fixed workload; review class "
        "regressions and maintenance before any production integration. "
        "Only unstable pairs extend to the frozen 27/63 ceilings; a stable "
        "interval crossing zero retains baseline without further repeats. "
        "Unresolved evidence retains baseline. Completion "
        "improvements are reported, not promoted from this timing-only study.",
        timing_design="Counterbalanced rotated arm orders and alternating "
        "backend order in independent process blocks. Every process pays at "
        "least three seconds validated warmup before one >=0.5 CPU-second "
        "sample, averaged over validated whole cohorts. Stable arms stop at "
        "nine pairs; only unstable candidates and their matched baseline "
        "continue. Parent holds exclusive "
        "machine lock through all children and extension rounds. "
        "Initialization, "
        "conversion, normalization, recovery and output/checkpoint validation "
        "remain inside each cohort. Cold/profiles are excluded.",
        reused_diagnostics="Original frozen kernel/stage diagnostics and "
        "separately instrumented baseline share/Amdahl evidence remain valid "
        "for their original workloads; no claim of identical shares on these "
        "new inputs. This follow-up measures complete ECM-only portfolios.",
    )
    write_new(bakeoff.PROTOCOL, protocol)
    paths = [
        Path(bakeoff.__file__),
        Path(__file__),
        bakeoff.PROTOCOL,
        bakeoff.CORPUS,
        study.FREEZE,
    ]
    write_new(
        bakeoff.FREEZE,
        dict(
            arms=list(kernels.ARMS),
            parent_revision=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=bakeoff.ROOT, text=True
            ).strip(),
            source_hashes={
                str(path.relative_to(bakeoff.ROOT)): common.digest(path)
                for path in paths
            },
        ),
    )
    bakeoff.verify_inputs()
    print(bakeoff.FREEZE)


if __name__ == "__main__":
    build()
