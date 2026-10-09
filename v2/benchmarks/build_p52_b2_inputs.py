"""Freeze A2/A3 controls and independent B2 training/held-out inputs."""

import hashlib
import json
import random
import subprocess
from pathlib import Path

from .build_phase_two_corpus import certified_prime, verify_certificates

INPUTS = Path(__file__).parent / "inputs"
REVISION = "94caf40"
MODULES = (
    "arithmetic",
    "constants",
    "utils",
    "prime_sieve",
    "budget",
    "preprocessing",
    "prac",
    "ecm",
    "pollard_rho",
    "factor",
    "schedules",
    "ecm_programs",
    "stage_jobs",
    "portfolio",
)


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as output:
        json.dump(value, output, indent=2)
        output.write("\n")


def build():
    root = Path(__file__).parents[2]
    sources = {
        name: subprocess.check_output(
            ["git", "show", f"{REVISION}:v2/{name}.py"], cwd=root, text=True
        )
        for name in MODULES
    }
    control = {
        "revision": subprocess.check_output(
            ["git", "rev-parse", REVISION], cwd=root, text=True
        ).strip(),
        "source": sources,
        "sha256": {
            name: hashlib.sha256(source.encode()).hexdigest()
            for name, source in sources.items()
        },
    }
    generator = random.Random(2026100952)
    certificates, fixtures = {}, []
    for split in ("training", "held_out"):
        for case, bits in (
            ("small", (16, 16)),
            ("medium", (26, 26)),
            ("uneven", (34, 100)),
            ("campaign", (133, 133)),
        ):
            for index in range(3 if case != "campaign" else 1):
                p, q = [
                    certified_prime(size, generator, certificates)
                    for size in bits
                ]
                fixtures.append(
                    dict(
                        id=f"{split}_{case}_{index}",
                        split=split,
                        case=case,
                        n=p * q,
                        factors=sorted([[p, 1], [q, 1]]),
                    )
                )
        # Distinct structured controls exercise powers and recursive factors.
        p = certified_prime(17, generator, certificates)
        fixtures.append(
            dict(
                id=f"{split}_structured",
                split=split,
                case="structured",
                n=p**3,
                factors=[[p, 3]],
            )
        )
    verify_certificates(certificates)
    corpus = dict(
        generation_seed=2026100952,
        seeds=[7, 19, 41],
        certificates=certificates,
        fixtures=fixtures,
    )
    control_path = INPUTS / "baselines/p52_b2_mainline.json"
    corpus_path = INPUTS / "corpora/p52_b2_corpus.json"
    write_new(control_path, control)
    write_new(corpus_path, corpus)
    protocol = dict(
        schema=1,
        revision=control["revision"],
        controls_sha256=hashlib.sha256(control_path.read_bytes()).hexdigest(),
        corpus_sha256=hashlib.sha256(corpus_path.read_bytes()).hexdigest(),
        backend="python-int",
        seeds=[7, 19, 41],
        memory_bytes=16 * 2**20,
        max_input_bits=329,
        segment_size=1024,
        program_bytes=8 * 2**20,
        work_limit=50_000_000,
        seconds=120,
        cpu_seconds=120,
        sampling=[[3, 9], [5, 31], [8, 63]],
        relative_iqr_limit=0.15,
        cases={
            "small": dict(tier=[50, 2000, 8], distances=[8, 16, 24]),
            "medium": dict(tier=[200, 20000, 16], distances=[32, 64, 96]),
            "uneven": dict(tier=[2000, 147396, 8], distances=[128, 384, 768]),
            "campaign": dict(
                tier=[11000, 1900000, 4], distances=[512, 1024, 2048]
            ),
            "structured": dict(tier=[50, 2000, 8], distances=[8, 16, 24]),
        },
        training_arms=[
            "streamed",
            "programs",
            "paired_0",
            "paired_1",
            "paired_2",
            "regenerated",
        ],
        selection="Fastest stable paired D per bound tier, with no training "
        "completion regression; ties choose smaller D. Structured "
        "uses small's selection. Regeneration is diagnostic only.",
        promotion="Fresh held-out complete factoring only: zero correctness "
        "failures, same limits, >=10% median reduction with a 95% "
        "bootstrap interval above zero or >=10 percentage points "
        "completion improvement with interval above zero; no other "
        "declared class loses >5 completion points. Unstable or "
        "inconclusive results retain defaults. Fixed nonsplitting "
        "campaigns cannot promote factoring defaults.",
        uncertainty="Timing bootstrap over repeated fixed-cohort samples; "
        "completion uncertainty clustered by fresh fixture, with "
        "all algorithm seeds kept together. Limited populations "
        "must be reported; no factor-size allocation claim.",
        timing_scope="Config, setup, context, programs, tables, products, "
        "recovery, recursion and checkpoint serialization included. "
        "Cold starts, profiles and work/throughput diagnostics "
        "are separate from whole factoring evidence.",
        machine_lock="/private/tmp/factor-performance.lock",
        status="Prepared; performance execution requires user go-ahead.",
    )
    write_new(INPUTS / "controls/p52_b2_protocol.json", protocol)


if __name__ == "__main__":
    build()
