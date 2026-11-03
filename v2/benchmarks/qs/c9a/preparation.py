"""Numerical structural preparation without collecting or factoring inputs."""

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

from .... import qs
from ....execution.budget import Budget
from ..c10a04.accounting import atomic, digest
from .configuration import BASE_COMMIT, SEEDS, configuration
from .control import MANIFEST
from .inputs import generate, validate
from .plan import build, limits
from .plan import validate as validate_plan

ROOT = Path(__file__).resolve().parents[4]
INPUTS = ROOT / "v2/benchmarks/inputs/c9a"


def structural(corpus):
    """Resolve multiplier, base and assignments without sieving."""
    resolved, work = {}, 0
    for row in corpus["fixtures"]:
        budget = Budget(work_limit=10**11, seconds=15, cpu_seconds=15)
        job = qs.SIQSJob(
            row["n"],
            seed=SEEDS[0],
            config=configuration(row["digits"]),
            budget=budget,
        )
        job._setup()
        if (
            job.divisor is not None
            or job.base is None
            or job.engine is not None
        ):
            raise ValueError(
                "structural input unexpectedly split or collected"
            )
        resolved[row["id"]] = dict(
            bound=job.base.bound,
            multiplier=job.multiplier,
            entries=len(job.base.entries),
            n_mod8=row["n"] % 8,
            kn_mod8=(row["n"] * job.multiplier) % 8,
            base_workspace_bytes=job.base.workspace_bytes,
            work=budget.used,
        )
        work += budget.used
    return dict(resolved=resolved, work=work)


def prepare(exclusions, *, band=None, output=INPUTS, seconds=590):
    """Generate once, retain proof inputs, then resolve literal arm domains."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if (output / "corpus.json").exists() or (
        output / "structural.json"
    ).exists():
        raise ValueError(
            "input identity already exists; never replace a corpus"
        )
    data = json.loads(Path(exclusions).read_text())
    corpus = generate(data["integers"], band=band, seconds=seconds)
    corpus["overlap_manifest_sha256"] = digest(Path(exclusions))
    atomic(output / "corpus.json", corpus)
    resolved = structural(corpus)
    atomic(output / "structural.json", resolved)
    return dict(
        work=corpus["work"] + resolved["work"], moduli=len(corpus["fixtures"])
    )


def source_manifest():
    head = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    names = subprocess.check_output(
        ["git", "ls-files", "v2"], cwd=ROOT, text=True
    ).splitlines()
    selected = [
        name
        for name in names
        if name.endswith(".py")
        or name.startswith("v2/benchmarks/inputs/c9a/")
        or name.startswith("v2/benchmarks/inputs/controls/c9a_")
    ]
    files = {name: digest(ROOT / name) for name in selected}
    return dict(
        commit=head,
        files=files,
        interpreter=str(Path(sys.executable).resolve()),
        interpreter_sha256=digest(Path(sys.executable).resolve()),
        implementation=platform.python_implementation(),
        version=platform.python_version(),
        control=json.loads(MANIFEST.read_text()),
        base_commit=BASE_COMMIT,
    )


def require_clean():
    dirty = subprocess.check_output(
        ["git", "status", "--porcelain", "--", "v2"], cwd=ROOT, text=True
    )
    if dirty.strip():
        raise RuntimeError("commit preparation before freezing execution")


def freeze(output):
    """Freeze an execution identity from committed preparation."""
    require_clean()
    corpus = json.loads((INPUTS / "corpus.json").read_text())
    validate(corpus)
    resolved = json.loads((INPUTS / "structural.json").read_text())["resolved"]
    plan = build(corpus, resolved)
    validate_plan(plan)
    data = dict(
        plan=plan,
        limits=limits(plan),
        source=source_manifest(),
        corpus_sha256=digest(INPUTS / "corpus.json"),
    )
    data["sha256"] = hashlib.sha256(
        json.dumps(data, sort_keys=True).encode()
    ).hexdigest()
    output = Path(output)
    if output.exists():
        raise ValueError("freeze already exists")
    atomic(output, data)
    return data


def load_freeze(path):
    data = json.loads(Path(path).read_text())
    expected = data.pop("sha256")
    if (
        hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()
        != expected
    ):
        raise ValueError("freeze contents changed")
    data["sha256"] = expected
    if source_manifest() != data["source"]:
        raise ValueError(
            "source, inputs or interpreter differ from the freeze"
        )
    return data
