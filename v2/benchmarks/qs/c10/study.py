"""Frozen, serialized C10 policy experiments on the settled C9 collector."""

import argparse
import fcntl
import hashlib
import json
import os
import platform
import random
import subprocess
import sys
import time
from contextlib import contextmanager
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

from ...infrastructure.performance.performance_audit import (
    fingerprint,
    verify_corpus,
)
from ..b1.b1_calibration import save
from ..c9 import c9_growth_confirmation as c9_confirmation
from ..c9 import c9_growth_study as c9

ROOT = Path(__file__).resolve().parents[4]
INPUTS = ROOT / "v2/benchmarks/inputs"
CONTROL = INPUTS / "controls/c10_training.json"
SELECTED = INPUTS / "controls/c10_selected.json"
CORPUS = INPUTS / "corpora/c10_fresh_corpus.json"
PROTOCOL = Path(__file__).with_name("protocol.md")
OUTPUT = ROOT / "v2/benchmarks/results/c10/experiment"
SEEDS = (7, 29, 47)
ARMS = ("control", "plan", "quality")
TRAIN_CAPS = {30: 10, 40: 60, 50: 600}
FRESH_CAPS = {30: 10, 40: 60, 50: 1200}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hashes():
    """Freeze runtime, all benchmark helpers and the measurement contract."""
    result = fingerprint(ROOT)
    for path in (PROTOCOL, *(ROOT / "v2/benchmarks").rglob("*.py")):
        if "results" not in path.relative_to(ROOT).parts:
            result[str(path.relative_to(ROOT))] = digest(path)
    return result


def configurations(control):
    """Two isolated policy changes; retain all C9 residual/storage fields."""
    return {
        "control": control,
        "quality": replace(control, assignment_policy="quality"),
        "plan": replace(
            control,
            collector=replace(control.collector, power_plan_bytes=2**21),
        ),
    }


def freeze(c9_revision, decision):
    if CONTROL.exists():
        raise ValueError("C10 training is already frozen")
    if len(c9_revision) != 40 or any(
        c not in "0123456789abcdef" for c in c9_revision
    ):
        raise ValueError("the final C9 commit is required")
    if not decision:
        raise ValueError("the settled C9 decision is required")
    source = INPUTS / "controls/c9_growth_selected.json"
    registry = json.loads(source.read_text())["configurations"]
    configs = {
        str(d): {
            a: asdict(c)
            for a, c in configurations(
                c9.decode(registry[str(d)]["challenger"])
            ).items()
        }
        for d in TRAIN_CAPS
    }
    corpora = {}
    for name in ("c1_fresh_corpus", "c9_growth_fresh_corpus"):
        corpora[name] = json.loads(
            (INPUTS / f"corpora/{name}.json").read_text()
        )
    prior = json.loads((INPUTS / "controls/c9_training.json").read_text())
    fixtures = [
        f for f in corpora["c1_fresh_corpus"]["fixtures"] if f["digits"] == 30
    ]
    growth = corpora["c9_growth_fresh_corpus"]["fixtures"]
    fixtures += [f for f in growth if f["digits"] in (30, 40)]
    fixtures += prior["fixtures"]
    fixtures += [f for f in growth if f["digits"] == 50][:1]
    if any(sum(f["digits"] == d for f in fixtures) != 3 for d in TRAIN_CAPS):
        raise ValueError("training requires three inputs per band")
    save(
        CONTROL,
        dict(
            schema=1,
            source=hashes(),
            c9_final_revision=c9_revision,
            c9_decision=decision,
            c9_configuration_sha256=digest(source),
            source_manifest_sha256=digest(
                INPUTS / "controls/c10_source_pins.json"
            ),
            configurations=configs,
            fixtures=fixtures,
            seeds=list(SEEDS),
            caps=TRAIN_CAPS,
            total_seconds=14400,
            prepaid_qa_seconds=900,
            training_seconds=5400,
            diagnostics_seconds=600,
        ),
    )


def load_control(path=CONTROL):
    value = json.loads(path.read_text())
    if value["source"] != hashes():
        raise ValueError("C10 source differs from its frozen identity")
    if "training_sha256" in value and value["training_sha256"] != digest(
        CONTROL
    ):
        raise ValueError("C10 selection belongs to different training")
    return value


def bind_manifest(output, phase, expected):
    """Refuse mixing a changed selection, corpus or runtime with prior rows."""
    path = output / f"{phase}-manifest.json"
    if path.exists() and json.loads(path.read_text()) != expected:
        raise ValueError("output belongs to another frozen experiment")
    if not path.exists():
        save(path, expected)


def expenses(output):
    return c9_confirmation.Expenses(
        output,
        {"spent_before_confirmation": {"wall": 900, "cpu": 900}},
        {"generation_seconds": 0},
    )


def phase_used(expense, phase):
    calls = [
        r for r in expense.data["calls"] if r["label"].startswith(phase + "-")
    ]
    return max(sum(r[k] for r in calls) for k in ("wall", "cpu"))


def charged(expense, phase, label, cap, function):
    """Apply the phase cap to warmups and failed calls as well as samples."""
    limit = {"train": 5400, "diagnostic": 600}.get(phase)
    if limit is not None and phase_used(expense, phase) + cap > limit:
        raise RuntimeError("C10 phase allowance cannot cover the next call")
    return expense.call(f"{phase}-{label}", cap, function)


@contextmanager
def machine_window():
    """Acquire without waiting; only the owner may release its lease."""
    with c9.LOCK.open("a+") as lease:
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        c9.check_quiet({os.getpid()})
        c9.OWNER.write_text(
            json.dumps(
                dict(
                    owner="C10",
                    pid=os.getpid(),
                    cwd=str(ROOT),
                    started_utc=time.strftime(
                        "%Y-%m-%dT%H:%M:%SZ", time.gmtime()
                    ),
                )
            )
        )
        try:
            yield
        finally:
            c9.OWNER.unlink(missing_ok=True)


def cold(fixture, arm, config, seconds, output):
    """One validated subprocess, reported separately from warm execution."""
    request = output.with_suffix(".request.json")
    save(
        request,
        dict(fixture=fixture, arm=arm, config=asdict(config), seconds=seconds),
    )
    started = time.monotonic()
    completed = subprocess.run(
        [
            sys.executable,
            "-B",
            "-m",
            __spec__.name,
            "cold-worker",
            "--request",
            str(request),
            "--output",
            str(output),
        ],
        text=True,
        capture_output=True,
        timeout=seconds + 30,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr)
    save(
        output.with_suffix(".startup.json"),
        dict(
            seconds=time.monotonic() - started,
            stderr=completed.stderr,
            classification="cold process diagnostic, not a warm sample",
        ),
    )


def warm(fixtures, configs, expense, output, phase, digits):
    """Repeat validated bounded work for >=3 wall and CPU seconds per arm."""
    for arm, config in configs.items():
        rows, wall, cpu = [], 0.0, 0.0
        while min(wall, cpu) < 3:
            row = charged(
                expense,
                phase,
                f"warm-{digits}-{arm}",
                5,
                lambda: c9.run_one(fixtures[0], 7, config, arm, 3.2),
            )
            rows.append(row)
            wall += row["seconds"]
            cpu += row["cpu_seconds"]
        save(
            output / f"{phase}-warm-{digits}-{arm}-{time.time_ns()}.json", rows
        )


def measure(output, phase, digits, lease_seconds):
    from .selection import instability

    control = load_control()
    expected = dict(
        control_sha256=digest(CONTROL),
        runtime=sys.version,
        platform=platform.platform(),
    )
    cap = (TRAIN_CAPS if phase == "train" else FRESH_CAPS)[digits]
    if phase == "train":
        fixtures = [f for f in control["fixtures"] if f["digits"] == digits]
        registry = control["configurations"][str(digits)]
    else:
        selected = load_control(SELECTED)
        corpus = json.loads(CORPUS.read_text())
        if corpus["selection_sha256"] != digest(SELECTED):
            raise ValueError("fresh corpus belongs to a different selection")
        verify_corpus(corpus)
        expected.update(
            selection_sha256=digest(SELECTED), corpus_sha256=digest(CORPUS)
        )
        fixtures = [f for f in corpus["fixtures"] if f["digits"] == digits]
        chosen = selected["choices"].get(str(digits))
        if chosen is None:
            return
        registry = {
            a: control["configurations"][str(digits)][a]
            for a in ("control", chosen)
        }
    configs = {a: c9.decode(v) for a, v in registry.items()}
    expense, started = expenses(output), time.monotonic()
    bind_manifest(output, phase, expected)
    finished = output / f"{phase}-finished-{digits}.json"
    if finished.exists():
        return
    # Cold probes have short diagnostic caps; measured groups use full caps.
    for arm, config in configs.items():
        path = output / f"{phase}-cold-{digits}-{arm}.json"
        if not path.exists():
            charged(
                expense,
                phase,
                f"cold-{digits}-{arm}",
                35,
                lambda: cold(fixtures[0], arm, config, 5, path),
            )
    warm(fixtures, configs, expense, output, phase, digits)
    rows = []
    for target in (9, 18):
        for index in range(len(rows) // len(configs), target):
            fixture, seed = fixtures[index % 3], SEEDS[(index // 3) % 3]
            order = list(configs)
            random.Random(101010 + digits * 100 + index).shuffle(order)
            paths = {
                a: output / f"{phase}-{digits}-{index:02d}-{a}.json"
                for a in order
            }
            missing = [a for a in order if not paths[a].exists()]
            if missing:
                reserve = len(missing) * (cap + 2)
                if (
                    time.monotonic() - started + reserve > lease_seconds
                    or expense.used() + reserve > 14400
                    or phase == "train"
                    and phase_used(expense, phase) + reserve > 5400
                ):
                    save(
                        output / f"{phase}-yield-{digits}.json",
                        dict(
                            next_group=index,
                            reserved_seconds=reserve,
                            aggregate_used=expense.used(),
                            phase_used=phase_used(expense, phase),
                        ),
                    )
                    return

                def group():
                    for arm in missing:
                        row = c9.run_one(fixture, seed, configs[arm], arm, cap)
                        row.update(
                            sample=index, measurement_process=os.getpid()
                        )
                        save(paths[arm], row)
                        print(
                            json.dumps(
                                {
                                    k: row[k]
                                    for k in (
                                        "id",
                                        "arm",
                                        "complete",
                                        "seconds",
                                        "reason",
                                    )
                                }
                            ),
                            flush=True,
                        )

                charged(
                    expense, phase, f"group-{digits}-{index}", reserve, group
                )
            rows.extend(json.loads(paths[a].read_text()) for a in order)
        reasons = instability(rows)
        save(
            output / f"{phase}-stability-{digits}-{target}.json",
            dict(
                reasons=reasons,
                samples=target,
                extension=target == 9 and bool(reasons),
            ),
        )
        if not reasons:
            break
    save(
        finished,
        dict(
            samples=len(rows) // len(configs),
            reasons=reasons,
            aggregate_used=expense.used(),
        ),
    )


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument(
        "action",
        choices=(
            "freeze",
            "train",
            "confirm",
            "select",
            "generate",
            "diagnose",
            "summarize",
            "cold-worker",
        ),
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--digits", type=int, choices=TRAIN_CAPS, default=30)
    parser.add_argument("--lease-seconds", type=int, default=3900)
    parser.add_argument("--c9-revision")
    parser.add_argument("--decision")
    parser.add_argument("--request", type=Path)
    args = parser.parse_args()
    c9.runtime()
    if args.action == "freeze":
        freeze(args.c9_revision or "", args.decision)
    elif args.action == "cold-worker":
        request = json.loads(args.request.read_text())
        original = c9.check_quiet
        # The parent owns the lock and is blocked waiting for this one child.
        with patch.object(
            c9,
            "check_quiet",
            lambda allowed: original(allowed | {os.getppid()}),
        ):
            save(
                args.output,
                c9.run_one(
                    request["fixture"],
                    7,
                    c9.decode(request["config"]),
                    request["arm"],
                    request["seconds"],
                ),
            )
    elif args.action in ("train", "confirm", "generate", "diagnose"):
        with machine_window():
            if args.action == "generate":
                from .selection import generate

                generate(args.output)
            elif args.action == "diagnose":
                from .diagnostics import run

                run(args.output)
            else:
                measure(
                    args.output, args.action, args.digits, args.lease_seconds
                )
    else:
        from .selection import select, summarize

        (select if args.action == "select" else summarize)(args.output)


if __name__ == "__main__":
    main()
