"""Frozen larger-band default comparison with durable finite leases."""

import argparse
import hashlib
import json
import os
import platform
import random
import resource
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path

from ...infrastructure.performance.performance_audit import (
    fingerprint,
    verify_corpus,
)
from ..b1.b1_calibration import save
from ..c9 import c9_growth_study as c9
from ..c10a05 import study as prior

ROOT = Path(__file__).resolve().parents[4]
INPUTS = ROOT / "v2/benchmarks/inputs"
TRAINING = INPUTS / "corpora/c10a05_plus_training.json"
CONTROL = INPUTS / "controls/c10a05_plus_training.json"
SELECTED = INPUTS / "controls/c10a05_plus_selected.json"
CORPUS = INPUTS / "corpora/c10a05_plus_fresh.json"
PROTOCOL = Path(__file__).with_name("protocol.md")
OUTPUT = ROOT / "v2/benchmarks/results/c10a05_plus/experiment"
ARMS = ("current", "c9_window", "c10_window", "collection")
SEEDS = (7, 29, 47)
CAPS = {"train": 600, "confirm": 720}
TOTAL = 36000
PHASE_LIMITS = {"train": 20000, "confirm": 12000}
PREPAID = 3500


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hashes():
    """Bind runtime, every live benchmark helper and this protocol."""
    result = fingerprint(ROOT)
    for path in (
        PROTOCOL,
        INPUTS / "controls/c10a05_sources.json",
        INPUTS / "controls/a7_r5_c10a05_defaults_sources.json",
        *(ROOT / "v2/benchmarks").rglob("*.py"),
    ):
        if "results" not in path.relative_to(ROOT).parts:
            result[str(path.relative_to(ROOT))] = digest(path)
    return result


def configurations():
    """Keep the C9/C10/pacing candidate arithmetic fields paired."""
    from ....qs import SIQSConfig

    prior_control = json.loads(
        (INPUTS / "controls/c10a05_training.json").read_text()
    )
    c10_window = c9.decode(prior_control["configurations"]["50"]["control"])
    c9_window = replace(
        c10_window,
        collector=replace(c10_window.collector, power_plan_bytes=0),
    )
    return {
        "current": SIQSConfig(),
        "c9_window": c9_window,
        "c10_window": c10_window,
        "collection": replace(c10_window, solve_policy="collection"),
    }


def freeze():
    if CONTROL.exists():
        raise ValueError("larger-band control already frozen")
    training = json.loads(TRAINING.read_text())
    verify_corpus(training)
    fixtures = training["fixtures"]
    if len(fixtures) != 3 or any(f["digits"] != 51 for f in fixtures):
        raise ValueError("51-digit training must have three products")
    if len({f["n"] for f in fixtures}) != 3:
        raise ValueError("training inputs must differ")
    config = {a: asdict(c) for a, c in configurations().items()}
    save(
        CONTROL,
        dict(
            schema=1,
            parent_commit="3d4bf42e2bacc9e5d558b60ac2b23b8ae97119b6",
            source=hashes(),
            training_sha256=digest(TRAINING),
            configurations=config,
            fixtures=fixtures,
            seeds=list(SEEDS),
            caps=CAPS,
            total_seconds=TOTAL,
            phase_seconds=PHASE_LIMITS,
            prepaid_seconds=PREPAID,
            selection="completions, then aggregate capped cost; tie current",
        ),
    )


def load(path=CONTROL):
    data = json.loads(path.read_text())
    if data["source"] != hashes():
        raise ValueError("larger-band source differs from freeze")
    expected = digest(CONTROL) if path == SELECTED else digest(TRAINING)
    if data["training_sha256"] != expected:
        raise ValueError("larger-band control or training corpus changed")
    return data


class Expenses:
    """Charge each arm durably, including interrupted reservations."""

    def __init__(self, output):
        self.path = output / "expenses.json"
        self.data = (
            json.loads(self.path.read_text())
            if self.path.exists()
            else dict(prior_wall=PREPAID, prior_cpu=PREPAID, calls=[])
        )
        if (
            self.data["prior_wall"] != PREPAID
            or self.data["prior_cpu"] != PREPAID
        ):
            raise ValueError("prepaid charge changed")

    def used(self):
        return max(
            self.data["prior_wall"]
            + sum(r["wall"] for r in self.data["calls"]),
            self.data["prior_cpu"] + sum(r["cpu"] for r in self.data["calls"]),
        )

    def phase_used(self, phase):
        calls = [r for r in self.data["calls"] if r["phase"] == phase]
        return max(sum(r[k] for r in calls) for k in ("wall", "cpu"))

    def write(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(".tmp")
        temporary.write_text(json.dumps(self.data, indent=2) + "\n")
        temporary.replace(self.path)

    def call(self, phase, label, cap, function):
        if self.used() + cap > TOTAL or (
            self.phase_used(phase) + cap > PHASE_LIMITS[phase]
        ):
            raise RuntimeError("larger-band allowance cannot cover next call")
        if (
            sum(
                p.stat().st_size
                for p in self.path.parent.rglob("*")
                if p.is_file()
            )
            > 256 * 2**20
        ):
            raise RuntimeError("larger-band capture disk cap")
        record = dict(
            phase=phase, label=label, wall=cap, cpu=cap, status="reserved"
        )
        self.data["calls"].append(record)
        self.write()
        started, cpu = time.monotonic(), time.process_time()
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        child_before = children.ru_utime + children.ru_stime
        try:
            result = function()
            record["status"] = "complete"
            return result
        except BaseException:
            record["status"] = "failed"
            raise
        finally:
            children = resource.getrusage(resource.RUSAGE_CHILDREN)
            wall = time.monotonic() - started
            spent_cpu = (
                time.process_time()
                - cpu
                + children.ru_utime
                + children.ru_stime
                - child_before
            )
            record.update(
                actual_wall=wall,
                actual_cpu=spent_cpu,
                wall=wall
                if record["status"] == "complete"
                else max(cap, wall),
                cpu=(
                    spent_cpu
                    if record["status"] == "complete"
                    else max(cap, spent_cpu)
                ),
            )
            self.write()


def read_rows(output, phase):
    return [
        json.loads(path.read_text())
        for path in sorted(output.glob(f"{phase}-51-[0-9][0-9]-*.json"))
    ]


def warm(expense, output, phase, fixture, configs):
    for arm, config in configs.items():
        wall = cpu = 0.0
        rows = []
        while min(wall, cpu) < 3:
            row = expense.call(
                phase,
                f"warm-{arm}",
                5,
                lambda: c9.run_one(fixture, 7, config, arm, 3.2),
            )
            rows.append(row)
            wall += row["seconds"]
            cpu += row["cpu_seconds"]
        save(output / f"{phase}-warm-{arm}-{time.time_ns()}.json", rows)


def measure(output, phase, lease_seconds):
    from ..c10a05.selection import instability

    control = load()
    expected = dict(
        control_sha256=digest(CONTROL),
        runtime=sys.version,
        platform=platform.platform(),
    )
    if phase == "train":
        fixtures = control["fixtures"]
        registry = control["configurations"]
    else:
        selected = load(SELECTED)
        corpus = json.loads(CORPUS.read_text())
        verify_corpus(corpus)
        if corpus["selection_sha256"] != digest(SELECTED):
            raise ValueError("fresh corpus belongs to other selection")
        expected.update(
            selection_sha256=digest(SELECTED), corpus_sha256=digest(CORPUS)
        )
        fixtures = corpus["fixtures"]
        chosen = selected["chosen"]
        if chosen == "current":
            return
        registry = {
            a: control["configurations"][a] for a in ("current", chosen)
        }
    if len(fixtures) != 3:
        raise ValueError("three independent inputs required")
    configs = {a: c9.decode(v) for a, v in registry.items()}
    expense, started = Expenses(output), time.monotonic()
    prior.bind_manifest(output, phase, expected)
    finished = output / f"{phase}-finished.json"
    if finished.exists():
        return
    for arm, config in configs.items():
        path = output / f"{phase}-cold-51-{arm}.json"
        if not path.exists():
            expense.call(
                phase,
                f"cold-{arm}",
                35,
                lambda: prior.cold(fixtures[0], arm, config, 5, path),
            )
    warm(expense, output, phase, fixtures[0], configs)
    for target in (9, 18):
        rows = read_rows(output, phase)
        for index in range(len(rows) // len(configs), target):
            fixture, seed = fixtures[index % 3], SEEDS[(index // 3) % 3]
            order = list(configs)
            random.Random(105105 + index).shuffle(order)
            paths = {
                arm: output / f"{phase}-51-{index:02d}-{arm}.json"
                for arm in order
            }
            missing = [arm for arm in order if not paths[arm].exists()]
            reserve = len(missing) * (CAPS[phase] + 2)
            if missing and (
                time.monotonic() - started + reserve > lease_seconds
                or expense.used() + reserve > TOTAL
                or expense.phase_used(phase) + reserve > PHASE_LIMITS[phase]
            ):
                save(
                    output / f"{phase}-yield-{time.time_ns()}.json",
                    dict(
                        next_group=index,
                        reserved_seconds=reserve,
                        aggregate_used=expense.used(),
                        phase_used=expense.phase_used(phase),
                    ),
                )
                return
            for arm in missing:
                row = expense.call(
                    phase,
                    f"sample-{index}-{arm}",
                    CAPS[phase] + 2,
                    lambda arm=arm: c9.run_one(
                        fixture, seed, configs[arm], arm, CAPS[phase]
                    ),
                )
                row.update(sample=index, measurement_process=os.getpid())
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
        rows = read_rows(output, phase)
        normalized = [
            dict(row, arm="control") if row["arm"] == "current" else row
            for row in rows
        ]
        reasons = instability(normalized)
        save(
            output / f"{phase}-stability-{target}-{time.time_ns()}.json",
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
            samples=len(read_rows(output, phase)) // len(configs),
            reasons=reasons,
            aggregate_used=expense.used(),
        ),
    )


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument(
        "action",
        choices=(
            "training-corpus",
            "freeze",
            "train",
            "select",
            "generate",
            "confirm",
            "summarize",
        ),
    )
    parser.add_argument("--lease-seconds", type=int, default=7200)
    args = parser.parse_args()
    c9.runtime()
    if args.action == "training-corpus":
        from .selection import training

        training()
    elif args.action == "freeze":
        freeze()
    elif args.action in ("train", "confirm"):
        with prior.machine_window():
            measure(OUTPUT, args.action, args.lease_seconds)
    else:
        from . import selection

        getattr(selection, args.action)(OUTPUT)


if __name__ == "__main__":
    main()
