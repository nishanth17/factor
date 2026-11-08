"""Frozen selected-kernel integration; no search or automatic adoption."""

import argparse
import fcntl
import gc
import json
import os
import platform
import random
import statistics
import subprocess
import sys
import time
from collections import Counter
from contextlib import contextmanager
from dataclasses import asdict, replace
from math import prod
from pathlib import Path

from v2.benchmarks.infrastructure.resource_control import recording_budget

from ....common import utils
from ....execution import budget as live_budget
from ....portfolio import _decode_siqs_config
from ....qs import defaults, siqs
from ...infrastructure.performance.performance_audit import fingerprint
from ..c9 import c9_growth_study as c9
from ..c10a04.accounting import Ledger, digest, output_lock, receipt
from . import analysis, inputs, runtime

ROOT = runtime.ROOT
CONTROL = inputs.INPUTS / "controls/a1a_integration_frozen.json"
WARM = inputs.INPUTS / "corpora/a1a_integration_warmup.json"
PROTOCOL = Path(__file__).with_name("protocol.md")
OUTPUT = ROOT / "v2/benchmarks/results/a1a_integration/frozen"
SELECTED = "density+presieve+scan_fused"
CAPS = {50: 480, 55: 1159}
WORK = 10**13
ARMS = ("control", "selected")
LIMITS = {
    "generation": dict(wall=120, cpu=120, work=0),
    "warm50": dict(wall=2800, cpu=2800, work=80 * WORK),
    "warm55": dict(wall=2800, cpu=2800, work=80 * WORK),
    "base50": dict(wall=36 * 485, cpu=36 * 485, work=36 * WORK),
    "base55": dict(wall=36 * 1164, cpu=36 * 1164, work=36 * WORK),
    "extension50": dict(wall=36 * 485, cpu=36 * 485, work=36 * WORK),
    "extension55": dict(wall=36 * 1164, cpu=36 * 1164, work=36 * WORK),
    "orchestration": dict(wall=600, cpu=600, work=0),
    "analysis": dict(wall=600, cpu=600, work=0),
    "disk_bytes": 128 * 2**20,
}
LIMITS["total"] = {
    key: sum(v[key] for v in LIMITS.values() if isinstance(v, dict))
    for key in ("wall", "cpu", "work")
}


def sources():
    """Bind all live runtime and transitive benchmark helpers, not results."""
    hashes = fingerprint(ROOT)
    for path in (ROOT / "v2/benchmarks").rglob("*.py"):
        if not {"results", "history", "upstream"} & set(path.parts):
            hashes[str(path.relative_to(ROOT))] = digest(path)
    for path in (PROTOCOL, runtime.BASELINE, WARM, inputs.EXCLUSIONS):
        hashes[str(path.relative_to(ROOT))] = digest(path)
    for path in (ROOT / "v2/tests").glob("test_a1a*.py"):
        hashes[str(path.relative_to(ROOT))] = digest(path)
    return hashes


def policy():
    configurations = {}
    for digits in CAPS:
        cfg = defaults.default_siqs_config(10 ** (digits - 1))
        selected = replace(
            cfg, collector=replace(cfg.collector, marking=SELECTED)
        )
        configurations[str(digits)] = {
            "control": asdict(cfg),
            "selected": asdict(selected),
        }
    return dict(
        schema=1,
        start_commit=runtime.START,
        study_commit="d61df6c53d8b52e851cb238e0285626b54ee3746",
        source=sources(),
        caps={str(k): v for k, v in CAPS.items()},
        work_per_call=WORK,
        limits=LIMITS,
        arms=list(ARMS),
        seed=7,
        selected=SELECTED,
        warm_blocks=1056,
        warm_cap=30,
        warm_max_calls=20,
        warm_min_seconds=3,
        warm_range_fraction=0.10,
        fresh_inputs_per_band=9,
        opposite_orders=2,
        drift_fraction=0.15,
        extension_max=1,
        alpha=0.005,
        bootstrap_draws=50000,
        configurations=configurations,
        python_version=sys.version,
        scope="native exact50/55; 1GiB owned memory",
    )


def require_committed(path):
    name = str(path.relative_to(ROOT))
    saved = subprocess.check_output(["git", "show", "HEAD:" + name], cwd=ROOT)
    if saved != path.read_bytes():
        raise ValueError("commit frozen source/control/corpus first: " + name)


@contextmanager
def machine_window():
    """Inspect coordination, then own the shared lease for this whole stage."""
    for name in (
        "factor-cpu-coordination.json",
        "factor-performance-owner.json",
    ):
        path = Path("/private/tmp") / name
        print(
            "coordination before lock:",
            name,
            path.read_text() if path.exists() else "absent",
            flush=True,
        )
    with c9.LOCK.open("a+") as lease:
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        c9.check_quiet({os.getpid()})
        previous = c9.OWNER.read_bytes() if c9.OWNER.exists() else None
        c9.OWNER.write_text(
            json.dumps(
                dict(
                    owner="A1a integration",
                    pid=os.getpid(),
                    cwd=str(ROOT),
                    control_sha256=digest(CONTROL),
                )
            )
        )
        try:
            yield
        finally:
            if previous is None:
                c9.OWNER.unlink(missing_ok=True)
            else:
                c9.OWNER.write_bytes(previous)


def load():
    control = json.loads(CONTROL.read_text())
    if control != policy():
        raise ValueError("integration source/configuration/protocol changed")
    runtime.checked_baseline()
    inputs.check_warmup(json.loads(WARM.read_text()))
    return control


def configuration(n, arm, control_runtime):
    """Change only marking at the exact tested band and owned allowance."""
    digits = len(str(abs(n)))
    if digits not in CAPS or arm not in ARMS:
        raise ValueError("outside frozen band/arm")
    # Decode the owned historical control into current classes. Inheriting
    # today's factory would change more than the frozen marking treatment.
    cfg = control_runtime[0].default_siqs_config(n)
    if arm == "selected":
        cfg = _decode_siqs_config(asdict(cfg))
        # Retain the previously integrated face/queue compatibility delta.
        cfg = replace(
            cfg,
            sign_policy="swap_high",
            collector=replace(cfg.collector, sparse_queue_bytes=32 * 2**20),
        )
    if cfg.backend != "python-int" or cfg.memory_bytes != 2**30:
        raise ValueError("outside frozen native/owned-memory envelope")
    if arm == "selected":
        cfg = replace(cfg, collector=replace(cfg.collector, marking=SELECTED))
    return cfg


def validate(row, fixture):
    """Check reconstruction, proper factors, certainty and finite resources."""
    if prod(row["factors"]) * prod(row["remaining"]) != fixture["n"]:
        raise ValueError("result does not reconstruct")
    if Counter(row["factors"]) - Counter(fixture["factors"]):
        raise ValueError("unexpected returned factor")
    if row["complete"] != (not row["remaining"]):
        raise ValueError("completion flag changed")
    if row["divisor"] is not None and not utils.valid_divisor(
        row["divisor"], fixture["n"]
    ):
        raise ValueError("improper returned divisor")
    # Fixture proofs validate primality independently. Expected API labels
    # still follow the classifier's proven range; proofs do not upgrade them.
    labels = [
        "proven_prime"
        if utils.deterministic_bases(p) is not None
        else "probable_prime"
        for p in row["factors"]
    ]
    if row["certainty"] != labels or "composite" in labels:
        raise ValueError("certainty label changed")
    if row["work"] > WORK or row["workspace_bytes"] > 2**30:
        raise ValueError("finite work/owned memory exceeded")
    if row["rss_bytes"] > 3 * 2**30:
        raise ValueError("finite process RSS exceeded")


def run_one(fixture, arm, seconds, control_runtime, *, blocks=None):
    """Time setup through exact output validation and graph teardown."""
    started, cpu = time.monotonic(), time.process_time()
    last_poll = started

    def poll():
        nonlocal last_poll
        if time.monotonic() - last_poll >= 1:
            c9.check_quiet({os.getpid()})
            if c9._rss_bytes() > 3 * 2**30:
                raise MemoryError("RSS allowance exceeded")
            last_poll = time.monotonic()
        return False

    qs_module, budget_module, utils_module = (
        control_runtime if arm == "control" else (siqs, live_budget, utils)
    )
    budget = recording_budget(
        budget_module.Budget,
        work_limit=WORK,
        seconds=seconds,
        cpu_seconds=seconds,
        cancelled=poll,
    )
    cfg = configuration(fixture["n"], arm, control_runtime)
    job = result = None
    factors, certainty, remaining = [], [], [fixture["n"]]
    divisor, reason, stats, base = None, "unstarted", {}, {}
    try:
        job = qs_module.SIQSJob(
            fixture["n"], seed=7, config=cfg, budget=budget
        )
        result = job.run(max_blocks=blocks)
        divisor, reason, stats = result.divisor, result.reason, result.stats
        if (divisor or 1) * result.cofactor != fixture["n"]:
            raise ValueError("SIQS result does not reconstruct")
        if job.base is not None:
            base = dict(
                count=len(job.base.entries),
                bound=job.base.bound,
                multiplier=job.multiplier,
            )
        if divisor is not None:
            if not utils_module.valid_divisor(divisor, fixture["n"]):
                raise ValueError("improper divisor")
            remaining = [divisor, fixture["n"] // divisor]
            for child in tuple(remaining):
                budget.consume(child.bit_length() ** 2)
                label = utils_module.classify_prime(
                    child, rng=random.Random(7)
                )
                if label is utils_module.Primality.COMPOSITE:
                    break
                factors.append(child)
                certainty.append(label.value)
                remaining.pop(0)
    except budget_module.BudgetExhaustedError:
        reason = budget.reason
    except MemoryError:
        reason = "storage_refusal"
    row = dict(
        id=fixture["id"],
        digits=fixture["digits"],
        seed=7,
        arm=arm,
        divisor=divisor,
        factors=factors,
        certainty=certainty,
        remaining=remaining,
        complete=not remaining,
        reason=reason,
        work=budget.used,
        stats=stats,
        base=base,
        configuration=asdict(cfg),
        workspace_bytes=stats.get("workspace_bytes", 0),
        rss_bytes=c9._rss_bytes(),
        cap_seconds=seconds,
    )
    validate(row, fixture)
    job = result = None
    gc.collect()
    row.update(
        seconds=time.monotonic() - started,
        cpu_seconds=time.process_time() - cpu,
    )
    return row


def semantic(row):
    return (
        row["divisor"],
        row["factors"],
        row["remaining"],
        row["certainty"],
        row["base"],
        [row["stats"].get(k) for k in ("blocks", "polynomials", "relations")],
    )


def warmup(ledger, digits, stage, control_runtime):
    fixture = next(
        f
        for f in json.loads(WARM.read_text())["fixtures"]
        if f["digits"] == digits
    )
    recent = {arm: [] for arm in ARMS}
    totals = {arm: [0.0, 0.0] for arm in ARMS}
    for attempt in range(20):
        # Both arms receive the same prefix count and bounded warmup schedule.
        for arm in ARMS if attempt % 2 == 0 else ARMS[::-1]:
            row = ledger.call(
                f"warm{digits}",
                f"warm-{stage}-{digits}-{attempt}-{arm}",
                35,
                WORK,
                dict(stage=stage, digits=digits, arm=arm, attempt=attempt),
                lambda arm=arm: run_one(
                    fixture, arm, 30, control_runtime, blocks=1056
                ),
            )
            if row["reason"] not in ("paused", "factor_found", "complete"):
                return False
            recent[arm].append(row)
            totals[arm][0] += row["seconds"]
            totals[arm][1] += row["cpu_seconds"]
        if semantic(recent["control"][-1]) != semantic(recent["selected"][-1]):
            raise ValueError("warmup semantic disagreement")
        stable = all(
            len(recent[arm]) >= 3
            and min(totals[arm]) >= 3
            and all(
                semantic(r) == semantic(recent[arm][-1])
                for r in recent[arm][-3:]
            )
            and all(
                (max(values) - min(values)) / statistics.median(values) <= 0.10
                for values in (
                    [r[field] for r in recent[arm][-3:]]
                    for field in ("seconds", "cpu_seconds")
                )
            )
            for arm in ARMS
        )
        if stable:
            return True
    return False


def measured_rows(ledger, digits):
    rows = []
    for entry in ledger.data["calls"]:
        if entry["phase"] in (f"base{digits}", f"extension{digits}"):
            wrapper = json.loads(
                (ledger.output / entry["receipt"]).read_text()
            )
            rows.append(wrapper["result"])
    return rows


def finish_stage(ledger, digits, stage, stable, started, cpu):
    """Record whole-stage clocks even when bounded calibration fails."""
    drift = analysis.stability(measured_rows(ledger, digits)) if stable else []
    entries = [
        e
        for e in ledger.data["calls"]
        if e["phase"] in (f"warm{digits}", f"{stage}{digits}")
        and e["identity"].get("stage", stage) == stage
    ]
    actual = [
        json.loads((ledger.output / e["receipt"]).read_text()) for e in entries
    ]
    wall, used_cpu = time.monotonic() - started, time.process_time() - cpu
    overhead_wall = max(0, wall - sum(r.get("actual_wall", 0) for r in actual))
    overhead_cpu = max(
        0, used_cpu - sum(r.get("actual_cpu", 0) for r in actual)
    )
    receipt(
        ledger.output / f"{stage}-finished-{digits}.json",
        dict(
            stable_warmup=stable,
            extension_required=bool(drift) and stage == "base",
            drift=drift,
            stage_wall=wall,
            stage_cpu=used_cpu,
            overhead_wall=overhead_wall,
            overhead_cpu=overhead_cpu,
            overhead_within_prepayment=max(overhead_wall, overhead_cpu) <= 60,
        ),
    )
    if max(overhead_wall, overhead_cpu) > 60:
        # This time is already spent. Record its finite excess even if an
        # allowance is exhausted, then stop; a room check cannot erase it.
        identity = dict(stage=stage, digits=digits, reason="overhead-overrun")
        path = ledger.output / f"{stage}-overrun-{digits}.json"
        row = dict(
            index=len(ledger.data["calls"]),
            phase="orchestration",
            receipt=path.name,
            identity=identity,
            status="reserved",
            wall=max(0, overhead_wall - 60),
            cpu=max(0, overhead_cpu - 60),
            work=0,
        )
        ledger.data["calls"].append(row)
        ledger.write()
        wrapper = dict(
            identity=identity,
            reservation=row["index"],
            result=dict(reason="overhead-overrun", work=0),
            charge_wall=row["wall"],
            charge_cpu=row["cpu"],
            charge_work=0,
        )
        receipt(path, wrapper)
        ledger.settle(row, wrapper)
        row.update(
            status="failed",
            error="orchestration prepayment exceeded",
            sha256=digest(path),
        )
        ledger.write()
        raise RuntimeError(
            "orchestration exceeded prepayment; audit before continuing"
        )


def execute(ledger, digits, stage, control_runtime, started, cpu):
    if any(e["status"] != "complete" for e in ledger.data["calls"]):
        raise RuntimeError("failed/interrupted measurement requires audit")
    if any(
        r.get("reason") == "interrupted" for r in measured_rows(ledger, digits)
    ):
        raise RuntimeError(
            "censored interrupted outcome retained; stop for audit"
        )
    marker = ledger.output / f"{stage}-finished-{digits}.json"
    if marker.exists():
        raise ValueError("stage already completed; no duplicate execution")
    if any(
        e["receipt"] == f"{stage}-overhead-{digits}.json"
        for e in ledger.data["calls"]
    ):
        raise RuntimeError(
            "unfinished stage retained; audit before any restart"
        )
    if stage == "extension":
        prior = json.loads(
            (ledger.output / f"base-finished-{digits}.json").read_text()
        )
        if (
            not prior["stable_warmup"]
            or not prior["extension_required"]
            or not prior["overhead_within_prepayment"]
        ):
            raise ValueError("extension is not justified by the frozen rule")
    ledger.call(
        "orchestration",
        f"{stage}-overhead-{digits}",
        60,
        0,
        dict(stage=stage, digits=digits),
        lambda: {"work": 0},
    )
    stable = warmup(ledger, digits, stage, control_runtime)
    if not stable:
        finish_stage(ledger, digits, stage, False, started, cpu)
        return
    corpus = json.loads(inputs.CORPUS.read_text())
    fixtures = [f for f in corpus["fixtures"] if f["digits"] == digits]
    offset = 0 if stage == "base" else 2
    for repeat in range(offset, offset + 2):
        for index, fixture in enumerate(fixtures):
            if not ledger.room(
                f"{stage}{digits}", 2 * (CAPS[digits] + 5), 2 * WORK
            ):
                raise RuntimeError("cannot admit whole matched pair")
            order = ARMS if (index + repeat) % 2 == 0 else ARMS[::-1]
            paired = []
            for arm in order:
                identity = dict(
                    sample=repeat * 9 + index,
                    repeat=repeat,
                    fixture=fixture["id"],
                    n=fixture["n"],
                    digits=digits,
                    seed=7,
                    arm=arm,
                    cap=CAPS[digits],
                )

                def call():
                    row = run_one(fixture, arm, CAPS[digits], control_runtime)
                    row.update(repeat=repeat, sample=identity["sample"])
                    return row

                paired.append(
                    ledger.call(
                        f"{stage}{digits}",
                        f"{stage}-{digits}-{index}-{repeat}-{arm}",
                        CAPS[digits] + 5,
                        WORK,
                        identity,
                        call,
                    )
                )
            if all(r["complete"] for r in paired) and semantic(
                paired[0]
            ) != semantic(paired[1]):
                raise ValueError("complete matched outcomes or paths disagree")
            print(
                f"completed {stage} {digits} input {index + 1}/9 "
                f"repeat {repeat - offset + 1}/2",
                flush=True,
            )
    finish_stage(ledger, digits, stage, True, started, cpu)


def main():
    started, cpu = time.monotonic(), time.process_time()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action",
        choices=("freeze", "verify", "generate", "run", "extend", "analyze"),
    )
    parser.add_argument("--digits", type=int, choices=(50, 55))
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    c9.runtime()
    if args.action == "freeze":
        if CONTROL.exists():
            raise ValueError("control already frozen")
        runtime.checked_baseline()
        inputs.check_warmup(json.loads(WARM.read_text()))
        receipt(CONTROL, policy())
        return
    control = load()
    if args.action == "verify":
        if inputs.CORPUS.exists():
            inputs.check(
                json.loads(inputs.CORPUS.read_text()), digest(CONTROL)
            )
        print("Frozen integration sources/control/corpus verify; no timing.")
        return
    require_committed(CONTROL)
    for name in control["source"]:
        require_committed(ROOT / name)
    with machine_window(), output_lock(args.output):
        ledger = Ledger(args.output, digest(CONTROL), LIMITS)
        receipt(
            args.output / "environment.json",
            dict(
                implementation=platform.python_implementation(),
                version=sys.version,
                executable=sys.executable,
                platform=platform.platform(),
                baseline_commit=runtime.START,
                control_sha256=digest(CONTROL),
            ),
        )
        if args.action == "generate":
            ledger.call(
                "generation",
                "fresh-generation",
                120,
                0,
                dict(binding=digest(CONTROL)),
                lambda: inputs.generate(digest(CONTROL)),
            )
            return
        require_committed(inputs.CORPUS)
        inputs.check(json.loads(inputs.CORPUS.read_text()), digest(CONTROL))
        if args.digits is None:
            raise ValueError("this stage needs --digits 50 or 55")
        if args.action in ("run", "extend"):
            execute(
                ledger,
                args.digits,
                "base" if args.action == "run" else "extension",
                runtime.load_control(args.output),
                started,
                cpu,
            )
        else:
            final = json.loads(
                (args.output / f"base-finished-{args.digits}.json").read_text()
            )
            if not final["overhead_within_prepayment"]:
                raise RuntimeError("orchestration overrun requires audit")
            if final["extension_required"]:
                final = json.loads(
                    (
                        args.output / f"extension-finished-{args.digits}.json"
                    ).read_text()
                )
                if not final["overhead_within_prepayment"]:
                    raise RuntimeError("orchestration overrun requires audit")

            def decide():
                if not final["stable_warmup"]:
                    return dict(
                        cost_gate_passed=False,
                        reasons=["bounded warmup failed; retain default"],
                    )
                return analysis.summarize(
                    measured_rows(ledger, args.digits), stable_warmup=True
                )

            report = ledger.call(
                "analysis",
                f"analysis-{args.digits}",
                300,
                0,
                dict(digits=args.digits),
                decide,
            )
            receipt(args.output / f"decision-{args.digits}.json", report)
            print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
