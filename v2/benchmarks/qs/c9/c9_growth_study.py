"""C9.2 allocator-grown storage with frozen collectors and charged studies."""

import argparse
import fcntl
import hashlib
import json
import os
import platform
import random
import time
from collections import Counter, deque
from contextlib import contextmanager
from dataclasses import asdict, replace
from math import prod
from pathlib import Path

from .... import qs
from ....common import utils
from ....execution.budget import Budget, BudgetExhaustedError
from ...ecm.p52.p52_realistic import check_quiet
from ...infrastructure.performance.performance_audit import fingerprint
from ..b1.b1_calibration import save
from ..c1.c1_implementation import arm_config
from ..phase_three.phase_three_reference import _rss_bytes

ROOT = Path(__file__).resolve().parents[4]
INPUTS = ROOT / "v2/benchmarks/inputs"
PROTOCOL = Path(__file__).with_name("c9_growth_protocol_v2.md")
CONTROL = INPUTS / "controls/c9_growth_training_v2.json"
JOINT = INPUTS / "controls/c9_growth_joint.json"
PRIOR = {"wall": 4634.658706497634, "cpu": 4588.290882917167}
MEMORY = 1024 * 2**20
ARMS = ("slp", "dlp", "grow", "grow_aligned")
LOCK = Path("/private/tmp/factor-performance.lock")
OWNER = Path("/private/tmp/factor-performance-owner.json")


def hashes():
    result = fingerprint(ROOT)
    for path in (PROTOCOL, *(ROOT / "v2/benchmarks").rglob("*.py")):
        if "results" in path.relative_to(ROOT).parts:
            continue
        result[str(path.relative_to(ROOT))] = hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
    return result


def runtime():
    if platform.python_implementation() != "PyPy" or (
        platform.python_version_tuple()[:2] != ("3", "11")
    ):
        raise RuntimeError("C9 requires PyPy implementing Python 3.11")


@contextmanager
def machine_window():
    """Own the shared machine lease, including heavy correctness checks."""
    with LOCK.open("a+") as lease:
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        check_quiet({os.getpid()})
        OWNER.write_text(
            json.dumps(
                dict(
                    owner="C9",
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
            check_quiet({os.getpid()})
        finally:
            OWNER.unlink(missing_ok=True)


def growing(config):
    """Remove count quotas; retain domains and explicit experiment memory."""
    options = asdict(config.collector)
    options.update(
        max_partials=None,
        max_atoms=None,
        max_relations=None,
        memory_bytes=None,
        split_call_limit=None,
    )
    return replace(
        config,
        memory_bytes=MEMORY,
        collector=qs.GrowingDoubleLargeSieveConfig(**options),
    )


def configuration(registry, arm):
    if arm in ("slp", "dlp"):
        config = decode(registry["slp_aligned" if arm == "slp" else "dlp"])
        return replace(
            config,
            memory_bytes=MEMORY,
            collector=replace(config.collector, memory_bytes=MEMORY),
        )
    config = growing(decode(registry["dlp"]))
    if arm == "grow_aligned":
        config = replace(
            config,
            collector=replace(
                config.collector,
                residual_bound=config.collector.large_prime_bound,
            ),
        )
    elif arm != "grow":
        raise ValueError("unknown growth arm")
    return config


def decode(value):
    value = dict(value)
    kind = (
        qs.GrowingDoubleLargeSieveConfig
        if value["collector"].get("storage_policy") == "grow"
        else qs.DoubleLargeSieveConfig
        if ("large_prime_bound" in value["collector"])
        else qs.SieveConfig
    )
    value["collector"] = kind(**value["collector"])
    return qs.SIQSConfig(**value)


def graph_diagnostics(collector):
    """Peel retained atom incidence; never confuse forest edges with cycles."""
    adjacency, edges, large = {}, [], Counter()
    for atom in collector._atoms.values():
        if atom.large_primes:
            left, right = atom.large_primes
        else:
            left, right = 1, atom.residual
        index = len(edges)
        edges.append((left, right))
        adjacency.setdefault(left, []).append(index)
        adjacency.setdefault(right, []).append(index)
        large.update(v for v in (left, right) if v != 1)
    seen, components, largest = set(), 0, 0
    for vertex in adjacency:
        if vertex in seen:
            continue
        stack, count = [vertex], 0
        seen.add(vertex)
        while stack:
            v = stack.pop()
            count += 1
            for index in adjacency[v]:
                a, b = edges[index]
                other = b if a == v else a
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        components += 1
        largest = max(largest, count)
    degree = {v: len(es) for v, es in adjacency.items()}
    queue = deque(v for v, d in degree.items() if d <= 1)
    removed = set()
    while queue:
        v = queue.popleft()
        for index in adjacency[v]:
            if index in removed:
                continue
            removed.add(index)
            for endpoint in edges[index]:
                degree[endpoint] -= 1
                if degree[endpoint] == 1:
                    queue.append(endpoint)
    return dict(
        retained_atoms=len(edges),
        distinct_endpoints=len(large),
        degree_histogram=dict(Counter(large.values())),
        components=components,
        largest_component=largest,
        retained_cycle_rank=len(edges) - len(adjacency) + components,
        core_edges=len(edges) - len(removed),
        core_vertices=sum(d > 1 for d in degree.values()),
        slp_above_dlp_endpoint=sum(
            a.residual > getattr(collector.config, "large_prime_bound", 10**12)
            for a in collector._atoms.values()
        ),
    )


def validate_result(row, fixture, memory_bytes):
    """Check reconstruction/certainty independently under this run's cap."""
    assert prod(row["factors"]) * prod(row["remaining"]) == fixture["n"]
    expected = Counter(fixture["factors"])
    assert not Counter(row["factors"]) - expected
    assert row["complete"] == (not row["remaining"])
    if row["complete"]:
        assert Counter(row["factors"]) == expected
    divisor = row["divisor"]
    if divisor is not None:
        assert 1 < divisor < fixture["n"] and fixture["n"] % divisor == 0
    assert len(row["certainty"]) == len(row["factors"])
    for factor, label in zip(row["factors"], row["certainty"]):
        wanted = (
            "proven_prime"
            if factor < 3317044064679887385961981
            else "probable_prime"
        )
        assert label == wanted
    assert row["work"] <= 10**13
    assert row["stats"].get("workspace_bytes", 0) <= memory_bytes


def run_one(fixture, seed, config, arm, seconds, *, prefixes=False):
    """Inclusive exact outcome validation; prefixes are attribution only."""
    check_quiet({os.getpid()})
    last_poll = time.monotonic()
    rss_refused = False

    def poll():
        nonlocal last_poll, rss_refused
        now = time.monotonic()
        if now - last_poll >= 1:
            check_quiet({os.getpid()})
            if _rss_bytes() > 3 * 2**30:
                rss_refused = True
                raise MemoryError("C9 process RSS cap")
            last_poll = now
        return False

    started, cpu = time.monotonic(), time.process_time()
    budget = Budget(
        work_limit=10**13, seconds=seconds, cpu_seconds=seconds, cancelled=poll
    )
    job = qs.SIQSJob(fixture["n"], seed=seed, config=config, budget=budget)
    observations, tiers = [], deque((60, 180))
    if prefixes:
        while True:
            result = job.run(max_blocks=256)
            if tiers and budget.wall_used >= tiers[0]:
                tier = tiers.popleft()
                observations.append(
                    dict(
                        tier=tier,
                        seconds=budget.wall_used,
                        stats=json.loads(json.dumps(result.stats)),
                        graph=graph_diagnostics(job.engine.collector),
                    )
                )
            if result.reason != "paused":
                break
    else:
        result = job.run()
    if (result.divisor or 1) * result.cofactor != fixture["n"]:
        raise AssertionError("result does not reconstruct")
    if result.divisor is not None and not utils.valid_divisor(
        result.divisor, fixture["n"]
    ):
        raise AssertionError("improper divisor")
    factors, labels, remaining = [], [], [fixture["n"]]
    reason = "rss_limit" if rss_refused else result.reason
    if result.divisor is not None:
        remaining = [result.divisor, result.cofactor]
        try:
            for child in tuple(remaining):
                budget.consume(child.bit_length() ** 2)
                label = utils.classify_prime(child)
                if label is utils.Primality.COMPOSITE:
                    break
                factors.append(child)
                labels.append(label.value)
                remaining.pop(0)
            budget.consume(0)
        except BudgetExhaustedError:
            reason = "classification_" + budget.reason
    row = dict(
        id=fixture["id"],
        digits=fixture["digits"],
        seed=seed,
        arm=arm,
        complete=not remaining,
        divisor=result.divisor,
        factors=factors,
        remaining=remaining,
        certainty=labels,
        reason=reason,
        work=budget.used,
        cap_seconds=seconds,
        stats=result.stats,
    )
    row.update(
        seconds=time.monotonic() - started,
        cpu_seconds=time.process_time() - cpu,
        rss_bytes=_rss_bytes(),
    )
    try:
        validate_result(row, fixture, config.memory_limit)
    except AssertionError:
        # Preserve rejected reports as failures, never usable timing evidence.
        failed = (
            ROOT
            / "v2/benchmarks/results/c9"
            / f"growth-validation-failed-{time.time_ns()}.json"
        )
        save(
            failed,
            dict(
                outcome=row,
                fixture=fixture["id"],
                memory_bytes=config.memory_limit,
            ),
        )
        raise
    # Inclusive timing ends only after independent output validation.
    row.update(
        seconds=time.monotonic() - started,
        cpu_seconds=time.process_time() - cpu,
    )
    row["capped_seconds"] = row["seconds"] if row["complete"] else seconds
    if prefixes:
        observations.append(
            dict(
                tier="terminal",
                seconds=budget.wall_used,
                stats=result.stats,
                graph=graph_diagnostics(job.engine.collector)
                if job.engine
                else {},
            )
        )
        row["observations"] = observations
        row["diagnostic_seconds"] = time.monotonic() - started
    return row


def freeze():
    if CONTROL.exists():
        raise ValueError("growth training already frozen")
    initial = json.loads((INPUTS / "controls/c9_training.json").read_text())
    save(
        CONTROL,
        dict(
            schema=2,
            source=hashes(),
            fixtures=initial["fixtures"],
            configurations={
                a: asdict(configuration(initial["configurations"], a))
                for a in ARMS
            },
            seeds=[7, 29],
            seconds=360,
            prefixes=[60, 180, 360],
            prior=PRIOR,
            total_study_seconds=14400,
            memory_bytes=MEMORY,
            checkpoint_bytes=2**20,
            rss_bytes=3 * 2**30,
            scope="Growth at fixed full-base C1 collector; alignment is "
            "a separately labelled admission change. All arms share RAM/time.",
        ),
    )


def expenses(output):
    from .c9_growth_confirmation import Expenses

    return Expenses(
        output, {"spent_before_confirmation": PRIOR}, {"generation_seconds": 0}
    )


def train(output, joint=False, lease_seconds=600):
    control = json.loads((JOINT if joint else CONTROL).read_text())
    if control["source"] != hashes():
        raise ValueError("growth source differs from freeze")
    manifest = output / ("joint-manifest.json" if joint else "manifest.json")
    if manifest.exists() and json.loads(manifest.read_text()) != control:
        raise ValueError("output has a different frozen manifest")
    save(manifest, control)
    expense = expenses(output)
    started, cpu = time.monotonic(), time.process_time()
    for index, fixture in enumerate(control["fixtures"]):
        arms = list(control["configurations"])
        random.Random(91020 + index + 20 * joint).shuffle(arms)
        for arm in arms:
            path = output / f"run-{index}-{arm}.json"
            if path.exists():
                continue
            if (
                max(time.monotonic() - started, time.process_time() - cpu)
                + 362
                > lease_seconds
            ):
                print(
                    "C9 lease yielded before the next reserved call",
                    flush=True,
                )
                return
            row = expense.call(
                path.stem,
                362,
                lambda: run_one(
                    fixture,
                    control["seeds"][index],
                    decode(control["configurations"][arm]),
                    arm,
                    360,
                    prefixes=True,
                ),
            )
            save(path, row)
            print(
                json.dumps(
                    {
                        k: row[k]
                        for k in ("id", "arm", "complete", "reason", "seconds")
                    }
                ),
                flush=True,
            )
    save(
        output / ("joint-finished.json" if joint else "finished.json"),
        dict(total_charged=expense.used()),
    )


def rank(rows, arms):
    from .c9_selection import ranked

    return ranked(rows, arms)


def joint_freeze(output):
    if JOINT.exists():
        raise ValueError("growth combinations already frozen")
    initial = json.loads(CONTROL.read_text())
    if initial["source"] != hashes():
        raise ValueError("initial source changed")
    rows = [
        json.loads(p.read_text()) for p in sorted(output.glob("run-*.json"))
    ]
    if len(rows) != 8 or not any(
        r["complete"]
        or r["stats"]
        .get("engine", {})
        .get("last_filter", {})
        .get("output_rows", 0)
        for r in rows
    ):
        raise ValueError(
            "growth advancement requires useful long-prefix structure"
        )
    chosen = rank(rows, ("grow", "grow_aligned"))[0]["arm"]
    old = json.loads((INPUTS / "controls/c9_training.json").read_text())
    half = growing(
        arm_config(decode(old["configurations"]["slp"]), "dlp_half")
    )
    if chosen == "grow_aligned":
        half = replace(
            half,
            collector=replace(
                half.collector, residual_bound=half.collector.large_prime_bound
            ),
        )
    slp = decode(initial["configurations"]["slp"])
    graph = growing(
        replace(
            slp,
            collector=qs.DoubleLargeSieveConfig(
                **asdict(slp.collector),
                large_prime_bound=2,
                large_product_bound=4,
            ),
        )
    )
    save(
        JOINT,
        dict(
            schema=2,
            source=hashes(),
            fixtures=initial["fixtures"],
            seeds=[7, 29],
            configurations={
                "half_grow": asdict(half),
                "graph_slp": asdict(graph),
            },
            selected_endpoint=chosen,
            seconds=360,
            rationale="C1 half-base bundle with the selected growth "
            "admission; reusable SLP uses the existing C1 graph_slp route. "
            "These are combinations, not isolated storage effects.",
        ),
    )


def select(output):
    from . import c9_growth_confirmation as confirmation

    if confirmation.SELECTED.exists() or confirmation.CORPUS.exists():
        raise ValueError("selection already frozen or fresh inputs revealed")
    initial, joint = (json.loads(p.read_text()) for p in (CONTROL, JOINT))
    if initial["source"] != hashes() or joint["source"] != hashes():
        raise ValueError("source changed after training freeze")
    paths = sorted(output.glob("run-*.json"))
    rows = [json.loads(p.read_text()) for p in paths]
    if len(rows) != 12:
        raise ValueError(
            "selection needs all eight individual and four joint calls"
        )
    slp_rank = rank(rows, ("slp", "graph_slp"))
    challenger_rank = rank(
        rows, ("grow", "grow_aligned", "half_grow", "graph_slp")
    )
    registry = initial["configurations"] | joint["configurations"]
    old = json.loads((INPUTS / "controls/c1_selected.json").read_text())[
        "configurations"
    ]
    configs = {}
    for digits in (30, 40):
        slp = replace(decode(old[str(digits)]["slp"]), memory_bytes=MEMORY)
        dlp = replace(
            decode(old[str(digits)]["dlp_half"]), memory_bytes=MEMORY
        )
        configs[str(digits)] = dict(
            slp=asdict(slp), dlp=asdict(dlp), challenger=asdict(growing(dlp))
        )
    configs["50"] = dict(
        slp=registry["slp"],
        dlp=registry["dlp"],
        challenger=registry[challenger_rank[0]["arm"]],
    )
    ledger = expenses(output).data
    prior = {
        kind: ledger["prior_" + kind] + sum(r[kind] for r in ledger["calls"])
        for kind in ("wall", "cpu")
    }
    save(
        confirmation.SELECTED,
        dict(
            schema=2,
            source=confirmation.hashes(),
            configurations=configs,
            choices={
                "slp_diagnostic": slp_rank,
                "challenger": challenger_rank,
            },
            spent_before_confirmation=prior,
            evidence=[
                dict(
                    path=str(p.relative_to(ROOT)),
                    sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                    outcome={
                        k: r[k]
                        for k in (
                            "id",
                            "arm",
                            "complete",
                            "reason",
                            "seconds",
                            "capped_seconds",
                        )
                    },
                )
                for p, r in zip(paths, rows)
            ],
            scope="Growth at 30/40-digit C1 DLP bundles; selected growing DLP "
            "challenger at 50. User-requested storage policy is separate from "
            "performance promotion. E1 and C10 remain open.",
            expense_convention="Original 14400 seconds; prior includes "
            "C9.1, interrupted-call reserve and 600 prepaid QA.",
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode", choices=("freeze", "train", "joint-freeze", "joint", "select")
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--lease-seconds", type=int, default=600)
    args = parser.parse_args()
    runtime()
    if args.mode == "freeze":
        freeze()
    elif args.mode in ("joint-freeze", "select"):
        if args.output is None:
            parser.error("--output required")
        (joint_freeze if args.mode == "joint-freeze" else select)(args.output)
    else:
        if args.output is None:
            parser.error("--output required")
        with machine_window():
            train(
                args.output,
                joint=args.mode == "joint",
                lease_seconds=args.lease_seconds,
            )


if __name__ == "__main__":
    main()
