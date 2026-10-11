"""Finite C9 endpoint/retention attribution with frozen collector controls."""

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
from pathlib import Path

from .... import qs
from ....common import utils
from ....execution.budget import Budget, BudgetExhaustedError
from ...ecm.p52.p52_realistic import check_quiet
from ...infrastructure.performance.performance_audit import fingerprint
from ..b1.b1_calibration import save
from ..c1.c1_feasibility import validate_control
from ..c1.c1_followup import load_followup
from ..c1.c1_implementation import arm_config
from ..phase_three.phase_three_reference import _rss_bytes

ROOT = Path(__file__).resolve().parents[4]
INPUTS = ROOT / "v2/benchmarks/inputs"
PROTOCOL = Path(__file__).with_name("c9_protocol.md")
CONTROL = INPUTS / "controls/c9_training.json"
ARMS = ("slp", "slp_aligned", "dlp", "aligned", "retention")
LOCK = Path("/private/tmp/factor-performance.lock")
OWNER = Path("/private/tmp/factor-performance-owner.json")


def hashes():
    result = fingerprint(ROOT)
    for path in (Path(__file__), PROTOCOL):
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


def configuration(base, arm):
    """Keep domains separate: alignment never narrows DLP candidate scoring."""
    if arm == "slp":
        return base
    if arm == "slp_aligned":
        return replace(
            base,
            collector=replace(
                base.collector,
                residual_bound=100 * base.base_bound,
            ),
        )
    config = arm_config(base, "dlp")
    if arm == "aligned":
        return replace(
            config,
            collector=replace(
                config.collector,
                residual_bound=config.collector.large_prime_bound,
            ),
        )
    if arm == "retention":
        return replace(
            config,
            collector=replace(
                config.collector,
                max_partials=16384,
            ),
        )
    if arm != "dlp":
        raise ValueError("unknown C9 arm")
    return config


def decode(value):
    value = dict(value)
    kind = (
        qs.DoubleLargeSieveConfig
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


def run_one(fixture, seed, config, arm, seconds, *, prefixes=False):
    """Inclusive exact outcome validation; prefixes are attribution only."""
    check_quiet({os.getpid()})
    last_poll = time.monotonic()

    def poll():
        nonlocal last_poll
        now = time.monotonic()
        if now - last_poll >= 1:
            check_quiet({os.getpid()})
            if _rss_bytes() > 2**30:
                raise MemoryError("C9 process RSS cap")
            last_poll = now
        return False

    started, cpu = time.monotonic(), time.process_time()
    budget = Budget(
        work_limit=10**13, seconds=seconds, cpu_seconds=seconds, cancelled=poll
    )
    job = qs.SIQSJob(fixture["n"], seed=seed, config=config, budget=budget)
    observations, tiers = [], deque((60, 120))
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
    reason = result.reason
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
    validate_control(row, fixture)
    row.update(
        seconds=time.monotonic() - started,
        cpu_seconds=time.process_time() - cpu,
        rss_bytes=_rss_bytes(),
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
    _, fixtures, configs = load_followup()
    save(
        CONTROL,
        dict(
            schema=1,
            parent="b3b3cfbea6105f08db0a8484ec09c8260d310e28",
            source=hashes(),
            fixtures=[f for f in fixtures if f["digits"] == 50],
            configurations={
                a: asdict(configuration(configs[50], a)) for a in ARMS
            },
            seeds=[7, 29],
            seconds=240,
            prefixes=[60, 120, 240],
            total_study_seconds=14400,
        ),
    )


def train(output):
    control = json.loads(CONTROL.read_text())
    if control["source"] != hashes():
        raise ValueError("training source differs from freeze")
    save(output / "manifest.json", control)
    start, cpu = time.monotonic(), time.process_time()
    # Training is instrumented attribution, not accepted timing. All calls
    # include setup, outcome validation and prefix inspection in their cap.
    for index, fixture in enumerate(control["fixtures"]):
        order = list(ARMS)
        random.Random(91010 + index).shuffle(order)
        for arm in order:
            if time.monotonic() - start + 241 > 2500:
                raise RuntimeError("training allowance exhausted")
            row = run_one(
                fixture,
                control["seeds"][index],
                decode(control["configurations"][arm]),
                arm,
                control["seconds"],
                prefixes=True,
            )
            save(output / f"run-{index}-{arm}.json", row)
            print(
                json.dumps(
                    {
                        k: row[k]
                        for k in (
                            "id",
                            "arm",
                            "complete",
                            "reason",
                            "seconds",
                        )
                    }
                ),
                flush=True,
            )
    save(
        output / "finished.json",
        dict(
            seconds=time.monotonic() - start,
            cpu_seconds=time.process_time() - cpu,
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("freeze", "train"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    runtime()
    if args.mode == "freeze":
        freeze()
    else:
        if args.output is None:
            parser.error("--output required")
        with machine_window():
            train(args.output)


if __name__ == "__main__":
    main()
