"""Prospective 60-digit comparison with fixed preliminary trend stops."""

import argparse
import json
import os
import platform
import random
import signal
import time
from collections import defaultdict
from pathlib import Path

from ..benchmarks.infrastructure.performance.performance_audit import (
    verify_corpus,
)
from ..benchmarks.qs.c10a04 import selection, study
from ..benchmarks.qs.c10a04.accounting import (
    Ledger,
    digest,
    output_lock,
    receipt,
)
from ..benchmarks.qs.c10a05_plus.selection import (
    GenerationRandom,
    existing_inputs,
)
from ..benchmarks.suites.build_phase_two_corpus import verify_certificates

ROOT = Path(__file__).resolve().parents[2]
CONTROL = ROOT / "v2/benchmarks/inputs/controls/c10a04_short.json"
CORPUS = ROOT / "v2/benchmarks/inputs/corpora/c10a04_short_fresh.json"
OUTPUT = ROOT / "v2/benchmarks/results/c10a04/short"
ARMS = ("control", "compact")
SEEDS = (7, 29, 47)
CAP = 1800
LIMITS = {
    "total": dict(wall=36000, cpu=36000, work=10**15),
    "setup": dict(wall=3600, cpu=3600, work=10**14),
    "confirm": dict(wall=34200, cpu=34200, work=9 * 10**14),
    "disk_bytes": 256 * 2**20,
}


def sources():
    """Add new orchestration identities without altering the old freeze."""
    return {
        str(path.relative_to(ROOT)): digest(path)
        for path in (
            Path(__file__),
            Path(__file__).with_name("__init__.py"),
            Path(__file__).with_suffix(".md"),
        )
    }


def load():
    study.c9.runtime()
    study.load()
    control = json.loads(CONTROL.read_text())
    study.require_committed(CONTROL)
    if (
        control["sources"] != sources()
        or control["original_control_sha256"] != digest(study.CONTROL)
        or control["arms"] != list(ARMS)
        or control["seeds"] != list(SEEDS)
        or control["cap_seconds"] != CAP
        or control["limits"] != LIMITS
        or control["checkpoints"] != [3, 6, 9]
        or control["early_saving_fraction"] != 0.20
        or control["deadline_unix"] != control["start_unix"] + 36000
    ):
        raise ValueError("short-study identity or rules changed")
    return control


def expenses(output, control):
    ledger = Ledger(output, digest(CONTROL), LIMITS)
    if not ledger.data["calls"]:
        # This nonrefundable reserve includes implementation and committed-only
        # QA before generation; later generation has its own charged call.
        prepaid = max(1800, time.time() - control["start_unix"])
        if prepaid > LIMITS["setup"]["wall"] - 125:
            raise RuntimeError("setup allowance exhausted before generation")
        ledger.data["calls"].append(
            dict(
                index=0,
                phase="setup",
                status="prepaid",
                receipt="setup-prepaid",
                wall=prepaid,
                cpu=prepaid,
                work=5 * 10**13,
            )
        )
        ledger.write()
    return ledger


def generate(output):
    control = load()
    if CORPUS.exists():
        raise RuntimeError("fresh corpus already generated; do not regenerate")
    ledger = expenses(output, control)

    def execute():
        # Exclude both previous products and their certified prime factors.
        seen, primes = existing_inputs(), set()

        def inspect(value):
            if isinstance(value, dict):
                for factor in value.get("factors", []):
                    if isinstance(factor, list):
                        factor = factor[0]
                    if type(factor) is int:
                        primes.add(factor)
                for child in value.values():
                    inspect(child)
            elif isinstance(value, list):
                for child in value:
                    inspect(child)

        for path in (ROOT / "v2/benchmarks/inputs/corpora").glob("*.json"):
            inspect(json.loads(path.read_text()))
        certificates = {}
        generator = GenerationRandom(control["generation_seed"])
        fixtures = selection.products(
            60, 3, generator, certificates, seen, primes, "c10a04_short_fresh"
        )
        verify_certificates(certificates)
        corpus = dict(
            schema=1,
            fixtures=fixtures,
            certificates=certificates,
            seed=control["generation_seed"],
            control_sha256=digest(CONTROL),
            sampling=(
                "certified balanced 60-digit products conditioned on a large "
                "p-1 factor; excludes prior versioned products and factors"
            ),
            generation_seconds=time.monotonic() - generator.started,
        )
        verify_corpus(corpus)
        receipt(CORPUS, corpus)
        return dict(work=10**13, corpus_sha256=digest(CORPUS))

    ledger.call(
        "setup",
        "fresh-generation",
        125,
        10**13,
        dict(control=digest(CONTROL), kind="generation"),
        execute,
    )


def schedule(corpus):
    """Cross every input with each seed before any outcome-dependent stop."""
    fixtures = corpus["fixtures"]
    if len(fixtures) != 3 or any(f["digits"] != 60 for f in fixtures):
        raise ValueError("expected three fresh 60-digit inputs")
    result = []
    for sample in range(9):
        order = list(ARMS)
        random.Random(104202613 + sample).shuffle(order)
        result.append(
            (sample, fixtures[sample % 3], SEEDS[sample // 3], order)
        )
    return result


def decision(rows):
    """Return descriptive trends, never a sequential significance claim."""
    groups = defaultdict(dict)
    for row in rows:
        if row["arm"] in groups[row["sample"]]:
            raise ValueError("duplicate cell")
        groups[row["sample"]][row["arm"]] = row
    if any(set(group) != set(ARMS) for group in groups.values()):
        raise ValueError("interim requires complete matched pairs")
    count = len(groups)
    if count not in (3, 6, 9):
        return "continue"

    pairs = [groups[index] for index in sorted(groups)]
    if not any(pair["compact"]["complete"] for pair in pairs):
        return "stop-no-compact-completions"
    if count == 6:
        by_input = defaultdict(list)
        for pair in pairs:
            by_input[pair["control"]["id"]].append(pair)
        if (
            len(by_input) == 3
            and all(len(group) == 2 for group in by_input.values())
            and all(pair["compact"]["complete"] for pair in pairs)
            and all(
                pair["compact"]["capped_seconds"]
                < pair["control"]["capped_seconds"]
                for pair in pairs
            )
            and all(
                sum(p["compact"]["capped_seconds"] for p in group)
                <= 0.8 * sum(p["control"]["capped_seconds"] for p in group)
                for group in by_input.values()
            )
        ):
            return "stop-preliminary-consistent-compact-trend"
    return "finished-nine-pairs" if count == 9 else "continue"


def read_rows(output, corpus):
    """Read and validate published measurements without ledger recovery."""
    expected = {
        (sample, arm): (fixture, seed)
        for sample, fixture, seed, order in schedule(corpus)
        for arm in order
    }
    rows = []
    for path in sorted(output.glob("pair-[0-9][0-9]-*.json")):
        wrapper = json.loads(path.read_text())
        row = wrapper["result"]
        fixture, seed = expected[(row["sample"], row["arm"])]
        study.validate(row, fixture)
        if (
            row["id"] != fixture["id"]
            or row["seed"] != seed
            or row["cap_seconds"] != CAP
            or wrapper["identity"]["control"] != digest(CONTROL)
            or wrapper["identity"]["corpus"] != digest(CORPUS)
        ):
            raise ValueError("measurement schedule or identity mismatch")
        rows.append(row)
    return rows


def summarize(output, corpus, reason):
    rows = read_rows(output, corpus)
    metrics = {}
    for arm in ARMS:
        group = [row for row in rows if row["arm"] == arm]
        metrics[arm] = dict(
            calls=len(group),
            complete=sum(row["complete"] for row in group),
            capped_seconds=sum(row["capped_seconds"] for row in group),
            actual_seconds=sum(row["seconds"] for row in group),
            work=sum(row["work"] for row in group),
        )
    per_input = {}
    for identity in sorted({row["id"] for row in rows}):
        per_input[identity] = {}
        for arm in ARMS:
            group = [
                r for r in rows if r["id"] == identity and r["arm"] == arm
            ]
            per_input[identity][arm] = dict(
                calls=len(group),
                complete=sum(r["complete"] for r in group),
                capped_seconds=sum(r["capped_seconds"] for r in group),
            )
    value = dict(
        per_input=per_input,
        reason=reason,
        interpretation="descriptive prospective validation; no promotion",
        metrics=metrics,
        control_sha256=digest(CONTROL),
        corpus_sha256=digest(CORPUS),
        receipts={
            p.name: digest(p) for p in sorted(output.glob("pair-*.json"))
        },
        utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    )
    receipt(output / "summary.json", value)
    return value


def interrupt_deadline(signum, frame):
    raise KeyboardInterrupt("short study absolute deadline")


def run(output):
    control = load()
    study.require_committed(CORPUS)
    corpus = json.loads(CORPUS.read_text())
    verify_corpus(corpus)
    if corpus["control_sha256"] != digest(CONTROL):
        raise ValueError("fresh corpus belongs to another plan")
    if (output / "summary.json").exists():
        raise RuntimeError("short study already stopped; do not restart")
    ledger = expenses(output, control)
    if any(row.get("error") for row in ledger.data["calls"]):
        raise RuntimeError("interruption/failure needs review; no auto retry")

    # Reserve final reporting time and never begin a pair that cannot fit.
    deadline = control["deadline_unix"] - 600
    if time.time() >= deadline:
        return summarize(output, corpus, "stop-total-budget")
    signal.signal(signal.SIGALRM, interrupt_deadline)
    signal.setitimer(signal.ITIMER_REAL, deadline - time.time())
    reason = "stop-total-budget"
    receipt(
        output / "measurement-process.json",
        dict(
            pid=os.getpid(),
            runtime=platform.python_version(),
            implementation=platform.python_implementation(),
            started_unix=time.time(),
            deadline_unix=deadline,
        ),
    )
    try:
        study.warm(ledger, "confirm", 60, corpus["fixtures"][0], ARMS)
        for sample, fixture, seed, order in schedule(corpus):
            reserve = len(order) * (CAP + 5)
            if time.time() + reserve > deadline or not ledger.room(
                "confirm", reserve, len(order) * study.WORK
            ):
                break
            for arm in order:
                identity = dict(
                    control=digest(CONTROL),
                    corpus=digest(CORPUS),
                    fixture=fixture["id"],
                    n=fixture["n"],
                    digits=60,
                    sample=sample,
                    seed=seed,
                    arm=arm,
                    cap=CAP,
                )

                def execute(arm=arm):
                    return dict(
                        study.run_one(fixture, seed, arm, CAP),
                        sample=sample,
                        measurement_process=os.getpid(),
                    )

                row = ledger.call(
                    "confirm",
                    f"pair-{sample:02d}-{arm}",
                    CAP + 5,
                    study.WORK,
                    identity,
                    execute,
                )
                print(
                    json.dumps(
                        {
                            k: row[k]
                            for k in (
                                "id",
                                "sample",
                                "arm",
                                "complete",
                                "seconds",
                                "reason",
                            )
                        }
                    ),
                    flush=True,
                )
            rows = read_rows(output, corpus)
            reason = decision(rows)
            if sample + 1 in (3, 6, 9):
                receipt(
                    output / f"checkpoint-{sample + 1}.json",
                    dict(
                        decision=reason,
                        pairs=sample + 1,
                        receipts={
                            p.name: digest(p)
                            for p in sorted(output.glob("pair-*.json"))
                        },
                    ),
                )
            if reason != "continue":
                break
        if reason == "continue":
            reason = "stop-total-budget"
    except KeyboardInterrupt:
        # The call has unwound and released its arithmetic state. Recovery now
        # records its full-cost censored cell; never while measurement is live.
        ledger.recover()
        reason = "stop-interrupted"
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    return summarize(output, corpus, reason)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("generate", "run"))
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    with study.machine_window(), output_lock(args.output):
        study.c9.OWNER.write_text(
            json.dumps(
                dict(
                    owner="C10a-04 short 60-digit study",
                    pid=os.getpid(),
                    cwd=str(ROOT),
                    started_utc=time.strftime(
                        "%Y-%m-%dT%H:%M:%SZ", time.gmtime()
                    ),
                )
            )
        )
        if args.action == "generate":
            generate(args.output)
        else:
            print(json.dumps(run(args.output)), flush=True)


if __name__ == "__main__":
    main()
