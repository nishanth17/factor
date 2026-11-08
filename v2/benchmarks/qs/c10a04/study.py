"""Coupled SIQS study with finite leases and complete-outcome receipts."""

import argparse
import json
import os
import platform
import random
import subprocess
import sys
import time
from collections import Counter
from contextlib import contextmanager
from dataclasses import asdict
from math import prod
from pathlib import Path

from v2.benchmarks.infrastructure.resource_control import (
    RecordingBudget as Budget,
)

from ....common import utils
from ....execution.budget import BudgetExhaustedError
from ....qs import SIQSJob
from ...infrastructure.performance.performance_audit import (
    fingerprint,
    verify_corpus,
)
from ..c9 import c9_growth_study as c9
from ..c10a05.selection import instability
from .accounting import Ledger, digest, output_lock, receipt
from .configuration import (
    ARMS,
    BANDS,
    COUNTS,
    ENDPOINT_CEILING,
    HELPER_BYTES,
    ConfigurationRefusedError,
    registry,
    resolve,
)

ROOT = Path(__file__).resolve().parents[4]
INPUTS = ROOT / "v2/benchmarks/inputs"
CONTROL = INPUTS / "controls/c10a04_frozen.json"
TRAINING = INPUTS / "corpora/c10a04_training.json"
SELECTED = INPUTS / "controls/c10a04_selected.json"
FRESH = INPUTS / "corpora/c10a04_fresh.json"
OUTPUT = ROOT / "v2/benchmarks/results/c10a04/experiment"
START = "e0dafa1ebce2307f1f2ac3d787af89dd54f5ecfc"
SEEDS = (7, 29, 47)
CAPS = {55: 1200, 60: 7200, 65: 21600}
WORK = 10**13
EXTENSION_CAPS = {d: 2 * c for d, c in CAPS.items()}
LIMITS = {
    "total": dict(wall=2950000, cpu=2950000, work=3 * 10**16),
    "stage1": dict(wall=3600, cpu=3600, work=10**14),
    "train": dict(wall=1630000, cpu=1630000, work=12 * 10**15),
    "confirm": dict(wall=1090000, cpu=1090000, work=8 * 10**15),
    "extension": dict(wall=200000, cpu=200000, work=2 * 10**15),
    "disk_bytes": 256 * 2**20,
}


def hashes():
    """Pin runtime and transitive benchmark helpers for the new study."""
    result = fingerprint(ROOT)
    paths = list((ROOT / "v2/benchmarks").rglob("*.py")) + [
        Path(__file__).with_name("protocol.md"),
        Path(__file__).with_name("research.md"),
        INPUTS / "controls/c10a04_sources.json",
    ]
    for path in paths:
        if "results" not in path.relative_to(ROOT).parts:
            result[str(path.relative_to(ROOT))] = digest(path)
    return result


def freeze():
    if CONTROL.exists():
        raise ValueError("control already frozen")
    corpus = json.loads(TRAINING.read_text())
    verify_corpus(corpus)
    for digits in BANDS:
        fixtures = [f for f in corpus["fixtures"] if f["digits"] == digits]
        if len(fixtures) != 3 or len({f["n"] for f in fixtures}) != 3:
            raise ValueError(
                "each acceptance band needs three independent inputs"
            )
    receipt(
        CONTROL,
        dict(
            schema=1,
            start_commit=START,
            source=hashes(),
            training_sha256=digest(TRAINING),
            configurations=registry(),
            eligible_counts={str(d): c for d, c in COUNTS.items()},
            endpoint_ceiling=ENDPOINT_CEILING,
            seeds=list(SEEDS),
            caps={str(d): c for d, c in CAPS.items()},
            extension_caps={str(d): c for d, c in EXTENSION_CAPS.items()},
            work_per_call=WORK,
            limits=LIMITS,
            arm_order_seed=1042026,
        ),
    )


def require_committed(path):
    """Require the exact committed identity before using fresh inputs."""
    relative = str(path.relative_to(ROOT))
    committed = subprocess.run(
        ["git", "show", f"HEAD:{relative}"],
        cwd=ROOT,
        capture_output=True,
        check=True,
    ).stdout
    if committed != path.read_bytes():
        raise ValueError("commit selection/corpus before proceeding")


def load():
    control = json.loads(CONTROL.read_text())
    if control["source"] != hashes() or control["training_sha256"] != digest(
        TRAINING
    ):
        raise ValueError("study sources or training corpus differ from freeze")
    if (
        control["configurations"] != registry()
        or control["eligible_counts"] != {str(d): c for d, c in COUNTS.items()}
        or control["endpoint_ceiling"] != ENDPOINT_CEILING
        or control["limits"] != LIMITS
        or control["caps"] != {str(d): c for d, c in CAPS.items()}
        or control["extension_caps"]
        != {str(d): c for d, c in EXTENSION_CAPS.items()}
        or control["seeds"] != list(SEEDS)
        or control["work_per_call"] != WORK
        or control["arm_order_seed"] != 1042026
        or control["start_commit"] != START
        or control["schema"] != 1
    ):
        raise ValueError("configuration or resource policy changed")
    return control


def validate(row, fixture):
    """Validate unresolved reconstruction, exact multiplicity and certainty."""
    if prod(row["factors"]) * prod(row["remaining"]) != fixture["n"]:
        raise ValueError("outcome fails exact reconstruction")
    if Counter(row["factors"]) - Counter(fixture["factors"]):
        raise ValueError("unexpected prime factor")
    if row["complete"] != (not row["remaining"]):
        raise ValueError("completion and unresolved cofactors disagree")
    if row["divisor"] is not None and not utils.valid_divisor(
        row["divisor"], fixture["n"]
    ):
        raise ValueError("improper divisor")
    labels = [
        "proven_prime" if p < 3317044064679887385961981 else "probable_prime"
        for p in row["factors"]
    ]
    if row["certainty"] != labels:
        raise ValueError("certainty changed")
    if row["work"] > WORK or row["stats"].get("workspace_bytes", 0) > 2**30:
        raise ValueError("outcome exceeds declared work or owned memory")


def run_one(
    fixture, seed, arm, seconds, *, max_blocks=None, allowed_processes=()
):
    """Include count resolution through final classification and validation."""
    c9.runtime()
    started, cpu = time.monotonic(), time.process_time()
    last_poll = started

    def poll():
        nonlocal last_poll
        now = time.monotonic()
        if now - last_poll >= 1:
            c9.check_quiet({os.getpid()} | set(allowed_processes))
            if c9._rss_bytes() > 3 * 2**30:
                raise MemoryError("study RSS allowance exhausted")
            last_poll = now
        return False

    budget = Budget(
        work_limit=WORK, seconds=seconds, cpu_seconds=seconds, cancelled=poll
    )
    config, divisor, job = None, None, None
    stats, base_info = {}, {}
    remaining, factors, certainty = [fixture["n"]], [], []
    reason = "unstarted"
    try:
        config, divisor = resolve(fixture["n"], fixture["digits"], arm, budget)
        if divisor is None:
            job = SIQSJob(
                fixture["n"], seed=seed, config=config, budget=budget
            )
            result = job.run(max_blocks=max_blocks)
            divisor, reason, stats = (
                result.divisor,
                result.reason,
                result.stats,
            )
            if (divisor or 1) * result.cofactor != fixture["n"]:
                raise ValueError("SIQS cofactor does not reconstruct")
            if job.base is not None:
                base_info = dict(
                    multiplier=job.base.multiplier,
                    discriminant_mod8=job.base.n_prime % 8,
                    count=len(job.base.entries),
                    endpoint=job.base.primes[-1],
                    exclusive_bound=job.base.bound,
                    singular_primes=[
                        p for p in job.base.primes if job.base.n_prime % p == 0
                    ],
                )
        else:
            reason = "configuration_gcd"
        if divisor is not None:
            if not utils.valid_divisor(divisor, fixture["n"]):
                raise ValueError("SIQS returned an improper divisor")
            remaining = [divisor, fixture["n"] // divisor]
            for child in tuple(remaining):
                budget.consume(child.bit_length() ** 2)
                label = utils.classify_prime(child)
                if label is utils.Primality.COMPOSITE:
                    break
                factors.append(child)
                certainty.append(label.value)
                remaining.pop(0)
    except BudgetExhaustedError:
        reason = budget.reason
    except MemoryError:
        reason = "storage_refusal"
    except ConfigurationRefusedError:
        reason = "configuration_refusal"

    row = dict(
        id=fixture["id"],
        digits=fixture["digits"],
        seed=seed,
        arm=arm,
        divisor=divisor,
        factors=factors,
        certainty=certainty,
        remaining=remaining,
        complete=not remaining,
        reason=reason,
        work=budget.used,
        cap_seconds=seconds,
        stats=stats,
        base=base_info,
        configuration=asdict(config) if config is not None else None,
        helper_reserve_bytes=HELPER_BYTES if arm != "control" else 0,
        rss_bytes=c9._rss_bytes(),
    )
    validate(row, fixture)
    row.update(
        seconds=time.monotonic() - started,
        cpu_seconds=time.process_time() - cpu,
    )
    row["capped_seconds"] = (
        row["seconds"] if row["complete"] else max(seconds, row["seconds"])
    )
    return row


def expenses(output):
    ledger = Ledger(output, digest(CONTROL), LIMITS)
    if not ledger.data["calls"]:
        # A conservative, nonrefundable Stage 1 charge includes all research,
        # corpus generation, failed screens, tests and committed-export QA.
        ledger.data["calls"].append(
            dict(
                index=0,
                phase="stage1",
                status="prepaid",
                receipt="stage1-prepaid",
                wall=3600,
                cpu=3600,
                work=10**14,
            )
        )
        ledger.write()
    return ledger


def rows(output, phase, digits):
    return [
        json.loads(p.read_text())["result"]
        for p in sorted(output.glob(f"{phase}-{digits}-[0-9][0-9]-*.json"))
    ]


@contextmanager
def machine_window():
    """Reuse the shared lease and identify this study's owner explicitly."""
    with c9.machine_window():
        c9.OWNER.write_text(
            json.dumps(
                dict(
                    owner="C10a-04",
                    pid=os.getpid(),
                    cwd=str(ROOT),
                    started_utc=time.strftime(
                        "%Y-%m-%dT%H:%M:%SZ", time.gmtime()
                    ),
                )
            )
        )
        yield


def cold(fixture, arm, output, name):
    request = output / (name + ".request.json")
    response = output / (name + ".response.json")
    receipt(request, dict(fixture=fixture, arm=arm))
    if response.exists():
        row = json.loads(response.read_text())
        validate(row, fixture)
        return dict(row, cold_seconds=None, recovered_startup=True)
    started = time.monotonic()
    result = subprocess.run(
        [
            sys.executable,
            "-B",
            "-m",
            __spec__.name,
            "cold-worker",
            "--request",
            str(request),
            "--output",
            str(response),
        ],
        capture_output=True,
        text=True,
        timeout=35,
    )
    if result.returncode:
        raise RuntimeError(result.stderr)
    row = json.loads(response.read_text())
    row.update(cold_seconds=time.monotonic() - started, stderr=result.stderr)
    return row


def warm(ledger, phase, digits, fixture, arms):
    for arm in arms:
        wall = cpu = 0.0
        for index in range(64):
            name = (
                f"{phase}-{digits}-warm-{arm}-{os.getpid()}-"
                f"{index}-{time.time_ns()}"
            )
            row = ledger.call(
                phase,
                name,
                5,
                WORK,
                dict(kind="warm", arm=arm),
                lambda: run_one(fixture, 7, arm, 3.2),
            )
            wall += row["seconds"]
            cpu += row["cpu_seconds"]
            if min(wall, cpu) >= 3:
                break
        else:
            raise RuntimeError("validated warmup allowance exhausted")


def measure(output, phase, digits, lease_seconds):
    control = load()
    corpus_path = TRAINING if phase == "train" else FRESH
    corpus = json.loads(corpus_path.read_text())
    verify_corpus(corpus)
    arms = list(ARMS)
    if phase == "confirm":
        require_committed(SELECTED)
        require_committed(FRESH)
        selected = json.loads(SELECTED.read_text())
        if selected["control_sha256"] != digest(CONTROL) or corpus[
            "selection_sha256"
        ] != digest(SELECTED):
            raise ValueError("confirmation identities disagree")
        chosen = selected["choices"][str(digits)]
        if chosen == "control":
            return
        arms = ["control", chosen]
    fixtures = [f for f in corpus["fixtures"] if f["digits"] == digits]
    if len(fixtures) != 3:
        raise ValueError("band requires three independent fixtures")
    expected = dict(
        control=digest(CONTROL),
        corpus=digest(corpus_path),
        runtime=sys.version,
        platform=platform.platform(),
    )
    receipt(output / f"{phase}-{digits}-manifest.json", expected)
    ledger = expenses(output)
    if (output / f"{phase}-{digits}-finished.json").exists():
        return
    leases = ledger.data.setdefault("leases", {})
    if leases.get(phase, 0) >= 100:
        raise RuntimeError("phase process-lease allowance exhausted")
    leases[phase] = leases.get(phase, 0) + 1
    ledger.write()
    started = time.monotonic()
    # Cold/startup and every process's warmup are charged but are not samples.
    for arm in arms:
        name = f"{phase}-{digits}-cold-{arm}"
        ledger.call(
            phase,
            name,
            40,
            WORK,
            expected,
            lambda arm=arm, name=name: cold(fixtures[0], arm, output, name),
        )
    warm(ledger, phase, digits, fixtures[0], arms)
    cap = control["caps"][str(digits)]
    for target in (9, 18):
        for index in range(target):
            fixture, seed = fixtures[index % 3], SEEDS[(index // 3) % 3]
            order = list(arms)
            random.Random(
                control["arm_order_seed"] + digits * 100 + index
            ).shuffle(order)
            missing = [
                a
                for a in order
                if not (
                    output / f"{phase}-{digits}-{index:02d}-{a}.json"
                ).exists()
            ]
            reserve = len(missing) * (cap + 5)
            if missing and (
                time.monotonic() - started + reserve > lease_seconds
                or not ledger.room(phase, reserve, len(missing) * WORK)
            ):
                receipt(
                    output / f"{phase}-{digits}-yield-{time.time_ns()}.json",
                    dict(
                        next_group=index, reserve=reserve, used=ledger.used()
                    ),
                )
                return
            for arm in order:
                identity = dict(
                    expected,
                    sample=index,
                    digits=digits,
                    fixture=fixture["id"],
                    n=fixture["n"],
                    seed=seed,
                    arm=arm,
                    cap=cap,
                )

                def execute(arm=arm):
                    return dict(
                        run_one(fixture, seed, arm, cap),
                        sample=index,
                        measurement_process=os.getpid(),
                    )

                row = ledger.call(
                    phase,
                    f"{phase}-{digits}-{index:02d}-{arm}",
                    cap + 5,
                    WORK,
                    identity,
                    execute,
                )
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
        # Resume may already have extension cells. Bind each immutable
        # marker only to its own prefix, never to later partial observations.
        prefix = [
            row
            for row in rows(output, phase, digits)
            if row["sample"] < target
        ]
        reasons = instability(prefix)
        receipt(
            output / f"{phase}-{digits}-stability-{target}.json",
            dict(samples=target, reasons=reasons),
        )
        if not reasons:
            break
    receipt(
        output / f"{phase}-{digits}-finished.json",
        dict(samples=target, reasons=reasons),
    )


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument(
        "action",
        choices=(
            "freeze",
            "screen",
            "grants",
            "train",
            "select",
            "generate",
            "confirm",
            "summarize",
            "cold-worker",
        ),
    )
    parser.add_argument("--digits", type=int, choices=BANDS, default=55)
    parser.add_argument("--lease-seconds", type=int, default=90000)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--request", type=Path)
    args = parser.parse_args()
    c9.runtime()
    if args.action == "cold-worker":
        request = json.loads(args.request.read_text())
        receipt(
            args.output,
            run_one(
                request["fixture"],
                7,
                request["arm"],
                3.2,
                allowed_processes=(os.getppid(),),
            ),
        )
    elif args.action == "freeze":
        freeze()
    elif args.action in ("train", "confirm", "screen", "grants"):
        with machine_window(), output_lock(args.output):
            if args.action in ("screen", "grants"):
                from . import selection

                if args.action == "screen":
                    selection.screen(args.output)
                else:
                    selection.grants(args.output, args.digits)
            else:
                measure(
                    args.output, args.action, args.digits, args.lease_seconds
                )
    else:
        from . import selection

        with output_lock(args.output):
            getattr(selection, args.action)(args.output)


if __name__ == "__main__":
    main()
