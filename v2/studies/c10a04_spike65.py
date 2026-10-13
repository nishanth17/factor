"""One prospective 65-digit matched pair with finite feasibility grants."""

import argparse
import json
import os
import signal
import time
from pathlib import Path

from ..benchmarks.infrastructure.performance.performance_audit import (
    fingerprint,
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
CONTROL = ROOT / "v2/benchmarks/inputs/controls/c10a04_spike65.json"
CORPUS = ROOT / "v2/benchmarks/inputs/corpora/c10a04_spike65_fresh.json"
OUTPUT = ROOT / "v2/benchmarks/results/c10a04/spike65"
ARMS = ("compact", "control")
CAP = 1800
LIMITS = {
    "total": dict(wall=4500, cpu=4500, work=10**15),
    "setup": dict(wall=800, cpu=800, work=10**14),
    "spike": dict(wall=3700, cpu=3700, work=9 * 10**14),
    "disk_bytes": 256 * 2**20,
}


def sources():
    result = fingerprint(ROOT)
    paths = list((ROOT / "v2/benchmarks").rglob("*.py")) + [
        Path(__file__),
        Path(__file__).with_suffix(".md"),
    ]
    for path in paths:
        if "results" not in path.relative_to(ROOT).parts:
            result[str(path.relative_to(ROOT))] = digest(path)
    return result


def load():
    study.c9.runtime()
    study.require_committed(CONTROL)
    control = json.loads(CONTROL.read_text())
    if control["sources"] != sources() or control["limits"] != LIMITS:
        raise ValueError("spike source identity or allowances changed")
    if (
        control["arms"] != list(ARMS)
        or control["cap_seconds"] != CAP
        or control["configurations"]
        != {arm: study.registry()["65"][arm] for arm in ARMS}
    ):
        raise ValueError("spike arms, configuration or grant changed")
    return control


def expenses():
    ledger = Ledger(OUTPUT, digest(CONTROL), LIMITS)
    if not ledger.data["calls"]:
        ledger.data["calls"].append(
            dict(
                index=0,
                phase="setup",
                status="prepaid",
                receipt="preparation",
                wall=600,
                cpu=600,
                work=5 * 10**13,
            )
        )
        ledger.write()
    return ledger


def generate(control):
    if CORPUS.exists():
        raise RuntimeError("corpus already generated; never regenerate")
    ledger = expenses()

    def execute():
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
        fixtures = selection.products(
            65,
            1,
            GenerationRandom(control["generation_seed"]),
            certificates,
            seen,
            primes,
            "c10a04_spike65_fresh",
        )
        verify_certificates(certificates)
        corpus = dict(
            schema=1,
            fixtures=fixtures,
            certificates=certificates,
            control_sha256=digest(CONTROL),
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


def run(control):
    study.require_committed(CORPUS)
    corpus = json.loads(CORPUS.read_text())
    verify_corpus(corpus)
    if corpus["control_sha256"] != digest(CONTROL):
        raise ValueError("corpus belongs to another plan")
    # A durable one-shot marker rejects duplicate or restarted measurement
    # before any mutable ledger recovery can inspect an accepted call.
    receipt(
        OUTPUT / "process.json", dict(pid=os.getpid(), started=time.time())
    )
    ledger = expenses()
    if any(c.get("error") for c in ledger.data["calls"]):
        raise RuntimeError("prior failure requires audit; no restart")
    deadline = control["deadline_unix"] - 600
    if time.time() + 2 * (CAP + 5) + 80 > deadline:
        raise RuntimeError("whole pair does not fit the absolute deadline")
    fixture = corpus["fixtures"][0]

    def interrupt(signum, frame):
        raise KeyboardInterrupt("spike measurement deadline")

    signal.signal(signal.SIGALRM, interrupt)
    signal.setitimer(signal.ITIMER_REAL, deadline - time.time())
    reason = "finished-one-pair"
    try:
        study.warm(ledger, "spike", 65, fixture, ARMS)
        for arm in ARMS:
            identity = dict(
                control=digest(CONTROL),
                corpus=digest(CORPUS),
                fixture=fixture["id"],
                seed=7,
                arm=arm,
                cap=CAP,
            )
            row = ledger.call(
                "spike",
                arm,
                CAP + 5,
                study.WORK,
                identity,
                lambda arm=arm: study.run_one(fixture, 7, arm, CAP),
            )
            print(
                json.dumps(
                    dict(
                        arm=arm,
                        complete=row["complete"],
                        reason=row["reason"],
                        seconds=row["seconds"],
                    )
                ),
                flush=True,
            )
    except KeyboardInterrupt:
        ledger.recover()
        reason = "stop-interrupted"
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    receipt(
        OUTPUT / "summary.json",
        dict(
            reason=reason,
            interpretation="one-input feasibility only",
            receipts={p.name: digest(p) for p in OUTPUT.glob("*.json")},
            finished=time.time(),
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("generate", "run"))
    args = parser.parse_args()
    control = load()
    with output_lock(OUTPUT):
        if args.action == "generate":
            generate(control)
        else:
            run(control)


if __name__ == "__main__":
    main()
