"""Explicitly approved continuation after a pre-measurement overlap failure."""

import json
import os
import signal
import subprocess
import time
from pathlib import Path

from . import c10a04_spike65 as spike

PLAN = spike.CONTROL.with_name("c10a04_spike65_resume.json")


def main():
    control = spike.load()
    spike.study.require_committed(PLAN)
    plan = json.loads(PLAN.read_text())
    if (
        plan["source_sha256"] != spike.digest(Path(__file__))
        or plan["original_control_sha256"] != spike.digest(spike.CONTROL)
        or plan["corpus_sha256"] != spike.digest(spike.CORPUS)
        or plan["deadline_unix"] != control["deadline_unix"]
    ):
        raise ValueError("continuation identity changed")
    spike.study.require_committed(spike.CORPUS)
    spike.study.c9.check_quiet({os.getpid()})
    with spike.output_lock(spike.OUTPUT):
        original = json.loads((spike.OUTPUT / "process.json").read_text())
        live = subprocess.run(
            ["ps", "-p", str(original["pid"]), "-o", "command="],
            capture_output=True,
            text=True,
        ).stdout.strip()
        if live:
            raise RuntimeError("original accepted process is still alive")
        path = spike.OUTPUT / "ledger.json"
        data = json.loads(path.read_text())
        if spike.digest(path) != plan["failed_ledger_sha256"]:
            raise ValueError(
                "failed ledger changed before approved continuation"
            )
        if len(data["calls"]) != 3 or (
            data["calls"][-1]["status"] != "failed"
            or "benchmark/test overlap" not in data["calls"][-1]["error"]
        ):
            raise ValueError(
                "continuation is only for the audited warmup failure"
            )
        if any(
            (spike.OUTPUT / (arm + ".json")).exists() for arm in spike.ARMS
        ):
            raise RuntimeError(
                "a measured cell already exists; never retry it"
            )
        deadline = control["deadline_unix"] - 600
        if time.time() + 2 * (spike.CAP + 5) + 80 > deadline:
            raise RuntimeError("whole pair cannot fit original remaining time")
        spike.receipt(
            spike.OUTPUT / "resume-process.json",
            dict(
                pid=os.getpid(),
                started=time.time(),
                plan_sha256=spike.digest(PLAN),
                failed_ledger_sha256=plan["failed_ledger_sha256"],
            ),
        )
        # Original process is dead; the failed row stays fully charged. A new
        # PID gives warmup distinct receipt IDs, never a selective timed retry.
        ledger = spike.expenses()
        fixture = json.loads(spike.CORPUS.read_text())["fixtures"][0]

        def interrupt(signum, frame):
            raise KeyboardInterrupt("original spike deadline")

        signal.signal(signal.SIGALRM, interrupt)
        signal.setitimer(signal.ITIMER_REAL, deadline - time.time())
        reason = "finished-one-pair-after-approved-warmup-retry"
        try:
            spike.study.warm(ledger, "spike", 65, fixture, spike.ARMS)
            for arm in spike.ARMS:
                row = ledger.call(
                    "spike",
                    arm,
                    spike.CAP + 5,
                    spike.study.WORK,
                    dict(
                        control=spike.digest(spike.CONTROL),
                        corpus=spike.digest(spike.CORPUS),
                        fixture=fixture["id"],
                        seed=7,
                        arm=arm,
                        cap=spike.CAP,
                    ),
                    lambda arm=arm: spike.study.run_one(
                        fixture, 7, arm, spike.CAP
                    ),
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
        spike.receipt(
            spike.OUTPUT / "summary.json",
            dict(
                reason=reason,
                interpretation="one-input feasibility only",
                receipts={
                    p.name: spike.digest(p)
                    for p in spike.OUTPUT.glob("*.json")
                },
                finished=time.time(),
            ),
        )


if __name__ == "__main__":
    main()
