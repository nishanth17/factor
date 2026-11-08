"""Finite native-only warmup diagnosis; never retries a factoring call."""

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

from v2.benchmarks.infrastructure.resource_control import (
    recording_budget,
)

from ....execution.budget import Budget
from ..c10a04.accounting import digest, output_lock
from . import host_state, pilots_v4, plan_v4

OUTPUT = plan_v4.ROOT / "v2/benchmarks/results/c10a01/canary-probe-v1"
MODULE = "v2.benchmarks.qs.c10a01.canary_probe"
ITERATIONS = 200_000
WORK = 110_000_000
LIMITS = dict(
    probe=dict(wall=40, cpu=40, work=220_000_000),
    analysis=dict(wall=10, cpu=10, work=1_000_000),
    orchestration=dict(wall=10, cpu=10, work=0),
    total=dict(wall=60, cpu=60, work=221_000_000),
    disk_bytes=2**20 * 4,
)


def source_identity():
    return ":".join(
        digest(path)
        for path in (
            Path(__file__),
            Path(host_state.__file__),
            Path(sys.executable).resolve(),
        )
    )


def worker():
    budget = recording_budget(
        Budget, work_limit=WORK, seconds=15, cpu_seconds=15
    )
    before = host_state.power_state()
    if before != dict(source="AC Power", low_power=0):
        raise RuntimeError("native probe requires the declared power mode")
    reference = None

    def once():
        nonlocal reference
        budget.consume(ITERATIONS)
        wall, cpu = time.monotonic(), time.process_time()
        value = hashlib.pbkdf2_hmac(
            "sha256", b"c10a01-host-canary", b"fixed-salt", ITERATIONS
        ).hex()
        sample = dict(
            wall=time.monotonic() - wall, cpu=time.process_time() - cpu
        )
        if reference is not None and reference != value:
            raise ArithmeticError("native output changed")
        reference = value
        return sample

    cold = [once() for _ in range(9)]
    started, cpu = time.monotonic(), time.process_time()
    warmed = 0
    while min(time.monotonic() - started, time.process_time() - cpu) < 3:
        if warmed == 500:
            raise RuntimeError("finite native warmup refused")
        once()
        warmed += 1
    warm_time = dict(
        wall=time.monotonic() - started, cpu=time.process_time() - cpu
    )
    warm = [once() for _ in range(9)]
    return dict(
        cold=cold,
        warm=warm,
        warmup=warm_time,
        warmed_calls=warmed,
        digest=reference,
        work=budget.used,
        before=before,
        after=host_state.power_state(),
        runtime=host_state.runtime_state(),
    )


def main():
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "worker"))
    args = parser.parse_args()
    if args.action == "worker":
        print(json.dumps(worker()))
        return
    with output_lock(OUTPUT):
        ledger = pilots_v4.BoundedLedger(OUTPUT, source_identity(), LIMITS)
        if ledger.data["calls"]:
            raise RuntimeError("probe already attempted")
        for index in range(2):

            def run():
                value = subprocess.run(
                    [sys.executable, "-m", MODULE, "worker"],
                    cwd=plan_v4.ROOT,
                    capture_output=True,
                    text=True,
                    timeout=20,
                )
                if value.returncode or len(value.stdout) > 32 * 1024:
                    raise RuntimeError(
                        "native probe worker failed: " + value.stderr[-512:]
                    )
                return json.loads(value.stdout)

            ledger.call(
                "probe",
                f"probe-{index}",
                20,
                WORK,
                dict(source=source_identity(), index=index),
                run,
            )
        for phase in ("analysis", "orchestration"):
            ledger.call(
                phase,
                phase,
                10,
                LIMITS[phase]["work"],
                dict(source=source_identity()),
                lambda: dict(work=0),
            )
        print(
            json.dumps(
                dict(
                    charged={
                        key: ledger.used(field=key)
                        for key in ("wall", "cpu", "work")
                    }
                )
            )
        )


if __name__ == "__main__":
    main()
