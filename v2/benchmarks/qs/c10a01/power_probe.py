"""Two bounded, source-pinned probes around an authorized power-mode change."""

import argparse
import json
import signal
import sys
from pathlib import Path

from ..c10a04.accounting import digest, output_lock
from . import execute_v4, host_state, pilots_v4, plan_v4, regression65

PLAN = Path(__file__).with_name("power_probe.json")
OUTPUT = plan_v4.ROOT / "v2/benchmarks/results/c10a01/power-probe"


def freeze():
    if PLAN.exists():
        raise RuntimeError("power probe already frozen")
    parent = regression65.check()
    row = next(r for r in parent["rows"] if r["arm"] == "control")
    sources = dict(parent["sources"])
    for path in (Path(__file__), Path(host_state.__file__)):
        sources[str(path.relative_to(plan_v4.ROOT))] = digest(path)
    limits = dict(
        warmup=dict(wall=70, cpu=70, work=580_000_000_000),
        analysis=dict(wall=5, cpu=5, work=1_000_000_000),
        orchestration=dict(wall=5, cpu=5, work=0),
        total=dict(wall=80, cpu=80, work=581_000_000_000),
        disk_bytes=8 * 2**20,
    )
    value = dict(
        identity="c10a01-power-probe-v1",
        sources=sources,
        runtime=parent["runtime"],
        row=row,
        limits=limits,
        call_seconds=35,
        call_work=290_000_000_000,
        authorization="User authorized temporary AC low-power disable "
        "and restoration on 3 November 2026.",
    )
    PLAN.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def check():
    value = json.loads(PLAN.read_text())
    for name, expected in value["sources"].items():
        if digest(plan_v4.ROOT / name) != expected:
            raise RuntimeError("power probe source drift: " + name)
    if (
        sys.implementation.name != "pypy"
        or digest(Path(sys.executable).resolve()) != value["runtime"]
    ):
        raise RuntimeError("PyPy binary drift")
    return value


def probe(row, label):
    expected = dict(source="AC Power", low_power=int(label == "low"))
    before = host_state.power_state()
    if before != expected:
        raise RuntimeError("probe requires declared active AC power mode")
    modules = execute_v4.runtime_modules(row)
    fixture = execute_v4.identity(row, json.loads(plan_v4.CORPUS.read_text()))
    result = pilots_v4.sentinel(row, fixture, modules)
    canary = host_state.native_canary()
    after = host_state.power_state()
    return dict(
        siqs=result,
        native=canary,
        before=before,
        after=after,
        runtime=host_state.runtime_state(),
        stable=before == after,
        work=result["work"] + canary["work"],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("freeze", "check", "low", "normal"))
    args = parser.parse_args()
    if args.action == "freeze":
        freeze()
        return
    value = check()
    if args.action == "check":
        print(json.dumps(value["limits"]))
        return
    with output_lock(OUTPUT):
        ledger = pilots_v4.BoundedLedger(OUTPUT, digest(PLAN), value["limits"])
        if any(
            c["receipt"] == args.action + ".json" for c in ledger.data["calls"]
        ):
            raise RuntimeError("probe already attempted; no implicit retry")
        execute_v4.OUTPUT = OUTPUT
        pilots_v4.host_admission(2)
        signal.signal(
            signal.SIGALRM,
            lambda *_: (_ for _ in ()).throw(TimeoutError("probe guard")),
        )
        signal.alarm(35)
        result = ledger.call(
            "warmup",
            args.action,
            35,
            290_000_000_000,
            dict(label=args.action, protocol=digest(PLAN)),
            lambda: probe(value["row"], args.action),
        )
        signal.alarm(0)
        print(json.dumps(result))


if __name__ == "__main__":
    main()
