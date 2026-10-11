# C9.2 longer-window diagnostic addendum

Declared 11 October 2026 while the frozen 360-second confirmation is running,
following the user's explicit request to allow longer 50-digit runs. This does
not change the existing confirmation source, configurations, outcomes or its
14,400-second aggregate cap. It authorizes one conditional diagnostic; it is
not another parameter search or a fresh performance claim.

## Trigger and selection

Finish or stop the original confirmation under its existing rules. In its
recorded execution order, select the first 50-digit input/seed with a validated
`wall_limit` or `cpu_limit` outcome, scanned positions greater than zero and
work below the existing 10^13 allowance. Ignore warmup and cold diagnostics.
Select by the first censored call in any arm, not by challenger success. If
there is no such call, do not run this follow-up. A short prefix with no
surviving matrix is not grounds for cancelling this pilot.

Reuse exactly the three frozen 50-digit configurations: calibrated SLP,
accepted fixed DLP and selected aligned growing DLP. All receive the same
1,200-second wall/CPU, 10^13-work, 1-GiB owned and 3-GiB process-RSS allowances.
Keep original factor-base, scoring, admission, split, family and dependency
policies. Use seed and input selected above; factors/certificates remain
validation-only. These inputs have already been revealed and this pilot is
explicitly diagnostic, not fresh confirmation.

## Finite allowance and execution

Reserve at most **4,500 additional active wall/CPU seconds** for this
user-requested follow-up, separately from the original study. This covers
three calls reserved at 1,202 seconds each, validated per-arm PyPy warmup,
inspection, validation and assessment. Unpublished/failed calls retain their
full reservation. Do not add inputs, seeds or replacement calls; do not extend
work, memory or a stopped resource limit. Observe the repository performance
lock and obtain a separate coordinated window after the current C9 and queued
C3 QA leases. Waiting and ordinary documentation are excluded as before.

The existing frozen `c9_growth_study.run_one` runner accepts the new duration
without a source change. Call it with `prefixes=True`; it records graph and
matrix progress near 60/180 seconds and at completion or the terminal cap.
The earlier uninstrumented 360-second observation remains a separate capture.
Do not compare instrumented prefix costs as accepted speed measurements.
The graph inspector's extra adjacency scratch and process-RSS limitation are
the same as the original training study and must be disclosed.

Runs start afresh because the confirmation captures are results, not resumable
checkpoints. Charge both the original attempt and every new call to the total
research expenditure; never present a restarted 1,200-second outcome as an
uncharged continuation. Preserve all unsuccessful/censored results. The old
matched comparison retains its own resource scope and acceptance verdict.

## Decision

Report complete factors, proper-divisor/reconstruction validation, elapsed
wall/CPU/work, collection/graph/filter/dependency progress, owned memory, RSS
and all resource stops. Distinguish delayed useful structure from continued
unproductive collection. Even a successful longer pilot establishes only a
feasibility witness on this selected input, not a new speed promotion, a
universal digit cutoff or general 60–100-digit capability. No automatic new
confirmation batch follows. Preserve the original C9.2 decision and leave
broader calibration to C10 and combined acceptance to E1.

## Selected diagnostic input

The rule selects `c9_growth_fresh_50_1`, seed 7: fixed DLP sample 0 is
the first qualifying censor, with 408,681,517 positions, 164,112,078,576
work units, 60,416 evictions and no surviving filtered matrix at 360 seconds.
The growing arm had already completed in 263.777 seconds. The condition is
collection/time censoring, not whether the challenger succeeds. The original
rule was committed in `d6f95fa` before the censored result. Exact selection,
arm order and unchanged resource envelope are recorded in
`inputs/controls/c9_long_window_protocol.json`.

## Reproduction runner

The following runner composes the already frozen, validated study functions.
Save the Python block as an ignored local `.py` file, run it from the repository
root with PyPy implementing Python 3.11, and use a new output directory. Its
source digest is retained in the diagnostic manifest. Fifty seconds of the
4,500-second allowance are prepaid for setup and assessment. It refuses to overwrite
an earlier attempt, reserves unsuccessful calls, validates every warmup and
outcome, and leaves the original confirmation directory untouched. The lease
must be separately coordinated before execution.

```python
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT))

from v2.benchmarks.qs.c9 import c9_growth_confirmation as confirmation
from v2.benchmarks.qs.c9 import c9_growth_study as study

study.runtime()
output = Path(sys.argv[1]).resolve()
protocol_path = study.INPUTS / "controls/c9_long_window_protocol.json"
protocol = json.loads(protocol_path.read_text())
if output.exists():
    raise RuntimeError("Preserve earlier pilot attempts; no replacement run")

with study.machine_window():
    started, cpu_started = time.monotonic(), time.process_time()
    selection, corpus = confirmation.load()
    assert hashlib.sha256(confirmation.SELECTED.read_bytes()).hexdigest() == (
        protocol["selection_sha256"]
    )
    fixture = next(
        case for case in corpus["fixtures"]
        if case["id"] == protocol["fixture"]
    )
    output.mkdir(parents=True)
    study.save(output / "manifest.json", dict(
        protocol=protocol,
        source=selection["source"],
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        runtime=sys.version,
    ))
    expense = confirmation.Expenses(
        output,
        {"spent_before_confirmation": {"wall": 50, "cpu": 50}},
        {"generation_seconds": 0},
    )

    def charged(label, reservation, function):
        if expense.used() + reservation > 4500:
            raise RuntimeError("Long diagnostic aggregate allowance")
        return expense.call(label, reservation, function)

    try:
        for arm in protocol["arms"]:
            config = study.decode(selection["configurations"]["50"][arm])
            warm, wall, cpu = [], 0.0, 0.0
            while min(wall, cpu) < 3:
                if len(warm) == 10:
                    raise RuntimeError("Validated warmup allowance")
                row = charged("warmup-" + arm, 5, lambda: study.run_one(
                    fixture, protocol["seed"], config, arm, 3.1
                ))
                warm.append(row)
                wall += row["seconds"]
                cpu += row["cpu_seconds"]
            study.save(output / ("warmup-" + arm + ".json"), warm)

            # Reserve the full call before starting. The existing runner
            # charges collection, splitting, verification and inspection.
            if time.monotonic() - started + 1202 > 3900:
                raise RuntimeError("Long diagnostic machine lease")
            row = charged(arm, 1202, lambda: study.run_one(
                fixture, protocol["seed"], config, arm, 1200, prefixes=True
            ))
            study.save(output / (arm + ".json"), row)
            print(json.dumps({key: row[key] for key in (
                "arm", "complete", "reason", "seconds", "work"
            )}), flush=True)
    finally:
        study.save(output / "receipt.json", dict(
            wall_seconds=time.monotonic() - started,
            cpu_seconds=time.process_time() - cpu_started,
            charged_seconds=expense.used(),
            interpretation="Selected-case diagnostic; no fresh speed claim",
        ))
```
