# C10a-07 consolidated face and queue comparison

Fresh-panel execution is in progress; nomination is complete. This is the
existing C10a-07 task, with one bounded
45/50/55/60-digit panel replacing execution of the two unstarted60-only screens.
Historical controls and accounts remain unchanged. No performance acceptance,
default change, merge, push or E1 promotion is claimed.

Execution resumed through a separately identified warmup harness repair after
the first12 fresh observations. All15 completed timed calls, including
nomination, are retained without reruns. The old warmup whitelist omitted the
native `stalled_yield` reason; it expected `stalled`. A finite audit showed the
45-digit prefix stops normally after80 blocks in0.120 seconds, rather than
exhausting time or memory. The repair recognizes that valid outcome during
warmup only. Timed native results, factoring code, corpus, order, caps and all
historical freezes remain unchanged.

## Fixed scope and identities

Worktree `/Users/nishanthmohan/Documents/code/factor-c10a07-panel`, branch
`codex/c10a07-panel`. Production source remains mainline e520ec9. New integration
source/protocol/tests: `9757105`; control commit: `8abbf09`. Freeze SHA-256:
`2d2564a0948411b7a3d144ff5178ed0d2d266d7109c00d1f5ab30c3679a1c5a2`.

The existing exact queue and face implementations are reused. A study wrapper
also supports unchanged nonstreaming reference jobs at45. Helper reservations
now fit **inside1 GiB total ownership**, preserving the factory-selected
marking, plan and pacing fields; this differs from the earlier fused-only queue
prefix experiment with additional helper memory. All arms retain proper-factor,
reconstruction, certainty and finite-budget checks. No factoring-library file
is changed.

Only exact PyPy3.11 at
`/opt/homebrew/Cellar/pypy3.11/7.3.23_2/libexec/pypybin/pypy3.11`, SHA-256
`3c45fc49fa01571d19bac1594495d7300c30589c79c0cad2c5bee8d59f9e5ef2`.
The runner acquires `/private/tmp/factor-performance.lock` itself. No agents or
simultaneous performance campaigns.

## Observations and advancement

1. Integration check:18 short128-block calls on existing50/55/60 inputs.
   Each queue/combined arm matched its corresponding control's exact retained
   state and coverage/recovery counters. Reference-mode arithmetic and resume
   are independently tested. This is not additional timing evidence.
2. Nomination: three complete-cost calls on the first exposed60-digit input,
   canonical/control, complement_high and swap_high. Rank the two faces by
   completion, capped wall, CPU and arm name. One exposed input cannot establish
   a gain or block fresh measurement merely because control wins. Commit the
   exact selection before generation.
3. Fresh panel: three new independent inputs at each45/50/55/60,12 total.
   Control/face/queue on every input, plus combined on all nine50/55/60 inputs.
   At45 the queue hook is bypassed, so combined is an exact face alias and its
   duplicate timed call is omitted. Maximum3+36+9=**48 measured calls**, not
  48 plus either old screen. The12 fresh moduli are independent of nomination.

The fresh inputs use the existing independently verified Pocklington-conditioned
generator, two balanced and one moderately uneven product per band. Exact digits,
achieved factor lengths, residue classes, distinct factors and overlap against
2237 earlier moduli/3344 factor entries are checked. This is not uniform RSA
sampling or the separately proposed second sampling family. The corpus and
exact arm order must be committed before fresh measurement.

Nomination selected `swap_high`, committed as `33bea9f` before generation.
All three calls completed exactly on the same exposed input:

| Arm | Wall seconds | CPU seconds |
| --- | ---: | ---: |
| Control | 1370.551 | 1364.678 |
| Complement | 1323.801 | 1323.479 |
| Swap | 1305.662 | 1305.139 |

The swapped face's4.7% wall/4.4% CPU difference is exploratory evidence from
one independent input, not accepted performance. Exact numerical results and
receipt hashes are committed in `c10a07_panel_nomination.json` in the execution
branch. The control preserves historical work195193348085 and scanned
positions1048674112; the historical cold timing is not a matched comparator.

The12 fresh inputs and order are independently validated and committed as
`f53a294`, before any confirmation observations. They have24 distinct prime
factors, two balanced and one moderately uneven case per band. Corpus SHA-256:
`c4e844fff20506fa852baa3a0de6302b040cf3ad069978ba3cba073fbcab40f3`.
Order SHA-256:
`095d7909053ab85a8510e30b8f1e35b44cb3cad083b9cfe47ae0990b17310ee2`.

Each new process requires >=3 validated wall AND CPU warmup seconds, with the
declared last-three stability gate. Warmup, startup, integration checks and
instrumentation are excluded from performance evidence. Measurements include
validation and GC cleanup. Fixed-control warm sentinels check process drift.

Per-call wall AND CPU caps:45=120s,50=300s,55=900s,60=1800s;1e13 work each.
Measured reservation maximum42720s; warm6000s; integration480s; generation180s.
Total49380 wall/CPU reservation seconds and1.72000001e15 work, separate analysis
120s/QA600s. Expected measured4–7 hours, recalibrated from clean observations.
No automatic cap increase or unstable-cell retry. All-arm unresolved first
blocks close the affected band inconclusive without blocking other bands.

Queue and combined primary comparisons use nine fresh active-band inputs;
45 queue cells check fallback. Face comparison uses the12-input panel. Bands
with three inputs remain descriptive. Input-cluster uncertainty, completion,
capped service cost versus complete cost, censoring, CPU and drift remain
separate. Missing prescribed observations or failed drift gates prevent an
acceptance recommendation. No percentage chaining or per-digit confirmed claim.

## Verification and live evidence

Repair source/tests/protocol commit `a46e482`, control commit `32e73eb`;
continuation control SHA-256:
`3ed32cf987d2bb91f533ffc8beacbe83227a6bf5095b9c97a242ea4ffbf21ae8`.
The isolated continuation deducts every parent/audit charge from the original
accounts and completes only unvisited cells. Its committed-only checkout ran
857 tests, nine skipped. The sole initial error was sandbox denial of read-only
`ps` in an existing test; that test passed with the required access. Full lint,
215 benchmark imports and the pinned loaders pass. Real45 warm readiness passes
all three arms with >=3 wall AND CPU seconds and the unchanged stability rule;
the measurement process rewarms separately.

The parent exited before its end-of-process drift checks. Those missing checks
cannot be restored after exit, so its timings remain exploratory. Continuation
sentinel starts are durable and ends are attempted on completion or failure.
Composite estimates cannot pass the original adoption gate using missing parent
drift. This limitation is retained rather than silently relaxing that gate.

First fresh input per active band, wall seconds (one independent input per row,
descriptive only):

| Digits | Control | Face | Queue | Face + queue |
| --- | ---: | ---: | ---: | ---: |
| 50 | 23.394 | 23.692 | 20.598 | 20.371 |
| 55 | 138.561 | 137.484 | 99.227 | 99.824 |
| 60 | 947.671 | 671.343 | 810.273 | 795.553 |

All twelve factorizations reconstruct exactly. Control/queue and face/combined
retained state and coverage counters match. The first60 face timing difference
is much larger than its0.8% reduction in scanned positions; a useful-yield
improvement or causal mechanism is not established. The combination is slower
than face alone on that input, so component percentages cannot be added.

The first45 matched block completes unresolved in every arm with the native
`stalled_yield` stop in under0.2 seconds. Only45 closes inconclusive under the
original futility rule; its other six cells remain unvisited. This is neither
time/memory exhaustion nor evidence about portfolio factoring. The second fresh
60-digit block is running. At this update18 total timed calls are saved.

The committed-only8abbf09 checkout passed **853 PyPy tests,9 skipped**, lint,
214 benchmark imports and frozen-loader verification. Full tests/lint took
95.166 seconds; generated QA stays ignored. All18 integration checks passed.
An independent receipt audit additionally verifies the actual nested collector
counters for all nine pairs; some automatic field comparisons addressed the
top-level rather than nested counter dictionary. Exact-state checks were already
effective. The frozen executable is unchanged. Per-block queue peak gauges sum
in raw aggregates; `result.queue_reserved_peak` is the checked helper peak.

Current measurement output:
`factor-c10a07-panel/v2/benchmarks/results/c10a07_panel/`.
Reissue the same stage command and output after interruption. Saved receipts
settle; missing timed receipts remain censored and charged; only unvisited cells
continue after rewarm. Arithmetic/provenance failure requires audit, not repair
inside this freeze. The handoff/protocol in the panel worktree provide exact
commands. Terminal numerical summary and durable evidence archive will be added
after execution. Existing recovery/source-gap archives remain intact.
