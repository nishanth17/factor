# 65-digit regression diagnosis, revision 1

The 3 November follow-up is explicitly authorized, including use of the
current model for diagnosis, design and fixes. It retains the original v4
balanced 65-digit pilot: control 1572.906 seconds, candidate 2310.505 seconds,
candidate identical-work CPU drift 2.504 times. No historical receipt or
source freeze is replaced. There is one exposed modulus, not fresh
confirmation. No defaults, integration or promotion is authorized.

The first question is whether unchanged work reproduces the deterioration,
and whether it is concentrated in a policy or affects process-wide execution.
Retained data show fewer candidate positions, work units and solver XORs,
but higher collection, preparation, filtering and elimination times. They
support investigation, not a retrospective causal attribution.

`regression65_v1.json` freezes exact configurations, inputs, current source
hashes and numerical allowances before execution. The owned control's QS,
common, budget and linear-algebra sources are unchanged in committed mainline
as of `99a2943aa9ba681091184007d6dd8aa7171969a2`; it remains explicitly
identified by its original `0076e25` snapshot. Both arms retain seed 7, their
literal policies, the 8500-prime base, 4096-position blocks, 131073-position
polynomial windows and matched Tier 2 ceilings (2 GiB owned, 5 GiB RSS).
Parent suballocations remain 32 MiB owned and 256 MiB RSS. Host admission is
rechecked before each call. Known factors are used only after the call.

First run four calls: control observer off/on, then candidate on/off. Each
sieves exactly 1024 polynomials unless it factors or reaches a native limit;
each has a 90-second wall AND CPU cap plus 5-second margin. On calls pause
every 128 polynomials for bounded GC and JIT-memory snapshots, then resume
the same in-memory job. Off calls have a single uninterrupted prefix. At
least three validated seconds on both clocks and nine identical-work
sentinels precede each measured call; nine sentinels follow it. At most 35
seconds per before/after call is reserved. Observation must preserve work,
atoms and deterministic counters, and its indicative cost must be acceptable
before the longer stage is launched. A single pair is not a speed claim.

The conditional longer stage comprises two calls, candidate then control,
each at most 8192 polynomials, with 1024-polynomial snapshots and a
600-second wall AND CPU cap plus 5-second margin. These cover 128 A
families when full prefixes finish. The longer stage requires an explicit
review and launch, not an automatic fallback after screen failure. Its
expected elapsed cost is about 7–15 minutes for both arms at retained rates;
the operational stage ceiling remains independent of that estimate.

| Reservation | Wall AND CPU seconds |
| --- | ---: |
| Four screen calls including margins | 380 |
| Two conditional longer calls including margins | 1210 |
| Twelve before/after warmup/sentinel calls | 420 |
| Analysis | 120 |
| Orchestration | 60 |
| Total | 2190 |

Work uses the inherited analytic `plan_v4.call_work(65, seconds)`: four times
the retained 8 billion sustained units per CPU second, plus separate bounded
setup, root-table/assignment construction and worst-case postprocessing.
The factor of four is operational headroom, not a throughput theorem.
Warmup uses the original finite 290 billion units per call. The JSON gives
every numerical per-call and cumulative work limit; sums are checked before
launch. Actual settled work and full interrupted reservations remain distinct.
No limit increases or selective retries belong to this protocol.

GC/JIT observations are diagnostics, not clean complete-factor timing.
`gc.get_stats(memory_pressure=False)` avoids a full heap pressure walk;
JIT memory is sampled with `get_stats_asmmemmgr`. `jit-summary` logging is
enabled only for observed calls. The installed PyPy 7.3.23 build aborts
inside `get_stats_snapshot` with a fatal RPython `NotImplementedError`;
the pre-experiment API probe is retained as a development failure and this
API is excluded. The initial hypotheses are host/core scheduling or clock
changes, JIT degradation, GC cost, and increased per-position policy work.
No JIT resets, GC tuning or host scheduling changes are silently applied.

Receipts are bounded to 512 KiB, child pipes to 64 KiB, and this stage to
16 MiB within the existing 512 MiB cumulative evidence ceiling. Recorder
refusal is distinct from factorer limits. Cold imports, warmup, measured
prefixes and post-sentinels remain separate. All outcomes, including native
limits, failed observation checks and unsuccessful hypotheses, are retained.

Only after this diagnosis may a separately frozen correction/ablation stage
or matched complete-factor pair be admitted. The complete pair is authorized
if needed, but is not run by this diagnostic protocol. A justified outcome
may be that no production change is supported by the evidence.

Commands from the task worktree, using PyPy implementing Python 3.11:

```sh
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.regression65 check
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.regression65 screen
# Only after reviewing screen receipts:
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.regression65 long
```

Runtime API references: [PyPy JIT hooks](https://doc.pypy.org/jit-hooks.html)
and [PyPy GC accounting](https://doc.pypy.org/gc_info.html). API availability
was checked against the actual installed build, separately from the current
online documentation.
