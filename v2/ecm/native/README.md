# Optional ordinary-GMP ECM packets and pthread slices

Python integers on PyPy3.11 remain the default. The isolated B4a branch admits
whole-classical scalar/pthread correctness on `fa2a4cc`: full QA runs 1,741 tests
with 20 explicit optional skips, lint and 308 ordinary benchmark imports. All
20 new whole-native test methods execute. The earlier `c-gmp` stage-one owner
remains an independent control; select `c-gmp-full` explicitly for the complete
C owner described below. Performance adoption, recalibration and integration
remain open after the cold/warm Tier 1 screen.

The Phase-2 experimental whole-ECM driver observes a childless macOS root
through one CFFI bridge using the same libproc usage/child queries as the
Python resource provider. Poll locations, cancellation, wall/CPU/RSS decisions
and publication ordering are retained. A live child, multiple roots or a
custom provider uses the original provider. Observations are never replaced
with a cancellation flag or an estimated memory quota. The provider is restored
after each poll, including errors, and its buffer is released before unload.
Committed-only QA on `ef061f3` passes 1,748 tests (20 optional skips), lint and
308 imports. Its first complete-cost screen establishes no improvement over
C1 and stops with three final A/A cells unrun. The child-list sizing query can
force fallback for childless roots; the actual-ID-query correction is a draft
awaiting native build/QA and remeasurement. No default or worker recipe changes.

Build explicitly from the repository root under the shared performance lock:

```sh
pypy3 -B -m v2.ecm.native.build --gmp-prefix /opt/homebrew/opt/gmp
pypy3 -B -m v2.ecm.native.build --gmp-prefix /opt/homebrew/opt/gmp --pthread
```

The explicit build also prepares an out-of-line CFFI ABI descriptor in the
ignored build directory. Its bounded bytes and filename are hash-checked
against the manifest before loading. Source, binary and GMP checks still run;
the first native call avoids parsing C declarations. Older manifests without
a descriptor retain inline binding. No compilation or preparation occurs on
default/ineligible calls. Declared missing/corrupt descriptors raise instead
of silently changing backend.

Both retain `-std=c11 -O2 -Wall -Wextra -Werror`. The scalar artifact is in
`build/`; pthread is separately in `build/pthread/`, with `-DB4_PTHREAD -pthread`.
Pthreads use libSystem on macOS; no OpenMP runtime or new dependency is needed.
The retained `--openmp-prefix` build is a separate legacy C capability, not the
current Python team owner and not validated by this continuation. A pthread
binary reports pthread capability and zero OpenMP capability. Explicit missing
backends raise named errors; the team never falls back to scalar or Python.
Factoring startup does not compile or install anything.

Manifests record compiler/version/flags/target, PyPy/CFFI, C sources including
`team.c`, GMP header/runtime and resolved library hashes, actual dylib linkage,
backend and binary hash. Canonical packet ABI1 is retained with additive team
entry points; source/capability validation refuses stale or mismatched builds.
The continuation's local `BUILD_PROVENANCE.json` also retains the compiler hash.
No binary is vendored.

Scalar selection is explicit, using `PortfolioConfig(ecm_kernel="c-gmp")` or
CLI `--ecm-kernel c-gmp`. `auto` may fall back only for unavailable dependencies
or ineligible widths; corruption/identity mismatch raises. Cheap exits do not
load the library. Scalar packets also work inside B12 spawned processes.
Python-level gmpy2 operators are a separate backend.

The packet accepts odd 3..1024-bit n, canonical unsigned big-endian values,
scalars below `2**4096`, and up to 16 verified records with at most 512 operations
and 16 point registers each. It preserves the fused Montgomery ladder, accepted
PRAC programs/eligibility, frontier order and exact Python recovery. Full n,
not the unknown factor size, determines arithmetic width. Public/checkpoint
values remain canonical. No interpreter-private integer layout, custom GMP
allocator, external factorer, new reducer or SIMD arithmetic is used.

GMP contexts, scratch, input/output transport and immutable verified programs
are reused. Their memory estimates are diagnostics, not quotas. H3 root wall,
aggregate CPU and observed current owned-process-tree RSS govern resources.
Process CPU already includes native threads; no thread CPU is added again.
GMP allocation failure may terminate the process; the last saved checkpoint
is its recovery boundary, retaining saved history rather than unsaved work.

## Explicit pthread ownership

```python
from v2.execution.budget import Budget
from v2.execution.ecm_native import NativeECMExecutor
from v2.execution.worker_contract import WorkerSettings

owner = NativeECMExecutor(
    seed=7,
    settings=WorkerSettings(workers=2),
    schedule="dynamic",
    chunk_size=1,
)
try:
    run = owner.run(
        NUMBER,
        config=CONFIG_WITH_EXPLICIT_C_GMP,
        budget=Budget(seconds=30, cpu_seconds=120),
    )
finally:
    owner.close()
```

The scheduled unit is **one independently owned continuation slice**, not
concurrent operations on one dependent curve point. Immutable program data may
be shared. The C entry checks distinct contexts; every job has private mutable
GMP state. Dynamic scheduling uses a bounded shared atomic chunk index; static
scheduling assigns cyclic chunks to fixed thread ranks. Supported workers are
1/2/4, `schedule="dynamic"/"static"`, and `chunk_size=1..32`. The ready capacity
is 1..32, default `2*workers`. Correctness covers chunks 1/3, uneven tails,
fewer jobs than threads, empty queues and repeated dispatch. No performance
choice of worker count, schedule or chunk is implied.

Contexts/programs and one fixed pthread team persist across useful calls.
Each generation publishes its queue only after the complete team exists and
drains before returning borrowed contexts to Python. Idle workers wait on a
condition variable; they join on close. Partial creation failure cancels and
joins the created subset before any arithmetic publication. Worker errors
cancel peers and drain the generation. A failed join invariant aborts rather
than return live borrowed state. Per-owner creation fault injection explicitly
rearms startup on an idle team; normal dispatch does not recreate live workers.
`NativeTeam.lifecycle` retains create/join/generation/worker-entry counters
after close for deterministic lifecycle diagnostics.

No native worker accesses Python or calls Python APIs/callbacks. Checked exports
and publication happen in the Python owner after join. All ready exports are
validated before any publication; stable queue order resolves simultaneous
successes. First proper factor prevents future slices, and discarded work stays
charged. `NativeTeam.cancel()` sets the C cooperative flag during active work;
root-budget cancellation/resource polls also run before/after native returns.
The per-owner C fault/start-gate hooks are inert unless explicitly injected by
diagnostics. They are not process-global runtime policy.

Close is idempotent, releases every context/program/transport and unloads the
pthread DSO once after its workers join and borrowing contexts are destroyed.
Active close and
use-after-close reject. The caller must keep the owner live until a concurrent
cancel call returns; close/cancel may not race. Parked workers retain their owner and DSO until close; no worker borrows a
context between generations. One active team per process is enforced by the
adapter.

C work targets five milliseconds and at most 128 actions per context, with
common wall/process-CPU deadlines and opcode/bit cancellation checks. The cheap
wall clock is checked at every action boundary; aggregate process CPU is checked
at slice entry and after each millisecond of observed wall progress. Root
wall/CPU/RSS and cancellation are still polled before/after foreign slices.
An exceptionally slow operation encounters the wall check at the next action
boundary. Startup, one GMP operation, join and OS descheduling can exceed the target;
it is not a hard real-time limit. H3 root resources include setup, conversion,
recovery and cleanup; no prepaid work/CPU/memory leases are restored.

In the retained `c-gmp` stage-one control, setup, prime/program metadata,
classical stage two, GCD/recovery and publication remain Python. Only its
accepted ladder/PRAC slices execute in C. The separate `c-gmp-full` owner below
executes all classical curve arithmetic/recovery in C; Python supplies metadata
and validates/publishes canonical results. Both refuse paired/wheel policy.
No universal modular-multiplication bottleneck is claimed.

Pthread owner checkpoints use `b4a-pthread-owner-h3-v3`. They retain stable
completed IDs, attempts and cumulative saved wall/CPU/work history after drain.
Unfinished curves restart at their assignment boundary, with no charge refund.
Source, settings, worker count, schedule, capacity and chunk changes reject on
resume. No pointers or partial native state serialize. Historical OpenMP owner
checkpoints do not migrate silently.

## Phase 1 delivery

Polling `7124730` and persistent-team `621da17` have deterministic injected-clock,
generation, partial startup, cancellation, drain and final-join coverage. Their
warm re-baseline establishes no consistent superiority or noise reduction;
empty-dispatch savings do not establish a whole-ECM speedup. `0df3dc3` then adds
the complete classical C owner. Full committed-only QA and cold/warm Tier 1
screening on `fa2a4cc` are complete. Full-native scalar loses these comparisons;
cold parallel gains are scoped and large warm drift makes most warm differences
ties. The earlier checkpoint-design stop is resolved, with its stopped delivery
preserved in `continuation_01/stopped_delivery/`. Exact reports, every failed
attempt and charges remain under `v2/audit/b4a/phase1_20261109/`.

## Validation and retained limits

Tests: `test_native_ecm`, `test_native_correctness`, `test_native_team`,
`test_native_pthread`, and Mach-timebase accounting regressions. The 11 fixed
fixtures run against Python/scalar and 1/2/4 pthreads with both schedules/chunks;
all committed boundaries, factors, failed searches and saturation/recovery agree.
The 60/80/100-digit fixtures have shallow B1=11/B2=47; additional tests cover
accepted production PRAC, ineligible ladder paths, independent affine checks,
carry/width/exceptional operands, actual native signal cancellation, root stops,
creation/worker/preparation failures and repeated saved-history resume.

Local evidence: `v2/audit/b4a/native_admission_research_20261109/`. A bounded warm
fixed-work calibration is retained, with no adopted speedup. Its CPU profiler
attempt failed PyPy timer handling and exhausted recorder output; profile files
are invalid. Exactness and clean timing results remain separate. No further
measurement or compiler tuning follows automatically. The previous two-curve
NEON loss remains scoped; no field-SIMD implementation is integrated.


## Whole-ECM C owner (Phase 1 continuation)

Select `PortfolioConfig(ecm_kernel="c-gmp-full")` explicitly, or CLI
`--ecm-kernel c-gmp-full`. `c-gmp` remains the retained stage-one control;
Python remains the automatic default. `NativeECMExecutor` dispatches whole
native continuations when its configuration selects `c-gmp-full`.

Suyama setup/inversion/nonunits, accepted ladder/PRAC stage one, normal GCD,
strict record/exceptional checked-ladder/prime-unit replay, classical baby/giant
preparation, prime-index traversal, cross terms, product windows and ordered
term replay execute in C using ordinary GMP mpz objects. Python supplies stable
sigma identities and prime/program metadata, checks canonical exports/proper
factors and performs stable parent publication. It remains the independent
arithmetic oracle. No external factorer or new residue/limb representation is
used. Setup and both stages use the same private context on the persistent team.

Native points, baby tables, cursors and scratch persist across bounded slices.
Immutable copied prime/program/index blocks are validated once per owner and
shared by its contexts. Ordinary field/scratch/table capacity is sized and
reused. The owner returns at meaningful checks, metadata refill, errors,
cancellation or the existing 128-action/five-ms cooperative target. Normal
completed-curve publication exports identity/outcome metadata without copying
native tables. Explicit diagnostics can export every complete canonical action.

Scalar checkpoints preserve the existing durable job representation. Canonical
state is committed separately from speculative ladder/record/giant/group work.
Close/checkpoint export retains complete baby/giant actions, the prime cursor,
product/ordered terms and replay positions. An interrupted compound action
exports its saved canonical input and restarts that action on resume; all
incurred charges remain recorded. Failed restore preparation leaves the
original saved job intact. No pointer or partially evaluated point serializes.
Pthread checkpoints keep the existing completed-assignment/history contract
and restart unfinished assignments. Source/kernel changes reject mismatched
saves. GMP allocator death recovers at the last externally saved checkpoint.

The explicit implementation supports odd 3–1024-bit moduli, classical Suyama
with the accepted PRAC/ladder policy, chunks/GCD windows 1–256, segments
1–32,768, uint64 bounds below UINT64_MAX and baby distance <=65,536.
Immutable block geometry supports 4,096 blocks / 64 MiB per owner.
Unsupported geometry or alternate formula/paired/wheel policies raise explicit
errors; native tests do not silently fall back. H3 root wall/CPU/current RSS
remain the resource guards. Geometry refusal is distinct from observed
production exhaustion.

`test_native_full` exercises actual scalar/pthread paths, every existing fixed
fixture/action, strict nonunits/exceptional retry, short tails, table boundaries,
root stops, preparation/worker/startup faults, signals, simultaneous factors,
active/repeated close and durable scalar/assignment-history resume. Detailed
receipts live in `v2/audit/b4a/phase1_20261109/continuation_01/`. The earlier
checkpoint-design stop was resolved by this implementation. Performance
adoption and cost/dispatch recalibration remain separate from correctness.
