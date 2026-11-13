# Optional ordinary-GMP ECM packets and pthread slices

The optional build is now selected lazily by normal eligible portfolio calls.
See the [automatic default scope](../../README.md#automatic-ecm-execution-scope)
for runtime, numerical, ownership and fallback conditions. Normal per-call
40–100-digit execution uses admitted whole-native scalar/pthread recipes;
below 40 and borrowed/retained 80 keep Python after inconclusive confirmation.
Explicit native overrides remain available. No build, CFFI import or native
library load occurs merely because an earlier-method call enables auto intent.
The dated implementation screens below retain their historical scope.

Pthread owners accept 1–16 threads; higher-level worker settings require
`allow_wide_workers=True` for counts outside 1, 2 and 4. The fixed handle and
argument arrays, startup-failure injection and cancellation barriers cover
this range. Capacity remains at most 32 contexts, with the same 128-action,
five-ms slices and root CPU/RSS enforcement. Supported counts are an opt-in
admission range; automatic selection uses the separately measured policy above.

The separate experimental build option `--extended-quantum` permits and selects
256 actions for whole-native scalar/pthread slices. The default remains 128.
Both retain the five-ms target, native clock/cancellation checks, root CPU/RSS
polling and publication ordering. An explicit smaller scalar `packet.actions`
override still takes precedence. Legacy generic packet callers keep 128 actions.

Python integers on PyPy3.11 remain the arithmetic backend. Historically, B4a admitted
whole-classical scalar/pthread correctness on `fa2a4cc`: full QA runs 1,741 tests
with 20 explicit optional skips, lint and 308 ordinary benchmark imports. All
20 new whole-native test methods execute. The earlier `c-gmp` stage-one owner
remains an independent control; `c-gmp-full` requests the complete C owner
explicitly. The later automatic activation supersedes the original opt-in-only
recommendation within its documented scope.

The Phase-2 experimental whole-ECM driver observes a childless macOS root
through one CFFI bridge using the same libproc usage/child queries as the
Python resource provider. Poll locations, cancellation, wall/CPU/RSS decisions
and publication ordering are retained. A live child, multiple roots or a
custom provider uses the original provider. Observations are never replaced
with a cancellation flag or an estimated memory quota. The provider is restored
after each poll, including errors, and its buffer is released before unload.
The initial sizing-query implementation is superseded by an actual child-ID
query with explicit fast-path admission. Focused correctness and matched
three-input complete-cost screens support the corrected driver and generated
bindings. Final ordinary and optional-limb committed-only builds on `b620d6f`
each pass 1,751 tests (20 optional skips), lint and 308 benchmark imports.
The final cold size matrix and separate warm anchors are complete; comparisons
retain Python, C1, uncertainty, historical failures and all charges locally.
That historical Phase-2 delivery did not change automatic selection. Later
activation uses fresh actual-production-plan calibration and confirmation.

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
ineligible or earlier-completed calls. Declared missing/corrupt descriptors
raise for explicit execution; optional automatic admission may select Python
before any curve action. Post-admission faults never silently rerun elsewhere.

Both retain `-std=c11 -O2 -Wall -Wextra -Werror`. The scalar artifact is in
`build/`; pthread is separately in `build/pthread/`, with `-DB4_PTHREAD -pthread`.
Pthreads use libSystem on macOS; no OpenMP runtime or new dependency is needed.
`--optimization O2|O3` and `--cpu-target NAME` are explicit experimental build
choices; generic O2 remains the default after the bounded screen. Unsupported
CPU names are compiler errors. These builds permit the compiler's usual
automatic vectorization; the old NEON diagnostic's disabling flags do not
apply. Separate LTO and PGO comparisons also retain ordinary O2 under complete
cost uncertainty. PGO uses disjoint training inputs and separately validated
scalar/pthread profiles. Six compiler configurations include baseline and the
training-only instrumented build; no LTO/PGO combination was tested. Exact
experimental compiler wrappers and profiles remain in the local phase audit.
`--limb-division` enables a separate 4/5/6-limb
GMP division prototype with ordinary-GMP fallback for other widths/operands.
It passes exactness tests but is disabled by default because complete-cost
evidence is inconclusive. No SME/NEON probe is a production backend.
`--coalesce-unit-yields` is a separate experimental driver option. It keeps
normal GCD checks inside the same bounded slice after a unit result, retaining
factor/nonunit/recovery, phase, metadata and atomic diagnostic returns. The
128-action/five-ms limits, cancellation and fresh root checks around every
foreign call and before publication remain. It remains disabled by default;
the scoped explicit delivery recipe enables it.
A separate whole-curve OpenMP experiment uses `b4_omp_step` with a real
parallel region and static/dynamic worksharing over the same `b4_slice` kernel.
It is selected only through `v2.benchmarks.ecm.b4a.omp_owner.OpenMPExecutor`;
production routing and the pthread owner remain unchanged. The legacy
`b4_step_until` capability remains a separate entry point.

Use the provided compatible LLVM/libomp explicitly, under the same build lock:

```sh
pypy3 -B -m v2.ecm.native.build --gmp-prefix /opt/homebrew/opt/gmp --coalesce-unit-yields --cc /opt/homebrew/opt/llvm/bin/clang --sysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk --openmp-prefix /opt/homebrew/opt/libomp
```

The explicit SDK overrides this host's stale LLVM default SDK path. OpenMP
builds use `-DB4_OPENMP -fopenmp -lomp`, record linkage/runtime hashes, and live
in `build/openmp/`. A pthread binary reports zero OpenMP capability and links
libSystem. Neither executor combines the two team implementations.

Requested OpenMP counts include the caller. Actual-size mismatch refuses
before kernel execution. Each implicit barrier drains all borrowers before
return; cancellation is owner-local. Contexts and metadata are freed at close,
while libomp's parked runtime threads and the loaded DSO have process lifetime.
Counters do not claim those runtime threads were owner-created or joined.
Use `OMP_WAIT_POLICY=PASSIVE KMP_BLOCKTIME=0 OMP_DYNAMIC=FALSE` for the common
blocking comparison. No compiler, runtime or SDK is installed automatically.

Separate pthread work/completion signaling was admitted and screened, then
removed from retained source because the apparent gain did not clear measured
A/A noise. Its opt-in implementation and exact evidence remain in history.

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
1..16 (with explicit wide-worker admission above),
`schedule="dynamic"/"static"`, and `chunk_size=1..32`. The ready capacity
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
Omitted selection follows the automatic policy above. `NativeECMExecutor` dispatches whole
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

## Remaining request delivery (9 November 2026)

The explicit recipe retains generic O2/ordinary GMP and exact root checks,
with coalesced unit yields, scalar 256 actions and pthread 128 actions.
Registration uses the same exact root transport. Independent GMP expression
fusion did not establish useful improvement and was removed. The later
automatic policy is distinct from 1–16-thread technical admission.

Build the selected scalar and pthread pair explicitly, under the shared lock:

```sh
pypy3 -B -m v2.ecm.native.build --gmp-prefix /opt/homebrew/opt/gmp --coalesce-unit-yields --extended-quantum
pypy3 -B -m v2.ecm.native.build --gmp-prefix /opt/homebrew/opt/gmp --pthread --coalesce-unit-yields
```

Use `NativeECMExecutor(settings=WorkerSettings(workers=4),
schedule="dynamic", chunk_size=1, capacity=32)` with an explicit
`c-gmp-full` configuration. This is a measured opt-in recipe, with existing
Python/B12 controls still applicable. Automatic engagement additionally uses
the production-plan and ownership evidence linked above.
Reuse the built dylib, checked ABI descriptor and manifest across calls;
rebuild when source, flags, toolchain or dependency identity changes.
A relocatable `.o` alone cannot be loaded by CFFI's dynamic-library loader.
Factoring calls never compile. CFFI is imported lazily by
`NativePacket._load` in `ecm/native_packet.py`; `ffi.dlopen` loads the checked
scalar or pthread dylib, and explicit owner close releases it after all users.

Fresh cold/warm 60/80/100-digit ECM and failed-search comparisons include
balanced 60/70 cases, startup/preparation and explicit cleanup. They retain
quota outcomes, A/A uncertainty and incomplete 80-digit warm confirmation;
no uniform-population/p90 gain or balanced completion pass is established.
Cheap paths never load CFFI/DSOs or start workers. Explicit native selection
still incurs cold Python-module/checkpoint setup; keep default Python exits.
H3/B12 child CPU units/history and combined-source engagement tables require
owner revalidation before broader integration.

Both committed-only configurations on `742cffb` pass 1,756 tests (14 optional skips), lint and 309 benchmark imports. The limb controller then hit a capture/delete cleanup race;
its failed receipt and full 900-CPU-second charge remain. A separately
validated observer rejects live access failures, safely handles proven
file disappearance, verifies saved evidence and completes owned cleanup.
The whole-curve OpenMP challenger is an explicitly selected experiment;
selected teams use pthreads and do not combine the two runtimes. Current
builds retain normal compiler auto-vectorization; no new SIMD kernel is
selected. Historical compiler comparisons, NEON and SME losses remain scoped.
Local delivery, complete cost/censoring tables and all resource accounting:
`v2/audit/b4a/remaining_native_ecm_20261109/DELIVERY.md`.

The 10 November coordination/backend phase retains pthread signaling and the
explicit recipe above. The real whole-curve OpenMP challenger is admitted;
common runtime/compiler results are inconclusive, with only a limited deep-cold
benefit in refinement. Short-active waiting adds CPU without an established
wall gain. Committed-only source `5194951` passes 1,761 tests (14 optional skips),
lint and 310 ordinary benchmark imports.
Current local costs, delivery and all separate charges:
`v2/audit/b4a/native_backends_20261110/DELIVERY.md`.


## Post-integration arithmetic and resource transport (11 November 2026)

The optional native implementation retains the same formulas and canonical
checkpoints. Stage-two difference products use signed `mpz_submul`; ordinary
products use one private preallocated numerator. On 64-bit, nail-free GMP and
4–6-limb positive moduli, the existing GMP `mpn_tdiv_qr` API reduces that product.
The bounded quotient buffer, disjoint limbs, zero/negative handling and wider
fallback preserve exact signed remainders. This uses no new modular domain or
custom reduction algorithm. The earlier `--limb-division` full-mpn product
experiment remains separate and unselected.

Builds expose `wide_gmp_remainder`; first native admission checks it against the
loaded library. The generated CFFI ABI descriptor, dylib, manifest and source/
dependency identities remain a unit. Cached owners reuse them; factoring calls
never compile. Native checkpoint identities change with these artifacts, while
canonical state and serialized formats retain their existing contract.

Root observation removes only the unused child-list sizing call. An actual
one-PID read still selects full-tree observation for any child, and zero with
errno refuses access. All CPU/RSS/wall/cancellation checks and cleanup remain.

Generic O2, scalar 256/team 128 and pthread 4/dynamic/chunk1/capacity32 remain
selected. Current-source 16/32-byte alignment did not improve deep complete
costs. Real OpenMP remains an alternative explicitly supplied build using
`-fopenmp -lomp`; scalar/pthread images do not link libomp. No new SIMD backend
or automatic selection rule is added.

The final exposed cold/warm controls, corrected complete owner/helper CPU,
recovery, cheap guards, exactness and CPU/heap/JIT attribution are documented in
`v2/audit/b4a/post_integration_tuning_20261110/HANDOFF.md`. Fresh 3F is frozen and
unexecuted. The measured 5.5% deep warm pthread/entry point effect has an
individual model interval; it is not a general factorization/default claim.

The current whole-owner capability is checked against the build manifest and
`b4_team_max_workers()`. The legacy stage-one `b4_team_step` control retains
its explicit1/2/4count contract. `NativeTeam.worker_activity()` reads idle-only
per-rank identities and task/curve counters; sampled explanatory thread clocks
are opt-in and excluded from clean comparisons. Process CPU already includes
all native threads; these diagnostic values are not added to H3 accounting.


The count screen validates whole-owner pthread/OpenMP counts through sixteen,
including partial startup, cancellation, root stops, draining and replay.
OpenMP `worker_activity()` reports logical ranks including the encountering
caller and their completed slice/curve totals. Dynamic short work can legally
leave a participating rank without a task; static admission fixtures separately
exercise every rank. The clean deep pthread screen retains eight: no higher
count clears its prospective complete-cost gates. Detailed exposed count costs,
resource/thermal limits and the required unexecuted 3F revisions are local in
`v2/audit/b4a/oversubscription_20261111/OVERSUBSCRIPTION_HANDOFF.md`.
