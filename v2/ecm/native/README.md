# Optional ordinary-GMP ECM packets and pthread slices

Python integers on PyPy3.11 remain the default. The isolated B4a continuation
admits scalar and pthread correctness on committed source `0f2fee3`: full QA
runs 1,712 tests with 20 explicit optional skips, plus lint and 308 ordinary
benchmark imports. Actual scalar/pthread arithmetic tests execute. Performance
adoption, engagement recalibration and mainline integration remain open.

Build explicitly from the repository root under the shared performance lock:

```sh
pypy3 -B -m v2.ecm.native.build --gmp-prefix /opt/homebrew/opt/gmp
pypy3 -B -m v2.ecm.native.build --gmp-prefix /opt/homebrew/opt/gmp --pthread
```

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

Setup, prime/program metadata, classical stage-two preparation/execution,
GCD/strict nonunit recovery and final publication **remain Python**. Only accepted
stage-one ladder/PRAC slices are dispatched to C. Paired/wheel continuation is
refused by this owner. Selected B12 metadata caches remain supported; no
whole-curve native coverage or universal modular-multiplication bottleneck is
claimed.

Pthread owner checkpoints use `b4a-pthread-owner-h3-v3`. They retain stable
completed IDs, attempts and cumulative saved wall/CPU/work history after drain.
Unfinished curves restart at their assignment boundary, with no charge refund.
Source, settings, worker count, schedule, capacity and chunk changes reject on
resume. No pointers or partial native state serialize. Historical OpenMP owner
checkpoints do not migrate silently.

## Phase 1 overhead continuation

Polling commit `7124730` and persistent-team commit `621da17` have committed-only
targeted correctness coverage, including deterministic injected clocks,
generation reuse, partial startup, cancellation, drain and final joins. The
same-cell warm re-baseline does not establish native superiority or consistently
lower A/A noise. The matched lifetime screen is inconclusive for complete
ECM; its empty-dispatch probe supports a narrower coordination-cost saving.
The full-ECM port stopped before implementation on an unresolved canonical
recovery/checkpoint export design. Setup, classical stage two and strict
recovery therefore remain Python. Exact reports, failed harness attempts and
charges are retained under `v2/audit/b4a/phase1_20261109/`.

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
