# Optional experimental ECM packets and OpenMP team

The default remains Python integers on PyPy3.11. This preparation implements an
optional ordinary-GMP C packet and a separate OpenMP owner. Build, executable
correctness checks and performance acceptance are pending. Do not interpret
source preparation as a demonstrated speedup or a completed optimization gate.

Build explicitly from the repository root, under the shared performance-window
rules. Factoring startup never compiles code or installs dependencies:

```sh
pypy3 -B -m v2.ecm.native.build --gmp-prefix /opt/homebrew/opt/gmp
pypy3 -B -m v2.ecm.native.build --gmp-prefix /opt/homebrew/opt/gmp \
  --openmp-prefix /opt/homebrew/opt/libomp --vector-report
```

The first build is independent of OpenMP. The second uses the supplied omp.h and
libomp with explicit Apple Clang flags and creates build/openmp separately.
Build manifests record compiler/flags/target, PyPy/CFFI, GMP version/library hash,
OpenMP header/runtime hashes, source, library hash and ABI1. No binary is vendored.
The initial build target is ARM64 macOS; installed files do not establish that
compilation/linking/runtime tests have passed. No system package changes occur.

Scalar use is explicit:

```sh
pypy3 -B -m v2.factor NUMBER --bounded --ecm-kernel c-gmp
```

The library option is PortfolioConfig(ecm_kernel="c-gmp"). "auto" falls back for
unavailable dependencies or ineligible modulus/scalar widths; corruption and
identity mismatch raise. Explicit c-gmp refuses unavailable/unsupported inputs.
Cheap exits do not load a library. Scalar packets also work in existing B12
spawned processes. Python-level gmpy2 operators are a different intervention.

The packet accepts odd3..1024-bit moduli, canonical unsigned big-endian values,
<=4096-bit combined scalars or<=16 verified accepted records (<=512 operations,
<=16 point registers each). It preserves formulas, B4's fused ladder, accepted
PRAC eligibility/programs, strict frontiers and Python exact recovery. Bounds,
curve families and prime coverage are unchanged. Public/checkpoint values remain
canonical. The adapter uses documented CFFI/GMP interfaces and explicit release;
no interpreter-private bigint casts or global allocator hooks are used.

Each scalar context reserves2MiB for bounded transport/program/native scratch,
in addition to the ordinary portfolio workspace. Inputs/descriptors/output
buffers and GMP scratch are reused. Canonical n/a24 identity bytes permit native
import reuse. Program packing/verification is cached within its funded owner.
Setup, file/dependency hashing, copies and unused preparation are charged;
mathematical/recovery work is funded by the caller before dispatch. GMP's default
allocation failure may terminate the process; it is not a catchable quota refusal.
The last external committed checkpoint and fully retained grant remain recovery
boundaries after process death.

OpenMP use is a separately selected library owner, with one active team per
process and no Python callbacks in native worker threads:

```python
from v2.execution.budget import Budget
from v2.execution.ecm_native import NativeECMExecutor
from v2.execution.portfolio_config import PortfolioConfig
from v2.execution.worker_contract import WorkerSettings
from v2.portfolio import factorize_bounded

budget = Budget(work_limit=2_000_000_000, seconds=30, cpu_seconds=120)
settings = WorkerSettings(
    workers=2, worker_memory=128 * 2**20, parent_memory=1088 * 2**20,
    memory_bytes=3 * 2**30, checkpoint_bytes=2 * 2**20,
    assignment_work=budget.work_limit, assignment_cpu=budget.cpu_seconds,
)
owner = NativeECMExecutor(seed=7, settings=settings, schedule="dynamic")
try:
    run = factorize_bounded(
        NUMBER, seed=7, config=PortfolioConfig(ecm_kernel="c-gmp"),
        budget=budget, ecm_executor=owner,
    )
finally:
    owner.close()
```

The initial owner requires per-curve caps covering the whole global grant and
refuses stricter ones. It does not silently ignore a requested per-curve CPU
policy. Process CPU already includes native thread CPU; no child CPU is added.
Setup/recovery/classical stage two remain in the Python parent. Accepted stage-one
PRAC/ladder continuations are parallelized; a whole-curve C port is not assumed.
The owner rejects paired/wheel continuation and supports the current classical
path, including selected B12 metadata caches. It owns stable B12 curve IDs and
private curve jobs, respecting tier order and a bounded ready set. First proper
factor stops future slices; stable queue order resolves simultaneously ready
results. Late/uncommitted work remains charged.

Choose workers1/2/4 and schedule="dynamic"/"static"/"guided". All use chunk1;
dynamic1 is the baseline hypothesis, not a measured optimum. Default ready-set
capacity is2*workers, explicit1..32 supported with simultaneous memory admission.
Before the first runtime load, the experiment launcher sets passive wait,
KMP_BLOCKTIME=0, OMP_DYNAMIC=FALSE, OMP_MAX_ACTIVE_LEVELS=1, OMP_THREAD_LIMIT=4,
KMP_STACKSIZE=4M and no hidden helper threads. Library users must control/account
for their runtime environment and other parallel work. No portable M4 P/E pinning
or process x thread hybrid is selected.

C checks a common wall/process-CPU deadline at point/opcode boundaries, with a
5ms cooperative work target and128 actions/context per call. Startup, one GMP
operation, barriers and OS descheduling can exceed this target; it is not a hard
real-time bound. Budget/signal publication checks run between native returns.
Every region drains before Python frees borrowed inputs. Runtime workers may
remain parked until process exit, so its DSO stays mapped; close clears curve
contexts/programs/transport, not libomp's process-wide pool. Runtime thread
creation failure may terminate; checked reduced-team admission publishes no
native result. The future owned-group watchdog covers process lifetime/failure.

The native owner checkpoint records completed curves, attempts and cumulative
charges after drain. Unfinished curves restart from their stable assignment
boundary on resume, without refunds; attempts and run count are finite. Changing
source, schedule, capacity or team count on resume rejects. No native pointers,
partially evaluated points or live thread state are serialized.

There is no new field-SIMD implementation. OpenMP simd/declare-simd and ARM NEON
were assessed separately from threading. Opaque mpz calls do not become vector
multiprecision arithmetic from a pragma; full60/80/100-digit moduli still need
full-width arithmetic. Metadata SIMD and compiler reports are separate from point
SIMD. No ornamental immintrin.h/arm_neon.h include or custom arithmetic is added.

Targeted tests are test_native_ecm.py and test_native_team.py. The B4a QA runner
builds both variants and then runs exact arithmetic, resource and scheduling
checks under a finite locked grant. Stage Q also runs the full suite and a fresh
committed-source archive check. Compilation, actual PyPy GIL/signal behavior,
OpenMP startup/teardown and generated-code diagnostics remain pending until that
window is released. The prospective runner permits no performance execution
without exact integrated source, passed Q, pins and a separately admitted freeze.
