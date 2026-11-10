# Native coordination and OpenMP decision — 10 November 2026

Base `bc8ea7d27893c3369c77ac0d63c551d6a5f51d37`; detailed source hashes,
read-depth inventory and preserved predecessor charges are local under
`v2/audit/b4a/native_backends_20261110/`. The previous six-item continuation
is complete. Keep ordinary GMP/generic O2, coalesced unit yields, scalar 256,
pthread 128 and four/dynamic/chunk 1/capacity 32 as the current explicit recipe.
Its censored confirmation, generic child-CPU dependency and earlier rejected
expression/compiler/SME/NEON results retain their scope.

## Source-backed decision

| Rank | Mechanism | Evidence and decision |
| --- | --- | --- |
| 1 | Separate work and completion signaling | Current `team.c` broadcasts one condition on entry and completion, waking both workers and caller. The retained deep 37 diagnostic uses 4,828 generations per batch, so this is repeated work, not an empty-dispatch-only opportunity. Implement one bounded candidate; compare complete costs against immutable control. Waiting sample shares are not removable fractions. |
| 2 | Whole-curve OpenMP | Provided LLVM/libomp 22.1.8 makes a real comparator feasible. Implement a scoped C parallel/for owner calling the same `b4_slice`; the legacy `step_until` entry is ineligible. Count the primary caller among compute threads and reject an actual-team mismatch before arithmetic. |
| 3 | Caller participation in pthreads | FLINT demonstrates a caller-plus-helper architecture. It may save a hop at small teams, but changes startup/failure accounting and is not established by idle samples. Reconsider only if complete costs/profiles after separate signals expose material coordination; otherwise retain the current caller role. |
| 4 | Wait policy or scheduling refinement | Common OpenMP blocking policy matches pthread sleep semantics. Test a small alternative only for an exposed wakeup/CPU/tail question. No repeated broad worker/compiler grid, weighted scheduling or QoS without imbalance evidence. |
| 5 | Plan/scratch/algorithm changes | Immutable blocks and private GMP scratch already exist. No new cache, layout, reducer, SIMD, chain, prime generation or continuation work is authorized by this decision. |

FLINT's pinned ECM caller distributes independent curves, shares immutable
prime/sigma data, allocates private scratch and executes one worker on the
caller. Its pool uses distinct work/completion conditions. YAFU's pinned ECM
caller dispatches curves and validates factors on synchronization; its pool
separates per-worker work signaling from the completion queue. These are
architectural evidence, not copied code or upstream speedup claims.
[FLINT caller](https://github.com/flintlib/flint/blob/0bb30ef32c7d136b29214fd5b40fb72f31f7a200/src/fmpz_factor/ecm.c),
[pool](https://github.com/flintlib/flint/blob/0bb30ef32c7d136b29214fd5b40fb72f31f7a200/src/thread_pool/init.c),
[YAFU caller](https://github.com/bbuhrow/yafu/blob/c7d90283972820f99b3b4920bfc1eff05c1f6342/factor/gmp-ecm/ecm.c),
[pool](https://github.com/bbuhrow/yafu/blob/c7d90283972820f99b3b4920bfc1eff05c1f6342/ytools/threadpool.c).

CADO's pinned batch cofactorization uses dynamic per-candidate scheduling for
variable method-list work and stores results per index before ordered output.
GMP-ECM's reviewed OpenMP callers parallelize fast-stage-two polynomial/residue
work with private modular scratch and explicit combination. They do not prove
that its ordinary ECM curve caller uses OpenMP, or authorize a continuation
change here.
[CADO](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/batch.cpp),
[GMP-ECM](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/pm1fs2.c).

OpenMP includes the encountering primary thread in the team and completes an
implicit barrier before returning. Actual team size can differ from a request.
LLVM's pinned fork/join, hot-team, dispatch and wait paths explain why retained
workers and active waiting must be measured; they do not predict an M4 winner.
The matching release source is `ca7933e47d3a3451d81e72ac174dcb5aa28b59d1`;
installed binaries/configuration are independently hashed, not asserted bitwise
identical to upstream. Review depth is selected caller/functions, not whole
project reentrancy certification. The MCS 1991 paper's abstract informs contention
risk only; full PDF access failed, so no unreviewed proof is claimed.
[OpenMP parallel semantics](https://openmp.org/spec-html/5.2/openmpse57.html),
[LLVM runtime](https://openmp.llvm.org/design/Runtimes.html),
[pinned fork/join](https://github.com/llvm/llvm-project/blob/ca7933e47d3a3451d81e72ac174dcb5aa28b59d1/openmp/runtime/src/kmp_runtime.cpp),
[primary synchronization paper](https://flint.cs.yale.edu/cs422/readings/papers/mcs91algorithms.pdf).

## Legal region and ownership

The unit is one bounded `b4_slice` on a distinct ready context. Share only
immutable plans, queues and deadlines; each job owns mutable GMP state and an
output slot. Validate nonaliasing before entry. Form the requested team, verify
its actual size, then execute static/dynamic loops with matched chunks. Reject
nesting or constrained teams without arithmetic. An owner-local atomic stop
flag is checked by kernels; do not rely on optional OpenMP cancel directives.
Drain all threads before export, resource checks or stable-order publication.

No Python APIs/callbacks run in workers. Preserve five-ms/action limits, root
wall/CPU/RSS polls, GCD/recovery, peer waste and saved assignment/history replay.
A benchmark-only OpenMP adapter may reuse the existing Python owner through
explicit backend hooks; its identity must reject cross-backend checkpoints.
Owner close frees private data; libomp workers may remain parked and the DSO
must stay retained until process exit. Record this memory/cleanup distinction.

## Execution boundary

Proceed automatically through the priced plan in
[`EXPERIMENT_PLAN.md`](../../benchmarks/ecm/b4a/EXPERIMENT_PLAN.md).
The new request supersedes the older standalone research-review stop. No broad
confirmation, production OpenMP integration or automatic promotion is implied.
