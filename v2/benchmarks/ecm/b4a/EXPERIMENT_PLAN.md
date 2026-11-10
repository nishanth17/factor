# Native coordination/backend experiment plan — 10 November 2026

The user pre-approves research through execution and final QA on
`codex/b4a-native-ecm`, base `bc8ea7d`. Local authority, immutable freezes,
receipts and new ledger: `v2/audit/b4a/native_backends_20261110/`.
Historical charges are separate and are never recycled into this grant.

| Stage | Question and decision | Expected wall / initial cap | Pool ceiling wall / CPU |
| --- | --- | --- | --- |
| Observer/build/admission | Does exact monitoring and each native candidate preserve arithmetic, lifecycle, resource and replay contracts? Refuse timing until admitted. | ~120s / 240s, revised by literal substage freezes | 2700 /10800 |
| 3I complete-cost screen | Do separate work/completion signals improve the retained two-anchor complete path? Keep simpler control under uncertainty. | ~200s /300s | Diagnostic/refinement1800 /7200 |
| 3O feasibility | Can same-kernel LLVM OpenMP run complete37-curve batches with exact teams and bounded CPU/RSS? Diagnose before the common grid. | ~30s /60s | Admission pool above |
| Common settings | Runtime effects at1/2/4 threads, static/dynamic, chunk 1, capacity 32, two anchors, three independent inputs. Retain strongest scalar and compiler controls. | ~700s /1000s, literal price after feasibility | 1800 /7200 |
| Named diagnostics/refinement | Profile overhead validated first; answer one residual wakeup/CPU/tail question, then at most one small justified refinement. No profile-per-configuration or repeated broad search. | ~150s /250s, literal price from results | Diagnostic/refinement pool above |
| Final committed-only QA | Full suite/lint/benchmark imports with retained native builds, actual signal/failure/replay coverage and safe cleanup. | ~350s /600s | 2700 /10800 |

These are prospective estimates, not execution grants; each issued substage
freezes exact source/build/input/cells/limits/stop rules first. Total ceilings
are9,000 wall /36,000 aggregateCPU seconds,8 GiB observed owned-tree RSS and 256 MiB
captures. Reserve cleanup within each pool; no automatic enlargement or selective
retry. A lock-busy attempt creates no reservation and does no workload.

## Measurement contracts

Use retained shallow `(2000,50000)` and deep `(11000,1100000)` anchors,
37 fixed curves/unchanged seeds and three independent representative moduli per
anchor. Three repetitions estimate short-call noise and remain input clusters.
Cold includes interpreter/import/JIT/preparation/validation/normal shutdown;
warm retains explicitly declared owners after at least three useful warmup
seconds, with separate close costs and start/end A/A. Record thermal telemetry
outside each timed cell. Preserve exact full-quota signatures and every stop.

Both runtime arms use the provided LLVM 22.1.8, generic O2, ordinary GMP, coalesced
yields and128 actions. Scalar remains256 actions. A bounded Apple/LLVM control
separates the necessary compiler change from runtime effects; no new O3/LTO/
PGO grid. Both backends keep the existing per-generation numerical semantics.
OpenMP counts its caller among requested compute threads; pthreads additionally
have a waiting caller. Count all process CPU/helpers/idle waits once, plus owned
RSS, initialization, output and cleanup. Do not infer physical core placement.

Common blocking coordination uses explicit `OMP_WAIT_POLICY=PASSIVE` and
`KMP_BLOCKTIME=0`; pthread workers sleep on conditions. This avoids retained
OpenMP spin from overlapping interleaved controls. Runtime-default/short-active
waiting is a separate conditional refinement, not silently omitted cost or a
claim that defaults are optimal. Record environment/runtime provenance.

Use 98⅓% paired-log-t input intervals and a measured A/A/history gate for each
initial comparison. A point estimate must clear noise as well as its interval;
ties retain the simpler valid mechanism. Results are exploratory component
comparisons, not simultaneous familywise, population/p90 or adoption evidence.
No fresh broad confirmation follows this grant. Do not time invalid candidates
or require the unrelated H3/B12 child-accounting repair for childless native work.

## Admission and preservation

Keep immutable baseline source and builds. Test actual arithmetic/canonical
oracles, private/aliased contexts, stable factor order, no worker Python calls,
partial startup/refusal, constrained OpenMP teams, cancellation/drain, close,
root CPU/RSS/wall refusal before publication and saved-history replay. Native
context/DSO ownership must survive failures and libomp's retained workers.

Retain useful changes only. Preserve all rejected/inconclusive candidates,
code/build identities, raw evidence, resource stops and charges. Detailed research
remains local; concise implementation decisions and this reusable protocol are
committed before implementation. Publish current configuration/limitations and
next-phase recommendation; no merge, push, routing or default changes.
