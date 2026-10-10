# A7 / P3.8-R5 reconciliation — 9 October 2026

A7's bounded interface/evidence acceptance is complete. SSS/SSSf remain
experimental, native serial remains the default, and E1 owns fresh comparisons
on the final integrated control. This is reconciliation of committed mainline
`76f06b0110941017b926caa51d341db6221f470b`, not a new collector or worker system.
No new factoring timings or performance promotions were made here.

## Acceptance matrix

“Verified” below means inspected implementation plus executable acceptance
coverage; historical measurements retain their original sources and scope.
Test paths are relative to `v2/tests/`.

| Contract | Implemented and verified | Evidence / disposition |
| --- | --- | --- |
| Explicit SSS/SSSf API | `SSSConfig`, `SSSJob`, recursive `PortfolioConfig(sss=...)`, and `--method sss|sssf`; SIQS/SSS alternatives are exclusive; library calls remain quiet. | `test_sss.py`, `test_sss_dispatch.py`: actual extraction, sign/powers, recursion, limits, CLI checkpoints and unchanged defaults. No automatic promotion. |
| Forced divisor and full exponents | Collision candidates carry `abs(F(x))/m`, where `m=M/q` (or M) divides the original F(x). Trees screen that quotient; admission divides the original polynomial value and recovers every base-prime valuation. | Enumeration oracle in `test_sss.py`; new `test_a7_reconciliation.py` independently recovers original valuations, including forced factors and repeated powers, in both modes. SSS currently requires A=1 and multiplier=1. |
| Loss policies and bounded scratch | SSSf's positive cutoff deliberately rejects non-small-base parts at or above the cutoff, except fully smooth parts. Zero disables this rejection, retaining two-stage detection. The second tree covers the full base. Candidate overflow refuses an unpublished assignment; trees, nodes, bits, retained atoms/rows/partials and checkpoint coexistence have finite caps. | `test_sss.py`, `test_performance_repairs.py`: candidate sets, distinct-prime/singular-root collisions, two-stage admissibility and cap refusal. Filtered seven-selection and unfiltered six-selection arms remain separately labelled. |
| Exact relations and provenance | Atomic verification checks the polynomial identity, sign, full factorization, square correction and a proven prime residual within the finite SLP bound. Combined rows retain unique atom IDs and checked residual squares. | `relations.py`, `extraction.py`; `test_qs.py`, `test_qs_pipeline.py`: malformed fields/provenance and composite residual rejection. C1 retains residual/DLP changes. |
| R3 ordering, identities and extraction | Mixed admission order is stable across later full-row admission; equal parity is not an equal row. Prepared/cache identities bind base, full row and pinned atoms. Dependencies are lifted and verified against original rows; extraction checks even totals, square congruence and both GCD signs. | `test_qs_r3.py`, `test_qs_pipeline.py`: exact duplicates, shifted indices, tampered rows, original-operator kernels, refused cache publication and proper divisors. Production cadence remains 1; 32 is scoped opt-in evidence. |
| Worker assignment and publication | Seeded family schedule is independent of worker count; chunk IDs encode family/Gray/block. Ordered parent admission verifies both arithmetic and assignment membership. Incomplete private prefixes replay under the same ID; completed pending batches retain admission cursors. | `test_qs_parallel.py`, `test_performance_repairs.py`: serial/thread/spawn equivalence, future/uncommitted-row rejection, cross-worker restart, exact fixed-work signatures. Race-dependent first-factor stopping is distinct from fixed-work equivalence. |
| Leases, refunds and cancellation | Parent reserves before submission; lease ceiling clips to available work. Admission and pending solving wait for in-flight refunds. Only reported unspent or undispatched work is refunded; consumed/discarded/cancelled work stays charged. Returns drain workers and release guards. | `test_performance_followup.py`, `test_qs_parallel.py`: restrictive grants, solver interruption, failure cleanup, pool reuse and early cancellation. No worker default change. |
| Aggregate resources | Owned-memory reservation includes parent, all workers, bounded coexisting batches and metadata; shared immutable bases are counted once. CPU includes parent and published spawned-process usage; final barrier accounts transfer/drain. | `test_qs_parallel.py`, `test_performance_followup.py`: CPU/wall/storage refusals and shared-base reservation. RSS sums process high-water marks, not simultaneous live RSS. Limits are cooperative, not OS enforcement. |
| Checkpoint identity/migration | Current SSS schema 3 binds target, seed, config, backend/build, base, mixed row order, assignment identity/cursor and solver prefix. Native versions 1/2 migrate with legacy order rules; current versions require mixed order. Parallel versions 1–4 are supported; version 1 cannot carry chunk progress. | `sss_checkpoint.py`, `parallel.py`; `test_sss_dispatch.py`, `test_qs_r3.py`, `test_performance_repairs.py`. New tests create genuine SSS v1/v2 snapshots using hash-checked committed baselines and restore both modes to v3. Backend/config/corrupt progress are rejected. New SSS solver fingerprints use the shared streamed hex-v1 encoding; absent encoding retains legacy decimal-v1, and unknown encodings are rejected. |
| Charged cumulative resume and terminal evidence | In-memory SSS resume extends the original Budget. Serialized restore retains prior work/wall/CPU before charged setup, exact store verification, assignment regeneration and solver replay. Skipped extraction trials are rechecked. Terminal splits and exhausted/refused schedules do not restart. | Both-mode interrupted assignment/solver/extractor coverage; new repeated restoration, exact-work refusal, cumulative deadline, terminal/refusal and original-budget tests. A too-small checkpoint can refuse while the in-memory job remains usable. |
| Independent output checks | Split/cofactor and resolved/unresolved results reconstruct; terminal certainty remains separate from proof-backed fixture knowledge. Upstream stdout is parsed as a candidate claim, then exact proper-split/product and certified-fixture checks apply. Relation-count text cannot establish a factorization. | `test_sss_dispatch.py`, `test_b1_calibration.py`; new invalid upstream stdout claims. Upstream does not expose raw relations to this harness; no independent upstream relation-verification claim is made. |

Polling checks are bounded by at most 64 atomic actions per adapter; external
stop detection can lag by 63 actions, the worker's 1 ms throttle, the parent's
10 ms wait, process startup before first publication and a native bigint action.
These are latency components, not a fixed millisecond overshoot guarantee.

## Accepted repair decisions and retained evidence

| Record / committed pin | Accepted | Retained, rejected or deferred |
| --- | --- | --- |
| [P3.6.1 first repair record](../../docs/studies.md#fresh-confirmation-and-decisions), [pre-repair snapshot](../../inputs/baselines/performance_audit_baseline.json) (31 modules) | Bounded polling, collision invariant reuse, sparse exponent recovery, preparation reuse, compact labels, finite reservations and smaller worker publication. Fresh four-input/two-seed 30-digit cohorts measured SIQS 23.324→6.840 s, SSS 6.181→1.412 s, filtered SSSf 9.460→2.020 s, all 72 timed attempts per arm complete. | Retain native serial/coarse defaults, strict P2 clocks/batching and full second SSSf base. Six-selection unfiltered SSSf beats the historical filtered arm, not SSS under that frozen promotion rule. Reject automatic residual-cap growth and disjoint-base trees. Larger probes remain diagnostic. These historical controls are not current calibrated algorithm rankings. |
| [Post-R3 follow-up](../../docs/studies.md#follow-up-restrictive-allowances-and-repeated-setup), [frozen control](../../inputs/baselines/performance_followup_baseline.json) (35 modules) | Clipped leases/refund-aware admission and solver continuation, disabled-stage reservation removal, exact resieve support bounds, local checked-base/coordination reuse and immutable identity caching. Serial-worker complete time fell 12.3%/12.7% small/medium and 32.2% at B=10000 chunks. Restrictive capacities improved 0/8→8/8 except medium 1M at 6/8. | No additional complete-native-SIQS or cold-process speedup was established. Thread fixed-work improvement still trails serial; unstable first-factor thread controls are excluded. Hensel caching, direct single-hit valuation and other prototypes did not justify adoption. Resieve 8/8 means collected windows, not completed factors. |
| [R3 record](../../docs/studies.md#p38-r3-stable-relations-matrix-and-cadence), [pre-R3 control](../../inputs/baselines/p38_r3_baseline.json) (35 modules) | Stable rows, complete identities, checked legacy solver replay and incremental pivot counts; cadence 32 has scoped 30-digit evidence. | Retain cadence 1 and dependency cache off; no default promotion for live-column compaction, pivot batching or merge histories. Standalone cadence 8 deferred. |
| [B1 calibration](../../docs/studies.md#b1-joint-qsmpqssiqs-calibration--9-october-2026) | Use selected 30-digit flyer/MPQS bundles and separately calibrated 40-digit controls; do not reuse the old narrow SIQS comparator as the final algorithm baseline. | Small independent populations and Pocklington sampling bias remain; inspected held-out inputs are no longer fresh. Wide QS/useful-yield and upper-band censored probes retain their original negative/diagnostic scope. |

The versioned [arm plan](../../inputs/controls/a7_r5_e1_arms.json) pins 148 current
runtime/adapter sources and 16 historical selection/protocol/corpus inputs.
Its `base_commit` identifies the starting mainline; `source_sha256` identifies
the prepared branch bytes. It is not a final combined E1 source freeze.
Original snapshots, source hashes, captures and verdicts were not rewritten.
Raw historical captures remain local; these concise public records and immutable
loadable sources are retained evidence, not newly regenerated timing evidence.

## Findings and accepted changes

- **Actual runner defect fixed:** B1 expected probable labels for every factor
  above 2**64, contrary to A10's accepted deterministic domains. Fresh validation
  now uses the accepted deterministic witness selector. Regression coverage
  checks a proven prime above 2**64 and retains probable status beyond the
  deterministic domain. Historical snapshots and legacy resume labels stay intact.
- **Actual checkpoint defect fixed:** SSS still hashed solver masks through
  decimal conversion and could fail on a wide valid lifted dependency at PyPy's
  decimal digit limit. Reuse SIQS's streamed hexadecimal fingerprint with an
  explicit encoding tag; absent tags retain legacy decimal replay. Current
  SSS schema 3 and supported native v1/v2 migration remain unchanged. Tests
  cover a valid 16,000-row lifted mask, actual legacy solver replay in both
  modes and unknown-encoding rejection.
- **Missing acceptance coverage supplied:** both-mode interrupted reconstruction,
  authentic SSS v1/v2 migration, cumulative reconstruction/refusal/deadlines,
  terminal states, original Budget identity, forced valuations and untrusted
  upstream output claims. No collector mathematics or production default changed.
- **Stale documentation corrected:** A7 no longer waits for “eventual” repair
  decisions or describes serialized SSS selection as absent. A11's usability gate
  was already accepted by the user; the committed outstanding status was stale.
  Its historical CLI measurement scope remains unchanged. Automatic handoff
  remains C3/G1; A11 is not reopened.
- **Genuine compatibility gap retained:** `ParallelConfig` supports the older
  finite reference-family schedule, at most eight A factors and half-width 8192.
  It cannot express B1's flyer stream/quota extension or 40-digit half-width
  65536. Truncating those settings would change the algorithmic schedule.
  A7 does not expand that API or claim calibrated worker superiority.
- **Missing performance evidence deferred:** no final integrated-source, fresh
  SSS/worker confirmation; no calibrated 40-digit SSS arm; no practical broad
  50–100-digit completion or array bottleneck evidence. P3.7 stays deferred.
  These are E1/conditional later gates, not defects requiring a new A7 engine.

## Comparable arms and the exact remaining comparisons

Run `pypy3 -B -m v2.benchmarks.qs.a7.a7_r5` to validate pins and construct the arms;
it performs no timing. `serial_configurations()` decodes exact B1 controls.
`run_serial(fixture, seed, band, arm)` reuses B1's complete-call budget and
classification adapter, selecting the existing SSSJob when appropriate.
Callers own the exclusive performance window, certified corpus and sampling.

| Study | Prepared arms / matched limits | Remaining E1 prerequisite |
| --- | --- | --- |
| 30-digit serial complete calls | Exact B1 QS/MPQS/SIQS selections versus repaired SSS, filtered SSSf and six-selection unfiltered SSSf. Seeds 7/29, 10**13 work, 5 s wall/CPU, 256 MiB owned, one core, native integers. SSS memory ceiling is explicitly normalized from its historical 128 MiB configuration. | Freeze any training-only SSS selection and final integrated source; generate untouched confirmation after selection. Compare useful dependencies, complete classified output, failure classes and total time. A bundle comparison is not an isolated algorithm-only effect. |
| 40-digit serial | Exact separately selected B1 QS/MPQS/SIQS bundles; seeds 7/29, 10**13 work, 30 s wall/CPU, 256 MiB, one core. | Train/freeze feasible SSS/SSSf parameters first, or explicitly defer this challenger class. No invented 40-digit SSS calibration or larger-band extrapolation. |
| Supported worker schedule | Accepted 256-position chunks for small/medium; serial, two/four threads and one/two/four spawned processes; same configuration and family IDs, seeds 7/29, 2*10**9 work, 5 s wall/aggregate CPU, 512 MiB owned and 1 GiB summed RSS high-water gate. Use `phase_three_parallel.run_one` and its native schedule control. | Keep fixed-work versus first-factor, reused versus cold pools separate; include transfer, wasted/cancelled work, drain, restart and aggregate resources. For an E1 calibrated-control claim, train a representable worker schedule or retain/defer workers with the incompatibility stated. This control is not B1 flyer. |
| Combined portfolio | Existing recursive portfolio with selected current-engine changes and resume checks. | B3 owns ECM chain integration; C1 owns residual/DLP; E1 settles selected integrations/deferrals and complete portfolio confirmation. C3/G1 retain allocation/handoff. |

Every pinned corpus has already been inspected and serves acceptance/training
only. E1 must freeze fresh input generation, certified factors, workload strata,
finite budgets, core counts, selection and regression limits before timing.
Use at least three seconds of validated PyPy 3.11 warmup and nine samples;
unstable cases extend to five seconds/15 samples, at most three attempts.
Stop at that finite study limit and label unresolved instability inconclusive.
Prespecify interleaved order and subject-level 95% uncertainty; preserve exhausted
outcomes and the roadmap's five-point class completion-regression ceiling.
Cold startup and instrumented profiles remain separate. No A7 bounded comparison
was needed: the two repair passes already settle its retained/rejected decisions.

Acquire `/private/tmp/factor-performance.lock`, check the owner/process inventory
and publish the performance owner before heavy checks or timings. A7 coordinated
its isolated file scope with B3; it edits no ECM, portfolio, stage-job, relation
or worker production modules. E1 must explicitly refresh the arm pins after
selected shared-source changes instead of bypassing pin failures.

## Research and licensing limits

This pass reused the [existing source/paper reconciliation](../qs_gnfs_research.md)
and accepted experiments; no unresolved mathematical question required another
survey, download or adaptation. References carried forward are [SSS §4,
version 2](https://arxiv.org/html/2301.10529v2#S4), [Bernstein Algorithm 2.1 /
Theorem 2.2](https://cr.yp.to/factorization/smoothparts-20040510.pdf), and the
[SSS upstream pin 8dbaf6d](https://github.com/sbaresearch/smoothsubsumsearch/tree/8dbaf6d39ab88a40380965d25ec2c363d7f27358).
The paper's 75–100-digit results are one-hour relation yields on two inputs
per size, not measured complete-factor speedups. The unchanged upstream arm
has its own SymPy/gmpy2 settings, monitoring and verified output parsing;
it is not a substitute for calibrated current SIQS. No explicit upstream root
license was found in the accepted pin review: clarify permission before copying
or further adaptation. No external code was adapted in A7.

## Verification

Full supported PyPy/GMP QA passes 503 tests and full lint. A committed-files-only
archive of `65f57ba` passes the same 503 tests/lint, 103 benchmark imports,
148 source pins, 16 input pins, three baseline loaders and seven certified
corpora. Only the external development runtime/toolchain is reused; ignored
repository source/evidence is unnecessary. The source-location-free AST review
limits existing function changes to the declared fixes/adapter and acceptance
tests. Raw logs stay in ignored `v2/benchmarks/results/a7-r5/`; commands and
counts are recorded in the [benchmark guide](../../docs/studies.md#a7--p38-r5-reconciliation--9-october-2026).
A7 is closed for its own bounded acceptance; final E1 confirmation remains open.
