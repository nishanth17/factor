# Factor status

Last updated: 4 November 2026 (America/Los_Angeles).
Integration base: `99a2943`; C10a-01 terminal source: `d37b809`.
A1b review record: `61899a9`; production defaults are unchanged.
QS roadmap review: `4d65d33`, consolidated into existing task scopes on
2 November; C10a-01 is now integrated with the selected policies enabled
for every fresh automatic SIQS configuration by explicit user direction.
The expanded 2 November ECM/QS literature and source review adds A1d
multi-polynomial buckets and conditional F9 genus-two/Kummer research to the
ECM/acceleration amendments. It changes no defaults and launches no runs.
Local supporting reviews: `v2/audit/research/ecm_qs_review_20261102.md` and
`v2/audit/research/ecm/roadmap_review_20261102.md`.
The broader portfolio follow-up adds B4a fused native ECM chunks and expands
E3/E4/E5/E6/A10a/A2a opportunity gates. Local supporting review:
`v2/audit/research/portfolio_review_20261102.md`. No experiment is launched.
The subsequent QS CFFI/SIMD pass refines B5a with queue-compatible buffers,
root/mark/scan packets and separate M4 NEON versus x86 AVX dispatch. Local
review: `v2/audit/research/qs_cffi_simd_20261102.md`.

**Target: complete-factorization wall time on the Apple M4 / 24 GiB / 10-core
machine; all cores are available.** Choose worker counts by measured latency
and aggregate memory, with total CPU reported as a secondary cost.

**Primary workload: uniformly sampled random integers at 60, 80 and 100 decimal
digits, with equal weight across the three size bands.** Report median and p90
complete-factorization time and completion rates for each band. Balanced
semiprime completion at 60 and 70 digits is a must-not-regress constraint.

**Native code is allowed as an optional backend; a working pure-Python/PyPy
fallback must remain available.** Select one justified implementation per
task; NumPy and CFFI are alternatives, not mandatory twin developments.
These are planning and acceptance objectives, not implemented default changes.
See [metric and resource rules](POLICY.md#scope-and-standing-directions).

Read this first, then the [live backlog](ROADMAP.md) and [stable policy](POLICY.md).
Full historical results and gates remain in the [history index](history/README.md).
Useful progress toward 60–80-digit factoring remains the immediate focus;
100 digits is also a measurement target, with completion feasibility still
unestablished. Existing experiments retain their frozen scopes and budgets.

## Current automatic defaults

These are current committed behavior, not claims that every routed band passed
an experimental gate. Sources: [SIQS selector](qs/defaults.py),
[portfolio routing](portfolio.py), [ECM selector](execution/ecm_defaults.py) and
[chain eligibility](ecm/chain_options.py). Explicit settings and saved policies
win. Native means Python integers on PyPy implementing Python 3.11.

The bounded portfolio runs exact preprocessing, rho, p−1 and ECM serially.
Omitted native configurations automatically enable SIQS at exactly 40 and
50–60 digits. Direct automatic SIQS calls use the same size selector. Explicit
`siqs=None` disables that fallback; reference `SIQSConfig()` keeps its own
fields. SSS and production workers remain opt-in/experimental.

| Decimal digits / scope | Owned memory allowance | Automatic SIQS selection | Adoption evidence and limits |
| --- | --- | --- | --- |
| Exactly 40, native | 1 GiB child; omitted outer portfolio adds 16 MiB | Growing DLP; B=5,000, 5 A factors, pool 32; 2 MiB power plan; parity_two marking; window pacing; coupled scored multiplier | [C10 plan adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10--source-derived-collection-strategies-closed-bounded-tranche). Growing versus fixed DLP was not separately confirmed at 40. Face/queue adoption below applies. |
| Exactly 50, native | Same | Growing DLP; B=30,000, 6 A factors, pool 32; 2 MiB plan; density+presieve+scan_fused+parity_two; collection pacing; coupled scored multiplier | [C9](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c9--source-derived-large-prime-yield-and-dlp-economics-complete-and-integrated), [C10](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10--source-derived-collection-strategies-closed-bounded-tranche), [C10a-05](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open), [A1a adoption](benchmarks/RESULTS.md#a1a-integration-and-extension). A1a's 50-digit integration gate passed. |
| Exactly 51, native | Same | Same concrete bundle as 50 | [51-digit completion adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open). A1a marking at 51 is an **unmeasured interpolation**. |
| 52–54, native | Same | Counted base 2,500 (endpoint ceiling 250,000), 7 A factors, pool 64; coupled odd/even scored multiplier (`multiplier=0`); 2 MiB plan; A1a marking plus parity_two; collection pacing | [C10a-04 adoption](benchmarks/RESULTS.md#c10a-04-larger-band-routing) and [A1a adoption](benchmarks/RESULTS.md#a1a-integration-and-extension): **unmeasured interpolated bands**. |
| Exactly 55, native | Same | Same counted bundle as 52–54 | Same [C10a-04 adoption](benchmarks/RESULTS.md#c10a-04-larger-band-routing) and [A1a adoption](benchmarks/RESULTS.md#a1a-integration-and-extension): **Adopted**; compact training unstable, A1a fresh integration Inconclusive due to drift. |
| 56–59, native | Same | Counted base 5,000, ceiling 250,000; 8 A factors, pool 64; coupled scored multiplier; 2 MiB plan; parity_two marking; collection pacing | [C10a-04 adoption](benchmarks/RESULTS.md#c10a-04-larger-band-routing): **unmeasured interpolated bands**. |
| Exactly 60, native | Same | Same counted bundle as 56–59 | [C10a-04 adoption](benchmarks/RESULTS.md#c10a-04-larger-band-routing): **Adopted**; preliminary fresh compact confirmation; original campaign incomplete. |
| Exactly 40 or 50–60, native, different explicit allowance | Caller grant; outer portfolio subtracts coexistence reserve before child selection | Same size/base/multiplier/growing policy; power plan off, parity_two marking, window pacing. Finite refusal remains possible. | [C9/C10 routing scope](qs/CONFIGURATION.md#explicit-and-automatic-selection); resource variants are **unmeasured**, not additional accepted gains. |
| Other sizes, or other backends, when SIQS is selected | Reference/caller allowance | Reference size parameters with the coupled policies below; no automatic queue. No size-based portfolio SIQS enablement. | [Retained policy](qs/defaults.py) and [C10a-07 adoption](benchmarks/RESULTS.md#c10a-07-face-and-queue). No calibrated size preset at 41–49 or 61+; the coupled-policy promotion below still applies. |

Every fresh automatic SIQS configuration at every supported size/backend
selects `sign_policy="swap_high"`, `polynomial_policy="parity2"`,
`multiplier_policy="parity2_wide"`, `multiplier=0`,
`assignment_policy="windowed_pair"`, `root_update_policy="batch"` and composed
`parity_two` marking. The C10a-01 all-size promotion is an explicit policy
decision before fresh confirmation; the [results and limits](benchmarks/RESULTS.md#c10a-01-coupled-polynomial-defaults)
retain the 65-digit failed canary and censored 70-digit control. Other sizes
keep generic reference parameters and grants with these selected policies;
there are no new calibrated 61+ presets or portfolio bands.
Native calibrated collectors above select a 32 MiB sparse-queue ceiling within
the existing owned allowance, with charged rebuild and exact fallback for
unsupported paths or capacity refusal. The combined C10a-07 policy is adopted.
Its 50/55/60 panel is exploratory; drift acceptance was unmet and the additional
18-call confirmation was not run. Other sizes,
backends and grants are unmeasured policy extrapolations. See the
[adoption and limitations](benchmarks/RESULTS.md#c10a-07-face-and-queue).

Omitted calibrated SIQS grants are 10^13 work units and 1,200 wall/CPU seconds
through 55 digits, 1,800 at 56–60. All calibrated bundles use half-width 65,536,
windowed paired assignments and a 64-polynomial family prefix. Explicit limits win;
these grants promise no completion. The [configuration guide](qs/CONFIGURATION.md#explicit-and-automatic-selection)
explains residual/endpoint/product controls and charged base-count resolution.

| ECM / shared portfolio scope | Owned memory / eligibility | Current automatic choice | Adoption evidence and limits |
| --- | --- | --- | --- |
| Exactly 30 digits | Existing funded cap; retain previous plan if coexistence fails | compact64: `(2000, 50000, 64)` | [C3 decision](benchmarks/RESULTS.md#c3-bounded-ecm-allocation-and-handoff-10-october-2026): **Adopted before fresh confirmation**; oversized comparison Inconclusive, source screen unstable. |
| Exactly 40 digits | Same; SIQS headroom above applies | escalating140: `(200, 7700, 10)`, `(2000, 81000, 30)`, `(10000, 554000, 100)` | Same [C3 adoption](benchmarks/RESULTS.md#c3-bounded-ecm-allocation-and-handoff-10-october-2026) and limitations. |
| Other sizes, or explicitly selected paired executors | Caller/retained cap | Retained fixed32 numerical schedule | [C3 scope](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#calibrate-automatic-selection-and-handoff-c3-then-g1); not a measured optimum for every size/backend. |
| Supported native schedule with B1=2,000, at least 8 curves, chunk 16 and B2 below 2^64 | Omitted eligible allowance resolves to 16 MiB; explicit cap must fund programs and chains without shrinking a child | Reused PRAC/batch16; unsupported/small-cap jobs use B4 ladder | [B3 adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#b3-production-integration-acceptance-9-october-2026). Wider allocation/population confirmation remains open. |
| Bounded p−1 / portfolio p−1 | Existing finite grants | Direct gap recurrence, p−1 chunk 64 | [A6 adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p53--improve-p1-powering-and-continuation-independently): stage evidence favorable; extended integer portfolio remains Inconclusive. |

## Recent completed review

**A1b is Done: retain integrated SIQS defaults and fixed-score opt-ins.**
All 30 planned prefix calls on nine exposed moduli completed; fixed/list lost
all matched pairs at 50/55/60 digits. This bundled scoring/scheduling screen
justified no complete-factor pilot or default change. Earlier harness failures,
missing historical clocks and the explicit accounting repair remain recorded.
See the [final decision](benchmarks/RESULTS.md#a1b-existing-fixed-score-review) and
[numerical evidence](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/a1b_final_results.json). Fixed-score queue support
is a separate deferred mechanism, not a continuation of this completed study.
Evidence and tests are integrated in mainline; task worktrees/scratch were
removed after archive verification. See the report's integration/cleanup record.

## C10a-01 delivered integration

The coupled arithmetic, finite multiplier/A selection, modular root updates,
independent exceptional-root/reconstruction oracles and versioned charged
resume are integrated. The all-size default decision follows the completed
60/62/65/70 exploratory feasibility work and the deep 65-digit host-state
investigation. No arithmetic regression was found; timing and confirmation
limits remain explicit in [results](benchmarks/RESULTS.md#c10a-01-coupled-polynomial-defaults).

Research, raw receipts, failures and sealed evidence are retained locally at
`v2/audit/archives/c10a01-20261103/`. Task worktree and temporary-state cleanup
is verified against the archive manifests. Source history and frozen controls
are retained; no E1 campaign or remote push accompanies this integration.

## Next eligible items

Use the [ranked queue](ROADMAP.md#ranked-next-work--expected-benefit-and-cost)
for the ordering under the wall-time objective and declared workloads above.
Task rows retain exact scope, dependencies and fixed-session model assignments.

1. **C3a/C11 current portfolio/ECM cost-yield and missing-policy screen — Sol / xhigh.**
   Start from completed C3 and existing deeper presets; select a concrete
   60–100-digit opportunity, completed-search credit or handoff delta.
2. **B12 → F4 opt-in independent ECM workers — Sol / xhigh.**
   Prepare stable RNG/lease/restart contracts. Isolated ECM feasibility does
   not wait for E1; measure 1/2/4 workers before wider M4 occupancy.
3. **A7a current-collector parallel QS collection — Sol / xhigh.**
   Reuse A7 workers, adapting growing-DLP and current collector identities
   before measuring 1/2/4-worker complete factoring at 60/70 digits. Existing
   workers reject growing DLP and cap windows below the calibrated serial
   width. Family-sharing alternatives are a later measured subquestion.
4. **C9a combined useful-yield candidate — Astra / xhigh.**
   Ranked next after parallel collection; research need not wait for a worker
   result. Rank intake/endpoints against current collection cost; use late-run
   replays for cycle/storage/cadence hypotheses.
5. **B4a fused native ECM chunks — Sol / xhigh, optional native track.**
   Price complete setup/conversion/recovery and scalar execution before SIMD
   or new reducers. A CFFI wrapper alone is no measured gain.
6. **A1d prime-major multi-polynomial buckets — Sol / xhigh.** Current queues
   are polynomial-local; amortize prime/root traversal across a funded Gray
   batch, retaining exact coverage and charging unused preparation.
7. **B3a larger-bound chains — Sol / xhigh, profile-gated; F3 polynomial
   continuation — Astra / xhigh.** Accepted chains stop at B1 = 2,000;
   neither deeper-bound savings nor a polynomial/backend crossover is measured.
8. **A1c/B5/B5a bulk kernels — Sol / xhigh.** Select a PyPy, NumPy or CFFI
   scan/root/mark packet with queue-compatible buffers and exact hit/score
   costs; retain the Python fallback. M4 NEON and x86 AVX are separate targets.
   M4 ARM SIMD is a source/host fact, not a measured gain or default migration.

The roadmap ranks C2/F2/F7, remaining A1a and matrix preparation behind these
initial opportunities, with explicit capacity/profile triggers. E3/A10a/E4/E5/
E6 can screen a named current cheap-path or overhead gap without restarting
E1/D7; integrated tuning retains its prerequisites. B14 remains
High for requested prime certification; accounting/correctness repairs take
precedence. A1b is Done, A1b-F1 stays Deferred, C10a-01 is integrated and E1
remains Deferred. Do not restart the stopped broad backend study or duplicate the completed
C10a-01 exploratory feasibility.

E6's prime-generation candidate is a source-bundled adaptation of the existing
C11 wheel-30 sieve through CFFI, with finite batches and a Python fallback.
It must not depend on a sibling checkout or an external prime-sieve library;
see the [native prime-provider gates](ROADMAP.md#e6--self-contained-native-prime-generation).
Local source review: `v2/audit/research/prime_provider_20261102.md`.

Separate worktrees isolate source changes. Preserve C10a-01 historical timing
freezes; later comparisons need a new committed source identity. Heavy QA and
measurements share one serialized performance window. A later separately
requested execution chat may use Sol / high where its task row permits it.

B12/F4 committed arithmetic source `f9473f0` passes 1,026 tests (9 optional
skips) and full lint. Observer source `27fb6c0` adds 11 passing targeted tests
with byte-identical arithmetic/execution sources. The opt-in v5 resource and
checkpoint contract adds bounded streamed/compact metadata reuse, parent
manifest reuse, lightweight worker imports and existing reduced stage-two
kernels; mathematical v1 identities and automatic defaults stay unchanged.

The 4/6/8/10/12/14/16-worker screen is followed by a narrowed reuse/parameter
study at 60/80/100 digits. Combined settings reduce selected fresh fixed-curve
CLI medians by 18–20% versus off settings on the same source; optimized serial
still wins cheap calls. Deeper-bound compact backing improves warm prefix
replay 16–17%. Warm count candidates vary10–16, while close rankings remain
uncertain. These are exposed/enriched component inputs; first-factor,
complete-call, balanced-tail and fresh independent-input confirmation remain
NOT RUN. PRAC's existing modulus guard excludes the 100-digit cases.

No engagement threshold or default adoption is established; B12/F4 stay open.
C9a's release was explicitly confirmed and runtime uses the shared lock.
The new phase retains all failed QA, recorder and refused-allocation charges;
the earlier feasibility group-exit acknowledgement remains unknown/full charged,
with normal worker joins reported. New successful lanes drain their owned groups.
See [reuse results](benchmarks/RESULTS.md#b12f4-bounded-metadata-and-import-reuse-4-november-2026)
and local `v2/audit/b12_f4/headroom_opt_v1/HANDOFF.md`. Earlier source-pinned
profiles/counts remain under `profiles-v1/` and `counts16_v1/`. E5a records the
separate serialization follow-up; it is not pursued in this headroom phase.

## Known gaps and unresolved conflicts

- **Deferred:** E1 combined portfolio confirmation. Individual accepted
  components do not establish a compounded win. H1 and its coverage review
  still precede GNFS activation.
- **Inconclusive/unmeasured:** the adoption limitations are attached to the
  tables above. No broad SIQS crossover or practical general 60–80-digit
  completion has been established. The historical 65-digit spike timed out
  under an older collector; C10a-01 now records exploratory completions at all four sizes, with
  the 65 rerun drift failure and censored 70 control. Consume those scoped
  results before a duplicate larger-size study; there is
  no promoted 61–80 route.
- **Open:** C10a-06 runner work; R3 general cadence/storage proofs; B14 proofs;
  final SSS/worker comparisons and upstream-license clarification. See the
  [live catalogue](ROADMAP.md#numbered-phase-catalogue-and-acceptance-gates).
- **Unresolved conflicts:** none of the checked stale claims requires choosing
  between unsupported conclusions. Historical contradictory instructions stay
  in the pinned Git revision and are superseded by dated closure/current code.
- **Unresolvable references:** none among the roadmap's “local audit material”
  labels. The [history reference register](history/ROADMAP_CLEANUP.md#local-audit-reference-register)
  names existing local paths and the committed M17 baseline. Historical source
  links use the pinned `v2/ROADMAP.md` Git revision and its original context.
