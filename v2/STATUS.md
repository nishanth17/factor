# Factor status

Last updated: 2 November 2026 (America/Los_Angeles).
Reviewed mainline/base: `b0813230a05396874c41c649fa01738484a1f55e`.
A1b review record: `61899a9`; production defaults are unchanged.
QS roadmap review: `4d65d33`, consolidated into existing task scopes on
2 November; implementation remains in the C10a-01 worktree.
The expanded 2 November ECM/QS literature and source review adds A1d
multi-polynomial buckets and conditional F9 genus-two/Kummer research to the
ECM/acceleration amendments. It changes no defaults and launches no runs.
Local supporting reviews: `v2/audit/research/ecm_qs_review_20261102.md` and
`v2/audit/research/ecm/roadmap_review_20261102.md`.

**Target metric: TBD.**

Read this first, then the [live backlog](ROADMAP.md) and [stable policy](POLICY.md).
Full historical results and gates remain in the [history index](history/README.md).
The current objective includes useful progress toward 60–80-digit factoring;
practical general completion there remains unestablished.

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
| Exactly 40, native | 1 GiB child; omitted outer portfolio adds 16 MiB | Growing DLP; B=5,000, 5 A factors, pool 32; 2 MiB power plan; sparse marking; window pacing; multiplier 1 | [C10 plan adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10--source-derived-collection-strategies-closed-bounded-tranche). Growing versus fixed DLP was not separately confirmed at 40. Face/queue adoption below applies. |
| Exactly 50, native | Same | Growing DLP; B=30,000, 6 A factors, pool 32; 2 MiB plan; density+presieve+scan_fused; collection pacing; multiplier 1 | [C9](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c9--source-derived-large-prime-yield-and-dlp-economics-complete-and-integrated), [C10](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10--source-derived-collection-strategies-closed-bounded-tranche), [C10a-05](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open), [A1a adoption](benchmarks/RESULTS.md#a1a-integration-and-extension). A1a's 50-digit integration gate passed. |
| Exactly 51, native | Same | Same concrete bundle as 50 | [51-digit completion adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open). A1a marking at 51 is an **unmeasured interpolation**. |
| 52–54, native | Same | Counted base 2,500 (endpoint ceiling 250,000), 7 A factors, pool 64; scored odd multiplier (`multiplier=0`); 2 MiB plan; A1a marking; collection pacing | [C10a-04 adoption](benchmarks/RESULTS.md#c10a-04-larger-band-routing) and [A1a adoption](benchmarks/RESULTS.md#a1a-integration-and-extension): **unmeasured interpolated bands**. |
| Exactly 55, native | Same | Same counted bundle as 52–54 | Same [C10a-04 adoption](benchmarks/RESULTS.md#c10a-04-larger-band-routing) and [A1a adoption](benchmarks/RESULTS.md#a1a-integration-and-extension): **Adopted**; compact training unstable, A1a fresh integration Inconclusive due to drift. |
| 56–59, native | Same | Counted base 5,000, ceiling 250,000; 8 A factors, pool 64; scored odd multiplier; 2 MiB plan; sparse marking; collection pacing | [C10a-04 adoption](benchmarks/RESULTS.md#c10a-04-larger-band-routing): **unmeasured interpolated bands**. |
| Exactly 60, native | Same | Same counted bundle as 56–59 | [C10a-04 adoption](benchmarks/RESULTS.md#c10a-04-larger-band-routing): **Adopted**; preliminary fresh compact confirmation; original campaign incomplete. |
| Exactly 40 or 50–60, native, different explicit allowance | Caller grant; outer portfolio subtracts coexistence reserve before child selection | Same size/base/multiplier/growing policy; power plan off, sparse marking, window pacing. Finite refusal remains possible. | [C9/C10 routing scope](qs/CONFIGURATION.md#explicit-and-automatic-selection); resource variants are **unmeasured**, not additional accepted gains. |
| Other sizes, or other backends, when SIQS is selected | Reference/caller allowance | Reference size policy; swapped signs below; no automatic queue. No size-based portfolio SIQS enablement. | [Retained policy](qs/defaults.py) and [C10a-07 adoption](benchmarks/RESULTS.md#c10a-07-face-and-queue). No measured promoted size default at 41–49 or 61+. |

Every fresh automatic SIQS configuration selects `sign_policy="swap_high"`.
Native calibrated collectors above select a 32 MiB sparse-queue ceiling within
the existing owned allowance, with charged rebuild and exact fallback for
unsupported paths or capacity refusal. The combined C10a-07 policy is adopted.
Its 50/55/60 panel is exploratory; drift acceptance was unmet and the additional
18-call confirmation was not run. Other sizes,
backends and grants are unmeasured policy extrapolations. See the
[adoption and limitations](benchmarks/RESULTS.md#c10a-07-face-and-queue).

Omitted calibrated SIQS grants are 10^13 work units and 1,200 wall/CPU seconds
through 55 digits, 1,800 at 56–60. All calibrated bundles use half-width 65,536,
flyer assignments and a 64-polynomial family prefix. Explicit limits win;
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

## Work in flight

These worktrees were present at review; presence establishes ownership, not an
experimental result. No other sessions were contacted.

| ID | One-line scope | Worktree / branch |
| --- | --- | --- |
| C10a-01 | Checked opt-in polynomial/multiplier/A-population/root-table bundle; v4 calibration and current 60/62/65/70 feasibility stay under their versioned freezes. Code is not integrated into mainline; see the [current scope](ROADMAP.md#c10a-01--coupled-polynomial-multiplier-a-population-and-root-setup). | `/Users/nishanthmohan/Documents/code/factor-c10a01` / `codex/c10a01-normalized-polynomials` |

## Next eligible items

Use the [ranked queue](ROADMAP.md#ranked-next-work--expected-benefit-and-cost)
for the current cross-engine ordering and its workload/uncertainty assumptions.
Task rows retain exact scope, dependencies and fixed-session model assignments.

1. **C3a/C11 current ECM cost/yield and missing-policy screen — Sol / xhigh.**
   Start from completed C3 and existing deeper presets; select a concrete
   60–99-digit opportunity, completed-search credit or handoff delta.
2. **B12 → F4 opt-in independent ECM workers — Sol / xhigh.**
   Prepare stable RNG/lease/restart contracts. Isolated ECM feasibility does
   not wait for E1; later combined acceptance retains it.
3. **C9a combined useful-yield candidate — Astra / xhigh.**
   Still the first QS task. Rank wider intake/endpoints against the current
   collection cost; use late-run replays for cycle/storage/cadence hypotheses.
4. **A1d prime-major multi-polynomial buckets — Sol / xhigh.** Current queues
   are polynomial-local; amortize prime/root traversal across a funded Gray
   batch, retaining exact coverage and charging unused preparation.
5. **B3a larger-bound chains — Sol / xhigh, profile-gated; F3 polynomial
   continuation — Astra / xhigh.** Accepted chains stop at B1 = 2,000;
   neither deeper-bound savings nor a polynomial/backend crossover is measured.
6. **A1c/B5/B5a bulk kernels — Sol / xhigh.** Compare PyPy, optional NumPy and
   separately selected C/CFFI tracks with exact hit/score and boundary costs.
   M4 ARM SIMD is a source/host fact, not a measured gain or default migration.

The roadmap ranks C2/F2/F7, remaining A1a and matrix preparation behind these
initial opportunities, with explicit capacity/profile triggers. B14 remains
High for requested prime certification; accounting/correctness repairs take
precedence. A1b is Done, A1b-F1 stays Deferred, C10a-01 stays in flight and E1
remains Deferred. Do not restart the stopped broad backend study or duplicate
C10a-01 feasibility under this review.

Separate worktrees isolate source changes. Keep current C10a-01 timing freezes
unchanged; parallel research/design and light edits can proceed, while heavy QA
and measurements share one serialized performance window. A later separately
requested execution chat may use Sol / high where its task row permits it.

## Known gaps and unresolved conflicts

- **Deferred:** E1 combined portfolio confirmation. Individual accepted
  components do not establish a compounded win. H1 and its coverage review
  still precede GNFS activation.
- **Inconclusive/unmeasured:** the adoption limitations are attached to the
  tables above. No broad SIQS crossover or practical general 60–80-digit
  completion has been established. The historical 65-digit spike timed out
  under an older collector; current C10a-01 60/62/65/70 feasibility is in flight.
  Consume its scoped results before a duplicate larger-size study; there is
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
