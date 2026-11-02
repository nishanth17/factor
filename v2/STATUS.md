# Factor status

Last updated: 1 November 2026 (America/Los_Angeles).
Reviewed mainline/base: `4a946b8d6294d51ac5d9c36f61880cc63d5b59fd`.

**Target metric: TBD — user to define.**

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
| Exactly 50, native | Same | Growing DLP; B=30,000, 6 A factors, pool 32; 2 MiB plan; density+presieve+scan_fused; collection pacing; multiplier 1 | [C9](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c9--source-derived-large-prime-yield-and-dlp-economics-complete-and-integrated), [C10](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10--source-derived-collection-strategies-closed-bounded-tranche), [C10a-05](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open), [A1a user-directed adoption](studies/a1a_default_promotion.md). A1a's 50-digit integration gate passed. |
| Exactly 51, native | Same | Same concrete bundle as 50 | [51-digit completion adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open). A1a marking at 51 is a **user-directed, unmeasured interpolation**. |
| 52–54, native | Same | Counted base 2,500 (endpoint ceiling 250,000), 7 A factors, pool 64; scored odd multiplier (`multiplier=0`); 2 MiB plan; A1a marking; collection pacing | [C10a-04 adoption](studies/c10a04_default_promotion.md) and [A1a adoption](studies/a1a_default_promotion.md): **user-directed, unmeasured interpolated bands**. |
| Exactly 55, native | Same | Same counted bundle as 52–54 | Same [C10a-04 adoption](studies/c10a04_default_promotion.md) and [A1a adoption](studies/a1a_default_promotion.md): **user-directed**; compact training unstable, A1a fresh integration Inconclusive due to drift. |
| 56–59, native | Same | Counted base 5,000, ceiling 250,000; 8 A factors, pool 64; scored odd multiplier; 2 MiB plan; sparse marking; collection pacing | [C10a-04 adoption](studies/c10a04_default_promotion.md): **user-directed, unmeasured interpolated bands**. |
| Exactly 60, native | Same | Same counted bundle as 56–59 | [C10a-04 adoption](studies/c10a04_default_promotion.md): **user-directed**, preliminary fresh compact confirmation; original campaign incomplete. |
| Exactly 40 or 50–60, native, different explicit allowance | Caller grant; outer portfolio subtracts coexistence reserve before child selection | Same size/base/multiplier/growing policy; power plan off, sparse marking, window pacing. Finite refusal remains possible. | [C9/C10 routing scope](qs/CONFIGURATION.md#explicit-and-automatic-selection); resource variants are **unmeasured**, not additional accepted gains. |
| Other sizes, or other backends, when SIQS is selected | Reference/caller allowance | Reference size policy; swapped signs below; no automatic queue. No size-based portfolio SIQS enablement. | [Retained policy](qs/defaults.py) and [C10a-07 adoption](studies/c10a07_default_promotion.md). No measured promoted size default at 41–49 or 61+. |

Every fresh automatic SIQS configuration selects `sign_policy="swap_high"`.
Native calibrated collectors above select a 32 MiB sparse-queue ceiling within
the existing owned allowance, with charged rebuild and exact fallback for
unsupported paths or capacity refusal. This is the **user-directed C10a-07
combined adoption**. Its 50/55/60 panel is exploratory; drift acceptance was
unmet and the additional 18-call confirmation was not run. Other sizes,
backends and grants are unmeasured policy extrapolations. See the
[adoption and limitations](studies/c10a07_default_promotion.md).

Omitted calibrated SIQS grants are 10^13 work units and 1,200 wall/CPU seconds
through 55 digits, 1,800 at 56–60. All calibrated bundles use half-width 65,536,
flyer assignments and a 64-polynomial family prefix. Explicit limits win;
these grants promise no completion. The [configuration guide](qs/CONFIGURATION.md#explicit-and-automatic-selection)
explains residual/endpoint/product controls and charged base-count resolution.

| ECM / shared portfolio scope | Owned memory / eligibility | Current automatic choice | Adoption evidence and limits |
| --- | --- | --- | --- |
| Exactly 30 digits | Existing funded cap; retain previous plan if coexistence fails | compact64: `(2000, 50000, 64)` | [C3 decision](benchmarks/RESULTS.md#c3-bounded-ecm-allocation-and-handoff-10-october-2026): **user-directed before fresh confirmation**; oversized comparison Inconclusive, source screen unstable. |
| Exactly 40 digits | Same; SIQS headroom above applies | escalating140: `(200, 7700, 10)`, `(2000, 81000, 30)`, `(10000, 554000, 100)` | Same [C3 adoption](benchmarks/RESULTS.md#c3-bounded-ecm-allocation-and-handoff-10-october-2026) and limitations. |
| Other sizes, or explicitly selected paired executors | Caller/retained cap | Retained fixed32 numerical schedule | [C3 scope](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#calibrate-automatic-selection-and-handoff-c3-then-g1); not a measured optimum for every size/backend. |
| Supported native schedule with B1=2,000, at least 8 curves, chunk 16 and B2 below 2^64 | Omitted eligible allowance resolves to 16 MiB; explicit cap must fund programs and chains without shrinking a child | Reused PRAC/batch16; unsupported/small-cap jobs use B4 ladder | [B3 user-directed adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#b3-production-integration-acceptance-9-october-2026). Wider allocation/population confirmation remains open. |
| Bounded p−1 / portfolio p−1 | Existing finite grants | Direct gap recurrence, p−1 chunk 64 | [A6 user-directed adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p53--improve-p1-powering-and-continuation-independently): stage evidence favorable; extended integer portfolio remains Inconclusive. |

## Work in flight

These worktrees were present at review; presence establishes ownership, not an
experimental result. No other sessions were contacted.

| ID | One-line scope | Worktree / branch |
| --- | --- | --- |
| C10a-01 | Coupled normalized odd/mod-8 arithmetic and multiplier hypothesis | `/Users/nishanthmohan/Documents/code/factor-c10a01` / `codex/c10a01-normalized-polynomials` |
| A1b | Reconcile remaining fixed-score options after the queue decision | `/Users/nishanthmohan/Documents/code/factor-a1b` / `codex/a1b-fixed-score-review` |

## Next eligible items

Retain the backlog's priority and predecessor order; these entries schedule no
measurements. [Execution rows](ROADMAP.md#open-and-conditional-execution-items)
carry model/effort, triggers and full gates.

1. A5 Williams p+1 binary baseline — Medium, independent coverage. A8 matrix
   control is Conditional on useful solve/capacity evidence.
2. B12 worker-contract preparation — Medium; B14 independent prime certificates
   — High for certified results. Both have settled foundation predecessors;
   B14 does not wait for GNFS or D7. Other B rows need their stated triggers.
3. C9a — ranked useful-yield/target attribution, then intake/cycle weight;
   C10a — remaining triggered collection questions; C11 — High, continuing
   provenance/calibration. Existing A1a follow-ups/A1b/A1c keep their scope.
4. Downstream D7/E3–E6, selected P6 challengers, G1/H1 and I1 wait for their
   recorded predecessors. E1 remains Deferred at user direction.

## Known gaps and unresolved conflicts

- **Deferred:** E1 combined portfolio confirmation. Individual accepted
  components do not establish a compounded win. H1 and its coverage review
  still precede GNFS activation.
- **Inconclusive/unmeasured:** the adoption limitations are attached to the
  tables above. No broad SIQS crossover or practical general 60–80-digit
  completion has been established. The 65-digit spike timed out both arms;
  [explicit 61/62 feasibility](ROADMAP.md#collection-subtranches-and-retained-follow-ups)
  remains Conditional, with no promoted 61–80 route.
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
