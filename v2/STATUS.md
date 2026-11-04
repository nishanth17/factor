# Factor status

Last updated: 4 November 2026 (America/Los_Angeles). Reviewed mainline:
`963d9cd8dd4f8a99c5ec8d093b0c28bad080ba2d`; acceleration-first planning revision
changes no implementation, defaults or active experimental grants.

Read this first, then the [live backlog](ROADMAP.md) and [stable policy](POLICY.md).
[Results](benchmarks/RESULTS.md) hold numerical evidence; the
[completed register](ROADMAP.md#completed-register) holds delivered milestones.

## Target metric and workload

Minimum complete-factorization wall time on the Apple M4 / 24 GiB / 10-core
machine. All cores are available; report total CPU, peak memory and startup/
cancellation cost. Use equally weighted uniformly sampled 60/80/100-digit
integers, with median/p90 time and completion per band. Balanced-semiprime
completion at 60/70 is a must-not-regress constraint. Feasibility at 100 is
unestablished. Optional native backends retain a working pure-Python/PyPy
fallback. These declared objectives and censoring rules are in
[policy](POLICY.md#scope-and-standing-directions).

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

## Work in flight

| ID | Current scope and state | Worktree / branch |
| --- | --- | --- |
| B12/F4 | **Open; ongoing.** ECM worker ownership, RNG/leases, diagnostics and latency evaluation. Completion and integration are not established by worktree presence. Preserve its active protocol and resource grants. | `/Users/nishanthmohan/Documents/code/factor-b12-f4` / `codex/b12-f4-ecm-workers` |

C9a and C10a-01 are delivered mainline scopes; they are not active worktrees.
No other ongoing work is inferred from an old preparation note.

## Recent completed review

- **C9a — Done, 4 November:** delivered and integrated; policies remain opt-in.
  Complete60/70 evidence is exploratory; one exposed input per size and no
  disjoint fresh confirmation establish no population speedup.
  [Completion and evidence](ROADMAP.md#c9a-completion).
- **C10a-01 — Done, 3 November:** integrated coupled policies are adopted for
  fresh automatic SIQS before fresh confirmation. The failed 65 canary,
  censored 70 control and unmeasured all-size gains remain explicit.
  [Completion and evidence](ROADMAP.md#c10a-01-completion).
- **A1b — Done, 2 November integration:** retain defaults after the fixed/list
  screen lost at 50/55/60; no complete-factor claim.
  [Final decision](benchmarks/RESULTS.md#a1b-existing-fixed-score-review).

<a id="c10a-01-delivered-integration"></a>

## Next eligible items

The [ranked queue](ROADMAP.md#ranked-next-work--expected-benefit-and-cost) now
prioritizes current acceleration and shorter later experiment turnaround:

1. **B12/F4:** finish the active ECM worker stream — Sol / xhigh.
2. **A7a:** current-collector parallel QS preparation — Sol / xhigh.
3. **B4a:** coarse ECM C/CFFI, then selected native SIMD/threads — Sol / xhigh;
   separately selected new arithmetic proofs retain Astra / xhigh.
4. **B5a:** QS C/CFFI and ARM NEON, with A1c/B5 controls — Sol / xhigh.
5. **E6:** bundled native prime batches/reusable-context opportunity screen —
   Sol / high screen, Sol / xhigh ABI and implementation.
6. **B15:** native helper/SIMD opportunity selection across E4/A10a/A2a/C5 —
   Sol / high evidence review; Sol / xhigh for boundary/protocol design. Deliver
   at most one justified candidate or a no-go; implementation retains its owner.
7. **C3a/C11:** policy/handoff screen on the settled worker cost identity —
   Sol / xhigh. Reuse accelerator receipts rather than repeating broad studies.

Then A1d, B3a, F3, broader C9b and the remaining ranked rows. F3 is the first
substantial Astra / xhigh algorithmic research choice; a concrete opportunity
can start before every conditional native port is settled.

A7a can prepare alongside B12 with pinned shared contracts; timed execution
uses the shared lock sequentially. Assess warm ECM→QS process reuse and parent
saturation. B4a/B5a need stable neighboring interfaces, not a prior worker win.
E6/B15 at ranks 5–6 are conditional screens. B15 is not a C3a prerequisite;
its selected helper implementation does not inherit rank 6 automatically.
Parallel latency gains do not imply lower CPU cost or higher batch throughput;
native-kernel gains and M4 SIMD remain unmeasured. x86-specific ports are outside
this priority increase. Small cost/A/A checks belong in each task; C9b is not
a gate on preparation. The local pressure test is
`v2/audit/research/acceleration_priority_20261104.md`.

This order is an expected-benefit judgment, not measured gains or campaign
authorization. Correctness/accounting repairs, C10a-06 and requested B14 proof
work take precedence. E1 remains Deferred and integrated tuning retains its
prerequisites. Closed C9a/C10a-01 evidence and active grants are unchanged.

## Known gaps and unresolved conflicts

- **Deferred:** E1 combined confirmation; GNFS until H1 and its coverage review.
- **Unmeasured/Inconclusive:** no general all-size speedup or broad crossover
  is established. Coupled SIQS policies apply at every size, but tuned size
  presets and automatic portfolio SIQS eligibility still end at 60. Other
  sizes/backends/grants retain the stated extrapolation and feasibility limits.
- **Open:** C9b timing variability and fresh confirmation; C10a-06 runner work;
  general R3/provenance and B14 proof scopes; current-worker/SSS comparison and
  upstream-license clarification. B12/F4 remains ongoing.
- **Resolved:** C9a's historical preparation/campaign instructions are superseded
  by scoped closure. Its uncalibrated 15% canary flags remain historical warnings;
  future judgments require measured variability and uncertainty.
- **Unresolved conflicts:** none found in the checked current behavior/closure
  records. No conflicting experimental verdict was upgraded.
- **Local references:** earlier C9a preparation files are preserved under
  `v2/audit/c9a/archive-20261104/base/c9a/`, with former worktree roots recorded
  in `v2/audit/c9a/archive-20261104/relocations.json`. Detailed chronology and prior
  quantitative records remain in [pinned history](ROADMAP.md#documentation-history-and-local-evidence).
  No missing audit reference was discarded.
