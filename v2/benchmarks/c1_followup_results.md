# C1 longer residual feasibility: evidence in progress

This report supplements the preserved [initial screen](c1_results.md).
It does not repeat its withdrawn final defer claim. The user-authorized
[follow-up protocol](c1_followup_protocol.md),
[observation addendum](c1_followup_resolution.md), and
[research reconciliation](c1_research.md) define the current experiment.
Production DLP and performance promotion have not yet passed their gates.

Collection runs source commit `204adc2` on PyPy implementing Python 3.11,
under the machine-wide performance lock. The driver loaded that source before
the subsequent benchmark-only cost-audit edits. Local captures are immutable
in `results/c1/followup/`, with input/source hashes and raw failed outcomes.
Each band uses two distinct certified training inputs; no inspected input is
fresh confirmation. SLP controls use B1's 40-digit calibration, a declared
50-digit interpolation, and its uncalibrated 60-digit diagnostic setting.

## Completed 40- and 50-digit cells

The following are diagnostic observations, not accepted timing ratios. DLP
products are bounded by 128 B², endpoints by 100B; the nested 64 B² cohort
uses the same positions. SLP retains its calibrated B² prime cap.

| Band/input | Blocks collected | Full / SLP / DLP residuals | Composite splits | Endpoint rejects | Failed splits | Terminal reason |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 40/0 | 3,762 | 300 / 13,939 / 7,432 | 7,481 | 49 | 0 | SLP factor |
| 40/1 | 4,719 | 264 / 15,937 / 8,687 | 8,753 | 66 | 0 | SLP factor |
| 50/0 | 34,896 | 678 / 44,783 / 20,076 | 20,212 | 136 | 0 | 65,536-record cap |
| 50/1 | 40,760 | 729 / 44,635 / 20,173 | 20,319 | 146 | 0 | 65,536-record cap |

The last classified candidate in capped cells is counted but cannot be
retained. Those final prefixes are censored. No split failure is reclassified
as a proven non-semiprime. Sampling retains 32,768 uniformly selected
positions per cell and explicitly ends near position 8.4 million; no
threshold-loss estimate is extrapolated to the later unsampled positions.

At 40 digits, policy 128 already yields independently verified complete
factors at retained prefixes 2,087 and 2,391 blocks, while same-position
SLP-only incidence has no dependency. These are directly checked witnesses,
not interpolated minimum collection times. Both original SLP controls finish,
in single diagnostic runs of approximately 4.7 and 6.1 seconds. Fewer
positions does not establish a faster complete factorization.

At 50 digits, input 1 has a verified policy-64 factor at 32,768 blocks while
SLP has no dependency. Retained refinement also finds a policy-128 factor at
23,552 blocks with no SLP factor. Input 0's coarse DLP analyses hit the local
exhaustive-validation deadline; unknown useful yield must not be recorded as
zero. The SLP controls respectively complete within their 120-second
allowance and retain the whole second input at timeout.

The independently reconstructed factorizations observed so far are:

- 40/0: 49,230,173,820,498,272,561 × 69,260,785,983,138,815,257.
- 40/1: 46,700,257,282,435,804,589 × 68,915,328,265,279,712,893.
- 50/0: 6,234,745,488,095,579,115,919,703 × 9,192,770,129,255,987,774,468,263
  (SLP evidence so far).
- 50/1: 6,629,519,674,062,239,063,796,497 × 8,334,069,097,930,144,147,346,597.

Runtime certainty labels and independent corpus certificates are retained
separately. Every extracted dependency reconstructs original exponents and
square corrections, checks the square congruence and both GCD signs, and
returns only proper divisors. Every unresolved outcome still reconstructs n.

## Cost-attribution issue and bounded repair

Policy-128 splitting plus certification costs approximately 0.879/0.975
instrumented CPU seconds in the entire 40-digit censuses and 3.196/3.258
seconds in the entire 50-digit censuses. This identifies affordable residual
splitting under the measured bounds; it excludes the extra candidate division,
production graph retention, filtering, provenance, matrix and extraction costs.

The literal frozen gate additionally charges a candidate for SLP and both
DLP policy analyses, repeated refinement, and exhaustive extraction after a
complete factor is already found. It will be reported unchanged. It is not a
measurement of one candidate's complete-factor pipeline. A separately frozen
[post-observation attribution audit](c1_cost_attribution.md), implementation
`df55a0b`, uses only retained records, at most 180 seconds and the unspent
part of the original 4,500-second study envelope. No new collection, parameter
search or bound expansion is authorized. Its post hoc status is explicit.

That audit charges all candidate graph/provenance/filter/matrix costs and
stops extraction after a verified complete factorization. Generated algebraic
kernel masks and independently checked square congruences are reported
separately. Whole-study exhaustive validation/search overhead is also retained.
An exploratory investment decision is not a performance promotion: a viable
production option still needs bounded ownership/eviction, independent complete
cycle oracles, charged cancellation/resume and fresh matched total-resource
confirmation after training freezes.

The 60-digit cells and attribution audit remain pending at this revision.

## Completed collection and investment decision

The six-cell collection finished in3,002.638 wall /2,970.117 CPU seconds.
Its literal frozen gate returns no qualifying policy; that verdict is preserved
in the original decision capture. The accounting audit then took24.869 wall /
24.773 CPU seconds. Even conservatively reserving all180 seconds originally
allowed for the earlier retained-record analysis, the complete study remains
within4,500 seconds. Twenty targeted PyPy3.11 tests passed in0.183 seconds.

The corrected attribution witnesses are:

| Input | Directly verified prefix | SLP constraints | DLP constraints | Charged split/certification + complete candidate processing CPU upper bound | Allowance |
| --- | ---: | ---: | ---: | ---: | ---: |
| 40/0 | 2,087 blocks | 303 | 633 | 1.474s | 15s |
| 40/1 | 2,391 blocks | 242 | 615 | 1.410s | 15s |
| 50/0 | 19,808 blocks | 735 | 1,630 | 5.132s | 60s |
| 50/1 | 23,552 blocks | 749 | 1,749 | 5.194s | 60s |

Every listed DLP prefix completes and the same-position SLP has no dependency.
The50/0 DLP factors match the independently certified factors listed above.
Forty-digit factors receive runtime PROVEN labels;50-digit factors receive
PROBABLE labels from the runtime, with separate corpus certificates. No label
was silently upgraded. Every generated algebraic dependency was checked against
original rows; every extracted square congruence was checked independently.
Unextracted algebraic masks are counted separately from verified squares.

Both60-digit runs remain unresolved after their longer collection allowances.
At106,548/97,258 blocks they retain12,179/13,699 DLP residuals. Policy128
increases LP-cancelled rows from534 to591 and628 to728 respectively, but all
rows are singleton-filtered away. At the terminal128 prefixes,45,543 of50,100
and51,201 of56,560 vertices have degree one. Split failures29/26 remain
unresolved; bounds and deadlines were not expanded. This diagnoses poor
useful density at these fixed parameters, not universal DLP infeasibility.

**Decision: proceed with bounded opt-in implementation at40/50.** This is an
explicitly post-observation investment decision after repairing cost attribution,
not a retrospective claim that the literal preregistered gate passed and not
an accepted speedup. The [implementation/confirmation protocol](c1_implementation_protocol.md)
freezes the next finite stage before production edits. Sixty-digit calibration,
80–100-digit claims and default promotion remain unpassed.
