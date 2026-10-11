# C9.2 allocator growth: protocol before new training

This supersedes C9.1's storage assumption at the user's explicit direction.
C9.1's protocol, captures, pins and stopped confirmation remain historical.
All seven C9.1 fresh inputs are now training/regression evidence, never fresh
confirmation for this implementation. No C9.1 performance promotion occurred.

## Implementation and question

Follow full SIQS implementations' growable relation stores: retain forest
edges and original atoms as needed, using PyPy's growable containers. There
is no implicit edge, atom, row, total-memory, checkpoint-byte or cumulative
split-call quota in `GrowingDoubleLargeSieveConfig` plus standalone
`SIQSConfig`. Explicit integer quotas remain supported. Legacy configurations
and schema 3/4 keep their fixed behavior; growing execution uses schema 5.
The outer finite portfolio still owns and resolves its child's memory share.
Per-residual effort, exact endpoint proofs and the 256-atom cycle policy stay
unchanged. This is a storage policy implementation, not a claim of optimality. Snapshot
encoding retains elapsed wall/CPU time even after a stopped computation budget;
restore charges verification/replay work, as with legacy emergency snapshots.

Compare complete factoring at the same resources. Removing a FIFO cutoff may
allow delayed graph structure, but preserving rare high endpoints can waste
memory and preparation time. Measure full-base growth first; then common
endpoint admission, the historical half-base bundle, and reusable SLP matching.
No new multiplier, polynomial, score, root, factor-base construction or matrix
algorithm is introduced. All-component cycles and original provenance remain.

## Finite study resources, distinct from production defaults

- Retain the original **14,400 active wall and CPU seconds** total. Charged
  before C9.2: **4169.706814 wall / 4124.095798 CPU seconds**, including stopped
  confirmation, a full interrupted-call reservation and 600 prepaid QA seconds.
  Research, ordinary edits and lock waiting are excluded. QA receipts reduce
  the prepaid allowance, never add free new compute; any excess is charged.
- Each 50-digit call: `10**13` work, **360 wall/CPU seconds**, **1 GiB owned
  memory**, **3 GiB process high-water RSS**, 1 MiB checkpoint configuration.
  The 1 GiB bound is an explicit common *experimental* allowance, not the new
  collector default. Legacy controls receive it too. Count quotas remain their
  declared policy difference. Record early quota refusals; never rescue a run.
- Cumulative attribution observations: **60 / 180 / 360 seconds**, or earlier
  terminal outcomes. Prefix inspection is charged. No early stop for an empty
  matrix. Track graph endpoint degrees/core/rank, original atoms, forest edges,
  emitted rows, filter survival, nullity when available, dependencies and factors.
- Eight initial training calls: the same two C1 50-digit training inputs,
  seeds 7/29 respectively; calibrated SLP, current accepted full-base DLP,
  growing full-base DLP, growing DLP with SLP endpoint aligned to DLP's bound.
  Alignment preserves the DLP scoring/product domain. Rank by completion,
  capped total cost, then declared order. One call per input is attribution,
  not accepted timing or an estimate of sampling noise.
- Advance only if a long-prefix arm completes or produces surviving filtered
  rows. Under a new explicit joint freeze, four calls total: historical C1
  half-base DLP plus the better growing admission; reusable SLP via the existing
  `graph_slp` construction with growth and the calibrated SLP score domain.
  C1 already supplied the reusable graph path; no new matching algorithm is
  claimed. Rank SLP and DLP candidates separately by complete capped costs.
- At most 360 seconds additional diagnostic/profile/checkpoint work, 120 for
  new certified corpus generation and the existing prepaid QA allowance.
  Checkpoint diagnostics may use 16 MiB checkpoint bytes inside the same 1 GiB
  total, with at most 64 MiB checkpoint files. They cannot replace timing calls.
  All captures together <=256 MiB. Actual RSS or resource failure stops the
  call/campaign; preserve censored outcomes. No 60–100-digit study follows.

## Confirmation and acceptance

Freeze selected source, configurations, runners and input generation before
new inputs are generated. Generate three independent balanced certified
50-digit and three 40-digit inputs, plus one 30-digit regression input, with
seed 202610110929. Reject all n already in versioned corpora. Pocklington
construction conditions primes on a large p−1 factor; this is not RSA sampling.
Known factors/certificates enter independent validation only.

Compare selected calibrated SLP, current C1 DLP and selected growing challenger.
At 30/40 digits the challenger changes storage on C1's half-base DLP bundle;
at 50 it uses the selected growth bundle. All receive 1 GiB owned memory,
3 GiB RSS guard, `10**13` work and identical wall/CPU caps: 5/30/360 seconds.
Observe **nine matched samples** per class (three inputs x seeds 7/29/47 at
40/50; one input x nine repeated samples at 30). A repeated seed is not an
independent input. Each arm/band gets >=3 wall and CPU seconds of validated
PyPy 3.11 warmup in the same measurement process. Cold processes are separate.
Interleaved deterministic randomized arm order limits drift.

Report paired capped-total-cost effects and 10,000-resample hierarchical
bootstrap 95% intervals over inputs then paired samples. One predefined
extension to 18 matched samples is requested for within-input sign reversals,
>10% repeated-cell first/last drift or relative IQR (when >=4 repeats). This
is a stability rule, not a minimum acceptable benefit. If the total study cap
cannot reserve the whole extension, keep the captures and declare timing
inconclusive. No tuning on fresh results. A band cannot be accepted below nine
matched samples. An interrupted process requires fresh same-process warmup.

Any sustained complete-run improvement may qualify: require a positive 95%
interval, correctness/resource contracts, acceptable complexity, and <=5
percentage points completion regression in each promoted class. Narrow scope
if necessary. User-requested adaptive storage can remain an explicit supported
policy even without a speed promotion; distinguish those decisions. Increasing
memory or duration alone is not a speedup. E1 combined acceptance, C10 broader
collector calibration, F6 higher arity and matrix redesign remain open.

All heavy work uses the repository machine-wide performance lock and owner
file. Every output, including unresolved cofactors, is independently checked.
Run full tests/lint and committed-files-only tests/imports. Preserve failed
attempts and reserve interrupted call costs. New controls/corpus are versioned;
raw captures and profiles stay ignored. No merge or push.

Machine leases may be interleaved with C3 in <=600-active-second segments.
Reserve the entire next call before starting it; yield if it cannot fit.
Resume immutable captures without selecting from partial samples. Confirmation
warms every arm again in each new measurement process and records that process
with each sample. Lock waiting does not consume the active study allowance.
