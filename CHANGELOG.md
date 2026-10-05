# Changelog

## Current development

- Repair P2 preprocessing/cursor overhead, sparse QS exponent recovery,
  skipped-prime score refinement, repeated prime-power lifting inverses,
  sparse matrix label storage and dependency traversal.
- Amortize SSS collision accounting while preserving seeded candidates and
  exact completed-assignment work; use disjoint SSSf smoothness bases.
- Reuse verified QS preparation only for identical immutable rows, bases and
  provenance atoms under a bounded collector-owned cache. Keep public and
  checkpoint verification exact, and amortize native collector/solver polling
  with forced stage-boundary checks and the original shared work ledger.
- Add bounded family/Gray/block worker assignments and checkpoint version 2,
  retaining legacy complete-family resume. Poll external limits at explicit
  intervals of at most 64 atomic actions while checking every work reservation.
  Preserve central exact verification, aggregate CPU and charged replay.

- Add an experimental coarse SIQS family executor for serial, thread and
  spawned PyPy process comparisons. Reverify bounded relation batches in a
  central store, preserve assignment identity across worker counts, charge
  parent-owned work leases and aggregate cooperative CPU, and drain cancelled
  work before returning. Add capped provenance checkpoints and charged replay
  of unfinished families. Fixed-work runs continue after direct residual splits.
- Retain serial collection after the P3.6 held-out experiment: worker arms do
  not pass the complete-factor promotion gate. Keep the standalone worker API
  experimental; medium-band batch-cap refusals remain explicit and resumable.

- Add an experimental bounded SSS/SSSf collector with seeded CRT/collision
  assignments, exact product/remainder-tree smoothness detection, complete
  exponent recovery and the shared QS relation/filter/extraction pipeline.
  Preserve candidate-prefix resume, finite storage and quiet calls.
- Add separate hash-checked upstream reproduction and trained held-out
  challenger runners. Optional comparison dependencies remain outside the
  factoring library; SSS has no automatic dispatcher promotion.
- Expose opt-in `--method sss` / `sssf` and `PortfolioConfig(sss=...)` under
  shared budgets and recursive result validation. Add capped SSS checkpoints
  with charged store/assignment/matrix reconstruction; retain version-2/3
  portfolio resume through version 4. Existing automatic defaults stay intact.
- Validate complete SSS/SSSf splits on the declared 30-digit held-out set.
  Keep larger capped runs diagnostic and retain the challenger pending a
  comparison with feasible trained SIQS; see the benchmark guide for limits.

- Preserve the original Python 2 implementation in v1; develop and test v2
  on PyPy implementing Python 3.11.
- Repair arithmetic boundaries, primality certainty, factor reconstruction,
  segmented sieve endpoints and saturated-batch recovery.
- Add bounded preprocessing and seeded rho/p−1/ECM portfolio execution with
  serializable checkpoints, shared allowances and explicit unfinished results.
- Add exact QS/MPQS polynomial and relation identities, bounded single-large-
  prime collection, provenance-aware filtering, GF(2) dependencies and modular
  factor extraction with in-memory resume.
- Remove redundant collector budget polling and proved-useless score marking.
  Keep experimental resieving and alternative scoring/filtering controls
  optional; see [measurement scope](v2/benchmarks/README.md).
- Keep code, tests, roadmap, research and required fixtures on GitHub.
  Exclude generated run captures, logs, profiles and local development records.

- Attempt factor extraction from retained QS relations before returning a
  storage-limit outcome; preserve resumable pending work and idempotent stops.
- Account for simultaneous QS collector, relation-preparation and matrix
  workspace before allocating private preparation storage.

- Add bounded experimental SIQS CRT/Gray families and incremental roots,
  with exact cached-root certification and compact family-only checkpoints.
- Complete bounded SIQS with shared verified relations, compact full-job
  checkpoints, charged root/matrix/prefix reconstruction and finite recovery.
- Add integer multiplier scoring, QS/MPQS/simple-family controls and optional
  ECM-to-SIQS dispatch under one allowance. Preserve version-2 portfolio resume.
- Validate disjoint certified inputs, including declared 30–80-digit capped
  comparisons. Keep ECM as the default: no broad SIQS crossover is established.

- Complete P3.4's bounded large-input evaluation with an independent certified
  30–80-digit and varied-shape corpus, trained frozen settings, cumulative
  checkpoint continuations, and declared negative outcomes.
- Add optional exact prime-power score marking, immutable base-column lookup,
  sparse matched-relation accounting and larger explicitly bounded capacities.
- Replace repeated global matrix rescans with exact deterministic incidence
  updates; the inspected 30-digit filtering control improves 5.920 → 3.425 s.
- Save wide binary dependency masks with compatible streamed hexadecimal
  digests and retain legacy checkpoint prefixes. A fresh balanced 50-digit
  SIQS run completes in 1,183.355 s; continuation time is recorded separately.

[P3.4](v2/audit/TODOS.md)'s bounded implementation/evaluation is complete;
ECM remains the default and SIQS stays opt-in. A general SIQS crossover and
practical 60–80-digit scaling remain unestablished. The earlier M30 completion
claim was corrected because 0.2-second probes and 24–26-bit fixtures were
insufficient; M31 supplies the declared longer experiment and its explicit
retain-baseline decision. See [measurements](v2/benchmarks/README.md#completed-larger-evaluation-and-filtering-repair-m31-4-october-2026).
P4.3 arithmetic assessment and bounded GNFS remain roadmap work.

PRAC optimization remains [P4.1](v2/audit/TODOS.md); production currently uses
the safe Montgomery ladder.

## P3.8 R1 isolated capacity implementation — 4 October 2026

- Add opt-in streamed SIQS nearest/flyer A assignments with bounded resident
  state, up to 32 A factors, and separate cumulative search and Gray quotas.
- Add exact external-square MPQS coefficients and preserve square corrections
  through atomic/combined relations, extraction and checked checkpoints.
- Add explicit monotone direct-job resume extension without resetting checked
  stores or consumed resources; distinguish exhausted assignment space.
- Add actual-base capacity diagnostics and sparse initial-incidence work
  charging while retaining conservative matrix storage bounds.
- Add independently certified sub-100-digit training/confirmation generation,
  frozen-control checks and matched full-call/cold/profile runners. Experimental
  defaults and overall R1 completion remain gated on confirmation/integration.
- Frozen R1 confirmation completes the fresh 30-digit control in all four
  arms (18/18 attempts each), while all balanced upper bands remain censored
  under one-second caps. Nearest/flyer/reference controls, a 3.0% legacy
  regression, cold lifecycles and separate profiles are recorded without
  changing defaults or claiming a broad crossover.
- Preserve the exact measured source snapshot and a bounded restoration tool;
  account for verification-cache growth in explicit quota extensions.
- Integrate R1 with the P2/P3.6.1 repairs under a single writer. The combined
  candidate passes 258 tests and full lint; broad upper-band feasibility and
  crossover work remain open, with no automatic parameter promotion.


## P3.8 R3 isolated preparation and provenance work — 4 October 2026

- Preserve full and matched relations in one stable admission sequence; bind
  prepared identities to complete immutable base, polynomial and atom payloads.
- Retain mixed row order in SIQS/SSS version-2 and parallel version-3
  checkpoints. Read preceding supported formats with their original grouped
  order and fully reverify every relation before charged solver replay.
- Count pivot nonzeros incrementally. Preserve exact original-row dependency
  lifting and conservative fill/provenance workspace reservations.
- Add opt-in `filter_row_growth` and `tested_dependencies` SIQS/QS settings,
  retaining every-change filtering and disabled tested-dependency caching as
  defaults. Pending work resumes immediately; terminal windows/storage limits
  force the remaining extraction. Skipped checked trials still count against
  the existing trivial-dependency allowance.
- Add independent certified fixtures, immutable controls and bounded matrix,
  compaction, history, cache and complete-factor experiments. Decisions and
  confirmation scope are recorded in the benchmark guide; matrix scaling and
  conditional packed-exponent/root work remain open.

- Confirm the scoped opt-in cadence-32 setting on fresh 30-digit fixtures:
  38.7% lower complete-cohort median than R3 cadence 1, 54/54 complete per arm.
  Keep cadence 8 and tested-dependency cache promotions deferred; neither this
  experiment nor passing recovery tests changes automatic dispatch defaults.

- Integrate only the R3 delta with the repair/R1 work. The combined candidate
  passes 273 PyPy tests, lint and all 44 benchmark imports from committed inputs;
  a stable matched bridge validates 36 complete outcomes in every arm.
