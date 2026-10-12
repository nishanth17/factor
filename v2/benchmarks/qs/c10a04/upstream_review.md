# C10a-04 current-upstream review

Audit observed 12 October 2026 UTC, following the user's challenge about
research depth. This supplements the Stage 1 freeze
`5dab76ec66c02049095ddf2f75e123f19d0d6627`. It does not change its
research/protocol hashes, candidate arms, inputs, allowances or decision rules.
The original survey traced real executable sources, but reused pinned
snapshots without checking every upstream head. Describing it as an exhaustive
current state-of-the-art survey would overstate the work.

The [audit manifest](../../inputs/controls/c10a04_upstream_review.json)
records observed heads, comparison URLs and file digests. Raw API responses,
diffs and source captures remain ignored under
`v2/benchmarks/results/c10a04/upstream-audit/`. Public native sources were
read, not compiled, benchmarked or incorporated into v2.

## Current heads and material changes

| Source | Observed revision | Finding |
| --- | --- | --- |
| YAFU | `8110dfbd8c6f9486d93b1a02de6eb7b180e55a80` | Still matches the freeze. |
| Alpertron | `93c5c8189cb2f149c8dae898eb11997c5f2e7980` | Still matches the freeze. |
| Yamaquasi | `181e282cbfcace30bd0d99c2b4bae181dee7ba69` | Three subsequent commits; reviewed diff contains Python conversion, Clippy and build/platform maintenance, with no change to the calibrated joint policy. |
| FLINT | `26b9bc9df973206db15616d0849a4c8be9d22059` | Material SIQS update missed by the frozen survey; detailed below. |
| CUDA-MPQS | `b1a9c4500da15a2bea56a59ecad7ef61cf054756` | Newly inspected parameter selection and autotuning paths. |
| CADO-NFS | `692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b` | Active factorization calibration path now traced; see the follow-up below. |
| msieve | Release 1.53, as frozen | Retained as an established implementation reference; this audit does not establish current development-head equivalence. |

[Yamaquasi comparison](https://github.com/remyoudompheng/yamaquasi/compare/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad...181e282cbfcace30bd0d99c2b4bae181dee7ba69)
and [FLINT comparison](https://github.com/flintlib/flint/compare/a4c9750d0d3d67bb01cf6d18c187591b313451c3...26b9bc9df973206db15616d0849a4c8be9d22059)
make the reviewed deltas reproducible.

## FLINT: a concrete missed update

The [11 October change](https://github.com/flintlib/flint/commit/b0f4d54bc915f5d2b298958143ea999ef95142d4)
retunes the base/interval table, changes sieve organization, adds double-partial
handling, replaces a problematic A search with bounded randomized selection
and corrective factors, and revises collection/recovery. These span several
owners: calibration here, population selection in C10a-02, kernels in A1a,
and residual/graph changes in C9a. They cannot be attributed to tuning alone.

The active table selects by input bit length. For the actual frozen fixtures,
its nominal base count and total interval normalize as follows:

| v2 fixture band / bits | FLINT nominal count | Total positions | Half-width |
| --- | ---: | ---: | ---: |
| 55 / 181–182 | 2500 | 81920 | 40960 |
| 60 / 200 | 4500 | 131072 | 65536 |
| 65 / 215–216 | 8000 | 327680 | 163840 |

These nominal counts have different special-prime bookkeeping from v2's
exact eligible count. The table reinforces scaling count and interval together;
it does not identify a PyPy optimum. The frozen challengers remain plausible
hypotheses, but omit intermediate and smaller interval choices.

## CUDA-MPQS: calibration and short-probe limits

The [October paper](https://arxiv.org/html/2610.07126v1), §VI-B, explicitly
reports that truncated probes can mispredict large-prime matching and that
autotuning may fail to repay its own cost. This supports measuring complete
factorizations with setup charged. Its reported record/performance results
were not independently reproduced here; GPU results do not establish v2
speed or parameter choices.

At the pinned [source](https://github.com/drjanosch42/cuda-mpqs/tree/b1a9c4500da15a2bea56a59ecad7ef61cf054756/src),
`determineParams` chooses an endpoint by size, defaults M to 262144 and
uses the classic sqrt(2N)/M target. `SieveParameterOptimizer::optimize`
optionally probes M, then searches endpoint/large-prime-bound pairs with cached,
budgeted evaluations. `estimateRuntime` uses truncated yield, ETA and a
birthday-model large-prime correction. This is a broader search than C10a-04.
Its native GPU cost model and projected winner are not transferred; the frozen
study selects using complete/censored outcomes instead.

## CADO-NFS: completed calibration trace

The initial screen was incomplete. The subsequent
[calibration trace](cado_calibration.md) now follows the pinned active
factorization profile lookup, task/CLI wiring, A enumeration, CRT/Gray norms,
base construction, residual recovery and filtering retry. Native execution
and performance reproduction remain outside this source review.

The concrete gap is smaller intervals coupled to full family reuse. CADO-NFS's
c60 profile has half-width 8192, eight factors and 128 variants; v2's frozen
60-digit challengers have half-width 65536/131072 and quota 64.
The current v2 quota is a Gray-order prefix, not a demonstrated diverse
subset. The frozen research's wording is corrected here while retaining
its hashed bytes. C10a-07 in the public roadmap now owns quota/traversal
research and experiments; C10a-04's identities remain unchanged.

## What the frozen experiment can establish

The existing protocol is suitable for a bounded first comparison of its three
specified bundles: matched inputs/seeds/grants, validated PyPy warmup, charged
setup/resume, unresolved outcomes, complete-factor selection and committed fresh
confirmation. Short-screen raw-row gains remain diagnostic only.

Its coverage has explicit limits:

- Three bundles cannot establish a global or state-of-the-art optimum.
  Fixed residual bounds also change their ratios to the enlarged base.
- Three independent inputs per band remain a small conditional sample.
  Nine seed observations are not nine independent moduli.
- Native implementations have different arithmetic representations, kernels,
  filtering/solvers and hardware cost structures. Source-derived trends are
  hypotheses; only matched v2 measurements can establish local improvement.
- The latest msieve development head was not audited. Native CADO-NFS
  performance and optional algorithm paths were not evaluated.

Stage 2 may run the existing freeze to answer that narrow question. Any claim
of comprehensive current-source calibration, or addition of new interval,
population or residual policies, needs a separately recorded research/design
revision and new freeze before measurements. This audit supplies no authority
to retune from accepted/fresh results. The Stage 1 completion claim concerns
the frozen first tranche and its executable handoff.
