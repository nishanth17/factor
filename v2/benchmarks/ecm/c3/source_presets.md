# Source-derived ECM transfer presets

These opt-in numerical schedules implement the requested upstream-default
comparison without changing SIQS parameters or protected admission. No upstream
code is copied. The current portfolio's curve arithmetic, stage two, cost model
and input prior differ; source counts are neither calibrated probabilities for
this implementation nor accepted performance claims.

| Preset | Independent tiers `(B1, B2, curves)` | Source meaning |
| --- | --- | --- |
| `yamaquasi_auto160` | `(200, 7700, 8)` | Automatic 65–160-bit input pretest |
| `yamaquasi_ecm64` | `(200, 7700, 10)`, `(2000, 81000, 30)`, `(10000, 554000, 100)` | ECM-only prefix through the 64-bit factor target, 140 curves total |
| `alpertron20` | `(2000, 200000, 25)`, `(11000, 1100000, 90)` | First two bound-ladder tiers, 115 curves total |
| `alpertron25` | Previous tiers, then `(50000, 5000000, 300)` | First three bound-ladder tiers, 415 curves total |
| `gmp_ecm20` | `(11000, 1900000, 74)` | Published expected-curve row for a 20-decimal-digit factor |
| `gmp_ecm25` | `(50000, 13000000, 214)` | Published default-polynomial expected-curve row for a 25-decimal-digit factor |

The suffixes describe factor targets except `auto160`, which names the source
input-bit band. They are labels, not automatic dispatch rules. GMP-ECM's 25-digit
count assumes its polynomial stage two; the v2 transfer does not acquire that
algorithm. Expected counts also are not near-certain success budgets. Alpertron
has a separate SIQS handoff table that can stop well before completing these
standalone ladder prefixes; the presets do not mislabel that table as a
415-curve automatic pretest.

## API and bounds

```python
from v2.execution.ecm_presets import ECM_PRESETS, with_ecm_preset

# config already specifies SIQS and an ECMAllocation with positive work,
# wall and CPU fallback floors.
config = with_ecm_preset(config, "alpertron20")
```

The helper replaces only `ecm_tiers`. It retains the identical SIQS and
allocation objects, including cumulative pretest ceilings and fallback floors,
and the outer memory cap. Normal configuration validation rejects a larger
workspace if it cannot coexist with SIQS storage. Execution retains shared
finite work/time budgets, charged setup, one-way handoff and schema-12 restore.
A reservation is an admission floor, not a promise of completion. Applications
must provide a finite total grant that funds the existing reservation.

This additive module is outside the frozen compact64/wide64 comparison path.
The completed calibration, selected table and existing comparison gates remain
unchanged. Transfer measurements will be separately frozen on revealed inputs
before any preset is selected; fresh confirmation follows candidate selection.
The table alone does not close C3 or justify a default change.

## Provenance

The prior source review pins Yamaquasi commit
`3f95f43682ed15d8c1ed206a9a702dd655d7c8ad`, Alpertron commit
`93c5c8189cb2f149c8dae898eb11997c5f2e7980` and GMP-ECM commit
`8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e`, with local capture hashes in
`inputs/controls/c3_research_sources.json` and
`inputs/controls/c3_round2_research_sources.json` (relative to `benchmarks/`).
Numerical rows were checked against the reviewed source:

- [Yamaquasi ECM](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/ecm.rs)
- [Alpertron ECM and separate SIQS handoff](https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/ecm.c)
- [GMP-ECM parameter table](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/README)

[The depth review](round2_depth_review.md) distinguishes these source meanings
and records the remaining 40-digit training survivors that motivate stronger
follow-up schedules.
