# Roadmap history

Source commit: `4a946b8d6294d51ac5d9c36f61880cc63d5b59fd`; source path: `v2/ROADMAP.md`; SHA-256: `0f1efd10c42d7a6f4c78f577888385c672f843ee78f5f05ba8e171a76169fba5`.

[STATUS](../STATUS.md) is the current entry point; [ROADMAP](../ROADMAP.md)
holds live work and [POLICY](../POLICY.md) holds stable rules.

The [pre-cleanup roadmap at `4a946b8`](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md) preserves detailed completed-item
narratives, every quantitative result, original decisions, limitations and
unexecuted gates. Git history is the historical record;
there is no duplicate full roadmap in the current tree. The completed register
and historical heading aliases link directly to that pinned revision. New dated
documents and current code govern live behavior. No historical freeze is rewritten.

Retrieve the exact original bytes locally with:

```sh
git show 4a946b8d6294d51ac5d9c36f61880cc63d5b59fd:v2/ROADMAP.md
```

## Original link context

Historical links open `v2/ROADMAP.md` at the pinned source commit. Its relative
URLs and heading fragments retain their original context. The local preservation
ledger checks the Git object's hash, historical anchors and destinations, and
links from the live documents. External research URLs are preserved references.

## Local audit reference register

These resolve every “(local audit material)” occurrence without publishing
local research. Source line numbers identify the pinned historical source. The
pre-removal links were checked in `328c823^:v2/audit/TODOS.md` where available;
current local files and anchors were checked at the paths below. Local paths
are named as paths, not presented as portable committed links. Repeated labels
share a row; all occurrences are accounted for in the local ledger.

| Historical source lines / labels | Existing destination |
| --- | --- |
| 13, 648, 1522 — Phases 3+ research pass / priorities / new research report | `/Users/nishanthmohan/code/factor/v2/audit/phases_three_plus_research.md` |
| 50, 4229 — C-to-Python transfer review / C transfer review | `/Users/nishanthmohan/code/factor/v2/audit/sieve_port_review.md` |
| 649, 2412 — Matrix research / matrix research report | `/Users/nishanthmohan/code/factor/v2/audit/gf2_matrix_research.md` |
| 777, 3898 — frozen baseline / M17 baseline | [versioned M17 baseline](../benchmarks/inputs/controls/phase_two_m17_frozen_baseline.json) |
| 797, 3894 — Phase 2 optimization report / research report | `/Users/nishanthmohan/code/factor/v2/audit/phase_two_optimization_research.md` |
| 941, 3277 — quadratic sieve research report / research report | `/Users/nishanthmohan/code/factor/v2/audit/quadratic_sieve_research.md` |
| 944 — source manifest | `/Users/nishanthmohan/code/factor/v2/audit/inputs/quadratic_sieve_research_sources.json` |
| 948 — literature/blog follow-up | `/Users/nishanthmohan/code/factor/v2/audit/quadratic_sieve_research.md#optimization-follow-up-m21` |
| 951 — supplement | `/Users/nishanthmohan/code/factor/v2/audit/inputs/m21_quadratic_sieve_optimization_sources.json` |
| 1541 — sub-100-digit reconciliation | `/Users/nishanthmohan/code/factor/v2/audit/phases_three_plus_research.md#below-100-digit-optimization-target` |
| 1635 — budget follow-up | `/Users/nishanthmohan/code/factor/v2/audit/phases_three_plus_research.md#budget-follow-up-large-job-capacity` |
| 2415 — Source provenance | `/Users/nishanthmohan/code/factor/v2/audit/inputs/m25_gf2_matrix_sources.json` |
| 4337 — validate.py | `/Users/nishanthmohan/code/factor/v2/audit/validate.py` |
| 4338 — extra_checks.py | `/Users/nishanthmohan/code/factor/v2/audit/extra_checks.py` |
| 4339 — end_to_end_checks.py | `/Users/nishanthmohan/code/factor/v2/audit/end_to_end_checks.py` |
| 4340 — prac_reference.py; test_prac_reference.py | `/Users/nishanthmohan/code/factor/v2/audit/prac_reference.py`; `/Users/nishanthmohan/code/factor/v2/audit/test_prac_reference.py` |
| 4341 — reduction_bench.py | `/Users/nishanthmohan/code/factor/v2/audit/reduction_bench.py` |

Unresolvable “local audit material” labels: none. Historical capture filenames
beside the scripts remain historical labels; generated captures stay local.
Do not treat a diagnostic script as an acceptance test or restart old accounts.

## Cleanup decisions and evidence

| Stale/conflicting entry | Resolution and source |
| --- | --- |
| Heading dated 3–5 October despite November entries | STATUS and ROADMAP use 1 November and the reviewed mainline commit. Original dates remain in the pinned Git revision. |
| Immediate allocation / lanes still start B1 | Mark B1's bounded tranche Done on 9 October; preserve larger feasibility and downstream owners. [Recorded B1 decision](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations). |
| Earlier C10a-07 “open”, “finish face comparison first” and unstarted screens | The completed 42-call panel superseded both unstarted 60-only screens. Closure supersedes the retain-default recommendation, with no upgraded experimental verdict. [Panel](../studies/c10a07_panel_results.md), [1 November adoption](../studies/c10a07_default_promotion.md), [changelog](../../CHANGELOG.md#c10a-07-combined-defaults-1-november-2026). |
| A1a six-call sparse-queue screen still unchecked | Panel replaced it; production queue adoption closes that delivered work under C10a-07. Keep quota-aware setup, lifted-root reuse and distinct conditional comparisons live. Same panel/adoption sources. |
| P3.4 says scored multipliers remain optional | Historical P3.4 decision retained. Current 52–60 defaults use scored odd multiplier (`multiplier=0`), confirmed by [selector](../qs/defaults.py), with the original [adoption limits](../studies/c10a04_default_promotion.md). |
| Broad “SIQS stays opt-in” / “handoff open” language | Scope current automatic routing to native 40 and 50–60 via [portfolio](../portfolio.py) and [selector](../qs/defaults.py); joint calibration remains G1/E1. Changelog wording is explicitly historical. |
| C6/B4/B3 ladder/default statements from separate dates | Keep isolated study verdicts; current supported native PRAC is the B3 route, with exact fallback from [eligibility code](../ecm/chain_options.py). |
| Phase 8 labelled wholly unimplemented despite A10/C3 deliveries | Show completed capability tranches separately from still-open integrated P8 work and deferred E1. [A10 acceptance](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p821--transfer-verified-v1-advantages-starting-with-primality-a10), [C3 decision](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#calibrate-automatic-selection-and-handoff-c3-then-g1). |
| Removed local-audit links | Resolve every marked occurrence through the register above. The original Git object remains unchanged. |
| Missing spaces around numbers | Normalize live summaries and affected guide/changelog prose. Keep math, identifiers, Git history and study reports unchanged. |

Unresolved conflicts: none among the checked entries. No source report,
experimental verdict, code, default, test, corpus or generated evidence was
changed. This follow-up changes documentation storage and links only. Publishing to
mainline is explicitly authorized; no new experiment is scheduled.
