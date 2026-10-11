# C9 mainline integration — 11 October 2026

The user authorized integration of C9, closure of its accepted tranche and
phasing the remaining literature transfers into C9a. C9 acceptance is pinned
at `868c266`; its subsequent source verification is `19bb25d`. Integration
consumes C3 mainline `d273cb8` and preserves both histories.

The accepted aligned growing-DLP configuration remains explicit and scoped
to the confirmed balanced 50-digit regime. C3's separate ECM defaults and
acceptance limitations remain unchanged. No new performance study, default
promotion, fresh corpus or combined E1 claim is introduced.

## Compatibility and evidence

C3's implicit schema-12 configuration decoder previously selected legacy DLP
when it encountered a growing collector. Decode `storage_policy == "grow"`
first so absent quotas remain valid and the finite resolved child memory
allowance survives restore. Independent regression cases cover an unadmitted
fallback and an active SIQS checkpoint, repeated implicit restore, complete
reconstruction, cumulative work/wall/CPU, one-way handoff and conflicting
policy rejection. Earlier growing-store arithmetic/kernel/ownership oracles
remain in the combined suite. No serialized layout or identity changes.

The separate `a7_r5_c3_c9_sources.json` manifest binds integrated live sources
and chains to the existing C3-defaults manifest. Explicit A7 arms and all
historical source manifests remain unchanged; this is preparation, not E1
acceptance. The accepted C9 selected control/corpus and their source hashes
are preserved. Its timing loader intentionally rejects the newer integrated
source; use the documented historical freeze for reproduction.

## Follow-up ownership

Roadmap C9 is crossed off for its accepted C9.2 tranche. C9a owns remaining
residual admission/ordered endpoints, yield-dependent target feedback, cycle
and duplicate attribution, and profile-gated compact storage/splitting/batch
or disk transfers. C10/C10a retains scoring/recovery and filter/solve scheduling;
B7/A8 retains matrix/provenance redesign, C1a anchor-only SLP, and F6 higher
arity. C9a must select and freeze its own bounded tranche and confirm complete
outcomes; hypotheses are not accepted gains. C10's frozen worktree is preserved.

## Verification and cleanup

Integration verification runs under `/private/tmp/factor-performance.lock`
with explicit ownership coordination. Results are recorded below after the
combined tests, lint and committed-files-only benchmark/loaders checks.
Useful raw C9 evidence and pinned-source captures stay ignored and are moved
to mainline's `v2/benchmarks/results/c9/` before removing the worktree.
Unrelated mainline edits are preserved outside the integration commits.

At integration revision `6ba8cdd634364d095bd17dfa876999894179db2a`,
`make -C v2 test` passes **690 tests**, and `make -C v2 lint` passes all
Ruff/format/pycodestyle checks. A `git archive` of that committed revision,
with only external PyPy dependencies, independently passes **690 tests**,
**140 live non-package benchmark imports**, snapshot/A7/B4/C6/C1 loaders and
the three independently certified C9 confirmation inputs. The historical C9
timing loader correctly refuses integrated source; accepted control and
corpus SHA-256 values remain unchanged. The temporary archive is removed.

The combined checks use 206.258 wall seconds and 216.667 child CPU seconds;
these are integration QA costs, separate from the historical C9 study's
closed allowance. No timing comparison or new study is conducted. The final
verification receipt update changes documentation only; runtime, tests,
benchmark runners and source manifests remain exactly as checked.
