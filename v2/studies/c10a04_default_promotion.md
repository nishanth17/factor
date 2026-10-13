# C10a-04 default promotion — 12 October 2026

The user explicitly requested automatic 50–60-digit routing, local mainline
integration and closure after reviewing the evidence. This is a policy
promotion at user direction, not completion of the original multi-band
experiment or a claim that all its statistical gates passed.

| Native input size | Fresh automatic configuration | Omitted wall/CPU grant |
| --- | --- | --- |
| 40 digits | Existing accepted bundle, unchanged | 1200 s |
| 50–51 digits | Existing accepted collection-pacing bundle, unchanged | 1200 s |
| 52–55 digits | Compact: 2500 eligible entries, seven A factors, pool 64 | 1200 s |
| 56–60 digits | Compact: 5000 eligible entries, eight A factors, pool 64 | 1800 s |

Both compact bands select an odd multiplier, retain half-width 65536 and
at most 64 Gray variants per family, and use existing growing-DLP collection,
collection pacing and the 2 MiB power plan at the 1 GiB owned allowance.
The actual base endpoint depends on the input and multiplier. Count resolution
has a finite prime-endpoint ceiling of 250000, checks temporary sieve memory
before allocation and shares the job's time/work budget. Replay repeats this
deterministic setup under remaining resources, without refunds.

Explicit budgets, memory grants, SIQS configurations, parameter flags and
`siqs=None` retain precedence. Other memory grants retain window pacing and
omit the power plan. Other backends, 41–49 digits and sizes above 60 retain
previous routing. Standalone `factorize()` is unchanged.

`CountedSIQSConfig` is a `SIQSConfig` subclass with an explicit positive
`base_count`; `base_bound` is its finite endpoint ceiling. Counted checkpoints
retain this policy and the existing multiplier/base identities in standalone,
portfolio and CLI restore. Historical `SIQSConfig` mappings acquire no field
or automatic policy. A changed count is not a resource-only resume extension.

## Evidence and limitations

The existing 50-digit collection-pacing improvement of 25.99% is already
accepted and enabled. It is preserved, not compounded with this study.
At 55 digits all 54 training calls completed: compact used 3094.613 seconds,
control 3450.669 and wide 3568.181 across 18 calls per arm. Compact's 10.32%
lower total cost failed the original persistent-instability rule; no fresh
55-digit confirmation ran.

The [fresh 60-digit follow-up](c10a04_short_results.md) stopped at its fixed
six-pair preliminary trend check: compact completed 6/6 versus 5/6, with
36.16% lower capped time and 31.72–40.08% saving per input. It covers three
independent conditioned balanced products, not six independent numbers or
nine accepted timing samples. Compact used more work and owned memory despite
less elapsed/CPU time. Routing 52–54 and 56–59 is an explicit extrapolation;
no measured gain against the prior automatic reference path is claimed there.

The original freeze `5dab76ec66c02049095ddf2f75e123f19d0d6627`, annotations,
receipts, interrupted reservation and short-study seal remain intact.
Integration changes current runtime hashes intentionally. Reproduce frozen
original and short-study loaders from committed state `225634b`; never update
their source manifests or restart their ledgers. A new chained A7 integration
manifest permits its unchanged explicit 30/40-digit arms on this runtime;
all prior manifests and the prepared arm plan remain byte-identical.

Tests cover automatic boundaries, unchanged accepted bundles, parity with the
measured 60-digit configurations, independent count arithmetic, proper GCDs,
finite work/memory refusal, CLI/portfolio serialization and charged replay.
Full PyPy tests, lint and committed-only benchmark imports are checked before
merge. Raw study and QA evidence remains ignored and is preserved before
worktree cleanup.

Close the delivered C10a-04 routing tranche at user direction. The original
campaign remains incomplete, E1 deferred and C10a-07 open. The separately
requested 65-digit feasibility spike uses 8500 entries and nine A factors,
with finite matched grants; extending defaults through 65 depends on its
result. A bounded 70-digit feasibility pilot remains a roadmap follow-up.
