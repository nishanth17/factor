# C10a-04 Stage 1 diagnosis and validation

This is a research/implementation freeze, **not** a completed training or fresh
confirmation campaign. Production defaults are unchanged; C10a-04 and E1
remain open. Raw captures are ignored under `v2/benchmarks/results/c10a04/`.
Required source/configuration controls and training inputs are versioned.

## Limited feasibility screen

One revealed balanced input per band and seed 7, all three arms at 45 seconds,
followed by separate 120-second grants at 60/65 digits. These are diagnostic
prefixes, not warmed matched acceptance observations or time-to-factor results.
The earlier prefixes were preserved; longer calls are new charged attempts,
not replacements or free continuation. Every output reconstructs exactly.

| Digits / grant | Control raw rows | Compact raw rows | Wide raw rows | Post-filter rows / complete calls |
| --- | ---: | ---: | ---: | --- |
| 55 / 45 s | 43 | 105 | 82 | 0 in every arm / 0 of 3 |
| 60 / 45 s | 5 | 57 | 58 | 0 in every arm / 0 of 3 |
| 65 / 45 s | 0 | 27 | 19 | 0 in every arm / 0 of 3 |
| 60 / 120 s | 14 | 167 | 157 | 0 in every arm / 0 of 3 |
| 65 / 120 s | 0 | 81 | 50 | 0 in every arm / 0 of 3 |

The larger bases preserve promising smooth/cycle opportunity, but no useful
matrix or dependency is established by this screen. At 120 seconds the 60-digit
compact/wide arms had 10/15 graph cycles; at 65 digits they had two each.
All 15 calls ended on wall allowance, with reachable assignment search left;
none reported search exhaustion or owned-storage refusal. New k/count/actual
endpoint were 3/2500/48491 at 55, 11/5000/104543 at 60 and 19/8500/188159 at 65.
Their D mod 8 classes were 5/1/1. Control counts were 1627/1638/1651, endpoint
29989/29989/29927. Singular multiplier primes 3/11/19 remain in the base but
are excluded from A selection. No even multiplier or normalized arm was run.

Screen inclusive saved-call wall/CPU totals were 1125.027/1063.013 s and
146214678025 charged work units. Its pre-freeze source manifest and exact new
helper sources are retained in `screen-v1/ledger.json` and `screen-v1-source/`.
Later accounting repairs did not change runtime, helper arithmetic, policy
bundles, inputs or these historical receipts; no screen contributes accepted
performance evidence. The final runner retains full time reservations even
on success, with actual time separately reported, and preserves lost measured
calls as censored outcomes. This is a new identity rather than a rewrite of
the frozen C10/C10a-05 runners.

A separate warmed 55-digit compact diagnostic, bounded at 600 seconds,
completed in 390.669 s inclusive wall and 389.201 s CPU. Both returned factors
reconstruct the certified input, with runtime `probable_prime` labels preserved.
It collected 2529 rows, retained 2384 after filtering, found 44 dependencies
and used 89,346,774,364 work units with 430,700,344 owned bytes. The 605-second
reservation is retained for wall and CPU; actual time is separately recorded.
Its `pilot/compact-55-600.json` receipt, warmup and `pilot-source/` snapshots
remain local. It establishes feasibility on one revealed number, not a matched
speedup or fresh acceptance; no ratio against the historical 421.872 s is used.

## Grant rationale

The accepted previous 55-digit diagnostic completed in 421.872 s under the
reused bundle, supporting a matched 1200-second training/confirmation grant
with headroom. At 60 digits the new richer-base prefixes retain about 160 raw
rows per 120 s against a 5000-entry base. A naive raw-row extrapolation is
already thousands of seconds, before considering graph acceleration, filtering
loss and matrix cost; 7200 seconds is an informative finite grant, not a
completion prediction. At 65 digits 50–81 raw rows per 120 s against 8500
entries motivates a substantially longer 21600-second grant. All raw rows can
still be singleton losses; fresh complete outcomes decide. The separate doubled
grants diagnose resource sensitivity without changing the acceptance envelope.
Their worst-case reservation is explicit in the protocol, not an instruction
to spend the entire allowance. No failure at these grants rejects SIQS at a
universal size or time boundary.

## Correctness and runner checks

Independent exhaustive residues supply 72 count/endpoint cases across four
odd input classes and six multipliers, including even k and singular primes.
The test oracle uses trial-division primes and brute-force squares, independent
of the production sieve/Euler/Tonelli paths. Tests also cover gcd splits,
gcd=N refusal, finite count/budget refusal, retained accepted collector fields,
proper factors, unresolved reconstruction and certainty labels. New joint
configuration checkpoint replay preserves exact base identity, prior wall and
charged work. Production arithmetic and checkpoint formats are unchanged.

Receipt tests cover identical republication, conflicts/tampering, a crash after
receipt publication, lost measured cells, fatal correctness failures, repeated
yield leases and interruptions between matched arms and within the extension.
Immutable nine-observation stability markers use only their own prefix on
resume, even when later extension cells already exist. A lost measurement
remains a full-cost censored cell; it cannot be retried selectively into a
successful acceptance outcome. An independent completion-bootstrap fixture
shows that a gain confined to one of three inputs retains a zero lower bound.

Initial focused-test fixture mistakes (the FactorBaseBuild field name and a
repeated test number) were corrected before freeze; both failed and successful
transcripts remain local. The twelve focused tests pass in the full suite. Full-suite, lint,
pilot, smoke and committed-export receipts are recorded in the final handoff.
A separate 27-cell one-block preflight covered all nine training products and
three arms: exact requested counts, reachable 100000-family streams, valid
reconstruction and the expected paused outcome. It charged 540 s wall/CPU
reservations and 328,709,181 work units; this is setup validation, not timing
evidence. The initial preflight asserted an incorrect name for the documented
paused outcome; that transcript and the original idempotently reused receipt
are retained without repeating the completed call.


The frozen commit `5dab76ec66c02049095ddf2f75e123f19d0d6627` passed all 757
tests and lint from a committed-files-only export; all 185 benchmark modules
imported and all 249 source identities/nine certified inputs verified. Frozen
CLI smoke/resume preserved cold receipts across two leases, charged fresh
warmup, yielded twice and started no measured samples. See the handoff for
exact identity digests, local receipt locations and remaining allowances.
The default experiment directory, selection and fresh corpus remain absent.
