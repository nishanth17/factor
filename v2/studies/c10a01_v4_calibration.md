# C10a-01 v4 short calibration and pilot handoff

**2 November 2026.** This is exploratory fixed-work calibration of the opt-in
coupled implementation. It contains no completed factorization. The later
complete-factor feasibility pilots are funded and source-pinned, but unrun.
The owned control is committed mainline `0076e256871722d52d849ab344152e0a7535c2ce`;
the challenger implementation is `4b9a6163898a0b485436665d31c2a274e51ef643`.
The primary execution source was `1f95fc4` and the diagnostic source was
`6eaaaec`. Current release hashes and exact literal configurations are in
[`nominee_v4.json`](../benchmarks/qs/c10a01/nominee_v4.json).

## Fixed-work observations

The primary comparison used two exposed, independently generated semiprimes
at each size. Every arm sieved 512 polynomials and 67,109,376 positions,
reached its fixed-work stop, and returned no factor. Timings below are wall
seconds, excluding separately accounted process warmup. Rows are retained
relations. The 70 balanced 30.340-second result was **only this prefix**, not
a complete 70-digit factorization.

| Digits | Balance | Control seconds / rows / k | Challenger seconds / rows / k | Wall ratio |
| --- | --- | --- | --- | ---: |
| 60 | balanced | 19.896 / 113 / 13 | 21.929 / 295 / 97 | 1.102 |
| 60 | 35:65 | 20.889 / 137 / 3 | 20.541 / 192 / 3 | 0.983 |
| 62 | balanced | 20.271 / 72 / 1 | 21.417 / 82 / 13 | 1.057 |
| 62 | 35:65 | 21.259 / 102 / 31 | 21.663 / 177 / 7 | 1.019 |
| 65 | balanced | 26.518 / 90 / 3 | 26.274 / 105 / 3 | 0.991 |
| 65 | 35:65 | 28.176 / 64 / 5 | 25.744 / 71 / 1 | 0.914 |
| 70 | balanced | 29.586 / 23 / 1 | 30.340 / 36 / 1 | 1.025 |
| 70 | 35:65 | 29.497 / 19 / 1 | 29.876 / 28 / 29 | 1.013 |

The challenger retained more rows on all eight inputs; elapsed differences
have mixed signs. The scored 70 balanced k=1 uses standard-form fallback, so
that cell cannot isolate the parity representation. No complete-cost or speedup
claim follows from these prefixes. The factor/modulus generators are conditioned
on size, residue and balance, not uniform RSA sampling.

## Component diagnostics and frozen selection

Nine counterbalanced root-construction repeats at each anchor measured the
same first eight A families and 64 Gray steps per family. Total measured CPU
seconds for standard / cached / batch were 9.724 / 8.597 / 8.559 at 62 and
17.527 / 15.494 / 15.666 at 70. Batch improved on standard in the paired
99% bootstrap intervals, but the batch versus cached intervals cross zero.
The predeclared no-common-winner rule retains **batch**.

At fixed k, windowed-pair versus paired-flyer yielded 82 versus 33 rows at
62, and 36 versus 25 at 70, each after 512 polynomials. These are deliberately
different A populations, with no trial-by-trial equality claim. The frozen
rule retains **windowed_pair**. Norm-preview versus wide-density setup used
2.567 versus 0.509 CPU seconds at 62 and 4.359 versus 0.821 at 70, selecting
the same k at both anchors. Offline density winners matched the norm policy's
selected k on all eight primary inputs. The frozen rule therefore selects
**parity2_wide**; its agreement on these inputs does not prove universal
equivalence.

The four forced-multiplier prefixes reached their declared positions: k=5
standard and parity at 60 (100,664,064 positions each), even k=2 standard
fallback at 62 (83,886,720), and eligible k=31 parity at 70 (67,109,376).
They are policy diagnostics, not scored-selector observations. Both observer
on/off pairs preserved work, positions and the recorded semantic digest.
That digest covers selected k, family cursor and atom IDs, not every relation
exponent. These short observer pairs cannot establish a timing effect.

The selected exploratory nominee is `root_update_policy="batch"`,
`assignment_policy="windowed_pair"`, `multiplier_policy="parity2_wide"`,
`polynomial_policy="parity2"`, and `marking_policy="parity_two"`. Defaults are
unchanged. The root cache, A population and multiplier policy are separately
implemented opt-ins; no additivity or global optimality is claimed.

## Charges, memory and rate audit

The corrected primary release has 16 completed ledger entries, 1,380 reserved
wall and CPU seconds, 182,037,200,010 settled work units, and 393.878 seconds
of inclusive measured wall. Diagnostics have 16 completed ledger entries
representing 22 planned cells, 1,530 reserved seconds, 232,825,562,254 work
units, and 455.039 seconds inclusive wall. All 38 cells therefore reserved
**2,910 wall and CPU seconds** of the new 3,400-second ceiling and settled
**414,862,762,264 work units**. The first execution release's ledger bug left
two valid raw receipts but zero settled work; its 180-second failed-release
reservation is retained separately. The original calibration had already
reserved 3,124 seconds. No failed receipt was rewritten.

The fastest sustained charged-work interval in the new single-cell prefixes
was 0.5412 billion units per CPU second. The fastest measured root group
member reached 1.370 billion, still below the retained 8-billion planning
rate; the old
matrix fixture measured 6.244 billion. No new prefix reached preparation or
solving. Thus the already frozen pilot guards keep 4 × 8 billion units per
CPU second, plus separately bounded setup and analytic postprocessing charges.
These work caps are finite safeguards, not runtime forecasts.

The maximum observed RSS among the new single-cell receipts was 133,152,768
bytes. Group-process RSS was not recorded. These warm-prefix RSS values cannot
estimate a mature relation store or matrix peak. Conditional owned scenarios
including parent and root cache are about 0.769 / 0.924 / 1.346 / 1.728 GiB at
60 / 62 / 65 / 70 under the stated population assumptions in
[`protocol_v4.md`](../benchmarks/qs/c10a01/protocol_v4.md). The matched
operational ceilings are 1 GiB owned / 3 GiB RSS at 60, 2 / 5 at 62 and 65,
and 4 / 8 at 70. They are capacities, not claims of required memory. Current
readiness checks saw 24 GiB physical RAM and about 18.04 billion bytes
available; the 70 tier requires about 12.88 billion including the 4 GiB host
reserve. RSS enforcement is cooperative and can overshoot in-flight.

| Digits | Structural support | Fixed-work result | Complete factoring | Main uncertainty |
| --- | --- | --- | --- | --- |
| 60 | Exact, reachable | Both inputs, both arms paused | Unrun v4 pilot | Mature cost and confirmation |
| 62 | Exact, reachable | Both inputs, both arms paused | Unrun v4 pilot | Relation/matrix growth |
| 65 | Exact, reachable | Both inputs, both arms paused | Unrun v4 pilot | Completion within 3,600 seconds |
| 70 | Exact, reachable; scored fallback at balanced k=1 | Both inputs, both arms paused; forced k=31 exercised parity | Unrun v4 pilot | Completion and useful parity frequency |

The 50/55 complete-cost and 51/52/55/56 boundary supplement is separately
designed but unrun. Its design turn used Astra high, below the user's Astra
xhigh condition for experiment design; it needs xhigh review before execution.

## Executable, unrun later-stage handoff

[`pilots_v4.py`](../benchmarks/qs/c10a01/pilots_v4.py) checks the frozen
nominee, both calibration ledgers when executing, the current-mainline source
snapshot, all inherited production hashes, the PyPy 3.11 binary, literal
configurations and host headroom. Each later cell has a fresh process, at
least three wall and CPU seconds of validated warmup, nine pre- and post-call
sentinels when available, a finite complete-factor grant, bounded worker
output and durable ledger reservations. A native censor stays in the record;
no automatic retry, larger grant or adoption test is built in.
Each later stage has a 16 MiB ledger/output allowance inside the shared
512 MiB study ceiling. The final local archive leaves about 76 MiB of that
ceiling available; later execution must recheck remaining space.

| Stage | Calls | Per-call wall and CPU caps | Stage wall and CPU ceiling | Cumulative work ceiling | Matched tiers |
| --- | ---: | --- | ---: | ---: | --- |
| A | 4, 60 and 62 | 1,800 / 2,400 seconds | 9,200 seconds | 3,983,055,278,297,769,690,620 | 1 / 2 |
| B | 2, 65 | 3,600 seconds | 7,600 seconds | 9,972,377,840,558,931,315,448 | 2 |
| C | 2, 70 | 7,200 seconds | 14,800 seconds | 31,435,080,283,294,435,139,152 | 3 |

Readiness-only commands, verified without factoring:

```sh
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.pilots_v4 check --stage A
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.pilots_v4 check --stage B
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.pilots_v4 check --stage C
```

Later execution commands, **not run in this session**:

```sh
/opt/homebrew/bin/pypy3 -u -m v2.benchmarks.qs.c10a01.pilots_v4 run --stage A
/opt/homebrew/bin/pypy3 -u -m v2.benchmarks.qs.c10a01.pilots_v4 run --stage B
/opt/homebrew/bin/pypy3 -u -m v2.benchmarks.qs.c10a01.pilots_v4 run --stage C
```

The full original-plus-new-plus-lower design ceiling is 42,300 wall and CPU
seconds. Including the failed release separately makes its historical maximum
42,480 seconds. Original actual reservations, the failed release, actual new
calibration reservations and the prospective lower/pilot ceilings total
**41,514 seconds**. The new calibration left 490 seconds unused. The 70 stage
is admitted on its
own arithmetic, configuration and host feasibility; it does not require 65
to complete. Retained 60-digit complete jobs took about 426–796 seconds under
an earlier combined policy. An old 65-digit 1,800-second failure preceded the
latest collector changes. Neither gives a defensible 70-digit scaling law.

No pilot or fresh confirmation has run. One exposed balanced modulus per size
would give feasibility evidence only. Adoption requires a separately frozen,
fresh, independent-input comparison with matched complete cost and uncertainty.
