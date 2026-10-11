# C10a-05 larger-band default-scope follow-up

## Decision and limits

**Adopt the existing C9+C10+collection-pacing bundle for native 51-digit
inputs at exactly 1 GiB owned SIQS memory. Retain the previous default for
52 digits and above.** The accepted 50-digit C10a-05 decision remains
unchanged. At 51 digits the current default stopped without a factor on all
nine fresh matched calls; the selected bundle completed all nine and passed
independent divisor, certainty and reconstruction validation. This is a
completion gain, **not a measured speedup over successful baseline calls**.
The source-derived 51-digit candidate arms held base, polynomial, interval,
collector admission, graph, storage, matrix, work, wall/CPU and memory fields
fixed; only the C10 power plan and then fixed solve pacing changed. The
comparison to the actual 51-digit default necessarily compares whole bundles:
its smaller base, finite 16-family search and 32 MiB owned allowance explain
the fast unsuccessful stop. Both arms had the same external work/time/RSS
allowances; their internal memory policies differed and are reported below.

The request for every input of 50 or more digits cannot be supported by a
finite bounded study. One inspected 55-digit product finished under the
unchanged bundle in 421.872 s, but a 60-digit product did **not** finish in
the existing 1,200 s omitted grant. Neither single-number pilot justifies a
default at 55, 60, or any unlimited higher band. The 1,200 s figure is a
finite automatic grant, not a universal hardness limit; explicit larger
budgets remain available. The unchanged 50-digit base/interval may need
separate calibration at larger sizes, outside this fixed-policy follow-up.

## Frozen comparison

The executable [protocol](protocol.md), committed runner at `4a197c0`,
[training corpus](../../inputs/corpora/c10a05_plus_training.json), and
[control](../../inputs/controls/c10a05_plus_training.json) were frozen before
accepted timing (`4c39e51`). [Selection](../../inputs/controls/c10a05_plus_selected.json)
at `76b2d4c` hashes all 36 training receipts and chooses collection pacing.
Only then were three [fresh certified products](../../inputs/corpora/c10a05_plus_fresh.json)
generated and committed (`6faa13d`). Each phase crossed three independently
generated 51-digit balanced semiprimes with seeds 7, 29 and 47; seeds are
repeated observations, not additional independent numbers. The generator's
prime certificates condition a large factor of p−1, so these are not an
unconditioned RSA-semiprime population. The first training number was already
inspected in feasibility; no inspected number entered fresh confirmation.

The 51/55/60 feasibility screen used one revealed number and seed per band.
Its first 51-digit product took 299.972 s with C9/window, 240.311 s with
C9+C10/window and 173.913 s with C9+C10/collection. At 55 digits the current
default stopped in 0.150 s with no rows; the candidate was censored at 240 s
with 483 rows, then completed a separate 1,200 s-grant probe in 421.872 s.
At 60 digits the current default likewise stopped in 0.160 s; the candidate
had 30 rows at 240 s and ended the longer 1,200.042 s call with 352 rows,
12 filters, one solve, a trivial dependency and no factor. Its 1,090.6 s
collection stage dominated; this was not a downstream-pacing bottleneck.
These pilots are excluded from acceptance uncertainty.

| Accepted 51-digit arm | Training complete | Training wall total | Fresh complete | Fresh wall total |
| --- | ---: | ---: | ---: | ---: |
| Actual current default | 0/9 | 1.100 s unsuccessful | 0/9 | 1.141 s unsuccessful |
| C9/window | 9/9 | 3,153.506 s | — | — |
| C9+C10/window | 9/9 | 2,495.983 s | — | — |
| C9+C10/collection | 9/9 | 1,792.294 s | 9/9 | 1,190.931 s |

On the nine training pairs, collection pacing reduced **complete** cost by
28.19% against the otherwise identical C10/window arm (95% input-cluster
bootstrap interval 27.43–29.07%). The C10 plan reduced complete cost by
20.85% against C9/window (interval 20.25–21.54%). These component contrasts
are training attribution; fresh confirmation compared the selected whole
bundle with the actual default. There was no predeclared instability trigger
in either phase, so neither phase extended to eighteen groups. Against the
fresh unsuccessful reference, the selected arm used 81.62% less **capped**
cost under the 720 s per-call allowance (95% interval 80.04–84.12%); that
ratio is a resource-and-completion statistic, not factoring speedup. The
reference actually returned in 1.141 s total with no factors.

The fresh selected arm spent 1,056.16 s collecting, 27.10 s preparing,
9.44 s filtering, 2.78 s eliminating and 0.19 s extracting. It used
646,915,427,800 charged work units, 478 filters, 52 solves and at most
360,230,952 owned bytes across its nine calls. The fresh current arm used
1,632,819,790 work units, 18 filters, nine solves and at most 13,102,904
owned bytes, but completed none. On training's matched C10/window and
collection arms, filter calls were 9,279 versus 479 and solves were 9,270
versus nine; preparation/filter/elimination totals were 754.30 versus
40.59 s. Collection times were 1,584.33 versus 1,593.81 s. These stage
figures explain the complete-cost result but were not the selection target.
All measured proper factors and reconstructed cofactors passed validation;
there were zero correctness failures.

Inclusive charged wall/CPU for the follow-up is 12,199.720/12,188.535 s,
including a conservative 3,500 s prepaid feasibility/generation/QA reserve.
The 36 training and 18 fresh accepted arm calls, cold captures and validated
PyPy warmups are recorded in the ignored `results/c10a05_plus/` tree. There
were no failed or uncharged reservations. Training yielded after six matched
groups at a whole-group lease boundary and resumed with the same 24 saved
receipts; no allowance was reset. The study stayed under its 20,000 s training,
12,000 s confirmation and 36,000 s aggregate caps. Generated receipts and
profiles remain ignored; committed controls, corpora and scripts are the
reproducibility inputs.

From the `6faa13d` checkout with PyPy implementing Python 3.11, the accepted
sequence is:

```text
PYTHONDONTWRITEBYTECODE=1 v2/.venv/bin/python -m v2.benchmarks.qs.c10a05_plus.study train --lease-seconds 7200
PYTHONDONTWRITEBYTECODE=1 v2/.venv/bin/python -m v2.benchmarks.qs.c10a05_plus.study confirm --lease-seconds 7200
PYTHONDONTWRITEBYTECODE=1 v2/.venv/bin/python -m v2.benchmarks.qs.c10a05_plus.study summarize
```

Reproduction must use that pinned checkout; later default integration
intentionally changes the study's strict source hash. The committed controls
and fresh corpus are immutable. The source/literature reconciliation, state
invariants, independent dependencies and charged-resume oracles remain in the
closed [C10a-05 research](../c10a05/research.md) and
[correctness](../c10a05/correctness.md) records. This follow-up changes only
the size-scoped default selector, not the collector, checkpoint schema or
explicit SIQS configurations. C9a feedback, other C10a items and E1 portfolio
acceptance remain open; no earlier percentage is added into a portfolio claim.
