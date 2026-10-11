# C9.2 amended 1,200-second confirmation

Frozen on 11 October 2026 at the user's explicit direction to raise the
50-digit comparison itself. This supersedes the unrun conditional diagnostic
in `c9_long_window_addendum.md`; that diagnostic consumed no measurement
allowance. It does not retune the selected collector or any arithmetic.

## Preserved evidence and source bridge

Stop the 360-second comparison with its raw records intact. Its 30/40-digit
bands each have 18 matched samples and no new growth promotion because the
prespecified stability gate remained unresolved. Its 50-digit band has only
two matched groups: aligned growth completed both, while calibrated SLP and
fixed DLP each completed one. Those observations are preliminary, censored
evidence; they are not nine-sample confirmation.

The original ledger charges **8,917.430828 active seconds**, including the
full 362-second unpublished reservation interrupted for C3's requested QA.
No cost is refunded. The prior 750-second QA reservation remains included.
The new stage has a separate **14,400-active-wall/CPU-second allowance**, with
300 seconds prepaid for its harness checks, source/corpus validation and final
assessment. Thus the amended aggregate allowance is **23,317.430828 seconds**,
not a reset that hides earlier expenditure. Research, ordinary editing and
waiting for coordinated machine ownership remain excluded. Stop at the cap;
there is no automatic further extension or new parameter search.

Production source, collector configurations and original training selection
are unchanged from the pre-confirmation runtime freeze. The new benchmark
entry point changes the 50-digit wall/CPU allowance, source/corpus identity,
generation seed and cold-worker entry point. Its source digest is included
in `inputs/controls/c9_long_selected.json`. This control retains the previous
selection's configuration and evidence hashes. Earlier source-frozen runners
remain reproducible from their original commits, rather than being silently
relabelled with the new allowance.

## Fresh inputs and matched design

Freeze the new source/control before generating **three independent balanced
50-digit inputs**, disjoint from all versioned corpora, using seed
`202610110619`. Certified primes are conditioned on a large known factor of
p−1; the population is not RSA-distributed. Expected factors and certificates
enter validation only. The earlier seven confirmation inputs are now revealed
regression evidence and are not reused here. No new 30/40-digit run is needed.

Use the same calibrated SLP, accepted fixed DLP and selected full-base aligned
growing DLP configurations. Each receives **1,200 wall/CPU seconds, 10^13 work,
1 GiB owned memory and 3 GiB process high-water RSS** per call. Successful calls
stop immediately. Keep explicit checkpoint, family, polynomial, scoring,
product, endpoint, splitter and dependency settings identical to the frozen
selection. A timeout is a censored outcome, not a digit-feasibility verdict.

Use three seeds (7, 29, 47) on each of the three fresh inputs, for a minimum of
nine matched groups, with the inherited deterministic randomized arm order.
Validate at least three seconds of both PyPy wall and CPU warmup for every arm
in every new measurement process. Keep one cold call per arm separate from
warmed comparisons. Validate every returned factor, unresolved reconstruction,
certainty label and resource envelope independently. Charge setup, failed
calls, splitting, proof checks, rebuilding, filtering and extraction. Existing
checkpoint/repeated-restore diagnostics remain separate and charged.

Use the inherited stability rules and hierarchical input/paired-sample
bootstrap. Sign conflicts or more than 10% repeated-cell drift/IQR require the
single predefined extension to 18 matched groups. That 10% is a noise trigger,
not an improvement floor. An extension must fit its full reserved allowance;
otherwise retain the baseline and label timing inconclusive. In particular,
27 additional 1,202-second reservations cannot fit this stage's allowance;
no exception is granted after observing the result. An aggregate stop before
nine matched groups also leaves confirmation incomplete.

Promotion requires a positive 95% benefit interval, completion regression no
larger than five percentage points, passed correctness/resource contracts and
acceptable costs. Report wall/CPU/work, complete outcomes, proper factors,
collection and splitting counters, graph occupancy/evictions, filtering and
kernel/extraction quality, owned memory and process RSS. The prior instrumented
60/180/360-second training observations remain attribution evidence; fresh
accepted timings omit the adjacency inspector. Preserve censored outcomes and
restrict any acceptance to this measured workload. C10 and E1 remain open.

## Reproduction

Use PyPy implementing Python 3.11 and a separately coordinated machine lease.
The corpus is now versioned; do not regenerate it or call it newly held out.

```sh
v2/.venv/bin/python -B -m v2.benchmarks.qs.c9.c9_long_confirmation confirm \
  --output "$PWD/v2/benchmarks/results/c9/long-reproduction" \
  --digits 50 --lease-seconds 3600
v2/.venv/bin/python -B -m v2.benchmarks.qs.c9.c9_long_confirmation summarize \
  --output "$PWD/v2/benchmarks/results/c9/long-reproduction"
```

Resume the same output directory after lease yields. Each process rewarms;
individual calls are admitted only when their complete reservation fits the
lease and aggregate. Record process identities when a matched group crosses a
lease boundary. Report the original study and this stage's expense ledgers
separately and their combined charged total. Neither incomplete groups nor
unsuccessful starts disappear from the resource account.
