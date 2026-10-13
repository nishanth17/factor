# C10a-04 65-digit feasibility spike

At the user's request, compare one once-generated fresh certified balanced
65-digit product at seed 7. Candidate is the existing 65-digit compact policy:
8500 eligible base entries, nine A factors, pool 64, selected odd multiplier,
half-width 65536 and a capped 64-variant Gray prefix. Control is the accepted
50/51-digit bundle explicitly applied at 65 digits, not the automatic reference
65-digit route. Both retain 1 GiB owned memory, 3 GiB RSS guard and 10^13 work.

Each arm receives 1800 wall/CPU seconds, including all setup and classification,
after at least three validated PyPy wall and CPU warmup seconds. Commit plan,
runner and runtime/source identities before generating; commit the corpus
before measurement. Exclude all previous versioned products and prime factors.
Generation is capped at 125 seconds and warmup at its finite ledger allowance.
Arm order is compact then control, fixed before outcomes. No repeated seeds,
extra inputs, retries, regenerated corpus or allowance resets are authorized.

A 90-minute absolute deadline includes preparation, with 600 seconds reserved
for reporting. Whole matched calls must fit remaining time and reservations.
An alarm unwinds accepted work at the measurement deadline; the runner then
retains its interrupted accounting. A failure stops for audit, never restart.
Read ledger JSON directly during measurement; never construct Ledger for a
live process. No CPU coordination, lock polling or scheduling with other chats.

This is a feasibility observation on one conditioned input. Completion can
justify planning a bounded larger study; a timeout says only that this input
and configuration did not finish under the grant. Preserve timeout cost,
work, stages, rows, dependencies, memory and identities. No statistical
speedup, stability or general 65-digit acceptance claim follows from one pair.
Default extension through 65 must be a separate explicit decision after the
result; any interpolation across 61–64 is unmeasured policy extrapolation.
70-digit feasibility remains open on the roadmap and is not part of this run.

Commands from the repository root, using PyPy Python 3.11:

```sh
PY=/Users/nishanthmohan/code/factor/v2/.venv/bin/python
$PY -B -m v2.studies.c10a04_spike65 generate
# Commit the once-generated corpus before run.
$PY -Bu -m v2.studies.c10a04_spike65 run
```
