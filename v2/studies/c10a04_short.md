# C10a-04 bounded 60-digit follow-up

The user replaced the remaining multi-band campaign with this separate,
prospective comparison. The original frozen experiment and every receipt,
including its interrupted 60-digit call and full reservation, remain unchanged.
The original Stage 2 campaign is incomplete. E1 remains deferred and C10a-07
open. No production defaults, algorithms, frozen mechanisms or sources change.

## Question and candidate selection

Existing exploratory results nominate compact against control at 60 digits.
All five completed 60-digit pairs favored compact; those results include
documented host sleep and are not fresh validation. The 55-digit band retained
control under its original stability rules. Wide and 65 digits are deferred.

Commit the candidate, this protocol, runner and control before generating
three fresh certified balanced 60-digit inputs once. Exclude prior versioned
products and factors. The generator conditions primes on a large p-1 factor;
results are scoped to this corpus distribution. Commit the fresh corpus before
measurement. Three inputs crossed with seeds 7, 29 and 47 give nine pairs,
not nine independent inputs. Use outcome-independent shuffled arm order.

## Budget and outcomes

The entire follow-up has a ten-hour absolute wall deadline starting with
preparation, and a 36000-second aggregate wall/CPU reservation cap. Preparation
and committed-only QA retain at least 1800 prepaid seconds; generation is
separately reserved at 125 seconds. Each measured call has matched 1800-second
wall/CPU caps, 10^13 work, 1 GiB owned memory and the original 3 GiB RSS guard.
Warm both arms for at least three validated wall and CPU seconds using the
original bounded warmup. Inclusive complete-factor outcomes include setup,
collection, filtering, matrix, extraction and independent reconstruction.
Timeouts remain outcomes; capped cost never treats a failed factor as cheap.

Reserve whole pairs before starting. Reserve 600 seconds before the absolute
deadline for final reporting. SIGALRM interrupts at the measurement deadline;
the existing accounting records a full-cost censored call after it unwinds.
The ten-hour budget may therefore stop before nine pairs. No extensions,
selective retries, fresh regeneration, resource resets or automatic restart.
A failure stops for audit. Raw receipts and attempts stay folder-ignored.

## Fixed early checks

Inspect only after 3, 6 and 9 complete matched pairs, covering all three inputs
at each checkpoint. Correctness failures stop immediately.

- After 3 or 6 pairs, zero compact completions stops as uninformative under
  this time grant; do not increase the grant.
- After 6 pairs, stop with a preliminary positive trend only if compact
  completes all six calls, every paired capped cost favors compact, and each
  input's aggregate capped cost is at most 80% of control over its two seeds.
- Otherwise finish at nine pairs or the total deadline and report all outcomes.
  Mixed direction or insufficient completion is inconclusive; no extra data.

These are pragmatic descriptive stopping rules, not a sequential hypothesis
test or 95% confidence guarantee. The 20% threshold is an early-stop signal
criterion, not a production promotion floor. An early stop below nine samples
is preliminary and does not meet the repository's ordinary accepted timing
gate. Even nine pairs cover only three input clusters, without broad size or
configuration optimality. Publish completion/capped cost per input and arm,
uncertainty and censoring; production promotion is a separate stage.

## Commands

Use PyPy implementing Python 3.11. New orchestration lives outside
v2/benchmarks so adding it does not change the original runner's recursively
hashed benchmark set. It reuses the original checked arithmetic and receipts.

```sh
PY=/Users/nishanthmohan/code/factor/v2/.venv/bin/python
$PY -B -m v2.studies.c10a04_short generate
# Commit v2/benchmarks/inputs/corpora/c10a04_short_fresh.json.
$PY -B -m v2.studies.c10a04_short run
```

Both commands acquire the shared performance flock. User-granted sole CPU
priority removes turn-taking; the flock still prevents accidental overlap.
Read ledger.json without constructing Ledger while a measured process lives.
