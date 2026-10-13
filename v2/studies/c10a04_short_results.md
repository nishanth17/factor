# C10a-04 bounded 60-digit result — 12 October 2026

Compact met the predeclared six-pair preliminary trend stop. It completed
6/6 fresh calls, versus control's 5/6, with lower capped cost in every pair.
Total capped time was 4588.288 seconds against 7187.725 seconds, a 36.16%
reduction. Each of the three inputs saved at least 20% over its two seeds.
This supports compact as the candidate for a separate, size-scoped promotion
review; it does not establish production acceptance or change defaults.

## Prospective design and scope

The [protocol](c10a04_short.md), runner and control were committed in
`405c11eced0a7b15eacb8ccde7367531026a08dc` before once-only fresh generation.
The certified fresh corpus was committed in
`7a738ca9d03a5a13ee52e95e238f8cd83b6bbff7` before measurement. Prior products
and prime factors were excluded independently. Candidate nomination used
five exploratory 60-digit pairs from the paused original campaign; their
sleep-overlapped timings contribute no fresh acceptance evidence here.

Three independent balanced 60-digit products were crossed with seeds 7 and
29 before stopping. Seed 47 and the last three planned pairs were not run,
as required by the positive stop. These are three input clusters, not six
independent products. Generation conditions primes on a large p−1 factor;
extrapolation to arbitrary products, other sizes or other mechanisms is open.
Arm order was shuffled independently of outcomes. Both arms received over
three validated wall and CPU seconds of PyPy warmup. Runtime was PyPy 7.3.23,
Python 3.11.15, using the existing Python-integer backend.

The coupled compact arm uses 5000 eligible base entries, selected multipliers
3/1/13 for inputs 0/1/2, eight A factors and pool size 64. Actual exclusive
base bounds were 104790/107844/103682. Control uses multiplier 1, exclusive
bound 30000, 1651/1534/1631 entries, six A factors and pool size 32. Both use
half-width 65536. This is a bundle comparison; it does not isolate which
parameter caused the difference.

## All observed outcomes

Times include configuration resolution, factoring stages and checked result
construction. Capped cost retains at least the full 1800-second grant for an
incomplete factorization; the timed-out call is included, without a retry.

| Input | Seed | Control seconds | Compact seconds | Control / compact complete |
| --- | --- | ---: | ---: | --- |
| 0 | 7 | 965.042 | 605.631 | yes / yes |
| 1 | 7 | 1145.509 | 998.387 | yes / yes |
| 2 | 7 | 862.820 | 631.353 | yes / yes |
| 0 | 29 | 1005.381 | 610.209 | yes / yes |
| 1 | 29 | 1800.001 | 1012.752 | wall limit / yes |
| 2 | 29 | 1408.972 | 729.956 | yes / yes |

| Input | Control capped seconds | Compact capped seconds | Saving |
| --- | ---: | ---: | ---: |
| 0 | 1970.423 | 1215.840 | 38.30% |
| 1 | 2945.511 | 2011.139 | 31.72% |
| 2 | 2271.792 | 1361.309 | 40.08% |
| Total | 7187.725 | 4588.288 | 36.16% |

Checkpoint 3 continued. Checkpoint 6 stopped with
`stop-preliminary-consistent-compact-trend`: compact completed all six,
every paired capped cost favored compact, and each input's aggregate cost
was at most 80% of control. An independent read-only audit recomputed those
conditions and matched both checkpoints and the final summary to every
receipt hash. All proper divisors, multiplicities and unresolved cofactors
reconstruct the certified inputs. Runtime `probable_prime` labels remain
distinct from the independent corpus certificates.

## Resource and stage evidence

Each measured call retained matched 1800-second wall/CPU, 10^13 work,
1 GiB owned-memory and 3 GiB process-RSS limits. Control's one timeout exceeded
1800 wall seconds by about 1.1 milliseconds and retained its full cost.
There were no command failures, interrupted calls, retries or pending
reservations in this follow-up. Measurement stopped normally at
21:16:07 PDT on 12 October, about 3 h 16 min after measurement began and
3 h 25 min after preparation began, well inside the ten-hour deadline.

| Metric across six calls per arm | Control | Compact |
| --- | ---: | ---: |
| Measured CPU seconds | 6918.191 | 4462.551 |
| Charged factoring work | 1104580096336 | 1232429487461 |
| Maximum recorded owned workspace, MiB | 385.225 | 584.692 |
| Additional count-helper reserve, MiB | 0 | 15.267 |
| Collection stage seconds | 6543.854 | 3976.129 |
| Preparation stage seconds | 9.860 | 96.864 |
| Filtering stage seconds | 2.988 | 57.346 |
| Elimination stage seconds | 5.814 | 6.892 |
| Extraction stage seconds | 0.047 | 0.108 |

Compact used 11.57% more charged factoring work and more owned memory,
despite lower elapsed and CPU time. Stage counters describe the executed
paths, including control's incomplete call; their sums omit outer setup,
roots and other overhead and are not an independent profile or causal
parameter attribution. The process high-water RSS observation reached
402.188 MiB. It is cumulative across the shared sequential process, so it
cannot be assigned independently to each arm. Recorded owned workspace plus
helper reserve stayed below 1 GiB; no storage refusal occurred.

The ledger retains 23595 wall and CPU seconds of aggregate reservations,
including 1925 setup and 21670 confirmation seconds. All twelve measured
calls retain their full 1805-second reservations; successful calls did not
refund time. Total charged work is 62338599908144, including conservative
setup and warmup. Saved output occupies 111774 bytes. All phase, total,
work and disk allowances pass. The original campaign's ledger and nominated
receipt identities still match the plan; its interrupted reservation was
not recovered or reset.

## Evidence identities and checks

Raw receipts, attempts, checkpoints, ledger and console remain locally
preserved under `v2/benchmarks/results/c10a04/short/` and `v2/local/`.
They are folder-ignored rather than published as repository contents.
The local final seal records every output-file hash, all per-call stages,
work, configuration, memory and certainty, and the independent decision.

| Artifact | SHA-256 |
| --- | --- |
| Frozen follow-up control | `46f520640270d046ae203c0815f2d12cbf3bd10c75fd3aa928cdd39a5faecf83` |
| Fresh corpus | `b31355b4b218b606ccb7029e6f61ea44ff3f0c031c92d4868b00cc65567c2a84` |
| Checkpoint 6 | `95c99f68b4fb01562ed6483c67fc5b18108d3bb3ec560a5e258c42de9dddfe32` |
| Final summary | `438bcbd2b555f09d1ae06e747aec1bb9360385d4fe44ee376b4f73ef7de54e33` |
| Final ledger | `0c4e3ca09938e4911975278819fd681482afb8ee097deb68c69117142255a64d` |
| Local final seal | `372b8a1388e3ded177ab7309e670f467e159c58379db81d873d81522eb93c380` |

Before measurement, a committed-files-only checkout passed
`make -C v2 test` (764 tests), `make -C v2 lint`, imports of all 155 live
benchmark modules plus the new study, and all 249 original frozen-source
checks. Historical input snapshots were excluded from the live module import
inventory. Fresh certificates, exclusions and committed corpus bytes were
independently verified. Final changes are documentation only: whitespace,
links, committed-source identities and sealed receipt hashes were rechecked;
no additional timing or broad test campaign was needed.

The protocol records generation and run commands. This completed run must
not be restarted. Its read-only audit command, from this worktree, is:

```sh
PYTHONPATH=. /Users/nishanthmohan/code/factor/v2/.venv/bin/python -B \
  v2/local/c10a04-short-monitor.py
```

## Verdict and promotion handoff

The preliminary fresh trend favors compact for this 60-digit distribution.
The user-authorized early stop provides the requested bounded answer without
collecting the last three pairs. It does not meet the ordinary nine-sample
accepted-timing gate, estimate sequential significance, establish global
optimality or confirm a general production default. No confidence interval
is claimed from three input clusters and an outcome-dependent stop.

Keep current defaults. A separate promotion stage must explicitly decide
whether this limited evidence is sufficient for a 60-digit opt-in or requires
a new prospective acceptance plan, then preserve the exact bundle, finite
budgets, reconstruction, certainty and resume contracts in any integration.
No further measurements are authorized automatically. The 55-digit original
stability verdict retains control; 65 digits and wide remain deferred.
The original multi-band Stage 2 stays incomplete, E1 deferred and C10a-07 open.
Nothing is merged or pushed.
