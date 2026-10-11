# C10 collection strategy experiments

Read [research and ranked gaps](research.md) and the
[premeasurement protocol](protocol.md). C10 remains open. Experimental source
and controls are versioned; raw captures live in ignored `results/c10/`.

Commands from the repository root use PyPy implementing Python 3.11:

```sh
pypy3 -B -m v2.benchmarks.qs.c10.study freeze \
  --c9-revision FINAL_C9_COMMIT --decision SETTLED_C9_DECISION
pypy3 -B -m v2.benchmarks.qs.c10.study train --digits 30
pypy3 -B -m v2.benchmarks.qs.c10.study train --digits 40
pypy3 -B -m v2.benchmarks.qs.c10.study train --digits 50
pypy3 -B -m v2.benchmarks.qs.c10.study diagnose
pypy3 -B -m v2.benchmarks.qs.c10.study select
pypy3 -B -m v2.benchmarks.qs.c10.study generate
pypy3 -B -m v2.benchmarks.qs.c10.study confirm --digits 30
pypy3 -B -m v2.benchmarks.qs.c10.study confirm --digits 40
pypy3 -B -m v2.benchmarks.qs.c10.study confirm --digits 50 --lease-seconds 3900
pypy3 -B -m v2.benchmarks.qs.c10.study summarize
```

Freeze only after the C9 handoff. Measurement commands fail immediately if
another owner holds the performance lock. Heavy QA uses that same lock.
Each new process warms all active arms again; captures are resumable and
completed groups are not rerun. Full group reservations protect matching near
aggregate/lease limits. A pending group may require another lease. An aggregate
cap is final for this tranche; rerunning cannot replenish it.

The training runner uses known C1/C9 inputs. Selection freezes before the
generator creates fresh certified inputs; no generation is needed for bands
without a selected challenger. Historical C9 source manifests remain immutable.
The separate A7 C10 integration pin updates loader compatibility, not historical
timing evidence. No C10 measurement is an E1 portfolio or worker comparison.
