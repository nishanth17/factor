# C3 round-two fitted candidate

All 2,052 frozen calibration assignments completed and reconstructed exactly;
no root curve was censored. Five coordinated processes charged 2,893.780
active seconds against the original 3,600-second calibration cap, including
validated reacquisition warmup and captures. This is instrumented training,
not accepted performance evidence. No fresh C3 round-two inputs exist.

The offline recurrence selects these ordinary finite tiers at call entry:

| Input band | Ladder | B1 | B2 | Independent curves | Conditional projected CPU |
| --- | --- | ---: | ---: | ---: | ---: |
| 30 digits | compact | 2,000 | 50,000 | 64 | 0.09461 s |
| 40 digits | wide | 2,000 | 147,396 | 64 | 2.82366 s |

Each band has nineteen independent subjects entering the root-ECM fit; seeds
are repeated observations of those subjects. Early preprocessing solves stay
in complete evaluation. Children inherit the configured schedule and share
the cumulative allowance; this table adds no per-child size dispatcher.

The other fitted ladder prefixes were wide32/deep20 in the smaller band
(projected 0.13226/0.13762 seconds), and compact64/deep20 in the larger band
(2.85278/2.98166 seconds). These projections omit some boundary bookkeeping,
use a matched no-ECM SIQS-cost proxy and have the declared synthetic-prior and
recursive-recovery limits. See [model interpretation](round2_model_limits.md).
They do not establish a speedup or justify changing production defaults.

The table and hashes of every calibration capture/ledger are committed in
[c3_round2_fitted_v1.json](../../inputs/controls/c3_round2_fitted_v1.json)
before any complete-call comparison. The primary calibration source/input
freeze remains [training v1](../../inputs/controls/c3_round2_training_v1.json);
its source bytes and assignments have not changed.

Before the first uninstrumented comparison, an additive measurement procedure
will warm both size bands and keep each matched control/fixed32/fitted group
inside one process/lease. This addresses an observed session CPU-rate change
without changing candidates, the 57x9x3 assignment, gates or the 3,600-second
comparison cap. The original comparison procedure has collected no data.
The fitted table must beat both controls under the predeclared complete-call
gates before a separately frozen fresh confirmation can begin.

Reproduce calibration/fit from the recorded source freeze:

    pypy3 -B -m v2.benchmarks.ecm.c3.round2_train calibrate RESULTS --lease-seconds 600
    pypy3 -B -m v2.benchmarks.ecm.c3.round2_train fit RESULTS

Use a new output tree for reproduction, retain all lease receipts, and do not
overwrite the committed selection. Required corpora, controls and immutable
baselines are versioned. Raw evidence remains under the ignored
results/c3/round2/calibration-v1/ tree. Full working-tree verification passed
621 PyPy tests and lint before the final calibration process.
