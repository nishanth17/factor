# Balanced complete-factorization measurement preparation

This extension measures complete calls after the earlier 30-second balanced
probes refused a QS checkpoint (60 digits) or had no automatic QS fallback
(70 digits). Increasing the clock alone does not address either condition.
It does not change production defaults or rerun E1.

Use committed-only source snapshots. S0 is production mainline
`0efaba1ee40f4fa9e68e0a1acc3bac830a5540fe`; AUTO is the candidate automatic
policy. The source-specific omitted ECM defaults are retained: S0 fixed32,
AUTO escalating140 at 60 digits; both fixed32 at 70 digits. This is a bundle
comparison, not identical-manifest attribution. Explicit serial control
ownership prevents a newer benchmark driver from silently enabling AUTO.

`completion.py --source-root ROOT --n N` produces literal configs for freezing
before execution. At 60 digits, R64 preserves calibrated QS numerical policies
and increases its checkpoint cap to 64 MiB and owned memory to 2 GiB (parent
2 GiB + 16 MiB). At 70 digits, Q70 explicitly supplies the retained C10a-01
10,500-prime coupled QS configuration, 64 MiB checkpoint, 4 GiB parent ownership.
Stock 70-digit automatic routing has no QS fallback. Thus S0[R64]/AUTO[R64]
and S0[Q70]/AUTO[Q70] are funded configurations, not stock-default claims.

The initial question is whether complete factoring finishes and whether ECM
savings survive fallback. Use the previously exposed certified
`c10a01_60_balanced` and `c10a01_70_balanced` moduli, algorithm seed 7919,
never passing factors or certificates to a child. Two sequential cold cells
per size, alternating arm order across sizes, precede further measurements.
Each cold cell is a fresh PyPy interpreter and includes interpreter/import,
pool construction, result construction, shutdown and controller-observed exit.
No warmed pool is relabeled cold. Primary cold CLI harness latency is distinct
from inner library latency; the production argument parser is not included.

| Lane | Cells | Per-cell wall | Aggregate CPU | Work |
| --- | ---: | ---: | ---: | ---: |
| QA/configuration | finite commands | 120 s total | 1,200 s total | n/a |
| 60-digit pilot | S0[R64], AUTO[R64] | 1,830 s (1,800 call) | 18,400 s | 10^18 |
| 70-digit pilot | AUTO[Q70], S0[Q70] | 10,840 s (10,800 call) | 108,400 s | 10^18 |

The pilot ceiling is 25,460 wall seconds and 254,800 aggregate CPU seconds,
including QA, controller, imports, validation and draining. These are maximum
reservations, not expected times. Retained direct-QS observations were roughly
205 seconds at 60 digits and 7,603 seconds at 70 digits; portfolio seeds and
whole-call costs differ. Expect minutes for the first pair and potentially
four or more hours for the second. The host has 24 GiB RAM: automatic ECM has
an 8 GiB aggregate owned cap; QS grants above include parent coexistence.
Require 12 GiB publicly reported available headroom before each heavy cell.
Owned estimates are not an RSS theorem. Report parent and child high-water
RSS separately; their sum is not simultaneous aggregate RSS. Evidence is
capped at 128 MiB; each stdout/stderr at 4 MiB. No profiling during clean runs.

Hold the existing shared performance-lock inode with nonblocking acquisition
throughout every runtime command. Busy means queued, without polling. Bound
owned process groups and verify their exit. A missing report retains its full
CPU reservation; incomplete calls, failures and charges remain in the ledger.
Stop on lifecycle/arithmetic/recorder faults. Do not automatically enlarge and
retry. A resource refusal requires a new diagnosed protocol identity.

One exposed input per size permits exploratory point comparisons only. No
population median/p90 or established speedup follows. If both 60-digit pilot
arms complete within 600 seconds each, reserve an optional follow-up before
execution: nine cold repetitions and nine warm timings per arm on that same
input (three separately retained contexts with three timings each), measured
warmup at least three validated seconds. Its literal caps derive from 1.5 times
the slower pilot, with startup/drain reserves, under an additional 43,200 wall /
432,000 aggregate CPU second maximum. If this is not informative, do not spend
it. A warm context pays its initial validated complete call outside timed warm
samples. Nine repeats address timing variability, not modulus independence.

Do not automatically repeat multi-hour 70-digit calls nine times. First assess
whether the policies actually engage different work and whether completion
cost permits an informative matched scope. Any subsequent population claim
needs independent balanced moduli, matched ceilings, uncertainty, canary drift
checks and a disjoint fresh confirmation freeze. The earlier frozen random and
ECM-enriched results remain separately labeled historical integration evidence.
