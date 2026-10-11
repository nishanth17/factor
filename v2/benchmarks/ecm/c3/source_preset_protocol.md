# Upstream ECM schedule transfer screen

This separate experiment compares four source-derived schedules through the
public protected-preset API against optimized fixed32 and the already committed
compact64/wide64 model. It does not alter the completed calibration or the
ongoing complete-call comparison. Numerical tiers are fixed in
`execution/ecm_presets.py`; source meanings and backend limitations are in
[source_presets.md](source_presets.md).

## Assignment and bounds

Use the first revealed round-two subject in each of these strata:

| Modulus digits | Classes |
| --- | --- |
| 30 | balanced, uneven_12, uneven_14 |
| 40 | balanced, uneven_12, uneven_16 |

The six arms are fixed32, fitted, yamaquasi_auto160, yamaquasi_ecm64,
alpertron20 and gmp_ecm20. This covers a cheap native automatic pretest and three
stronger standalone ECM prefixes. The larger Alpertron25 and GMP-ECM25 presets
remain explicit API options; they are not silently added to this assignment.

All arms use the original C1 SIQS settings and identical cumulative protected
allocation, 288 MiB outer memory, 10^13 work and 5/30 seconds wall/CPU. Each
transferred preset is constructed from the protected32 configuration using
`with_ecm_preset`, which preserves the exact SIQS and allocation objects.
Capture and compare the full SIQS settings, outer cap and allocation for every
matched group; reject mismatches. Larger curve workspace must fit the same cap.
No factors, hidden labels or observed outcomes configure dispatch.

Start with nine samples, using seeds 17,43,89,113,151,181,223,269,307. A declared
uncertainty extension adds seeds 347,389,431,479,523,569,617,661,709, then at most
751,797,839,887,929,977,1021,1069,1117. Rotate arm order by sample/input and
reverse alternate samples. Complete six-arm groups remain in one acquired
lease and prepay all service grants plus one second per arm: 36 seconds for
30 digits, 186 for 40. The first block is exactly 324 calls; full extensions
reach 648 or 972. No calls or favorable incomplete groups may be added.

Every measurement process warms each arm in both bands for at least three
validated CPU and wall seconds on the corresponding balanced subject, and
warms forced SIQS separately. Each warmup task prepays 60 seconds and checks
its 25-second continuation threshold before another full call. A failed
warmup or invalid factorization makes this phase inconclusive.

## Finite ledger and interpretation

Reuse the frozen whole-group capture, hash and prepaid ledger implementation
inside an isolated process context that changes only its benchmark manifest
and arm enumeration. That context always restores the earlier globals. Pin the
reused runner and original manifests along with this runner, tests, presets,
source table and procedure in `c3_source_presets_v1.json` before measurement.
All pinned files must be committed. The original training identity remains a
nested dependency; the group comparison identity is this separate manifest.

The entire screen has a fixed 3600-active-second cap including every acquired
warmup, measurement, capture and assessment. Individual measurement leases are
60..1200 seconds, with two seconds reserved for cleanup and 183 seconds retained
for the next assessment. Each assessment prepays a hard 180-second grant.
Interrupted grants stay charged. A partial matched group, interrupted report,
correctness error or changed capture permanently ends the study inconclusively.
Only hash-checked complete groups can resume; only a hash-recorded assessment
can authorize the next sample block. Coordinate every shared machine handoff
with C9. A phase that cannot fund its next full grant receives no extension of
its time allowance.

Assess only complete assignments at 9, 18 or 27 samples. Unresolved output pays
at least the full service cap. Use the existing equal-band/class/input weighting
and 10,000 paired bootstrap draws with seed 2026101109 against both controls.
A promising candidate has positive point CPU improvement against both,
weighted wall ratios <=1.05 and no class completion loss greater than five
percentage points. A stable promising candidate additionally has positive
95% lower endpoints and interval widths <=20 percentage points against both.

If any promising candidate is unstable, extend the entire matched assignment
from 9 to 18, then at most 27 samples, provided the unchanged phase cap funds it.
Otherwise stop. At the final stop select the stable promising candidate with
lowest weighted CPU cost, if any; ties follow declared preset order. No
selection follows an incomplete assignment or unstable final candidate.
This is revealed training screening, with one independent subject per stratum;
bootstrap seed repetition cannot establish broad input-population uncertainty.
Neither multiplicity across four candidates nor absent easy-factor strata is
resolved by this screen. A selected schedule still needs a separately frozen,
disjoint, broadly stratified fresh confirmation before any acceptance claim.

Generated captures and detailed receipts remain ignored. Publish concise
complete counts, validated output, cost/coverage differences and limitations.

    pypy3 -B -m v2.benchmarks.ecm.c3.source_preset_study measure CAPTURES --lease-seconds 1200
    pypy3 -B -m v2.benchmarks.ecm.c3.source_preset_study report CAPTURES
