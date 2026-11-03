# C10a-01 primary execution revision 2

The first execution source (`execution_v4_sources.json`, commit 2bbc654) ran
three declared primary reservations. Both balanced 60-digit cells completed
67,109,376 positions and returned exact results. Their worker receipts show
8,654,079,332 and 15,954,850,282 total charged-work units, respectively,
but the generic ledger settled each as zero because the worker returned
`total_work` without the `work` key expected by `Ledger.call`. The next uneven
60-digit challenger reservation was interrupted while stopping this flawed
stage and retains its full 60-second and 2,607,786,404,459,813,954,784-work
charge. This is a harness accounting failure, not native algorithm exhaustion.
None of the three entries is overwritten or silently accepted as clean v4
calibration evidence.

Revision 2 binds a separate source manifest and `calibration-v4-r2` output.
Every declared primary cell is rerun in its original order; there is no
outcome-based selection. The new worker total is checked against prewarm,
measurement and postwarm components, then is supplied as `work` to the ledger.
A focused regression settles the exact 100+300+200 example as 600, which the
old code would have settled as zero. All time/work, memory, input and position
limits in `experiment_v4.json` remain literal. The old receipts and ledger
hashes are checked before this stage starts.

The failed release reserved 180 wall AND CPU seconds, including the interrupted
cell, and the interrupted cell retains its 2,607,786,404,459,813,954,784-work
reservation. The two settled receipt work amounts (24,608,929,614 combined)
are reported in addition to their historical zero ledger fields. Revision 2
receives a separately identified prospective 3,400-second primary/diagnostic
envelope. The maximum combined whole study envelope therefore becomes
42,480 wall AND CPU seconds (11 h 48 min); replacing the original 3,600-second
allocation by its actual 3,124 reserved seconds gives 42,004 seconds.
These are reservations, not expected runtime. The later stage grants remain
unchanged and are not running in this primary stage.

No extra source or arithmetic policy changes accompany the revision. The
original paired observations remain visible as invalid-accounting evidence.
The new primary pair may be compared only within revision 2. The remaining
diagnostic groups require their order-aware adapters and their own committed
source release before measurement. A lower-band supplement needs Astra xhigh
review of its design before use: session turn metadata records Astra high for
that supplement, whereas the main 60–70 plan and implementation used Astra
xhigh. The current execution runner uses Sol under the user's clarification
that Astra xhigh was required for research, implementation and experiment
design, not running.
