# C9.2 pre-confirmation correctness and restart bridge

This addendum is prepared before any C9.2 fresh corpus generation. It preserves
the original training protocol and captures. The selected configuration and
ranking are not changed by the repairs below. The final selected control pins
the corrected source before generation; the versioned bridge records both
source hashes and the original selected-control hash.

Review found that uncapped checkpoint encoding with an explicit total-memory
cap must count currently retained preparation and solver storage, in addition
to original relations. The repair counts those live containers, not a historical
peak after they have been released. Restore likewise reserves the actual
encoded input's existing eight-representation allowance while rebuilding.
That reservation is released only after checked reconstruction; the returned
job retains its original public `checkpoint_bytes=None` policy and execution
identity. An explicit caller cap can reject oversized input before JSON decode.
Repeated restore, unchanged child allowance and cap refusal are regression cases.
Legacy explicit-checkpoint paths and serialized relation semantics are unchanged.

The joint training restart initially failed because the exclusive manifest
writer was called even after verifying an existing identical manifest. No
experiment call began. The identical first-segment manifest was preserved,
and the launch was conservatively charged two seconds without moving the
shared-window deadline. The runner repair skips existing identical manifests
and completed-phase receipts; conflicting manifests still fail. Repeated
metadata-only restart has a regression test.

The source audit compares ASTs and requires all functions outside these four
to remain identical: checkpoint `_reserve_growing_checkpoint`,
`_checkpoint_restore_config`, `restore_job`, and training `train`. Every timing
configuration has a positive explicit checkpoint-byte allowance; the changed
dynamic reservation branches are inactive there. Explicit-size diagnostics
also retain their original route. No sieve, residual admission, splitting,
graph, filter, solver, extraction or selection function changes in this bridge.

Original individual/joint training can be reproduced at `9707dab` using its
versioned controls. Its segmented manifest restart needs the documented
preserve-and-recreate step; the final runner repairs that metadata behavior.
Fresh confirmation uses the later committed selected-control source identity.
There is no retrospective relabelling of old captures as new-source timings,
and no new candidate or measurement-extension rule is introduced.

The full working suite hit a 65-second subprocess timeout in the slower
session; its focused new tests had passed. The failed attempt is preserved
and conservatively charged 123 wall/CPU seconds because the wrapper did not
collect its terminated descendant's CPU. A process inventory confirmed no
remaining test child. Final QA uses process-group termination and a full
remaining-allowance deadline. Before fresh generation, **150 seconds are
transferred from the unspent original 14,400-second cap to final QA**, making
its total prepaid reserve 750 seconds. Unused QA reserve is not refunded.
The aggregate cap, candidate grid and timing extension rules do not change.
Full lint and all 20 focused growth/runner/A7 tests pass after the repair;
complete working and committed-only suites remain mandatory final gates.
