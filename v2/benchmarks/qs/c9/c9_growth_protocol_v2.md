# C9.2 validation repair and new training freeze

The [original C9.2 protocol](c9_growth_protocol.md), configuration
`inputs/controls/c9_growth_training.json` and source commit `cae314a` remain
immutable. This correction repairs the validator's resource envelope
and failure capture. It also rejects the new growing class in SSS/coarse-worker
adapters that lack its ownership/checkpoint integration; serial arithmetic
is unchanged. All training arms, admission/scoring, durations, fresh
generation and promotion criteria remain as stated in that original
protocol. C3 explicitly approved a <=1200-active-second restarted training
lease after its current lease, preceded by <=90 seconds QA. Other leases
retain the 600-second default unless separately coordinated; all study totals
are unchanged. No fresh C9.2 inputs have been generated.

The first legacy DLP training call completed in 102.900501 seconds, with
102.951893 inclusive charged wall seconds and 102.195085 CPU seconds. The
next growing call reached C1's reused validator, which still asserted a
256 MiB owned-memory maximum. C9.2 had already declared 1 GiB for every arm.
This was a harness mismatch, not a reason to raise the frozen allowance.
The failed call's report was not published, so it supplies no usable timing
or completion result. Reserve its full 362 wall/CPU seconds.

The restarted study begins at **4634.658706 wall / 4588.290883 CPU** charged
seconds, repeats the entire eight-call initial grid, and retains the
original 14,400-wall/CPU-second aggregate cap. The two old inputs are training;
no selection is made from the failed attempt. Raw original logs and ledger
remain in the ignored `results/c9/growth-training/` directory. New captures
use `results/c9/growth-training-v2/`.

The new validator checks exact reconstruction, proper divisors, expected
factors, certainty labels, work and the actual frozen configuration memory
allowance. Unit tests require 700 MiB to pass under 1 GiB and fail under
256 MiB, and reject corrupt arithmetic/labels. Future validation failures
persist their unaccepted reports before raising.

Charged prior wall: 4634.658706497634. Charged prior CPU: 4588.290882917167. Both include the original
600-second QA reserve; no additional allowance or refund is introduced.

Confirmation retains the calibrated legacy SLP control from C9.1 (3M at
50 digits), so reusable-SLP growth can be compared directly against it.
The single challenger is selected across **grow, grow_aligned, half_grow,
graph_slp**, by completion then capped total cost then that listed order.
This replaces the original split SLP/DLP selection rule before the restarted
grid: reusable SLP is a candidate mechanism, and can win without DLP. No
C9.2 reusable-SLP outcomes have been observed. No extra arms or search budget
are introduced. The original protocol and failed freeze remain preserved.
