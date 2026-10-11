# C10a-05 state transitions and invariants

Keyword-only `SIQSConfig.solve_policy` defaults to `window` (legacy). `collection` applies
only to shared-store SIQS. Direct QS/MPQS and SIQS with unshared stores keep
per-window finalization. No collector arithmetic, polynomial order, residual
admission, matrix representation or filtering algorithm changes.

| State/event | Transition |
| --- | --- |
| Ordinary block or shared polynomial end | Finish pending solver/extractor first. Inspect new retained rows for zero parity with charged work. Otherwise apply existing row-growth and post-filter excess gates. Window completion permits the next polynomial, not an unconditional solve. |
| New zero-parity row | Bypass cadence; prepare/verify the original store, filter, solve and extract normally. Never return a factor from unchecked parity alone. The cursor advances past a square only after completed downstream work. |
| More families/width recovery available | Continue with the unchanged retained originals and cadence counts. Root/setup work remains charged. |
| Finite family/assignment search, stalled/trivial stop | Persist `pending_final_reason` before final processing. Bypass cadence/excess once at the actual collection stop. Commit `finished_reason` only after processing or an explicit matrix-memory refusal. |
| Full relation/atom/memory store | Existing `_finish_storage` forces retained rows once. Budget refusal leaves its completion flag false and pending solver/extractor intact. |
| Setup/transition memory refusal with retained rows | Try the same bounded terminal drain; an impossible matrix reservation is an explicit terminal memory stop. |
| Work/wall/CPU/cancellation refusal | No final solve outside the budget. Return unresolved cofactor, keeping the pending terminal reason and resumable downstream prefix. |
| Repeated run / serialized restore | Completed stops do not repeat work. Pending terminal work resumes before collection. Restore reverifies originals, rebuilds roots and replays matrix/solver/extractor prefixes against cumulative allowances. |

An interrupted preparation/filter pass has no published matrix; resume rebuilds
it and charges the rebuild. Existing solver and extractor progress retains its
current prepared store. Collection cannot append rows while downstream work is
pending. Atomic identities, original exponents, signs, square corrections and
combined-row ownership remain unchanged. Singleton pruning affects the current
matrix only; no original row/atom is discarded.

`zero_checked_count` is an append-only row cursor. A refused row inspection
does not advance it. At a detected square it stays pending until downstream
completion; successful extraction attempts move it to the retained row count.
Unshared polynomial changes reset it. The finite family, polynomial, stall,
trivial-dependency, work/time and owned-memory bounds continue to apply.

Schema 6 stores `solve_policy`, `pending_final_reason` and the zero-row cursor,
plus the existing store and checked progress. Its collector type is decoded
from the saved configuration (ordinary, fixed DLP or growing DLP). Schemas
1–5 retain their window semantics; new window-policy snapshots retain the
legacy schema/configuration shape. A monotone resource extension cannot change
the scheduling policy. Existing terminal-stop extensions remain supported.
Nested portfolio checkpoints omit the default policy field to preserve their
legacy configuration identity. An empty growing-store restore now rounds its
temporary encoding cap up to the existing 4096-byte API minimum; both legacy
and schema-6 snapshots exercise this correction without format changes.

The scheduling decision changes which prefix is solved and can change work,
relation counts and stopping factors. It promises exact results and charged
continuation, not identical work to legacy execution. Public callers continue
to validate proper divisors, distinguish probable/proven primes, and reconstruct
all factors and unresolved cofactors.
