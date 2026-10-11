# C9.2 constraints and transfer audit

The [pinned comparison](c9_research.md) is the full mechanism account. Its
source archives were rechecked at the same revisions before implementing
allocator growth. “Other implementations have no limits” is incorrect.
The relevant distinction is a growable full-SIQS relation history versus
admission, cofactoring effort, path, representation and stopping policies.

| Mechanism / limit | Actual source behavior and v2 transfer |
| --- | --- |
| Relation-store growth | YAFU full SIQS doubles vertex arrays (`filter.c:1177`) and separately doubles in-memory relation records (`tdiv.c:1117`). Its default 70-digit in-memory cutoff selects memory versus disk; 32,768 is initial allocation, not a relation ceiling. msieve doubles vertex arrays and grows loaded relation arrays by 1.5x; originals live in the savefile. FLINT doubles relation/factor storage and grows lookup entries by 1.4x. Java and Yamaquasi use growable containers. **Implement:** PyPy's own list/dict growth, retain partial history, no implicit FIFO or record quota. Explicit caller limits still work. Native packed byte sizes are not Python reservations. |
| Not every YAFU path grows | `siqs_test.c:179–207` dispatches inputs below 115 bits to `tinysiqs`; this path has `relation_list[4*MAX_RELATIONS_TINY]`, with `MAX_RELATIONS_TINY=512`. `cofactorize_siqs.c:1334–1378` drops unmatched partials when a fixed-depth hash bucket is full and stops if the full/partial regions meet. The same cofactorizer is used in GPU residual splitting. This is a counterexample to a global “no storage bounds” claim, **not** the full SIQS path relevant to the 50-digit study. |
| Cycle length | Both YAFU and msieve full graph extraction have two 100-element traceback arrays; YAFU's `filter.c:3776ff` rejects overly long paths. v2's explicit 256-original-atom cap also loses yield. Preserve it for this tranche: retained evidence records zero such rejections; report the counter. Removing it requires separate arithmetic/scratch/provenance work, not a copied native threshold. |
| Admission versus scoring | Common SLP/DLP endpoint limits in YAFU/msieve avoid retaining primes that cannot join DLP edges. Yamaquasi's nominal score product is not its hard intake product: recovered residuals can be accepted to endpoint product L². v2's high SLP allowance is not mathematically wrong; uncapped history can recover its extra matches. Test alignment separately after growth, keeping DLP score/product coverage fixed. C10 owns broader candidate recovery. |
| Split effort and certainty | Upstreams use different probable-prime screens and native rho/SQUFOF/ECM regimes. v2's sampled losses show no failed splits, and sieve/root work dominates. Retain proven endpoints and current per-residual effort. Remove the unrelated implicit cumulative split-call quota in growing mode; each attempt remains charged. No residual-certainty shortcut is justified. |
| Reusable SLP matching | For k originals sharing an endpoint, a star gives k−1 independent LP cancellations; disjoint pairs give floor(k/2). v2's ordinary SLP matcher uses disjoint pairs, but C1's existing `graph_slp` route already supports the star. Evaluate that route with growth and the same SLP candidate domain. Do not call it an entirely missing v2 algorithm. |
| Final matrix size and oversampling | Yamaquasi `siqs.rs:160–162` truncates final relations to factor-base size plus `MIN_KERNEL_SIZE`; msieve targets factor-base size plus 96, with an optional relation stop. These are collection/matrix policies, not evidence that every retained partial has a fixed storage budget. v2's 32-row growth/excess and forced polynomial-window solves are heuristic. Keep them frozen here; a cadence change may save repeated preparation or delay useful dependencies. |
| Other retained restrictions | Exact residual size supports the proven endpoint contract; factor-base/input/position representation guards and finite caller work/time remain. The 2 MiB verification cache is evictable derived data, not partial-relation ownership. Streaming family/Gray work allowances and scoring choices remain C10; portfolio allocation remains C3. No claim that all settings are calibrated or optimal. |

`GrowingDoubleLargeSieveConfig` makes absent storage quotas explicit in the
API and schema 5. Python containers grow on insertion, so adding a second
manually doubled “capacity” array would duplicate runtime allocation and
potentially increase copying. Conservative owned-memory accounting continues;
`None` means no user quota, not a huge preallocation. No disk-spill algorithm
is claimed. Under an explicitly bounded outer portfolio, the child receives
that parent's remaining space; it cannot allocate the parent's space twice.

Allocation failure during a staged plan leaves the published store intact.
An allocator exception or interruption in the final multi-container publication
can leave that *live object* unusable: it is marked failed and cannot collect,
solve or publish another checkpoint. Recovery uses the previous complete
checked checkpoint and charges replay. This does not promise recovery from
an operating-system kill or successful checkpoint creation after out-of-memory.
Cancellation/work/memory refusals before publication remain resumable.

Investment order is (1) storage growth, (2) separately measured endpoint
admission, (3) existing reusable SLP and historical half-base combinations.
Splitter replacement and solve-cadence redesign are deferred by current loss
attribution. This ranking can be rejected by complete outcomes; larger graphs,
more cycles or smaller matrices alone cannot accept a candidate.

## Subsequent literature verification

The [11 October transfer verification](c9_literature_verification.md) checks
current source revisions and records implemented, evaluated, deferred and
previously unmapped mechanisms. Asymmetric ordered LP bounds were missing
from the earlier explicit comparison; recovered-product admission, target
feedback and cycle-weight optimization remain untested. Only the named C9.2
tranche is accepted. This addendum changes no frozen source/configuration or
historical result and does not assert exhaustive adaptation or optimality.
