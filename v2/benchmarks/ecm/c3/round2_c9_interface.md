# C3 / C9 growing-storage integration contract

The current C3 experiment consumes frozen C1 collectors. It does not include
C9's growing-storage implementation or claim combined-source acceptance.
C9 introduced that separate class at 8f3b866; later C9 study/checkpoint revisions
remain in its own worktree. The calibration source is pinned by
c3_round2_training_v1.json. Do not cherry-pick C9 changes into that freeze.

Three integration points need combined verification after the studies:

1. A growing SIQS child can have memory_bytes=None. PortfolioConfig must
   resolve it against the finite outer allowance, subtracting concurrent
   workspace and the existing metadata reserve. Serialize the resulting
   finite child configuration so implicit restore repeats the same ownership.
   C9 prepared this narrow constructor change in coordination with C3.
2. Before that resolution, ECM default-options calculation can encounter the
   same None child memory. Disable automatic extra ECM workspace in that case;
   retain explicit caller choices and validate coexistence normally. C9 also
   prepared this narrow change.
3. C3 schema 12 restores the full configuration when config is omitted. Its
   existing collector decoder selects legacy DLP solely from the presence of
   large_prime_bound. A growing collector also has that field, plus
   storage_policy='grow' and optional count limits. Select
   GrowingDoubleLargeSieveConfig for that storage policy before the legacy
   DLP/SLP branches. Preserve every serialized field; unknown policies must
   still fail validation. C9's branch predates this C3 schema-12 decoder, so
   its standalone resume tests do not establish this combined path.

For combined acceptance, use an actual growing collector with explicit finite
outer memory and C3 allocation. Check fresh construction, active SIQS pause,
implicit config restore, explicit config restore, cumulative work/CPU/wall,
resolved memory identity, retained partial relations and exact reconstruction.
Also retain legacy SLP/DLP restore and explicit finite-child behavior. Existing
C9 guards keep growing collectors out of unsupported SSS/coarse-worker paths.
These tests need both runtime implementations, not a fabricated collector shim.

The third point was identified by reading both branches during C9's shared
machine window and communicated to its task. It is a future integration
requirement, not a change to the current numerical policy or an invitation to
rewrite either source freeze. No merge or push is performed by this study.
