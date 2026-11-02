# C10a-01 v4 diagnostic execution release

This release follows the 16 completed revision-2 primary prefixes. Their ledger,
analysis and individual receipt hashes are mandatory for execution; a source-only
checkout can verify the code without pretending ignored raw evidence is present.
No primary outcome is changed or retried. The earlier failed release remains
preserved with its 180-second reservation.

The 22 diagnostic calls retain the original numerical caps. Root construction
and multiplier setup are counterbalanced in four group workers: two root groups
contain three policy cells each, and two ranking groups contain two policy
cells each. Each group reserves all constituent cells and their warm contexts
before its first sample; the group receipt carries every individual sample,
pre/post warmup, active wall/CPU and charged work. The group ledger charge is
the sum of its actual components. A, forced-multiplier and observer calls have
individual receipts. This gives 16 durable group/individual receipts for 22
planned cells. No sample order depends on an observed result.

For root groups, the production SIQS setup selects the fixed multiplier,
constructs one exact factor base and chooses the first eight A assignments
once. Each policy builds fresh production CoordinateFamily objects on that
identical shared base and assignment tuple, visiting all 64 specified Gray
steps in each family. Shared setup is timed/charged once and reported separately
from inclusive family build/traversal. A complete round compares final family
fingerprints across the three policies. Earlier direct-integer/CRT oracles
remain the independent transition proof; timed fingerprints are a passivity
check. Ranking cells make nine fresh production setup calls each, including
selector preview work. Root and ranking samples are required to finish all
nine rounds inside their original per-policy caps.

Single-cell diagnostics retain the full runner's exact factor reconstruction,
position counters, selected multiplier, bounded observer and RSS checks. The
70-digit A comparison runs its two policies in reverse order from 62 as the
protocol requests. Observer pairs compare semantic digests, charged work and
positions after their second receipt; a difference stops the stage without
silently deleting either observation. The old primary source remains available
at commit 1f95fc4; this release owns its changed diagnostic worker sources.

Primary reservations consumed 1,380 wall AND CPU seconds and 182,037,200,010
charged work units. Remaining diagnostic cells reserve 1,530 wall AND CPU
seconds, including their warm contexts. Analysis/orchestration and the two
unused warm slots leave 490 seconds within the original 3,400-second revision-2
calibration ceiling. New diagnostic receipts share the original 64 MiB study
output bound and the global 512 MiB ceiling. No late pilot or lower-band call
is admitted by this release.
