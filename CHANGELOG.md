# Changelog

## Current development

- Preserve the original Python 2 implementation in v1; develop and test v2
  on PyPy implementing Python 3.11.
- Repair arithmetic boundaries, primality certainty, factor reconstruction,
  segmented sieve endpoints and saturated-batch recovery.
- Add bounded preprocessing and seeded rho/p−1/ECM portfolio execution with
  serializable checkpoints, shared allowances and explicit unfinished results.
- Add exact QS/MPQS polynomial and relation identities, bounded single-large-
  prime collection, provenance-aware filtering, GF(2) dependencies and modular
  factor extraction with in-memory resume.
- Remove redundant collector budget polling and proved-useless score marking.
  Keep experimental resieving and alternative scoring/filtering controls
  optional; see [measurement scope](v2/benchmarks/README.md).
- Keep code, tests, roadmap, research and required fixtures on GitHub.
  Exclude generated run captures, logs, profiles and local development records.

Next: SIQS families, incremental roots, serialized checkpoints and bounded
dispatcher integration under [P3.4](v2/audit/TODOS.md).

PRAC optimization remains [P4.1](v2/audit/TODOS.md); production currently uses
the safe Montgomery ladder.
