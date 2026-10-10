# B3 production chain integration: source and proof boundary

Base: committed mainline `76f06b0110941017b926caa51d341db6221f470b`.
The [C6 optimized study](c6_optimization.md) already supplies independently
verified PRAC/Lucas records, coordinate-factor certificates and strict recovery.
Its 6.73% reused native saving and 5.68% fresh regression are historical
research observations, not production break-even estimates. No search,
register allocation, reducer, kernel variant, curve family or pairing is added.

## Primary research and licensing

- [Bernstein–Cottaar–Lange](https://antsmath.org/ANTSXVI/papers/BernsteinCottaarLange.pdf),
  sections 1.1 and 2: differential scalar identities are separate from a
  weighted runtime advantage; the minimum-length claims concern the CF family.
  Re-read for the composition boundary. No search is repeated.
- [GMP-ECM ecm.c, 8ea5e214](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/ecm.c):
  inspected the versioned C6 source, including its false-infinity example
  `33554520197234177`, sigma `2046841451`, B1=373, and bounded Lucas decoder.
  Remote pinned-source retrieval returned cache misses during B3; the exact
  local source and its LGPL-3.0-or-later header remain authoritative evidence.
  C6's required catalog is reused; no new C source is adapted or linked.
- [CADO bytecode, 692ecb7e](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/bytecode.c):
  reuse C6's recorded separation of verification, encoding and execution.
  Its inspected LGPL-2.1 implementation is a design reference. B3 does not
  copy its source; remote retrieval of this pin also failed.
- [GMP-ECM maintainer discussion](https://gmplib.org/list-archives/gmp-discuss/2008-May/003158.html)
  records a historical Lucas generation defect. [Clift's author notes](https://additionchains.com/Lucas.html)
  describe a corrected overflow/pruning defect. These reinforce verification
  discipline; neither supplies an arithmetic guarantee or performance claim.

Upstream notices/licenses remain in
[the C6 input directory](inputs/upstream/c6_gmp_ecm/NOTICE.md). Production
reuses this repository's independently written C6 verifiers/interpreters and
B4 residue-equivalent kernels. There is no new upstream source adaptation.
The record verifier, strict interpreter and checked ladder preserve the C6
function bodies; frontier verification changes only local module references.
The native kernels preserve C6/B4's bodies. The tuple execution loop preserves
C6's instructions and masks, binding the selected kernels once per action.

## Integration proof obligation

No new elliptic-curve formula or optimality theorem is required. The new
obligation is **compositional coordinate-factor coverage and atomic commit**:

1. The independent integer interpreter proves each record's exact scalar and
   each known sum/difference. A production chunk traverses the same inclusive
   prime-power schedule in increasing order. Composition therefore multiplies
   by the product of those powers, including repeated-prime multiplicities.
2. C6's forward-bitset certificate proves that each intermediate X/Z factor
   divides a guarded coordinate or an output coordinate. A record's output
   becomes the next record's input. Induction across at most 16 records makes
   the accumulated guards plus final X/Z cover every discarded coordinate,
   including record-boundary points. Polynomial divisibility remains true
   modulo composite and prime-square moduli; no field-only cancellation is
   used. Reduced kernels compute the same residue polynomials.
3. GCD one on that aggregate proves all covered coordinates are units modulo
   n. In particular a zero or ambiguous known-difference coordinate cannot
   pass certification. Nonunit or saturated aggregates publish no fast point:
   replay begins at the original chunk point with strict X/Z checks before
   overwrite. The unchanged strict interpreter has one finite checked ladder
   retry; a failed prime power has at most its exact exponent's unit retries.
   Every returned divisor then passes the existing production proper-divisor
   validator. Saturation is never a split or a prime classification.
4. If the strict prefix cannot certify the entire chunk, retain the original
   job point and powers and enter the existing durable prime-unit replay.
   Its actions keep their normal per-unit reservations. Fast execution and
   its immediate worst-case strict/unit replay are reserved together before
   arithmetic, even when the recovery allowance goes unused.
5. Only certified chunk points/cursors enter checkpoints. Schema 10 records
   catalog, bound, backend/version, selected kernel, interpreter/recovery and
   batch identity. Resume independently checks the consumed prime prefix,
   pending powers, chunk counters and canonical unit-Z point, then charges
   those checks. It reconstructs an empty run-owned cache and pays every miss.
   Reapplying a whole M(B1) to a completed stage is prohibited by the cursor;
   increased-B1 ECM continuation remains outside B3.

## Finite ownership and cost model

Plans retain immutable code/masks and decoded tuples, never curve points.
Each miss reads/hashes/parses exactly 202,461 catalog bytes, then charges four
verification/decode work units per instruction plus register/record costs and
B1 units for schedule construction. Failed preparation retains consumed work
but publishes no partial plan. A hit costs one lookup unit; execution reserves
point/guard/GCD units plus the strict and repeated-prime worst case. These new
work units differ from scalar-bit ladder charges; comparisons use equal total
allowances and disclose completion and actual elapsed/CPU costs.

The chain cap reserves 4 MiB construction/recovery scratch plus at most 4 MiB
per retained plan; minimum 8 MiB. Conservative plan ownership is
`4096 + sum(4096 + 512 * instruction_count)` over 333 records. Scratch covers
raw/parsed catalog coexistence, verification bitsets, decoding, one saved
point, at most 16 working points, guard products and strict arithmetic
for the existing 4,096-bit input cap. Before a miss, LRU eviction reserves a
complete replacement alongside all remaining retained plans. Refused or
cancelled misses cannot publish incomplete records. Process/JIT RSS is a
separate measurement. Catalog hashes and metadata have bounded sizes.

The independent schedule store has its own existing cap. Both caps are added
to PortfolioConfig's simultaneous workspace reserve, including checkpoint
output and context storage. A recursive portfolio shares one store across
curves and cofactors; each resumed invocation starts a fresh store. Setup,
verification, evictions, output checks and unsuccessful curves are timed.

The [frozen protocol](inputs/controls/b3_protocol_v3.json) fixes two certified
12-input cohorts, repeated seeds, controls, scopes, bounds and equal resource
caps before timing. The routing trial is explicit `reuse`, B1=2,000,
chunk16, tier at least eight curves. Fresh/smaller/unsupported tiers retain
B4. Default `off` retains all old schemas and arithmetic. Low-level eviction
stress alternates 2000/1999/2000 and is outside enabled production routing.
The acceptance decision and measurements are published only after the gate.

Pre-timing preflight correction: explicitly release the evicted local plan
reference before replacement allocation. The original `b3_protocol.json`
remains byte-identical and had no timing captures. `b3_protocol_v2.json`
updates source hashes only, retaining all inputs, parameters, selection and
uncertainty rules; no result informed this correction.

The final pre-timing v3 freeze also replaces two equivalent complex slice
endpoints with named indices to satisfy both Ruff and pycodestyle. v1/v2
remain unchanged; no timing capture or candidate selection preceded v3.
