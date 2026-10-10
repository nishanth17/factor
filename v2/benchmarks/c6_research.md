# C6 / P4.5: bounded precomputed Lucas chains

The user activated C6 on 9 October 2026. The control is mainline
`bcf5f3d1e57304694b48ba6e7ef8b4ea2ffd0db0`; production remains unchanged.
The accepted A4 campaigns already reject checked PRAC promotion. This study
asks whether precomputation and compact execution change that decision.

## Sources and mathematical scope

| Source inspected | What it establishes | C6 use |
| --- | --- | --- |
| [GMP-ECM LucasChainGenerator, pinned 8ea5e214](https://github.com/sethtroisi/gmp-ecm/tree/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/LucasChainGenerator) | McLaughlin/Zimmermann generator produces prime-chain codes described as optimal or near-optimal. Search enumerates increasing lengths with divisibility and reachability pruning; ties favor more doublings, then shorter codes. | Reproduce the upstream search at B1=2,000; preserve code and decoder inputs, independently verify all results. C6 certifies scalar action, not global minimality. |
| [GMP-ECM ecm.c at the same pin](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/ecm.c) | A 64-bit code is decoded to offsets; stage one retains 16 projective points in a rolling buffer. The PRAC comments document a false-infinity example. | Extract decoder verbatim with assertions; compare the ring and last-use allocation on identical chains. Keep checks and recovery in the Python executor. |
| [CADO bytecode, pinned 692ecb7e](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/bytecode.c) | PRAC rule bytecode separates generation, weighted costs, compression and integer checking. Its cache keys include the cost-model pointer. | Inspiration for separating chain choice from dispatch. No CADO code is copied. Our four-byte register instructions deliberately use one executor for PRAC and Lucas. |
| [Bernstein, Cottaar, Lange, *Searching for differential addition chains*](https://doi.org/10.1007/s40993-024-00604-8) | *Research in Number Theory* 11, article 45, published 27 March 2025; ePrint 2024/1044 is the earlier preprint. Minimum length means the continued-fraction subclass. Improved pruning and the left-length meet-in-the-middle variant are exact there; the left-interval variant can miss minima. Runtime is distinct from weighted M/S/constant costs. | Recognize existing CF records and evaluate Algorithm 1 with three live points. New offline search is conditional on the frozen compact-executor gate. No global-optimum claim. |
| [Kruppa, 2010 thesis](https://docnum.univ-lorraine.fr/public/SCD_T_2010_0054_KRUPPA.pdf), §4.5.1 | Precomputing PRAC rule sequences amortizes repeated stage-one schedules; interpreter compression is an implementation decision. | Explains why both fresh construction and reuse need measurement. |
| [GMP-ECM 6.2 maintainer announcement](https://gmplib.org/list-archives/gmp-discuss/2008-May/003158.html) | Technical discussion records a Lucas-chain generation bug affecting P+1/ECM. | Historical reminder to validate actions, not merely count instructions. |
| [Kushagra Singh's FLINT implementation blog](https://iamkush.me/implementing-the-ecm-stage-i-in-flint/) | Search located the author's account of ECM stage one and PRAC. Full page retrieval failed in this environment. | Discovery only; no correctness or PyPy performance claim relies on it. |

The paper's compressed CF state is `(a,b,c)`, initially `(1,2,3)` with
`c=a+b`. Each bit keeps either `a` or `b`, retains `c`, and computes its
sum with the retained entry using the other entry as known difference.
This allows three persistent points. The dispatch diagnostic admits only
upstream records with the same additions and initial doubling; a matching
integer sequence containing extra doublings is excluded. General Lucas
chains need not have
that shape; the GMP-ECM ring supports older differences and extra doublings.
PRAC can also use subtraction, so its record verifier checks both permitted
sum/difference directions. ECM x-coordinate chains are not a drop-in
implementation of Williams p+1 Lucas-sequence arithmetic.

## License and reproducibility

The exact GMP-ECM generator, decoder source and notices are versioned in
[inputs/upstream/c6_gmp_ecm](inputs/upstream/c6_gmp_ecm/NOTICE.md), under
LGPL-3.0-or-later with COPYING and COPYING.LIB. C6 invokes separate executables;
production does not link them. The only generator build change reduces the
six-million-entry allocation capacity to 4,096; search logic is unchanged.
One generator thread is used. CADO's LGPL-2.1 COPYING and source were inspected,
but no code was adapted. The Python allocator/verifier is independently written.

No licenses are inferred from a repository's popularity. This repository has
no top-level license granting blanket rights over all its files; upstream
notices apply to the vendored inputs and their extracted decoder.

## Independently verified contracts

* Every upstream offset record is decoded by upstream C, then translated to
  SSA, checked against its claimed integer values, and verified by A4's integer
  interpreter. A separate mutable-register interpreter checks the compact
  record, including overwritten and uninitialized slots and final scalar.
* Tests compare all 303 prime records with full-coordinate affine arithmetic,
  exercise small nonsingular curves, prime powers, CRT composite moduli,
  prime squares, proper nonunits, saturation, and the published regression
  `n=33554520197234177`, `sigma=2046841451`, B1=373. Projective validity is
  checked separately from equality; `(0,0)` cannot pass vacuously.
* Compact execution checks X and Z before discarding a point. A4 retains
  earlier X coordinates for exceptional recovery; a compact ring cannot
  discard those factors. The extra X GCD is charged, including during the
  one permitted checked-ladder retry. Infinity/order-two differences trigger
  recovery before ambiguous differential addition.
* A saturated prime power replays only from that power's starting point,
  at most its exponent times. A proper divisor is returned immediately;
  unresolved inputs remain explicit. B3's stage jobs/checkpoints are untouched.

## Frozen scope and finite allowances

`inputs/controls/c6_protocol.json` pins source hashes, certified inputs,
selection rules, bounds and seeds before accepted timing. Required inputs
are versioned; raw captures, failed runs and process transcripts stay in
`results/c6/` (unpacked upstream scratch is in local `audit/results/c6/`).

The scope is B1=2,000 / B2=147,396, eight curves per attempt and twenty seconds
wall/CPU per attempt, on ten existing certified inputs: balanced and 10-digit
small-factor cases at each input size 40/50/60/70/80 digits. Balanced target
factors are 20/25/30/35/40 digits. This is a finite feasibility comparison,
not adequate coverage for balanced 80-digit factoring. Seeds 41001/48920
are fixed before timing; the next two corpus seeds are reserved for any
qualifying fresh confirmation. Each timing sample repeats the complete
fixed cohort, so repeated samples are not new independent factoring trials.

Prime generation is offline, one thread, B1=2,000, with sixty-second CPU/wall
limits, a 512 MiB process-RSS watchdog and 16 MiB per-output-file limit.
The decoder is an isolated process consuming only those generated codes.
It uses a 64-element upstream array and assertions; Python rejects records
above 64 elements. This is not a public decoder for arbitrary hostile codes.
The parent also imposes the subprocess timeout and validates every output.

Python records have 32-bit scalars, at most 512 instructions and 16 retained
register points. Execution additionally retains one original point for retry
and one result temporary. CF execution has three persistent working points,
plus original/recovery and arithmetic temporaries. A bound-owned program
has at most 512 unique records and 1 MiB of bytecode. Catalog reads are capped
at 1 MiB. Construction temporarily owns the <=303-entry Lucas SSA catalog;
no curve points are cached, and the A4 compiler cache is cleared after each
program build. No full-lcm scalar is searched; the ladder/oracle may evaluate
its ordinary exact scalar as before.

Separate PyPy 3.11 workers isolate int and mpz traces. Each warmed group uses
at least three seconds of validated warmup and nine repeated cohort samples;
relative IQR above 15% extends to 5 seconds/18 samples, then 8/27. Instability
remaining at the cap cannot earn promotion. Cold processes and upstream C
construction are reported separately. The shared flock and process inventory
exclude B4/A6 heavy checks and timing overlap.

The original protocol recorded a 5% compact-Lucas/compact-PRAC search margin,
a 10% distance from the ladder, and a 10% promotion threshold. These fields
remain historical evidence. During the stage-fresh capture, the user explicitly
removed the universal 10% floor and requested the revised roadmap policy,
read from the B4 worktree. The versioned
`inputs/controls/c6_policy_amendment.json` records that change without modifying
sources, arms, seeds or resource limits. Small sustained gains can qualify;
positive uncertainty and fresh held-out confirmation remain necessary.
The search decision is also reassessed for a credible remaining opportunity,
not rejected solely for missing an arbitrary percentage threshold.

## Reproduction and interpretation

Run from this worktree with the GMP-enabled PyPy 3.11 environment. The runners
acquire `/private/tmp/factor-performance.lock` and refuse competing benchmark
or test processes. Reserve the whole sequence with other study owners.
Use new output/scratch paths; preserve earlier captures instead of overwriting.

```sh
v2/.venv/bin/python -B -m v2.benchmarks.build_c6_inputs --scratch v2/benchmarks/results/c6/generation-new
v2/.venv/bin/python -B -u -m v2.benchmarks.c6_study --scope stage_reuse --cold --output v2/benchmarks/results/c6/reuse-new.json
v2/.venv/bin/python -B -u -m v2.benchmarks.c6_study --scope stage_fresh --output v2/benchmarks/results/c6/fresh-new.json
v2/.venv/bin/python -B -u -m v2.benchmarks.c6_study --scope campaign --output v2/benchmarks/results/c6/campaign-new.json
v2/.venv/bin/python -B -m v2.benchmarks.c6_costs --upstream-build v2/benchmarks/results/c6/generation-new --output v2/benchmarks/results/c6/costs-new.json
v2/.venv/bin/python -B -m v2.benchmarks.c6_report v2/benchmarks/results/c6/reuse-new.json v2/benchmarks/results/c6/fresh-new.json v2/benchmarks/results/c6/campaign-new.json --output v2/benchmarks/results/c6/summary-new.json
make -C v2 test PYTHON="$PWD/v2/.venv/bin/python"
make -C v2 lint
```

Generation rewrites the required catalog deterministically; check its hash
against the frozen protocol before running. For multiple cost runs, use a
fresh parent directory: each keeps nine `generator-repeat-*` subdirectories.
Do not invoke `--freeze` to bless a changed control during a comparison.

Stage samples sum ten one-curve attempts (five balanced input sizes and two
seeds). Campaign samples sum twenty complete attempts (both input shapes,
five sizes and two seeds), each allowed eight curves. These are actual
standalone ECM attempts, including setup, failed curves and stage two, not
complete recursive portfolio factorizations. A successful split still leaves
its cofactor explicit; failures retain the original composite. Production
primality labels, stage budgets and checkpoint formats are untouched.

Timed execution includes conversion, record construction/verification where
applicable, dispatch, X/Z checks and recovery. Independent certified-corpus
and affine-oracle validation is performed outside warmed execution timers;
fresh-process totals include it. The revised policy's output-validation cost
would need inclusion in any future promotion capture. No candidate receives
promotion on these diagnostic execution-only totals. Confidence intervals
from `c6_report` describe repeat-timing uncertainty on this fixed cohort,
not independent success-probability estimates or held-out confirmation.

## Decision boundaries for B3

The reusable output is a verified, bounded research candidate: immutable
records, a common executor, explicit nonunit/retry outcomes, and a reproducible
catalog. B3 would still have to provide work reservations, chunk/replay
semantics, program identity in checkpoints, bounded ownership across resumed
jobs, and complete portfolio validation. None of those gates is closed here.

Comparing Lucas with compact PRAC isolates chain choice under the same
interpreter and safety checks. Comparing the rolling and last-use layouts
isolates retained-point layout. Comparing either with checked A4 also changes
verification placement and factor-retention cost: A4 can inspect retained X
coordinates during recovery; the compact executor checks X before overwrite.
Removing those GCDs to win a timing would need a new factor-preservation proof.

A shorter integer chain is not a scalar-action proof or a runtime result.
Even an exact optimum in the continued-fraction family need not minimize
weighted arithmetic cost, and a weighted minimum need not minimize Python
execution time. Any future search extension must declare its family, pruning
argument, bounds and cost model before selection. This study does not search
or claim optimality for the full stage-one lcm.

## Accepted evidence and bounded stop

All 411 PyPy/GMP tests and full lint passed before measurement. The final
captures verify frozen source hashes. The [benchmark summary](README.md#c6--p45-bounded-lucas-study--9-october-2026)
records all five arms and both backends. Local evidence is
`results/c6/{stage-reuse,stage-fresh,campaign,costs,summary,operation-counts}.json`.
The measured implementation is `6d54a14`, based on `bcf5f3d`; later changes
only document policy/results and add exploratory per-class report output.

The 303 prime-power actions use 3,688 additions / 629 doublings for PRAC and
3,437 / 864 for Lucas. The illustrative 6A+5D cost is 25,273 versus 24,942,
a 1.31% reduction. This model counts squares like multiplications and the
curve-parameter multiplication in doubling; it omits dispatch, allocation,
GCDs and recovery. The ordinary lcm scalar has 2,878 bits. Integer counts
are diagnostics, not measured instruction costs or an optimality certificate.

Lucas reused stage time is 0.260573 s per ten native-int stages and 0.456377 s
per ten GMP stages. Compact PRAC takes 0.262837 / 0.483140 s: Lucas reduces
these point estimates by 0.86% / 5.54%, but the ladder takes only
0.091371 / 0.328431 s. Ring layout takes 0.263638 / 0.460233 s and supplies
no material advantage. A substantially shorter checked chain would still
need to overcome this observed gap; the near-optimal upstream records and
small weighted-cost change supply no evidence for such an opportunity here.

Exploratory campaign class breakdowns (ten attempts each) are:

| Backend / class | Ladder seconds | Lucas seconds | Lucas/ladder ratio, 95% interval | Splits ladder/Lucas |
| --- | ---: | ---: | --- | --- |
| int, balanced | 1.539843 | 3.013910 | 1.957 [1.921, 2.169] | 0/0 |
| int, 10-digit factor | 0.104120 | 0.324317 | 3.115 [2.947, 3.419] | 10/10 |
| GMP, balanced | 6.011574 | 7.259992 | 1.208 [1.201, 1.242] | 0/0 |
| GMP, 10-digit factor | 0.356613 | 0.495821 | 1.390 [1.306, 1.491] | 10/10 |

Every arm consumes 90 curves in this cohort. Lucas performs 80 stage-two
calls on balanced cases and five on small-factor cases; five small-factor
splits occur in stage one. No campaign prime-power replay or timeout occurs;
exceptional recovery is exercised independently by adversarial tests.
Per-class medians need not sum to the median of the complete cohort.
No class was selected after timing for promotion.

### Construction, retained storage and amortization

Nine identical upstream code files have SHA-256
`7b9e696ffe7d78e9d866fea4e511ce86083ff8abe735247ec5112b7f1aecb9ab`.
They hold 299 64-bit codes (2,392 bytes), plus four decoder small-prime cases.
The versioned decoded catalog is 59,130 bytes. Nine fresh C-process runs
measure generation at 310.428 ms median, decoding at 217.221 ms and parent
integer verification at 9.231 ms. The median of the nine total costs is
526.836 ms (range 473.554–772.383 ms); this is mixed cold-process construction,
not warmed PyPy speedup evidence. The separately timed first cold build
costs 1.081 s generator compilation plus 0.107 s decoder compilation.

Validated warmed program construction (>=3 seconds, nine samples unless
extended) is 21.528 ms checked PRAC, 24.381 ms compact PRAC, 8.780 ms Lucas,
and 9.974 ms ring Lucas (18 samples). Catalog load/decode/verify is 5.707 ms.
The owned 317-record PRAC program has 4,394 instructions / 17,576 compact
bytes / at most five register points. Lucas has 4,378 instructions /
17,512 bytes / at most four last-use registers, or sixteen rolling slots.
Recovery-only unit records explain the difference from stage-action counts.
Byte counts exclude Python object overhead; the combined cost-diagnostic
process peaks at 102,907,904 bytes RSS, not a measurement of cache size.

An accounting model for N stage executions sharing one generated catalog
and one loaded program adds `(generation + load/verify)/N` to reused stage
execution, plus one-time compilation/N when applicable. With the measured
medians, the Lucas setup term is about 535.6 ms/N (0.536 ms at N=1,000).
This arithmetic amortization example is not a fitted throughput prediction.
As N grows, setup approaches zero but the measured execution remains slower
than the ladder: **there is no observed break-even reuse count against it**.
Fresh-stage and complete-campaign captures charge actual construction in
place; generation and cold startup are reported separately rather than hidden.

### Three-point interpreter and search stop

Only 45 upstream chains match the CF family with identical arithmetic
(one initial doubling and the same subsequent differential additions).
For these 45 records on field 1,009, the compact/ring/three-point diagnostic
medians are 115.542/113.083/116.750 microseconds for int and
4.058/4.080/3.989 milliseconds for GMP. GMP compact/ring extend to 18 samples;
all final relative IQRs are below 11%. The GMP three-point point estimate is
1.69% below compact, but there is no independent confirmation or workload
stage evidence for that small effect. The field, parameter, point and
eligible subset are deliberately fixed; this cannot be extrapolated to
40–80-digit composites or to all 303 records.

**Decision: accept the verified research candidate, retain the production
ladder, and stop expansion into new offline CF search in this tranche.**
The original numerical search gate was not met. Reassessment under the user's
no-universal-floor policy reaches the same conclusion: all full-stage and
complete-campaign candidates lose, and the same-executor comparison gives
no credible remaining chain-quality opportunity. This is a bounded empirical
stop, not a proof that every possible differential chain/interpreter is slow.
Reserved confirmation seeds remain unused because no candidate qualifies;
there was no held-out tuning. New search and B3 integration gates remain open.

For B3, retain these records/verifiers as an opt-in correctness reference and
reuse the test corpus. Do not route compact PRAC, Lucas or the CF diagnostic
into production defaults. Reopen performance work only with a concrete cost
hypothesis and fresh frozen full-stage/campaign evidence, including output
validation and the integration obligations described above. Do not combine
these numbers with B4 kernels or A6 changes: the arithmetic control stays pinned.

## Committed-only acceptance receipt

An archive of commit `860bab4` (with only the external PyPy development
runtime linked in) passes `make -C v2 test` with 411 tests, `make -C v2 lint`
with all 148 Python files formatted, and imports all 72 benchmark modules.
Frozen C6 hashes, the 303-record catalog and the independently certified
campaign corpus load without any raw evidence or ignored audit inputs.
Logs and the validation receipt remain in `results/c6/committed-*`.
The shared timing/check window was released to B4 after these checks.
Subsequent changes only add this documentation receipt; no measured source
or required input changed. Mainline stays at `bcf5f3d`; no merge was performed.
