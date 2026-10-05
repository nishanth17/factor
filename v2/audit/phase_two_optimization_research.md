# Phase 2 optimization follow-up and future portfolio work

Date: 3 October 2026. Milestone: M23. Scope: source, literature and
implementer-blog review; documentation only. The
[source manifest](m23_phase_two_optimization_sources.json) records retrieved
bytes, hashes, immutable implementation revisions and inspection scope.
The verification records
content/link checks and preservation of existing code and evidence.

The useful next work is calibration of the existing bounded portfolio,
proof-backed preprocessing shortcuts, and measured reductions in setup and
recovery overhead. Phase 2's core is accepted. Brent batching, streamed
schedules, chunk recovery, checkpoints and M13's local arithmetic loops
already exist. They are controls for new experiments, not missing features.
Native implementation timings do not establish which configuration wins on
PyPy Python 3.11.

At the user's request, remaining portfolio optimization is consolidated in
[future Phase 8](TODOS.md#phase-8--revisit-phase-2-portfolio-optimization).
Historical P2 IDs and open experiment gates remain intact. The separately
completed M22/P3.2 collector leaves P3.3 filtering/extraction next; Phase 8
does not delay the committed SIQS/GNFS
workstream. Individual experiments may start earlier when profiling supplies
a reason, but portfolio integration needs working SIQS and the bounded small
GNFS baseline.

## Reconciliation with completed work

The [M17 frozen baseline](phase_two_m17_frozen_baseline.json),
fresh confirmation,
small-band controls, and
[public changelog](../../CHANGELOG.md) are the authoritative local evidence. M17 accepted
the core with 97 tests; the separate M20 QS acceptance later recorded 120
tests. M22's collector acceptance records
134 tests and passes P3.2, with measured small-fixture regressions rather
than performance promotion. None closes an unmeasured Phase 2 tuning gate.

| Original item | Already implemented or decided | Remaining work and owner |
| --- | --- | --- |
| P2.1 | Exact squares/higher powers, multiplicities, classification reuse, bit-based removal of twos, bounded optional Fermat | Held-out trial cutoff and preprocessing-cost selection; P8.2 |
| P2.2 | Seeded serial dispatcher, shared work/wall/CPU allowances, cancellation and reconstructible exhaustion; large-cofactor rho; extreme legacy ECM jump removed | Rho batches/restarts and joint ECM B1/B2/curve allocation; P8.3–P8.4. SIQS integration remains P3.4; GNFS crossover remains P7.7 |
| P2.3 | Packed base primes, private bounded marking buffers, half-open prime streams, on-demand stage 2 and optional capped cache | Cold/warm cache and setup amortization if the workload warrants it; P8.6 |
| P2.4 | Stage-1 chunk actions, fine saturated-chunk replay, stage-2 term replay, versioned JSON checkpoints and corrupt-state rejection | GCD/checkpoint/recovery granularity and serialization costs; P8.5 |
| P2.5 | Frozen 882-input corpus, independent certificates hidden from algorithms, seeded isolated runner with both output modes and resource accounting | Fresh tuning/held-out design, broad repeated evaluation and actual competitor execution; P8.1. Publication remains P6.4 |
| P2.6 | Reusable contexts and rolling-offset arm validate; M12 rolling factoring regressed, so rolling stays off | Revisit only after a changed bound/backend/profile supplies a reason; P8.6 |
| P2.7 | Wheel-6/bytearray control and validated wheel-30, pre-sieve, bitset and packed-output arms; alternatives lacked broad factoring evidence | Conditional consumption/end-to-end evaluation, not another automatic wheel rewrite; P8.6 |
| P2.8 | M14 feasibility passed: 48 configurations, 432 samples, cold/reused workers and first-factor cancellation; serial retained | Broader process promotion remains P6.3; integrate any accepted winner in P8.7 |

M13's loop improvement is already banked: under M17's declared 50 ms
evaluation caps, balanced 20-digit complete-factorization completion rose
from 70.8% for the M12 control to 91.1% for the bounded candidate. This is
historical matched evidence, not an M23 measurement or a new parameter win.
The larger-band screen is one
repetition under five seeds; its limited completions do not establish a
universal digit limit or justify tuning on those outcomes as held-out data.

Retain trial bound 25,000, rho batch 64, chunk size 16, stage-2 GCD batch 128,
wheel-6/bytearray, cache/rolling off and serial execution. ECM uses
B1/B2 = 2,000/147,396; the production allowance is 32 curves, whereas M17
evaluated two. The distinct p−1 bounds are 2,000/200,000. Fermat remains off
by default. M17's input cap was 512 bits versus the library's 4,096-bit cap;
these are configurations, not measured coverage promises. The bounded API
remains opt-in.

## Implementation references and transfer limits

| Reference inspected | Immutable revision | Useful comparison |
| --- | --- | --- |
| [GMP-ECM mirror][gmp-readme] | `8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e` (2026-03-09) | Factor-size probability models, joint stage bounds, incremental prime powers, stage-2 blocking and saved progress |
| [YAFU autofactor][yafu-auto] | `963dbe9c45283cc06e9a71d830b0676bc1b0d343` (2026-09-24) | Explicit pretesting depth and locally tuned QS/NFS time estimates |
| [FLINT factor driver][flint-factor] | `00cc19b350b5302876b4fa8c877df5a676c2061c` (2026-10-02) | Exact power handling and bounded small-factor work before QS; word-sized rho specialization |
| [SymPy factor helpers][sympy-factor] | `06fdc356176f6ae82c865a365701f64c5ee67aa9` (2026-10-03) | Perfect-power search reduced by known trial progress; staged factorization and explicit method controls |
| [primefac][primefac] | `2401341d8689ea2ab517ef311f4f91ecd33d14d0` (2023-06-23) | Python Brent product batches/replay and differing factoring/primality contracts |

GMP-ECM's pin is a public GitHub mirror, not a certification of the latest
official distribution. Reuse the existing
[competitor pins](../benchmarks/inputs/phase_two_competitors.json); they remain
unexecuted. This is a design comparison across mature implementations, not a
current fastest-implementation ranking. Imported code would also require
license review; no upstream code was executed or copied into Factor.

## Findings worth folding into the plan

### Preprocessing with fewer exact root attempts

[Bernstein's 1998 paper][bernstein] establishes perfect-power detection and
describes practical residue screening in its final section. Testing prime
exponents is already sufficient and already implemented here. The new
candidate is to avoid unnecessary exponents and root work, not to add power
detection again.

SymPy's `_perfect_power(n, next_p)` uses a proven lower bound on possible
prime factors to constrain the exponent search. If completed trial division
proves every prime factor of a remaining cofactor is at least L, a
nontrivial power n = b**k has b >= L and therefore L**k <= n. Compute that
bound with exact integer comparisons. Preserve the proof of completed trial
progress across splits/resume; partial trial division cannot justify the
same bound. SymPy also uses floating estimates in its helper; those are not
suitable for Factor's exact-arithmetic contract. [Implementation][sympy-factor]

For a prime k, a small prime q with k dividing q−1 can reject a kth power
when `pow(n % q, (q - 1) // k, q)` is neither 0 nor 1. Passing this filter
is not proof: retain exact root/exponent verification. Charge table creation,
modular powering and cache storage. Test adversarial inputs that pass all
filters as well as powers and neighboring values. [Paper][bernstein]

A [GMP developer's 2008 account][gmp-power-notes] describes p-adic roots and
bit-length screening before a costly final exponentiation. It is useful
implementation evidence, but its native large-input timings do not transfer
to our domain. Keep a full root-algorithm replacement conditional on a root
bottleneck; cheaper exact exponent pruning/residue tests come first.

The optional Fermat path currently recomputes `isqrt(n)` and `a*a-n` per
candidate. An engineering experiment can initialize once and update
`D(a + 1) = D(a) + 2*a + 1`, with exact square-residue rejection before
`isqrt(D)`. Keep its finite allowance and default off unless declared
close-factor workloads support a full-run benefit. These are code-derived
hypotheses; no Fermat gain was measured here. Owner: P8.2.

### Rho tuning while preserving bounded recovery

[Brent's 1980 paper][brent] supplies the cycle-finding/batched-factorization
foundation. Factor's `_rho_step` already doubles cycle lengths, accumulates
products, saves the batch start and performs bounded individual recovery.
M13 already moves arithmetic into locals between committed boundaries.

The useful remaining comparison is GCD batch size jointly with walk length,
restart allowance and actual saturation cost. Reuse P2.2's 32/64/128/256
grid as candidates, holding assigned seeds, total resources and accounting
fixed. Larger batches save GCDs but delay detection and may increase replay
and cancellation latency. Primefac confirms the familiar batching design,
but its random batch size and unbounded search/replay are unsuitable for our
API. [Python source][primefac]

[Algorithmica's implementer chapter][algorithmica] discusses batched GCDs
and fixed-width Montgomery arithmetic; it also reports false negatives,
overflow concerns and the need for replay when a product saturates. Its
word-sized C++ task and error-tolerant benchmark differ from arbitrary-size,
validated complete factorization. Its batch 1,024, reciprocal divisibility
tricks, native throughput and parallel scaling estimates are not defaults
to import. FLINT's word-sized rho similarly specializes reduction and
batching; its constant 100 is another candidate, not an optimum for PyPy.
[Native source][flint-rho]

Start with policy calibration against the current Brent control. Reducer
implementation remains P4.4 and process search remains P6.3. A product-tree
recovery variant is eligible only if saturation replay dominates and its
nodes/work fit caps; rho would also need retained/reconstructed differences.
Owner: P8.3, with shared recovery evaluation in P8.5.

### Joint stage bounds and earlier handoff to the relation engines

GMP-ECM chooses B1/B2 and curve counts by target factor size and expected
time, with different costs for its continuation variants. Its smoothness
model accounts for curve parametrization; random integers alone are an
incomplete group-order model. Use this to structure a PyPy training grid,
then fit measured stage costs and completion. Do not transplant native
tables or claim to know the hidden smallest factor. [README][gmp-readme],
[probability implementation][gmp-prob]

YAFU maintains pretesting plans and completed ECM depth, then estimates
QS/NFS cost from tune data. FLINT's current driver also spends limited
rho/ECM work on small factors before QS. The transferable idea is a
measured handoff policy. Their digit/bit thresholds and native timing fits
are not Factor crossovers. [YAFU policy][yafu-auto],
[plan documentation][yafu-doc], [FLINT driver][flint-factor]

Compare bounded rho, p−1 and ECM allocations against the cost of reaching
SIQS/GNFS, using only observable input/configuration and recorded progress.
Stratify evaluation by known factor sizes without exposing factors to the
algorithm. Compare marginal completion per CPU-second on the remaining
cofactors, with full recursion/setup costs. Repeated failures condition the
remaining input population; do not assume its observed success rate stays
constant. Policy identity and prior expenditure must survive checkpoints.

Repeating the same p−1 base and bounds is duplicate work; unlike fresh ECM
curves it does not sample a new group order. Enlarging B1 can instead apply
the missing prime-power exponent ratios: old primes may need higher powers,
so skipping every old prime is wrong. This is visible in GMP-ECM's
`B1done` handling. Same-curve extension should compete with fresh curves,
not receive assumed benefit. [README][gmp-readme], [stage 1][gmp-pm1]

P5.3 owns incremental powering/continuation implementations; P8.4 owns their
budget calibration. P3.4/P7.7 own relation-engine dispatch and crossover.
Neither a native heuristic nor a feasibility screen closes those gates.

### Separate four kinds of granularity

GMP-ECM's stage-2 implementation checks cancellation around substantial
blocks and exposes memory-dependent continuation choices. Its historical
[ECM survey][ecm-survey] describes blocking to trade setup/throughput against
storage. Factor already has streamed continuation and chunk replay; adopting
a polynomial continuation still belongs to P6.2. [Current source][gmp-stage2]

Measure these controls separately before changing them:

1. GCD polling and saturated-product recovery.
2. Cooperative budget/cancellation checks and atomic state commits.
3. Optional external durable checkpoint writes by the caller.
4. Final JSON serialization and resume verification.

The present library maintains in-memory progress and emits a checkpoint at
return; it does not write a file every arithmetic iteration. A checkpoint
frequency claim must specify which control changed. Larger atomic actions
can reduce overhead while increasing cancellation/deadline overshoot;
native big-integer calls remain indivisible. Keep strict work reservations,
finite replay, consumed CPU/wall allowances, canonical state validation and
reconstruction. Measure pack/unpack and cold resume without removing checks.
Version any changed state or accounting semantics. Owner: P8.5.

### Setup, schedules and classification are conditional opportunities

In `factorize_bounded`, a context sized for the largest configured bound is
constructed before cofactor classification. Profile the wasted cold setup
on primes, powers and early rho successes. A lazy context or staged growth
is a concrete candidate, but must retain upfront configuration/memory checks,
one-time work charges, exact streams and resume identity. Sharing bounded
integer prime/power/gap schedules is safe only independently of modulus;
point tables and powered residues remain curve/modulus-specific. Owner:
P8.6; stronger powering/table algorithms remain P5.2–P5.3.

Revisit existing cache/rolling/wheel arms only when the larger B2, backend or
reuse pattern changes the measured bottleneck. Include prime consumption,
decoding, setup, cache misses, memory and full factoring; reuse the
[C sieve transfer review](sieve_port_review.md). M12's losses remain valid
evidence for their tested workload. A warm schedule utility win does not
overturn them.

SymPy's [primality implementation][sympy-prime] combines deterministic
ranges with probable-prime BPSW for larger values; Factor currently honors
its documented Miller–Rabin rounds and certainty. If classification cost
dominates, a Lucas/BPSW-based extra composite rejection may be an experiment,
with independent pseudoprime fixtures and charged state. Do not silently
replace the configured rounds, label probable primes proven, or reuse
classification across a changed cofactor. Owner: P8.2; a certainty-contract
change would require its own API decision.

### Curve improvements and other algorithms keep their existing owners

[Barbulescu, Bos, Bouvier, Kleinjung and Montgomery][ecm-galois] prove
divisibility properties explaining ECM-friendly families. The 2012 preprint
appeared in the 2013 ANTS-X proceedings. This supports P6.1's whole-engine
family experiment, including setup/conversion and success, rather than a
new duplicate Phase 8 curve implementation.

The [Zimmermann/Dodson survey erratum][ecm-errata] reports that testing many
curves against one prime gave misleading average valuations; differing
prime congruence classes changed the result. Use many independently selected
prime factors and seeds, with residue-class controls where relevant. A
favorable single input or lower operation count is insufficient evidence
for curve/default promotion. Its 2023 indexing correction also matters if
P6.2 later ports the survey's d1/d2 continuation; do not copy the uncorrected
formula.

[Hart's 2017 implementation blog][hart] explains the role of small-factor
methods versus QS in a portfolio. Its historical coverage ranges are not
current PyPy dispatch thresholds. FLINT's bounded SQUFOF branch is interesting
for word-sized cofactors, but adding another engine is lower priority than
calibrating the present rho path and completing SIQS/GNFS. No new SQUFOF
implementation gate is added without a demonstrated coverage/cost gap.

PRAC, fused kernels, GMP, custom reducers, p+1, paired stage 2, curve families,
polynomial continuations and worker execution remain P4.1–P4.4,
P5.1–P5.3 and P6.1–P6.3. P8.7 integrates only individually accepted winners
and records retain/defer/reject decisions; it does not require every optional
method to be implemented.

## Future experiment sequence and acceptance

1. **P8.1: establish the new control and untouched evaluation data.** Preserve
   M17 captures. Freeze the current working SIQS/small-GNFS configuration,
   profile separately, and use fresh training/held-out inputs. M15/M17 inputs
   have now informed this review and cannot be recycled as new tuning data
   followed by an independent-validation claim. Execute pinned competitors
   only as explicit, matched arms; retain infeasible/unavailable outcomes.
2. **P8.2–P8.4: test one policy or preprocessing change at a time.** Prefer
   proof-backed root pruning and rho/ECM allocation. Joint tuning follows
   individual cost/success evidence. Do not use hidden factors in dispatch.
3. **P8.5–P8.6: reduce demonstrated overhead.** Include recovery, cancellation,
   serialization, context allocation and output consumption. More elaborate
   roots, product trees, caches and wheels need a bottleneck first.
4. **P8.7: integrate and confirm on untouched inputs.** Reuse the TODO
   promotion policy: zero correctness failures, identical declared resource
   limits, uncertainty-aware full-run completion/time benefit, and bounded
   regressions in other classes. Retaining the control is a valid decision.

Run eventual experiments only on supported PyPy Python 3.11. Warm relevant
paths for several seconds, verify stability and record JIT/runtime versions;
keep cold startup and warm workers distinct. Use repeated inputs/seeds,
medians/spread, censored exhaustion, total CPU and aggregate RSS. Check both
factor-one and complete-factorization outputs and exact reconstruction.
Owned-workspace caps, consumer output and process/JIT RSS are distinct.
When an accounting unit changes, version/disclose it and compare actual
wall/CPU limits rather than pretending the numerical work caps are equivalent.

M23 changes no runtime defaults, closes no additional implementation or
performance gate, and runs no competitor, tests or factoring benchmarks.
The workspace's QS collector/build files changed during this review and
M22 accepted P3.2 separately; its evidence and observed hashes are recorded
separately from this documentation-only task.
The thirteen M17 Phase 2 production files still match the frozen control.
The ePrint byte fetch failed; the published ECM-family PDF was retrieved and
inspected instead. Task-owned downloads/extracted text are removed after
hash verification; only provenance metadata is retained.

[gmp-readme]: https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/README
[gmp-prob]: https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/rho.c
[gmp-pm1]: https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/pm1.c
[gmp-stage2]: https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/stage2.c
[yafu-auto]: https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/autofactor.c
[yafu-doc]: https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/docfile.txt
[flint-factor]: https://github.com/flintlib/flint/blob/00cc19b350b5302876b4fa8c877df5a676c2061c/src/fmpz_factor/factor_no_trial.c
[flint-rho]: https://github.com/flintlib/flint/blob/00cc19b350b5302876b4fa8c877df5a676c2061c/src/ulong_extras/factor_pollard_brent.c
[sympy-factor]: https://github.com/sympy/sympy/blob/06fdc356176f6ae82c865a365701f64c5ee67aa9/sympy/ntheory/factor_.py
[sympy-prime]: https://github.com/sympy/sympy/blob/06fdc356176f6ae82c865a365701f64c5ee67aa9/sympy/ntheory/primetest.py
[primefac]: https://github.com/lucasaugustus/primefac/blob/2401341d8689ea2ab517ef311f4f91ecd33d14d0/primefac.py
[brent]: https://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf
[bernstein]: https://cr.yp.to/papers/powers-ams.pdf
[ecm-survey]: https://members.loria.fr/PZimmermann/papers/ecm.pdf
[ecm-errata]: https://members.loria.fr/PZimmermann/papers/#ecm
[ecm-galois]: https://msp.org/obs/2013/1-1/obs-v1-n1-p04-s.pdf
[algorithmica]: https://en.algorithmica.org/hpc/algorithms/factorization/
[gmp-power-notes]: https://gmplib.org/list-archives/gmp-devel/2008-May/000798.html
[hart]: https://wbhart.blogspot.com/2017/02/integer-factorisation-in-flint.html
