# C10a-01 implementation and experimental handoff

The exact parity-normalized SIQS representation and its coupled multiplier
score are implemented **opt-in**. Standard arithmetic and defaults remain the
control. Independent arithmetic, roots, extraction and resume checks pass.
Short calibration is complete; complete-factor pilots at 60, 62, 65 and 70
are designed but **not run**. No speedup or default adoption is established.

The [arithmetic contract](c10a01_arithmetic.md) distinguishes existing division
by A from the new parity population. For odd D=kN≡1 (mod 4), odd B and odd a,
store A=2a and F=2P with P=a*x²+B*x+(B²−D)/(4a). Then
U=2ax+B and U²−D=4aP=AF. Even D and D≡3 (mod 4) use standard fallback.
The new score includes finite even squarefree multipliers with explicit
standard-form semantics. Direct P storage, native multiplier tables, paired
kernels, hypercube batching and new storage/solver formats are deferred.

The retained and current source review follows official msieve, YAFU, FLINT,
Yamaquasi, PARI/GP, CADO's actual SIQS build, Alpertron, SIQS.NET and JavaMath.
Primary literature includes Alford–Pomerance, Contini, Boender–te Riele,
Carrier–Wagstaff and Bradford–Monagan–Percival, with Percival's implementation
material and monograph support. Detailed operative callers, exact revisions,
hashes, mechanisms, alternatives, inaccessible material and limitations are
retained locally in `v2/audit/c10a01/research.md`. This is not exhaustive SOTA
coverage or a native-reference performance comparison.

## Short observations

Every scored prefix visited 16777344 positions: 128 polynomials, two A
families, 4096-position blocks and the one-position interval tails. Two exposed
independent inputs per size are not nine samples or fresh confirmation.
All calls completed fixed work; none completed factoring. Columns below list
balanced / approximately 35:65 results. Wall timings are indicative only.

| Digits | Control seconds | Challenger seconds | Retained rows, control | Retained rows, challenger | Structural support / complete factoring |
| --- | --- | --- | --- | --- | --- |
| 60 | 4.840 / 4.808 | 4.904 / 4.930 | 24 / 32 | 28 / 46 | Reachable and exercised / not run |
| 62 | 5.293 / 5.336 | 5.252 / 5.436 | 30 / 28 | 16 / 30 | Reachable and exercised / not run |
| 65 | 6.520 / 6.231 | 6.664 / 6.240 | 28 / 12 | 28 / 19 | Reachable and exercised / not run |
| 70 | 7.282 / 7.081 | 7.211 / 7.128 | 4 / 9 | 4 / 7 | Reachable and exercised / not run |

The scored balanced 70-digit cell chooses k=1 in both arms, so its standard
fallback has identical arithmetic/counters. It cannot test normalization.
The separately frozen forced-policy diagnostic selects k=31, the highest
ranked eligible candidate, before any complete-factor outcome. It completed
512 polynomials / 67109376 positions in 28.576 seconds, retaining 30 rows.
Stage C explicitly uses this forced coupled candidate; it is not a claim
about the automatic scored selector.

A fixed-k=5 residue-5 comparison at 60 completed 768 polynomials per arm:
27.297 / 28.681 seconds and 151 / 151 retained rows. This one observation
provides no speedup evidence. An even-k=2 fallback at 62 completed 640
polynomials in 25.247 seconds. All four forced calls exercised construction,
roots, scores, exact refinement and retained relations. Four observer on/off
pairs preserved exact state and work; ratios 0.972–0.984 are a short overhead
screen, not evidence that recording speeds execution. Each arm/configuration
received at least three validated wall AND CPU seconds of process-local warmup.

A-product envelopes were reachable, but selected products were not uniformly
near target: the balanced 62 challenger visited a family at 0.266 times its
target; forced 70 visited one at 0.284 times target, with the other seven near
1. Norm samples did not establish a general reduction. These population
uncertainties remain in the handoff; no retired assignment/face/quota grid
was repeated or parameter retuning performed.

## Calibration consumption and retained failure

Measured completed functions used **291.007 wall / 290.769 CPU seconds**.
This excludes unmeasured orchestration and one failed report-aggregation call;
those retain conservative charges. Total charged allocation is **3124 wall
AND CPU seconds**, leaving 476 of the 3600-second ceiling unused. Charged
work is **138632524390**, including the failed report's full 10 billion units.
No factoring observation failed, exhausted a cap or was retried.

| Phase | Calls | Charged wall AND CPU seconds | Measured wall / CPU seconds |
| --- | ---: | ---: | --- |
| generation | 1 | 600 | 0.097 / 0.097 |
| warmup | 18 | 630 | 63.913 / 63.894 |
| prefix | 16 | 980 | 95.162 / 95.108 |
| analysis | 31 | 114 | 1.125 / 1.042 |
| observer | 8 | 160 | 13.614 / 13.609 |
| forced | 4 | 200 | 109.807 / 109.733 |
| memory | 2 | 100 | 7.289 / 7.288 |
| orchestration | 1 | 340 | conservative allowance |

The aggregation error counted its own pending reservation as an incomplete
prior call. Its failed ledger entry remains, and corrected analysis uses a new
identity. This is distinct from native memory/work exhaustion. Recorder output
stayed bounded; numerical observation is now closed. Research, implementation
and ordinary correctness QA were separate from this experimental allocation.

## Resources and later stages

Literal configurations hold interval, block, LP domains, plans, queues,
backend, pacing and finite storage policies fixed within each pair. Memory
never silently selects another factory default. The outside assignment and
bounded harness reserves are suballocated from the total owned tier.

| Digits, both arms | Base / A factors | Rows / atoms / partial caps | Dense nonzero capacity | Owned / RSS tier | Central conditional owned scenario |
| --- | --- | --- | ---: | --- | ---: |
| 60 | 5000 / 8 | 12000 / 65536 / 32768 | 60012000 | 1 / 3 GiB | 0.733 GiB |
| 62 | 6000 / 9 | 16000 / 98304 / 49152 | 96016000 | 2 / 5 GiB | 0.888 GiB |
| 65 | 8500 / 9 | 24000 / 131072 / 65536 | 204024000 | 2 / 5 GiB | 1.308 GiB |
| 70 | 10500 / 10 | 32000 / 196608 / 98304 | 336032000 | 4 / 8 GiB | 1.688 GiB |

The central scenarios assume about 5F atoms, 4F forest edges, 1.2F rows,
20 exponents/atom and two source atoms/combined row; they are not RAM theorems
or completion forecasts. Exact conditional ownership formulas, shared-store
lifetimes, dense-fill reservations and upper scenarios are in the local
resource audit. Independent caps need not fit simultaneously: native admission
refuses excess ownership. Keep the 1 GiB comparison at 60. No automatic larger
tier retry. Measured prefixes used about 79–104 MiB owned and 100–114 MiB RSS;
the 70 capacity fixture reached 297 MiB RSS with sparse synthetic rows.
Host admission verified 24 GiB actual RAM and about 16.56 GiB available.
RSS checks are cooperative. Large emergency checkpoint serialization has an
existing preallocation limitation and is excluded from these pilots; small
charged resume remains verified.

| Stage | Calls | Per-call wall AND CPU caps | Whole-stage ceiling | Matched tiers |
| --- | ---: | --- | ---: | --- |
| A: 60 / 62 | 4 | 1800 / 2400 seconds | 9200 seconds | 1 / 2 |
| B: 65 | 2 | 3600 seconds | 7600 seconds | 2 |
| C: 70 | 2 | 7200 seconds | 14800 seconds | 3 |

[The executable freeze](../benchmarks/qs/c10a01/pilots.json) contains exact
per-call and cumulative work limits. Collection's maximum observed prefix
rate was 0.4445 billion units/CPU second; the 60 matrix fixture reached about
6.244 billion. The model rounds up to 8 billion, uses fourfold rate headroom,
10 billion setup units and a separate conservative analytic bound on all
postprocessing passes. These very large finite guards are not runtime
predictions. Full derivation is in local `v2/audit/c10a01/work.md`.

Including calibration, the worst-case envelope is **35200 wall AND CPU
seconds (9 hours 46 minutes 40 seconds)**. The latest combined 60-digit
426–796-second jobs motivate the anchor; the old 65-digit 1800-second timeout
predates that collector. Above 60 there is no calibrated completion-time law.
70 has its own arithmetic/resource admission and does not depend on a
65-digit completion. No automatic repeated pairs or enlarged retries.

## Verified handoff and commands

Base mainline: `4a946b8d6294d51ac5d9c36f61880cc63d5b59fd`.
Opt-in arithmetic: `e51b3f5`; scored-prefix freeze: `b666bde`;
forced/capacity freeze: `ce7d66f`. Required exposed inputs and recursive
Pocklington certificates are committed. Factor calls receive N only. Sampling
is size/residue/balance-conditioned, not uniform RSA sampling. The achieved
factor-bit lengths are 99/99, 69/129, 102/103, 71/134, 107/108, 75/140,
116/116 and 81/151 in size/balance order.

The task worktree is `/Users/nishanthmohan/Documents/code/factor-c10a01` on
`codex/c10a01-normalized-polynomials`. Source manifests bind the later runner
and every runtime Python source. Historical raw receipts preserve their own
source identities. The later runner uses bounded streaming state hashes to
avoid a whole-store JSON copy; it does not alter production arithmetic.

Readiness only (performs no factoring):

```sh
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.pilots check --stage A
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.pilots check --stage B
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.pilots check --stage C
```

Later execution commands, **not run during this session**:

```sh
/opt/homebrew/bin/pypy3 -u -m v2.benchmarks.qs.c10a01.pilots run --stage A
/opt/homebrew/bin/pypy3 -u -m v2.benchmarks.qs.c10a01.pilots run --stage B
/opt/homebrew/bin/pypy3 -u -m v2.benchmarks.qs.c10a01.pilots run --stage C
```

Each stage takes the shared performance lock, verifies source/input/host
readiness, reserves every call, warms both clocks and records nine short
sentinel samples before and after each measured job. A 15% drift screen is
recorded; failure prevents performance claims but does not erase outcomes.
Arithmetic/provenance/harness failures stop for audit. Native time/work/memory
and storage exhaustion remain censored outcomes; complete every valid matched
cell. Ratios are descriptive only when both arms complete. Otherwise report
finite-budget completion and capped service cost. No pilot yields an adoption
claim. Future confirmation needs a new selection identity before revealing
fresh independent inputs, adequate input-level uncertainty and matched checks.

Detailed research, conditional resource/work proofs, raw failures and receipts
remain under ignored `v2/audit/c10a01/` and
`v2/benchmarks/results/c10a01/`. No merge, push, default promotion or E1 work
was performed.

Validation: `make -C v2 test` passes 841 tests (9 optional-backend skips),
and `make -C v2 lint` passes. A committed-only export passes the same suite,
lint and all 204 benchmark/study module imports without local audit, results
or upstream files. Focused checks include 50/55-digit opt-in setup, routing
boundaries, observer passivity, original-row kernels and pilot censoring rules.
All three readiness-only commands pass source, exact PyPy binary/version,
input, resource and independent budget checks. Pilot execution remains zero.
