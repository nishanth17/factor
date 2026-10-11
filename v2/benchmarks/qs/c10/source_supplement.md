# C10 source expansion: normalization, diversity and feedback

Reviewed 2026-10-11, before C10 measurements. These findings refine the
[ranked gaps](research.md), without changing the first tranche's two isolated
arms. The [manifest](../../inputs/controls/c10_source_pins.json) pins the
additional files. Only source reading and paper inspection were performed;
neither upstream implementation was compiled or run.

## Literature leading to executable sources

Kechlibar's [implementation paper](https://www.karlin.mff.cuni.cz/~krypto/Implementace_MPQS_SIQS_files/main_file.pdf),
§2.3.2 and §2.3.8–10, treats the normalized square identity, the endpoint-based
A target, and disjoint core/corrective pools. It distinguishes unique A values
from their quality and discusses early-combination product bias. Its numerical
integral-minimizing approximation is an author calculation, not a constant
adopted here. The existing integer oracle establishes the transfer identity
independently. Printed pages 18–19 were also visually inspected to check the
normalization and interval convention. PDF SHA-256 is
`14c268f31e6cb7d3d18efb61ec59fbac8ba3759d8cdb153b09ee6f28f75379e7`.

The [publisher's 2010/2015 implementation](https://www.karlin.mff.cuni.cz/~krypto/Implementace_MPQS_SIQS.html)
provides the executable comparison below. Its source archive hash is
`61cad09211add0eaf3e3b012b17d39247da4c20b15665d040e21091b2d7f4e90`.
The publisher specifies LGPL without a version; the inspected source headers
do not resolve that version. No licensed implementation text is incorporated
into v2.

Percival's [2005 author note](https://www.daemonology.net/blog/2005-11-30-quadratic-sieve-constant.html)
identifies per-prime initialization and no-hit prime visits as obstacles to
shorter intervals. This motivates checking sparse-prime amortization in actual
code. Its projected speedups and crossover sizes are historical predictions,
not v2 measurements.

The [CUDA-MPQS paper, v1](https://arxiv.org/html/2610.07126v1), dated 2026-10-05,
describes polynomial batching, accumulator saturation and expensive truncated
autotuning. It reports that tuning overhead can exceed the saving within one
run and that short probes underobserve large-prime matching. These are useful
mechanism and failure-mode leads. GPU factorization records, speedups and
implementation-specific dependency failures were not reproduced and are not
evidence about v2. Claims below come from separately inspected source; the
paper and tagged source differ in relevant details.

## Prague SIQS: active and historical paths

The source paths here are relative to the pinned archive. `ks/ks/mpqs.cpp`
constructs `QSParameters`, whose defaults select `ESiqs`. In
`QuadraticSieve::PerformSieving` (`quadratic_sieve.cpp:2419`), that mode goes
through `LineSiever::RunSieve`; the later loop is the older implementation.
Both call the shared polynomial generator and candidate admission code.

| Investigation | Executable finding | v2 classification / consequence |
| --- | --- | --- |
| N01/N08: mod-8 representation | `CalculateKRoot:3296` normally chooses k=n mod 8, so odd n gives D=1 mod 8. `CreateNewSiqsPolynomial:5829` makes b odd. The final `ComputeC:6215` divides by **4A**; an earlier temporary division and commented-out line do not define the final coefficient. | Normalized representation is active, unlike v2's standard F. This corroborates the oracle contract, not the desirability of restricting multiplier choice to this one residue. The optional extended-k branch is a separate heuristic, not Knuth–Schroeppel scoring. |
| N02: actual product versus log approximation | `CreateNewAForSiqs:5991` constructs a unique odd-index core, then binary-searches quantized log sums of even-index triples. Its prose says closest product, but the executable choice compares log-distance surrogates. | v2's corrective flyer and new quality batches compare exact integer products. Do not describe the Prague choice as exact product minimization. Copying a large triple table would add construction and storage costs not justified by current evidence. |
| N02: diversity | The same function explicitly records that it does not check a two-entry difference from the previous A. Disjoint pools ensure identifiable cores, not low overlap between families. | Supports keeping the product-error and overlap censuses separate. A bounded diversity selector remains an uncalibrated future policy; no mathematical implication from uniqueness to useful dependency diversity is assumed. |
| N03: active feedback caller | `LineSiever::RunSieve:1354` invokes `GetOptimizer().AnalyzeAndOptimize` every 500 lines while collection is unfinished, then writes back the byte initialization complement. `SaveCandidate:8199` increments accepted/rejected counters after recovery and residual handling. | This is active per-job feedback, not merely a dormant helper. It optimizes accepted/rejected candidate ratio, with no elapsed-cost or post-filter dependency objective. v2's exact block/value thresholds exist; this feedback policy does not. |
| N03: policy and lifecycle caveats | `QSParameters::InitDefaults` sets policy 0; `-to=` can override it. `ThresholdOptimizer::AnalyzeAndOptimize` requires 10 accepted or 100 rejected. The zero-rejected branch adds 3 before the policy switch and does not clear counters. Ordinary branches choose a ratio-based step and reset counters. | The nominal off setting does not guard every mutation in this source. This is a source-level observation, not a runtime reproduction. A v2 controller needs an unconditional disabled path, bounded counters/steps and explicit checkpoint state. Do not port these branches or their logging into the quiet API. |

The historical `PerformSieving` loop also contains threshold adjustment code;
it is not the default `ESiqs` caller. The active line-siever caller is essential
to the classification above. Residual admission semantics remain C9-owned.

## CUDA-MPQS v1.0.8: transferable mechanisms and limits

Pinned tag `v1.0.8` resolves to
[`99c61967585b97a71b28fbaf19690af83ac88c17`](https://github.com/drjanosch42/cuda-mpqs/tree/99c61967585b97a71b28fbaf19690af83ac88c17).
All inspected files match the Git tree's blob hashes. Source SPDX identifies
LGPL-3.0-only; `LICENSE.md` adds the NVIDIA CUDA Toolkit linking exception.
GPU architecture, fixed-width integers and CUDA memory geometry are native
constraints, not PyPy prescriptions.

| Investigation | Pinned executable path | v2 classification / consequence |
| --- | --- | --- |
| N02: A scheduling | `prime_algorithms.cu::init_a_factors` searches actual products in a target-to-100/95-target window using every other eligible index, starting at index 150. `advance_a_factors` exchanges neighboring primes by Gray steps, then slides the window. | Actual-product windowing exists, but this is neither a nearest-product ordering nor a minimum-overlap guarantee. v2's core permutation/flyer has different population semantics. Index 150 and this native size range are not portable calibration. |
| N09 / A1a: sparse-prime amortization | `device_sieving_controller.cpp` selects batch geometry; `kernel.cu` initializes family data, fills buckets across polynomial/subinterval pairs and dispatches sieve/scan/recovery kernels. | Cross-polynomial sparse planning is distinct from v2's per-polynomial capped power plan. A1a should measure remaining sparse setup after N04 and account for all simultaneously live roots/buckets before choosing this larger kernel change. CUDA atomics, shared-memory occupancy and packed layout are native-only. |
| N06: score saturation | `DeviceSievingController::initiate:117–206` selects a wide accumulator and then permits saturating bytes if a sampled target plus safety margin passes its gate. The comment at 168 explicitly says the three-point sample and one A bit length are **not a proven per-location/per-polynomial upper bound**. | The paper's description should not be upgraded into a general exact-bound theorem. v2 already saturates scores and thresholds together (`sieve_collector.py::_clip_threshold` and `_add`), allowing extra candidates conservatively. That differs from claiming bit-identical candidate selection. Preserve the existing exact-value checks and exception tests. |
| N03/N05: optional per-job autotuning | `orchestrator.cpp:311` constructs `AutotuneController` only for autotune-only mode or enabled full/sieve-only modes. `run:91` skips inputs below 70 digits. Its stages tune kernels, estimate runtime and optionally search F/L. | This path is optional and outside C10's current 30–50-digit population even in its own implementation. Persistent history and joint parameter search are separate policies, not proof that a small per-job pilot pays for itself. |
| N03: objective and observer bias | `runtime_estimator.cpp:128–190` projects SLP matches with a birthday model, caps the estimated contribution at 35%, then adds heuristic matrix/linear-algebra costs. `autotune.cpp:703` sets 40 search probes and a remaining wall allowance; `sieve_optimizer.cpp::cachedEstimate` caches repeated configurations. | A useful dependency forecast requires evidence beyond recovered candidates. The paper reports budget overrun; this pin also prevents rerunning stage 3 with a fresh cache. Do not treat either the comment or design cap as a verified hard total bound. v2 pilots must precharge all work, include tuning overhead and confirm full outcomes. |

The paper's large-prime dependency failures are specific to that pipeline.
They do not justify new C9 admission limits, an SLP ban for small inputs, or a
matrix change in C10. The current tranche retains settled C9 residual and
storage controls. Revisit N03 only with a bounded lost-candidate/useful-row
audit; revisit N09 only with a measured remaining sparse-setup cost.
