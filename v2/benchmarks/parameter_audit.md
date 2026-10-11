# Cross-engine parameter audit and calibration ledger

First pass: 10 October 2026. Owner: **C11** in [the roadmap](../ROADMAP.md).
QS/MPQS/SIQS and DLP retain their existing **C10/C9** owners and
[parameter audit](qs_parameter_audit.md). This document covers preprocessing,
primality, rho, p−1, ECM, SSS-specific search, shared execution and portfolio
choices. It records current choices and missing justification; it does not
claim that the parameter space has been optimized or change defaults.

Paths below refer to the active restructured working tree. The static inventory
covered 54 production Python files, including QS to identify ownership, and
recorded 4,025 numeric literal occurrences locally. Targeted source review
produced the grouped ledger below. The literal count includes exact
formula coefficients, residue tables, counters and serialization fields; it
is not a count of tunable parameters. Research downloads, source hashes and
the full mechanical inventory remain ignored under
`v2/audit/parameter-reconciliation/`. The grouped ledger below is the public
review artifact. C11 stays open until every performance-bearing default,
embedded threshold, cap and duplicate entry point has a disposition.

## C3 allocation disposition (10 October 2026)

The [C3 mechanism review](ecm/c3/research.md) and
[fresh decision](ecm/c3/acceptance.md) add explicit cumulative pretest ceilings,
finite campaigns and protected fallback admission. These caller-set amounts
are **L** finite policies, not success estimates. Default omission retains
existing bounds/counts and budgets. Schema 12 and the allocation version are
**I** identities; simultaneous relation/context/output caps remain owned-memory
reservations, separate from RSS.

Six jointly specified bound/count bundles reuse the first three v1 bound pairs
and a source-supported cheap pretest. Training selects quick8 for its smaller
bundle, but fresh cost is 51.56% higher with uneven regressions and instability;
retain numerical defaults. The larger selection retains control. No empirical
factor-size posterior, universal digit threshold, cross-engine work conversion
or deeper-campaign optimum follows. C3's 500M/15B-work and 1/10-second floors
are explicit service hypotheses; four of five smaller 2M-work probes refuse
unfunded admission. The unchanged default is not silently enlarged.

G1/E1 and continuing C11 calibration retain broader ownership. C9/C10 collector,
residual, splitter and graph settings remain frozen in this tranche.

## Where PRAC, Lucas and CF live

| Role | Current path / entry point |
| --- | --- |
| Bounded PRAC generator | [`ecm/prac.py`](../ecm/prac.py), `get_chain`, `_cached_chain`, `_prac_chain` |
| Lucas reference generator wrapper | [`build_c6_inputs.py`](ecm/c6/build_c6_inputs.py), `prepare` / `generate`; pinned C in [`inputs/upstream/c6_gmp_ecm/`](inputs/upstream/c6_gmp_ecm/) |
| CF reference search wrapper | [`build_c6_cf_inputs.py`](ecm/c6/build_c6_cf_inputs.py), `generate`; pinned search in [`inputs/upstream/c6_dacbench/`](inputs/upstream/c6_dacbench/) |
| Certified catalog preparation | [`build_c6_fast_inputs.py`](ecm/c6/build_c6_fast_inputs.py); existing C6 record/certificate conversion |
| Production PRAC / legacy explicit GMP Lucas | [`ecm/chains.py`](../ecm/chains.py), `ChainPlan`, `Action`, `execute` |
| Production optional Lucas / CF | [`ecm/chain_options.py`](../ecm/chain_options.py), `ChainPlan`, `ChainPlans` |
| Catalogs | [`c6_fast_records.json`](inputs/controls/c6_fast_records.json): ladder/PRAC/Lucas; [`b3_cf_records.json`](inputs/controls/b3_cf_records.json): production CF |

The production route consumes completed C6 research. It never invokes the
Lucas/CF searches during factoring. Catalog SHA-256 and exact byte lengths
are content identities/integrity checks, not arithmetic proofs or tuning
parameters. Scalar/frontier verification establishes the arithmetic guarantee.
The catalog digest joins backend, executor and recovery identities in plan and
checkpoint policy; a changed catalog needs explicit compatibility/migration.
Hashing/verification is charged on a miss, and a cache hit has its own charge.
A digest is not an authenticity signature. A manually maintained version could
identify a catalog, but would not itself detect accidental byte changes.

## What a justification must establish

- **P — proof/invariant:** exact formula, theorem endpoint, residue coverage or
  necessary precondition. Preserve the proof; do not sweep it as a speed knob.
- **R — representation:** packed word/index width, sentinel, record count or
  serialized format. Derive it from the format; changing it needs migration.
- **L — finite policy/cap:** a user allowance, cancellation bound or conservative
  owned-memory reservation. Finiteness is justified; the numerical allowance
  need not be economically optimal. Memory coefficients need an object/payload
  and simultaneous-lifetime derivation, distinct from measured process RSS.
- **M — bounded measurement:** supported by a named local experiment on its
  declared backend, inputs and reuse pattern. It is not a global optimum.
- **U — uncalibrated heuristic/baseline:** inherited, reference or convenience
  choice without adequate local selection evidence. Retain it as a control;
  make the missing experiment explicit.
- **I — identity:** digest, schema or version. Correctness/compatibility metadata,
  rather than an optimization knob.

Several rows have more than one role. For stochastic complete factoring,
"optimal" requires a workload distribution, certainty target, resource cap,
backend and objective. Upstream recommendations or a mathematical chain-length
minimum do not prove minimum PyPy time. The required outcome is an honest,
band-specific measured choice, sensitivity and uncertainty, or a documented
deferral—not an unsupported global-optimality claim.

## Primary implementations, papers and transfer limits

| Reference inspected | Exact pin / license | Applicable comparison |
| --- | --- | --- |
| GMP-ECM README, `pm1.c`, vendored `ecm.c` / Lucas generator | [`8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e`][gmp]; inspected C is LGPL-3.0-or-later; [preserved notices](inputs/upstream/c6_gmp_ecm/NOTICE.md) | Target-factor-size bounds/curve effort, p−1 reuse, stage-two memory/time policy and PRAC costs. Its native kernels and polynomial stage two differ from ours. |
| CADO-NFS `sieve/ecm/bytecode.c` | [`692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b`][cado]; root `COPYING` supplies LGPL-2.1; verify file-specific terms before any adaptation | Operation-specific cost structures, ten compiled PRAC multipliers with a conditional 18-choice table, compression, integer checking and an explicit cache-free function. Its global growing cache does not satisfy our finite per-job ownership by itself. |
| FLINT `fmpz_factor/factor_smooth.c`, `factor.c`, `pollard_brent.c` | [`a4c9750d0d3d67bb01cf6d18c187591b313451c3`][flint]; inspected file headers LGPL-3.0-or-later | Trial cutoff by prime count, factor-bit ECM table, finite rho tries/iterations. Source says the ECM tuning is rough in larger bands. |
| Yamaquasi `pollard_rho.rs`, `pollard_pm1.rs`, `params.rs`, `lib.rs`, ECM notes | [`3f95f43682ed15d8c1ed206a9a702dd655d7c8ad`][yama]; BSD-3-Clause | Size-tiered rho/pretest policy, stage-two table costs and automatic versus explicit ECM. Its native Montgomery/Edwards/FFT implementation is not ours. |
| SymPy `factor_.py`, `primetest.py` | [`2f22a5f81e2f4124380be3739a092e9ff20128de`][sympy]; inspected root license BSD terms, with per-file notices | Bounded early trial effort, perfect-power/Fermat checks, escalating factor searches and certainty. Its rho entry point is not the same Brent executor. |
| primesieve README / license | [`c8cfc3ed9065e9d42910a1a106a1bca2f39b7147`][primesieve]; BSD-2-Clause | Cache-aware segment sizing, wheel/bucket context and bytes-per-prime accounting. KiB packed native storage is not a count of Python candidates. |
| Brent, *An improved Monte Carlo factorization algorithm* (1980), §§4,7,8 | [Author's paper][brent] | Batch-GCD tradeoff and cycle schedule under a stated abstract cost/random-map model; not a proof of our batch/retry limits. |
| Bernstein–Cottaar–Lange, *Searching for differential addition chains*, §§1.1,2 | [ANTS paper][chains-paper] | CF-family length optimality and operation-cost caveat; length is not runtime. |
| Sorenson–Webster, *Strong Pseudoprimes to Twelve Prime Bases* | [Primary paper][mr-paper]; [accepted A10 reconciliation](primality/a10/a10_primality.py) | Strict deterministic endpoints; evidence must justify the exact base set and inequality. |
| Menezes–van Oorschot–Vanstone, *Handbook of Applied Cryptography*, §4.2.3, Fact 4.25 | [Authors' chapter][hac-primality] | Per-composite randomized MR bound with independently sampled witnesses; fixed reproducible seeds do not establish those probability assumptions. |
| Hittmeir, *Smooth Subsum Search*, v2, §§3–4; accepted upstream review | [Paper][sss-paper], [upstream `8dbaf6d`][sss-code] | Search/filter choices and relation-yield experiments. The accepted review found no explicit root license: no new copying/adaptation is authorized by availability. |
| GMP 6.3 number-theoretic API | [Official manual][gmp-primality] | `reps` combines BPSW and MR in current GMP, so it is not numerically equivalent to our MR witness count. |

Primary raw source retrieval succeeded for the implementation pins above;
local captures retain URL and SHA-256. Some browser views of pinned GitHub
URLs failed, so those views alone are not evidence; the retrieved source bytes
and previously accepted vendored files were inspected. No upstream code is
adapted in this audit. License approval and attribution remain required for
any later adaptation. "Leading" here denotes relevant reference implementations,
not a new measured ranking or a claim that each is fastest on every workload.

## Current source map

| Ledger area | Active implementation |
| --- | --- |
| Shared defaults and legacy dispatch | [`constants.py`](../constants.py), [`factor.py`](../factor.py) |
| Trial, powers, Fermat and primality | [`common/preprocessing.py`](../common/preprocessing.py), [`common/utils.py`](../common/utils.py), [`common/arithmetic.py`](../common/arithmetic.py) |
| Prime sieves and contexts/caches | [`common/prime_sieve.py`](../common/prime_sieve.py), [`execution/schedules.py`](../execution/schedules.py) |
| Rho | [`rho/brent.py`](../rho/brent.py) |
| p−1 | [`pm1/core.py`](../pm1/core.py), [`pm1/bounded.py`](../pm1/bounded.py), [`pm1/tuning.py`](../pm1/tuning.py), [`pm1/gaps.py`](../pm1/gaps.py) |
| ECM formulas, schedules and stage two | [`ecm/core.py`](../ecm/core.py), [`ecm/programs.py`](../ecm/programs.py), [`ecm/paired.py`](../ecm/paired.py), [`ecm/wheel.py`](../ecm/wheel.py) |
| Portfolio, bounded jobs and allowances | [`portfolio.py`](../portfolio.py), [`execution/stage_jobs.py`](../execution/stage_jobs.py), [`execution/budget.py`](../execution/budget.py), [`execution/work_budget.py`](../execution/work_budget.py) |
| SSS-specific search and workers | [`qs/sss.py`](../qs/sss.py), [`qs/smooth_batch.py`](../qs/smooth_batch.py), [`qs/parallel.py`](../qs/parallel.py); shared collector choices remain in the [C9/C10 ledger](qs_parameter_audit.md) |

## Preprocessing, primality and dispatch

| v2 choice and location | Role / present justification | Primary comparison and remaining work |
| --- | --- | --- |
| `constants.TRIAL_BOUND=25_000`; `PortfolioConfig.trial_chunk=64` | U/L. Trial bound is inherited from v1; chunk bounds progress/polling. No universal optimum. | SymPy's early pass uses a bound `2**15` with a separate 600-miss early stop. FLINT uses a count of primes (up to 3,512), not a prime-value bound. Normalize actual prime work/yield; E3 owns the experiment. |
| `SIZE_THRESHOLD_RHO=10**20` in legacy `factor.py` | U. Inherited total-cofactor cutoff, about 67 bits; bounded portfolio instead has explicit rho effort. | Yamaquasi's automatic rho path has explicit size tiers and skips >64-bit inputs. This reflects its other native engines, not a transferable cutoff. Compare legacy/bounded paths and handoff marginal value through E4/G1. |
| `fermat_steps=0` in bounded config/CLI | U/L. Disabled search is a finite policy, not evidence that near-square preprocessing never pays. | SymPy tries a short Fermat pass; compare close-factor and ordinary inputs, charging failed probes. E3 owns activation and effort. |
| Power prefilter `_POWER_MODULI`: `{2:5,3:7,5:11,7:29,11:23,13:53,17:103,19:191,23:47,29:59,31:311,37:149}` | P/U. Every pair has prime `q` and `k | q-1`; rejection is a necessary-condition proof, followed by an exact integer root. The particular selection/extent is not time-optimal. | Exact-root/perfect-power pretests in FLINT/SymPy are the control designs. Profile screen cost/rejection by exponent before changing which valid moduli to use. Preserve the independent implication proof. |
| `SMALL_PRIMES=(2,...,37)`; shortcut `n < 41**2` | P/U. The shortcut follows complete divisibility checks through 37; its boundary is not independent. The amount of trial screening can be tuned jointly. | Compare screen cost versus deterministic/probable tests; do not change `41**2` independently of the screened set. |
| Deterministic base sets and strict limits `2**64`, `318665857834031151167461`, `3317044064679887385961981`; short tiers `(31,73)`, `(2,7,61)` | P/M. Accepted A10 theorem/range choices, not guessed thresholds. | Preserve strict inequalities and exact bases from the primary results. Reorder/use a cheaper proved tier only with proof plus full-call cost evidence. Certainty labels cannot be relaxed by tuning. |
| `PRIMALITY_ROUNDS=30` in native/GMP wrapper and portfolio | L/U. A confidence policy rather than a time optimum. The per-composite MR bound is at most `4**(-30)=2**(-60)` under independent uniform witnesses ([HAC Fact 4.25][hac-primality]); it is not a proven-prime label, cryptographic RNG claim or unconditional whole-run bound. | GMP 6.3 `reps=30` means BPSW plus `reps-24` MR tests, not 30 MR witnesses. C11 must declare per-call/cumulative certainty targets and seed assumptions before testing cheaper policies; preserve A10/legacy resume. |
| Backend `python-int` default, optional `gmpy2-mpz` | M/U. Accepted A2/B4/C6 comparisons are scoped; no universal bit-size crossover. | Native FLINT/GMP arithmetic costs cannot prescribe PyPy dispatch. Keep separate tracks and compare coarse conversion-inclusive complete calls; no new arithmetic backend/kernel is part of parameter calibration. |

## Prime sieves and shared schedules

| v2 choice and location | Role / present justification | Primary comparison and remaining work |
| --- | --- | --- |
| `SMALL_THRESHOLD=60`, `UNDER_60`, Atkin residue/DFG tables | P/R for exact tables; U for dispatcher cutoff. | A larger tiny-table cutoff needs complete prime coverage. Residue coefficients encode exact congruences, not tunable speed weights. Atkin is retained as a checked alternative, not selected by an unmeasured threshold. |
| `ERAT_THRESHOLD=3_500_000` | U. Inherited v1 Eratosthenes/segmented switch. | primesieve uses cache-aware segmented native storage. Tune the Python crossover on actual schedule endpoints and short/high-offset intervals; E6 owns this. |
| `LOWER_SEGMENT_SIZE=65_536` candidate positions | U/L. Current standalone segmented-sieve default, not 64 KiB of all live storage. | Normalize candidate count, byte flags, marking slices and base-prime headers against native cache sizing. Measure full schedule generation and factoring, not count-only throughput. |
| `UPPER_SEGMENT_SIZE=2_097_152` | U, dormant. No live production consumer found in this pass. | Do not include a dead knob in a timing grid. C11 should record compatibility/dead status; removal is a separate cleanup with API review. |
| `SieveContext`: 1 MiB default, segment 4,096; bounded ECM/p−1/portfolio segment 1,024 | U/L. Distinct callers and finite scratch policies, not one optimized size. | Compare source-supported cache scaling, short terminal segments, setup amortization and allocations under E6. Preserve exact half-open endpoints. |
| `ScheduleCache`: 65,536 bytes, 8 entries; portfolio `schedule_cache_bytes=0`, `rolling=False` | M/L/U. Prior reuse/cache/rolling studies retain current choices on their workloads; cache identity/lifetime is exact. | Reopen only if measured workload/reuse changes. Charge build/miss/eviction and resume; primesieve's persistent native structures are not free cross-call ownership for v2. |
| Sieve reserves `4096 + 18*base_count + 3*segment_size`; cached schedule `1024 + 16*count` | L. Conservative simultaneous-owned-memory model, not measured optimal constants or process-RSS caps. | Derive each live container/copy and packed payload; reconcile allocator/JIT RSS separately. C11 may tighten a coefficient only after the ownership argument and refusal/resume tests survive. |

## Rho

| v2 choice and location | Role / present justification | Primary comparison and remaining work |
| --- | --- | --- |
| `RHO_ATTEMPTS=16`, `RHO_EVALUATIONS=50_000`; bounded config instead `4 × 5_000` | U/L. Finite standalone/portfolio effort policies, not equivalent searches or proven optima. | Yamaquasi uses size-tiered iterations; FLINT exposes tries and iterations. E4 should compare total-work-normalized restart/long-walk bundles and marginal completed-factor yield. |
| `RHO_BATCH_SIZE=64` | U/M. Accepted Brent implementation with prior loop evidence, no demonstrated universal batch optimum. | Brent §7 trades gcd frequency against accumulated-product cost; Yamaquasi checks a product every 128 iterations after its initial phase. Different loops/count units forbid copying 128 as a prescription. |
| `RHO_RECOVERY_LIMIT=128`; shared recovery value also appears in p−1 config | L/U. Bounded replay/abandonment policy; not independently calibrated for every engine. | Measure saturated/mixed-factor recovery and wasted long-walk work. Separate engine meanings before changing a shared constant. Every found divisor remains validated. |
| Cycle-length doubling (`<<=1`) | Algorithm/U. Brent's convenient near-optimal schedule under his stated model, not a proof for all PyPy costs. | Preserve it as control. A changed cycle schedule needs its own termination/replay proof and marginal evidence; do not silently fold an algorithm change into a batch sweep. |
| Seed, initial point and polynomial offset sampling | L/U. Reproducible finite assignment policy, not an optimal deterministic polynomial sequence. | Freeze assigned seeds and compare retries without cherry-picking lucky factors; distinguish workload variation from repeated timing variation. E4 owns any changed policy. |

## p−1

| v2 choice and location | Role / present justification | Primary comparison and remaining work |
| --- | --- | --- |
| `PM1_B1=2_000`, `PM1_B2=200_000` | U/L. Fixed 100× ratio and modest control; A6 did not optimize the entire bound surface. | GMP-ECM chooses method-specific B2, while Yamaquasi has stage-two cost tables. Tune smooth-factor work/coverage and portfolio opportunity cost; C3 owns bound/allocation policy with A6 as retained executor. |
| Standalone `PM1_ATTEMPTS=3`, initial base 2 then `base+attempt`; bounded `pm1_attempts=1` | U/L. Different finite entry-point policies. | GMP-ECM explicitly cautions that repeating the same p−1 bounds with another base usually adds little; saturation/base exceptions can still matter. Compare useful recovery versus redundant smoothness work, not ECM-like independent-curve assumptions. |
| `PM1Config.chunk_size=16`; fresh `factorize_pm1_bounded`/portfolio use recurrence + chunk 64 | M/L. User-promoted A6 complete-stage result; explicit config preserves legacy executor. | A6 measured a bounded choice, not a proof that 64 is best at all B1/input sizes. Reuse its source/confirmation rather than rerun it indiscriminately; recalibrate only for a new workload or cost gap. |
| Shared `GCD_BATCH_SIZE=128` | U/L. Stage-two product/recovery granularity, not a universally measured optimum. | Compare gcd/product/replay cost with both current bounded configurations and GMP-ECM's different polynomial continuation. C11 needs a separate p−1 and ECM disposition. |
| Recurrence cache 64 even powers, exceptional gap cache 64; `gap_entries=64` (cap 256) | M/L/U. A6 supports recurrence in its bounded experiments; capacity 64 is not a proven economic optimum. | Profile observed gap coverage, setup, fallback exponentiation and empty-cache rebuilding. Retain exact residue/B1 identity and charged reconstruction. |
| `chunk_bits=0`; nonzero range 32–4,096; wheel `0` with allowed 30/210 | M/L/R. A6 experimental alternatives, finite exponent cap and supported wheel encodings. | Existing accepted/rejected evidence is the control. Source wheel geometry is an algorithm/representation choice; do not promote on scalar-bit counts or lower reservation charges alone. |
| `segment_size=1024`, `memory_bytes=8 MiB`, `max_input_bits=4096`; chunk/gcd validation ceiling 256 | L/U. Finite API/cap policies. | No upstream memory/iteration limit proves these optimal. Derive simultaneous workspace, then measure setup/recovery/continuation with matched total resources and exact legacy resume. |

## ECM, chain selection and stage two

| v2 choice and location | Role / present justification | Primary comparison and remaining work |
| --- | --- | --- |
| `ECM_B1=2_000`, `ECM_B2=147_396`; `compute_bounds(n)` currently returns fixed values | U/L. Traceable to v1's first bound-table entry; v2 deliberately retained a modest baseline. It is not a tuned size policy. | GMP-ECM's target-factor-size table and FLINT's factor-bit table use different probabilities and native costs. C3 must fit finite pretest/campaign tiers to measured marginal factor yield; input digit count alone cannot identify the smaller factor. |
| `MAX_CURVES_ECM=32`; portfolio repeats literal 32 in its tier | U/L. Finite cap; default bridge tested 32, not all alternatives. | FLINT explicitly targets roughly 1/3 success with its table; GMP-ECM publishes expected effort by factor size. Duplicate choices need an entry-point consistency disposition, then C3/G1 allocation evidence. |
| Sigma lower bound 6; `MAX_RANDOM_ECM=2**63`; bounded `randrange(6,2**63)` versus standalone inclusive `randint` | P for Suyama preconditions; U/R/L for sampling envelope. Inherited reproducible seed range, not an optimal curve distribution. | Inspect exact accepted parametrization/setup assumptions. Preserve saved seeds/endpoints; curve-family changes remain F2. Do not silently unify inclusive/exclusive sampling because that changes reproducible assignments. |
| Suyama coefficients 5,4,16; `(A+2)/4`; ladder/differential formulas | P. Exact algebraic definitions/coordinate invariants. | These are not knobs to tune. New formulas/curve families require B4/F2 proof and comparison, separate from constants calibration. |
| Stage-one chunk 16; default native reduced PRAC; optional Lucas/CF; automatic GMP ladder | M/L. C6/B3 evidence and explicit user default direction; first-use and per-call regressions remain documented. | CF-family minimum length does not imply minimum runtime (Bernstein–Cottaar–Lange §1.1). Reuse accepted records and kernels; compare preparation/reuse and full calls separately per backend. No repeated C6 search is needed for C11. |
| Chain route B1 exactly 2,000; >=8 planned curves; cofactor `[10**39,10**80)`; chunk exactly 16 | M/L, evidence-support boundaries. Eight curves is observed amortization in a declared stage cohort, not a guarantee a job will execute eight curves before finding a factor. | Charge actual misses/early successes/eviction/late resume. Wider routes require independent coverage and confirmation; they do not follow from this audit. |
| PRAC `add_cost=6`, `double_cost=5` | U/M. Matches pinned GMP-ECM `ADD=6.0`, `DUP=5.0`; the local docstring counts `4M+2S` and `3M+2S`, treating S like M and omitting Python dispatch/reductions. This is an abstract model, not measured PyPy cycles. | CADO supplies operation-specific cost structures. The chain paper separates M/S/constant/addition costs. Use existing B4/C6 cost evidence first; any reweighting needs independently verified records plus construction and whole-job evidence. No new search in this audit. |
| Ten PRAC ratios over `10**17`; neighboring splits ±1 | U/M/R. All ten decimal values match pinned GMP-ECM's reciprocal continued-fraction choices: all ones, or a 2 at the indicated position. Exact integer rounding avoids float dependence; ±1 extends a bounded heuristic neighborhood. No proof says ten ratios or that neighborhood is optimal. | GMP-ECM tries at most ten choices, limited by modulus limb count; v2 examines up to 30 splits independent of modulus size. CADO compiles ten multipliers and includes a conditional 18-choice table with an additional family; normalize its reciprocal representation before comparison. Those are candidate rationales, not permission to redo completed C6 search. Document construction cost and a new workload trigger before any extension. |
| PRAC nine-rule coefficients/guards: `4*d <= 5*e`, divisibility by 2/3/6, `d <= 4*e` | P/U. Exact integer forms of Montgomery's Table 4/GMP-ECM reductions; the local invariant preserves `d*a+e*b=odd`, known differences, gcd and a decreasing positive `d+e`. These are not independent arbitrary batch constants. | The accepted integer verifier proves every emitted record. A changed reduction/guard/order needs a fresh termination/invariant argument and independent scalar/coordinate-factor coverage; a cheaper abstract cost alone does not establish valid arithmetic or fastest execution. |
| Offline Lucas B1=2,000, one thread, capacity4,096 and decoder64 elements; CF depth18/node10M; subprocess60 s CPU/wall, 512 MiB RSS watchdog, 16 MiB file cap | L/R/M. Accepted C6's bounded research envelope. Lucas changes only the original six-million-entry allocation capacity; CF keeps upstream branch order with the documented exact-integer patch. | These are search/storage safeguards, not mathematically optimal search parameters. Preserve the frozen records and independent scalar/CF-family-minimum checks. A larger bound or changed search is a separate research gate, not required for C11 or production factoring. |
| PRAC scalar width 32 bits, 512 steps, 512 LRU records; cost key width <=32 | L/R. Search/cache safeguards; above the supported scalar width the ladder remains valid. | No mathematical theorem makes 32/512 optimal. First profile whether this public generator path matters; production catalogs already avoid runtime search. Keep data-only caches, cancellation and bounded construction. |
| Catalog bound 2,000; 333 scalar records; 512 steps, 16 point slots, 4-byte operations, double sentinel255 | R/L. Current verified format/catalog coverage; record counts and bytes derive from frozen inputs. | Limits/encoding are not economic recommendations. A changed representation needs independent scalar/frontier verification, memory/recovery derivation and checkpoint compatibility. |
| Plan 4 MiB + simultaneous scratch4 MiB; automatic total16 MiB; packed programs >=512 KiB | L/U. Conservative ownership and separately user-approved capacity, not an empirically optimal split. | Derive per-record/container lifetimes, never count cache verification/reconstruction as free. A global/cross-call cache needs a separate finite lease API; native GMP-ECM lifetime assumptions do not supply it. |
| Stage-two `D=min(isqrt(B2),(B1-1)//2)`; even-baby/giant table | Algorithm/U. Valid B4-stage geometry and bounded memory, not demonstrated optimal D. | Yamaquasi's `stage2_params` accounts for totient/polynomial-degree tradeoffs; GMP-ECM changes block count to trade time for memory. Those engines differ. C3/P5.2 owns measured D/table/bound decisions, B2 owns pairing. |
| Optional pairing distance/wheel defaults `None`; packed programs/cache caps; cache hits/misses | M/L. Existing A3/B2 retained decisions and exact range constraints; no universal win. | Prior construction and pairing regressions stay controls. Reopen on a changed reuse workload, count unmatched primes, all construction and recovery, and preserve coordinate-factor coverage. |
| Catalog SHA/byte length, kernel/family/backend/schema identities | I/R. Exact content and compatible execution accounting. | Do not tune away verification to claim arithmetic speed. If reducing identity/checkpoint overhead is material, measure it while preserving exact compatibility and independently verified records. |

## SSS-specific search and parallel execution

C10 owns shared QS relation/scoring/factor-base choices; C11 must not duplicate
those experiments. The following SSS-specific parameters need an explicit
calibration owner even though their implementation lives under `qs/`.

| v2 choice and location | Role / present justification | Primary comparison and remaining work |
| --- | --- | --- |
| `SSSConfig.base_bound=1000`; `small_bound=0` chooses first fifth of base | U/M. Reference split; accepted SSS repairs do not establish a universal optimum. | Hittmeir §§3–4 and upstream settings are hypotheses; normalize prime count versus prime magnitude and current collector cost. E1 owns confirmation against calibrated SIQS. |
| `selection_size=6`, `collision_min=3`, `search_rounds=256` | U/M/L. Existing six/seven-selection arms and collision correctness; finite rounds. | Separate forced divisibility, candidate yield and actual dependencies. Freeze a small training grid only after profiling; retain both unfiltered and explicitly lossy modes. |
| `max_candidates=4096`, tree bit cap262144 and node cap8192; hard ceilings1048576 bits/8192 nodes | L/U. Finite candidate/product-tree capacity; exact saturation/refusal paths. | Bernstein's [smooth-parts method][smooth-parts] supplies the exact batch algorithm, not these capacities. Measure truncation/refusal, full-base processing, simultaneous tree storage and charged recovery. |
| `filter_divisor=10` chooses first tenth for SSSf first stage; `filter_bound=0` disables lossy cutoff | U/M. Repaired two-stage admissibility and measured retained decisions; first/second base fractions differ. | Paper relation-yield evidence is not a full-factor default. A positive cutoff deliberately loses candidates: report that loss, never present it as an exact rejection proof. |
| SSS memory64 MiB, checkpoint256 KiB; shared collector residual10000/atoms4096/rows2048/partials512 | L/U. Finite reference capacities; collector knobs belong to C9/C10. | Measure useful retained yield/copy/rebuild cost under one envelope. No new root-license permission was found for upstream adaptation. |
| `ParallelConfig`: assignment10M work, max batch atoms1024, batch width0, poll64 | M/L/U. Coarse worker lease/publication policy; accepted repairs are scoped. | A7/E1 own workers and confirmed bundle comparisons. Measure first-factor latency, wasted completed work, serialization and aggregate CPU; native thread counts are not Python worker optima. |
| Parallel parent64 MiB / worker32 MiB / total512 MiB / checkpoint4 MiB; worker family defaults | L/U. Conservative finite coexistence and small reference workload. | Derive peak parent/workers/result copies; C10 owns family/collector tuning, E1 worker feasibility. Increasing a limit alone does not establish useful scaling. |

## Shared budgets, reservations and ownership

| v2 choice and location | Role / present justification | Remaining obligation |
| --- | --- | --- |
| `Budget`: 2M work, 30 s wall, 30 s CPU; CLI also30 s | L/U. Explicit finite service policy, not an optimal engine budget. | G1 must calibrate marginal handoff under matched total resources. Different work currencies can censor different searches; preserve coverage and completion alongside time. Never claim speed by finishing less unsuccessful work. |
| Portfolio trace256, max input4096 bits, max64 tiers, max256 arithmetic chunks, trial64 | L/R/U. Bounded state/result/cancellation policy, no demonstrated economic optimum. | Document whether changing each value affects RNG, event evidence, checkpoint size, atomic latency or search effort; migrate/pin any changed saved policy. |
| `PollingBudget.interval=64`, hard maximum64 | L/M. Existing bounded cooperative-polling contract, at most63 additional bounded actions after an external stop before polling. | This is an action-count bound, not a wall-time guarantee. Profile atomic action sizes and first-factor/cancellation latency; do not relax the contract to obtain an unexplained timing gain. |
| Portfolio memory8 MiB for unsupported/off/GMP and automatic16 MiB for eligible native chains | L/U. Finite user policy; explicit caps never raised; resume restores saved cap. | Compare total owned memory versus RSS separately. No upstream native allocator setting proves these the best budgets. |
| Workspace coefficients: coordinate `128+bits//8`; base16384; segment256×; input256×; trace1024×; chunks256×; output/context8192 | L. Conservative accounting assumptions, not arithmetic identities or measured exact PyPy heap sizes. | Provide an explicit live-object/payload/copy derivation for each coefficient and simultaneous high-water case. Tightening reserves needs proof/cap/refusal/cancel/resume checks, not only fewer estimated bytes. |
| Chain miss verification units, point/GCD units, worst-case strict/unit replay and no recovery refund | L. Exact cumulative bookkeeping and finite upper-bound reservations; units are not CPU cycles. | Preserve whole-chunk atomic reservation. Benchmark matched real elapsed/CPU and common search coverage; retuning a work model is a separately versioned API/accounting change. |
| Checkpoint byte caps, repeated reconstruction and identity hashes | I/L/R. Compatible finite serialization/rebuild policy. | Charge validation, empty caches and reconstruction. A cheaper path must preserve proper-divisor checks, every unresolved cofactor, certainty, canonical points, RNG and cumulative allowances. |

## Ownership and bounded experiment queue

C11 is the umbrella ledger/calibration gate, not a new competing factoring
implementation. It coordinates these existing owners:

| Area | Executor owner | C11 deliverable |
| --- | --- | --- |
| Trial/powers/Fermat/primality | E3; A10/B14 for certainty/proofs | Confidence/coverage-preserving candidate rationale and marginal full-call cost |
| Rho batch/walk/restart | E4 | Separate standalone/portfolio settings and total-effort-normalized selection |
| ECM/p−1 bounds, curves, stage-two geometry and handoff | C3, P5.2/P5.3, G1 | Factor-size/reuse-aware policy, including unsuccessful searches and preparation |
| Sieve/segment/cache/setup | E6 | Endpoint/reuse/context-based cost and conservative storage justification |
| QS/MPQS/SIQS / DLP | C10 / C9 | Cross-reference the existing ledger and selected bands; do not duplicate |
| SSS-specific search / workers | E1 using A7 contracts | Train missing feasible arms, then confirm versus calibrated controls |
| Combined selected settings | E1/D7/G1 | Frozen integrated control and untouched complete-factor confirmation |

The first tranche should resolve cheap high-value provenance and observed cost
questions: duplicate entry-point defaults, dormant knobs, p−1 repeated-bound
utility, ECM bound/curve allocation, and setup/cap assumptions. It should not
restart completed C6 chain search or rejected kernel/reducer work. Record the
chosen priority using a separate profile and pilot yield evidence, not guessed
importance. Where the existing workload shows no limiting cost, defer the knob
with a concrete trigger rather than sweep it automatically.

For each bounded tranche:

1. Record every candidate's v2 symbol/current value, caller, units, role,
   upstream revision/license, model assumptions, dependency and missing evidence.
   Include embedded literals, validation limits, derived formulas and overrides,
   not just `constants.py`. An invariant, dormant field or retained user cap is
   a valid explicit disposition; an unexplained numeric choice is not.
2. Specify the objective (marginal proper-factor yield, complete time,
   completion at fixed resources, memory or cancellation latency) and input/
   smaller-factor/reuse bands. Derive a small grid from source policies, measured
   bottlenecks and neighbors of the retained control. Publish why each candidate
   is included; bound the search, simultaneous storage and total experiment time.
   Pilot failures have a stop/deferral rule. Candidate grid values are hypotheses,
   not new defaults or already proven optima.
3. Freeze source, catalogs, training/untouched-confirmation inputs, seeds,
   bounds, curves, total allowances, selection rule and early-stop policy before
   accepted timings. Do not learn smaller factors from the confirmation answers
   and feed them to dispatch. Separate one-knob attribution from selected-bundle
   evaluation, and keep integer/GMP tracks separate.
4. Use PyPy implementing Python 3.11, >=3 s validated warmup and >=9 samples,
   extending instability by a prespecified rule. Separate cold startup and
   instrumented profiles. Serialize accepted timing/heavy checks with the
   machine-wide performance lock. Run no new sweep just to repeat settled work.
5. Charge setup, catalogs/proofs, cache miss/hit/eviction, unsuccessful search,
   recovery, recursive classification, serialization and rebuilding. Validate
   every proper divisor and full reconstruction including unresolved cofactors.
   Report completion, search coverage, censoring, CPU, RSS/owned memory,
   chronological drift and seed/workload variation. Time samples alone do not
   establish population uncertainty or mathematical optimality.
6. Apply the current revised roadmap policy, not a universal 10% floor. Confirm
   selected bundles on fresh untouched inputs and retain defaults where evidence
   is inconclusive. Publish sensitivity/bracketing, rejected/deferred choices,
   limitations and reproducible commands. Do not select a new optimum from the
   same confirmation cohort. Any certainty, work, memory or resume change needs
   its own proof/compatibility acceptance and relevant tests/lint/committed-only
   imports before production promotion.

No timing experiments or production changes accompany this first-pass audit.
C9/C10, C3/E3/E4/E6/E1/G1 remain open; adding C11 does not mark their experiments
complete. Maintain the ledger whenever engine/backend or feasible workload
changes, and close only named, confirmed tranches.

[gmp]: https://github.com/sethtroisi/gmp-ecm/tree/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e
[cado]: https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/bytecode.c
[flint]: https://github.com/flintlib/flint/blob/a4c9750d0d3d67bb01cf6d18c187591b313451c3/src/fmpz_factor/factor_smooth.c
[yama]: https://github.com/remyoudompheng/yamaquasi/tree/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad
[sympy]: https://github.com/sympy/sympy/blob/2f22a5f81e2f4124380be3739a092e9ff20128de/sympy/ntheory/factor_.py
[primesieve]: https://github.com/kimwalisch/primesieve/blob/c8cfc3ed9065e9d42910a1a106a1bca2f39b7147/README.md
[brent]: https://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf
[chains-paper]: https://antsmath.org/ANTSXVI/papers/BernsteinCottaarLange.pdf
[mr-paper]: https://arxiv.org/abs/1509.00864
[hac-primality]: https://cacr.uwaterloo.ca/hac/about/chap4.pdf
[sss-paper]: https://arxiv.org/html/2301.10529v2#S4
[sss-code]: https://github.com/sbaresearch/smoothsubsumsearch/tree/8dbaf6d39ab88a40380965d25ec2c363d7f27358
[smooth-parts]: https://cr.yp.to/factorization/smoothparts-20040510.pdf
[gmp-primality]: https://gmplib.org/manual/Number-Theoretic-Functions
