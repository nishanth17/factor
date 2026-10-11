# C9.2 adaptive storage and useful large-prime yield

Status: **C9.2 bounded tranche accepted, 11 October 2026.** Adopt the explicit
full-base, endpoint-aligned growing-DLP bundle for the measured balanced
50-digit regime. Nine matched comparisons reduce complete-call time by
**38.10% versus calibrated SLP** and **53.08% versus fixed DLP**, with 9/9
validated completions in every arm. Retain prior 30/40-digit recommendations
and runtime defaults; C10 and E1 remain open.

## Contribution and scope

The [pinned research](c9_research.md) follows candidate recovery, exact
residuals, splitting, relation storage, graph cycles, filtering and extraction
in YAFU, msieve, Yamaquasi, JavaMath and FLINT. The
[constraint audit](c9_growth_audit.md) distinguishes initial allocation from
limits on admission, effort, paths and representation. Full SIQS implementations
grow their relation history; YAFU tiny-QS is a concrete fixed-storage exception.
The comparison records licenses, arithmetic backends, certainty assumptions,
caller branches and paper sections. No upstream code or native constants were
copied, and no upstream native timing predicts PyPy speed.

The implemented `GrowingDoubleLargeSieveConfig` uses PyPy's growable containers,
retains partial history and removes implicit FIFO, atom, row, checkpoint-byte,
total-memory and cumulative split-call quotas from the standalone growing
mode. Explicit quotas remain available. A finite outer portfolio resolves the
child's memory allowance from its own remaining reservation. This adds a
supported storage policy; legacy configurations, CLI defaults and schema 3/4
keep their behavior. Growing SIQS has schema 5 and rejects unsupported SSS and
coarse-worker adapters. The [API guide](../../../README.md#c9-allocator-grown-dlp-storage)
explains usage and recovery.

This is not a claim of globally optimal factoring. Exact endpoint proofs,
per-residual splitting effort and the 256-original-atom cycle policy remain.
The latter rejected no observed training cycle. Finite work/time, configured
family/stall limits and explicit caller caps still stop work. No disk spill,
higher-arity relation representation, matrix redesign, new scoring, polynomial
or multiplier algorithm is introduced. C10 owns broader collector calibration;
C3 owns allocation/handoff; E1 owns combined portfolio acceptance; F6 remains
separate.

## Preserved historical evidence

C1's selected balanced 40-digit bundle reduced confirmed complete-run cost by
**31.73%** (95% interval 29.21–33.99%). Its fresh full-base 50-digit SLP/DLP
comparison completed **0/18 per arm within 120 seconds**, with DLP median
evictions **21,504 / 28,672** on its two inputs. The separate half-base
checkpoint/resume witness completed in **141.063 cumulative seconds**, using
512 MiB, an 8-MiB checkpoint and a 300-second allowance. That changed both
configuration and allowances; it does not isolate duration or retention.

[C9.1](c9_results.md) records the initial policy study and stopped confirmation.
The user's storage-growth direction superseded its implementation assumption;
its partial fresh captures supply no promotion. Its seven inputs became
training/regression evidence for C9.2. Historical source, scope and verdict
remain intact. In particular, the 16K-edge refusal at 256 MiB rejected that
reservation bundle, not adaptive growth.

## Why storage and admission were selected

The retained C1 census and C9.1 replay cover rejected score/refinement domains
as well as survivors. Bernoulli position samples fully divide rejected values;
their 32,768-position coverage ends near 8.4M positions, earlier than the full
65,536-record censuses. This limitation is disclosed rather than treating the
samples as representative of arbitrary later prefixes.

On those fixed record prefixes, high SLP endpoints contributed only 94/99
independent repeat combinations despite 35,017/34,765 records. Under FIFO,
endpoint alignment increased retained LP cycles from 1,136/1,193 to
2,427/2,494 and yielded verified factor-base kernels and proper factors. The
unaligned filtered matrices had zero nullity. Three trivial dependencies
preceded the proper factor on one aligned input. This distinguishes raw
records, LP cycles, matrix survival, kernel dimension and useful extraction.

The ranking was: (1) preserve collision history through growth; (2) test a
common SLP/DLP endpoint ceiling separately at identical scoring/product domain;
(3) combine the better admission with C1's historical half-base bundle and
compare existing reusable SLP matching. A star provides k−1 independent LP
cancellations from k equal partials versus floor(k/2) disjoint pairs, but only
complete factoring can establish economic value. v2 already had the graph-SLP
route; C9 does not claim to have invented it.

Splitter replacement, certainty shortcuts and solve-cadence redesign were
lower-value investments: retained splits had no failures, while profiles
identified sieve/root work as dominant. Forced polynomial-window solves can
still incur repeated preparation; this is a documented future hypothesis,
not evidence for an unmeasured cadence change. Some earlier prefixes had no
surviving matrix before later prefixes completed, so zero short-prefix rows
were never a termination criterion.

## Frozen design and accounting

The [original C9.2 protocol](c9_growth_protocol.md) and
[validation correction](c9_growth_protocol_v2.md) predeclare a **14,400-active-
wall/CPU-second aggregate** cap, including unsuccessful studies, stopped
confirmation, diagnostics, generation and an initial 600 prepaid QA seconds. The pre-confirmation correctness bridge
transfers another 150 seconds from the same aggregate allowance to final QA;
no overall cap increased in that original stage. The later user-requested
1,200-second amendment below has its own declared allowance. Research,
ordinary edits and waiting for the machine lock are excluded. Interrupted or
unpublished failed calls retain their full reserved cost. No study-cap extension
was made for growth.

The measured runtime is PyPy 7.3.23 implementing Python 3.11.15 on an Apple
M4 with 10 logical CPUs and 24 GiB RAM, macOS 26.6.2. Runs use the repository
PyPy environment, not a CPython fallback. Collector arithmetic is the existing
exact-integer path. Main-process CPU reported per call differs from the
inclusive expense ledger, which also charges monitoring/validation subprocess
costs. Process high-water RSS persists across calls in a measurement process;
it is not an independently reset heap peak for each arm.

The original C9.2 50-digit calls receive identical **360-second wall/CPU, 10¹³-work,
1-GiB owned-memory and 3-GiB process high-water RSS** limits. The memory limit
is a matched experimental allowance, not a new production default. Prefixes
are 60/180/360 seconds or earlier terminal outcomes. Independent validation
checks proper divisors, reconstruction, certificate-backed expected factors,
runtime certainty labels and the actual frozen resource envelope. The old
C1 validator's 256-MiB assertion caused one initial C9.2 harness failure; its
unpublished outcome remains unusable and its full 362-second reservation is
charged. The repaired initial grid was repeated in full.

The two old C1 50-digit inputs, seeds 7/29, are attribution/training only.
Individual configurations froze in `2a8eac5`; joint configurations in `9707dab`.
Source hashes and configurations are versioned under `inputs/controls/`.
All timings and heavy verification use the machine-wide lock/owner procedure,
with explicit C3 handoffs. Prefix inspectors build extra adjacency/counter
scratch outside the job-owned reservation, from the finite retained store.
The process RSS guard polls during collection; the terminal inspector follows
the recorded per-call RSS read, so that row alone is not a strict peak bound
for its later inspection scratch. These are instrumented attribution calls,
not accepted timing or total-memory comparisons. Fresh timing omits the
inspector; its job-owned and process RSS measurements remain distinct.

## Individual growth attribution

Full base means numeric bound 30,000, half-interval 65,536, six A factors,
endpoint 3M and product 115.2B. Legacy DLP's SLP ceiling is 900M; alignment
changes it to 3M while preserving the DLP score/product domain. Calibrated
SLP uses 3M with its narrower score domain. Only `grow` versus legacy DLP
isolates storage/count-quota policy at the same collector and admission.

| Arm | Input 0 seconds | Input 1 seconds | Complete | Owned peak MiB, 0 / 1 |
| --- | ---: | ---: | ---: | ---: |
| Current fixed DLP | 100.863 | 112.799 | 2/2 | 193.74 / 194.85 |
| Growing DLP | 54.825 | 64.352 | 2/2 | 595.31 / 590.82 |
| Growing DLP, aligned endpoints | 47.475 | 54.697 | 2/2 | 317.61 / 317.97 |
| Calibrated SLP | 72.088 | 84.721 | 2/2 | 164.88 / 166.12 |

These single-seed, instrumented calls select candidates; they are not accepted
speed estimates. Fixed DLP evicts 78,848/75,776 edges. Growth evicts zero and
retains 37,432/36,894 original atoms, exceeding the old 32,768 atom quota.
Admitted original atoms, graph edges and combined rows are different counts.
Growth's filtered matrices have 1,522/1,560 rows, compared with 1,237/1,329
for fixed DLP. Both complete on both inputs, with fewer collection positions
under growth and more retained memory. On input 0, growth produces 1,234
large-prime combinations from 78,512,727 positions, versus 674 combinations
from 186,779,025 positions for fixed DLP. Both reach proper factors; this is
a useful-yield improvement along the complete trajectory, not just occupancy.

Aligned growth retains 17,615/17,463 original atoms and produces 1,529/1,561
filtered rows. Its extra admission loss reduces storage/preparation while
preserving useful structure on these two inputs. It is a deliberate policy,
not a proof that high endpoints can never help. Split failures, duplicate
rejections and cycle-length rejections are zero throughout this grid.

The growing unaligned input-0 call takes 41.40 seconds in collection versus
80.93 for fixed DLP; aligned growth takes 34.30 seconds. Preparation/filtering/
elimination/extraction are all included. Process high-water RSS ranges from
191 to 279 MiB in the shared initial process and is cumulative, not an isolated
per-arm heap measurement. Owned reservations are conservative and include
simultaneous storage/scratch; they need not equal RSS.

Retained graph development at the terminal observations is recorded below.
The incidence rank includes smooth atoms as zero-LP columns; it is not a
count of DLP-only cycles. Degree counts exclude endpoint 1. Final solver
kernels are diagnostic; each completed call produced one proper factor.

| Arm / input | Distinct LP endpoints | Degree-one endpoints | 2-core edges | Atom-incidence nullity | Last solver kernels |
| --- | ---: | ---: | ---: | ---: | ---: |
| dlp / 0 | 9196 | 7975 | 2254 | 1587 | 1 |
| dlp / 1 | 10349 | 9065 | 2258 | 1610 | 2 |
| grow / 0 | 37258 | 31865 | 3997 | 1620 | 2 |
| grow / 1 | 36679 | 31221 | 4176 | 1654 | 7 |
| grow_aligned / 0 | 17430 | 12021 | 4034 | 1622 | 3 |
| grow_aligned / 1 | 17224 | 11745 | 4188 | 1653 | 5 |

## Combination selection and limitations

| Growing combination | Input 0 seconds | Input 1 seconds | Complete | Capped sum |
| --- | ---: | ---: | ---: | ---: |
| Historical half-base plus aligned growth | 46.425 | 117.317 | 2/2 | 163.742 |
| Reusable graph-SLP at calibrated 3M | 65.611 | 159.555 | 2/2 | 225.166 |

The predeclared completion/capped-cost ranking selects **full-base aligned
DLP growth** (102.172 seconds), followed by full-base unaligned growth
(119.177), half-base aligned growth (163.742), and graph-SLP (225.166).
Calibrated legacy SLP remains the control (156.809), not a post-hoc
per-input mixture. Selection was committed in `c550c5c` before fresh generation.

These are screening ranks, not precise mechanism speed estimates. The two
input-1 combination calls ran in a later process/session with substantially
more CPU seconds per work unit than earlier calls. Cold/JIT, session and
mechanism effects are not separated. No additional training was launched to
seek a preferred ranking. The half-base input-1 arithmetic trajectory retains
9,929 atoms, 914 rows and a 901-by-899 filtered matrix; its lower work and
storage do not establish a faster complete call in this capture. Graph-SLP
scans fewer positions than ordinary SLP on that input (158.2M versus 167.9M)
and produces 907 versus 880 matches, yet takes longer. More matching yield
alone therefore cannot accept it, and these unreplicated timings do not prove
an inherent reusable-matching regression. C10 retains broader calibrated
base/interval and SLP/DLP bundle selection.

## Profiles and charged repeated restore

Two separate 50-second profiles are censored attribution, not accepted timings.
Fixed DLP spends 31.69 cumulative profile seconds in sieving, including 22.98
in prime-power roots; aligned growth has 32.87/23.56. Residual handling is
1.47/1.33, classification 0.63/0.57 and splitting 0.31/0.53 seconds. These are
nested times and must not be added. They support retaining the current
splitter/certainty contract and routing sieve/root optimization to C10.

Aligned full-base growth pauses near 10 and 20 seconds and restores twice
under the same 240-second/10¹³-work allowance. It releases old stores before
rebuilding and independently rechecks original atoms, graph paths and square
corrections. Checkpoints are 378,773/765,860 bytes; replay charges
10,920,247/19,212,103 work units. Post-restore cumulative wall is
10.631/20.854 seconds. The correct proper factor and complete reconstruction
arrive at **100.172 wall / 94.638 main-process CPU seconds**. The separate
1-GiB owned/16-MiB checkpoint envelope reserves simultaneous encoding/rebuild
space; process high-water RSS is 203.52 MiB. This is a correctness diagnostic,
not a speed comparison with the 1-MiB checkpoint timing arms.

The [pre-confirmation bridge](c9_growth_freeze_bridge.md) records a subsequent
uncapped-checkpoint reservation repair and metadata-only training restart fix.
It preserves the selected configuration and original capture identities. AST
review confines source changes to four named checkpoint/training functions;
normal factoring and the explicit checkpoint-size experiment branches are
unchanged. Corrected source froze in `fbfd662` before the seven independently
certified fresh inputs in `f62432a`. The bridge also records the failed full-QA
timeout, conservative charge and final verification reserve transfer.

## Fresh confirmation and acceptance

The original 30/40-digit bands each completed 18 matched samples, with correct
complete factors from every arm. Paired signs and repeated-cell drift remained
unstable after the one permitted extension, so neither receives a new growth
speed promotion. Those inputs comprise one independent balanced 30-digit and
three independent balanced 40-digit integers; repeated seeds do not increase
that count.

The smaller-band totals below include 18 complete calls per arm. Growth does
not reach a relevant legacy storage ceiling in these controls; identical
arithmetic work in DLP and growth confirms that their small timing difference
is not a demonstrated yield improvement.

| Digits | SLP total seconds | Fixed DLP total seconds | Growth total seconds | Growth benefit versus DLP, 95% interval |
| --- | ---: | ---: | ---: | --- |
| 30 | 7.996 | 7.082 | 6.912 | 2.40%, −4.79% to 9.26% |
| 40 | 271.066 | 186.059 | 185.097 | 0.52%, −1.67% to 4.14% |

The 40-digit growth-versus-SLP estimate is 31.72% (24.80–35.66%), consistent
with the pre-existing C1 bundle rather than evidence of a new growth gain.
There were 19/11 remaining stability reasons at 30/40 digits. Retain the
previous C1 recommendation and record these bands as regression controls,
with no additional timing promotion. Owned peaks are 45.04/114.98 MiB for
DLP and growth, versus 42.18/114.99 MiB for SLP. Process high-water RSS is
cumulative and reaches approximately 139/223 MiB in the two bands.

The original 360-second 50-digit study stopped at the user's explicit request
to raise the comparison allowance. Its two completed matched groups remain
preliminary observations, not acceptance. A third group's first call was
interrupted at C3's requested QA handoff; its full 362-second reservation is
charged, and no unpublished result is used.

| Original 50-digit input / seed | SLP wall seconds | Fixed DLP wall seconds | Aligned growth wall seconds |
| --- | ---: | ---: | ---: |
| 0 / 7 | 129.555, complete | 186.715, complete | 90.832, complete |
| 1 / 7 | 360.001, censored | 360.001, censored | 263.777, complete |

On input 1, fixed DLP scans 408,681,517 positions and admits 69,290 atoms, but
60,416 evictions leave only 419 combinations and a final 949-row input matrix
that filters to zero rows. This is a time-censored observation, not a proof
that longer collection cannot complete. The separate proposed 1,200-second
pilot was superseded before execution and consumed no measurement allowance.

The [amended confirmation protocol](c9_long_confirmation_protocol.md) changes
the comparison itself to **1,200 seconds per arm**, preserving the selected
collector configurations and 1-GiB owned/3-GiB RSS/10¹³-work envelope. The
existing production source is unchanged; the new benchmark entry point also
routes cold subprocesses to the amended identity. Source and control froze in
`9fb5aea` before three new independent certified 50-digit inputs in `5c7b686`.
No old confirmation input is reused in this new stage. Expected factors and
certificates remain validation-only. Prime generation conditions p−1 on a
large known factor; this is not an RSA sample population.

The original study is charged **8,917.431 seconds**, including unsuccessful
and unpublished reservations and its prepaid QA. The amended stage has a
separate finite **14,400-second** allowance, including 300 prepaid QA/assessment
seconds; the combined declared bound is **23,317.431 seconds**. This records the
user-directed change explicitly without hiding previous search cost. No
further parameter search or automatic study follows. The new stage finishes
at **6,607.145 charged wall / 6,586.598 charged CPU seconds**, including cold
calls, repeated warmup, generation and the prepaid QA/assessment allowance.
Together the two stages charge **15,524.576 seconds**, below the amended bound.
The new bootstrap analysis takes 0.732 wall / 0.664 CPU seconds; final source,
corpus and acceptance consistency checks take 0.233 wall / 0.167 CPU seconds,
within that prepaid allowance. Machine-window elapsed time also includes
ordinary editing/waiting and is not substituted for the active expense ledger.

### Completed 1,200-second comparison

Every arm completes **9/9**, produces **9 proper factors**, and independently
passes reconstruction, certainty-label and resource validation. Each of three
fresh inputs has seeds 7/29/47; no partial group enters acceptance. The frozen
stability check has no reasons for extension. All nine paired directions favor
growth against both controls. The hierarchical input/paired bootstrap gives:

- **38.10%** complete-call reduction versus calibrated SLP, 95% interval
  **35.23–39.81%**.
- **53.08%** reduction versus fixed DLP, 95% interval **52.36–53.84%**.

Completion regression is zero percentage points, within the five-point gate.
There is no minimum benefit threshold beyond the frozen uncertainty gate.
These intervals describe this small, generation-conditioned population; they
are not a guarantee for arbitrary 50-digit integers. The selected bundle
combines storage growth and endpoint alignment. Only the earlier individual
training ablation isolates storage policy alone.

| Arm | Complete / proper factors | Total wall = capped seconds | Median wall seconds | Inclusive CPU seconds | Charged work, billions | Peak owned MiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Calibrated SLP | 9 / 9 | 1,931.013 | 202.524 | 1,926.086 | 449.884 | 166.89 |
| Fixed DLP | 9 / 9 | 2,547.200 | 273.486 | 2,541.613 | 1,824.300 | 198.15 |
| Aligned growing DLP | 9 / 9 | 1,195.263 | 127.464 | 1,192.632 | 901.150 | 325.99 |

Inclusive CPU includes the monitoring children recorded by the expense ledger.
Raw main-process CPU totals are 1,828.672 / 2,413.262 / 1,132.141 seconds in
that order. Per-call expense wall including harness overhead totals
1,931.707 / 2,547.927 / 1,195.924 seconds, so the small reporting overhead does
not drive the result. All arms have the same 1-GiB owned allowance. Growth
uses more peak owned memory and about twice SLP's charged work units, while
reducing CPU/wall cost. Work units are algorithmic charges, not CPU cycles;
this result does not establish dominance in every resource dimension.

Process high-water RSS reaches **239.27 MiB** for each arm in the shared
measurement processes. This cumulative statistic is not an isolated per-arm
heap measurement. The first matched group runs in process 42333; the remaining
eight run in 46024. Each process performs at least three seconds of both
validated wall and CPU warmup per arm. No matched group crosses that boundary.

| Fresh input | Mean SLP seconds | Mean fixed-DLP seconds | Mean growth seconds | Growth benefit vs SLP / DLP |
| --- | ---: | ---: | ---: | ---: |
| 0, three seeds | 199.629 | 272.128 | 128.661 | 35.55% / 52.72% |
| 1, three seeds | 253.185 | 325.463 | 152.351 | 39.83% / 53.19% |
| 2, three seeds | 190.857 | 251.476 | 117.409 | 38.48% / 53.31% |

All fresh calls finish in 114.248–334.279 seconds, before the common 1,200-second
cap. This new cohort does not resolve the old 360-second censors or isolate
the effect of extending duration. It also does not establish that 360 seconds
is sufficient for other inputs. The superseded selected-censor pilot remains
unrun.

### Useful yield and downstream cost

The following are per-call medians. Admitted atoms are cumulative admissions,
not simultaneous retained-edge occupancy. LP combinations, original atoms,
filtered rows and factor-base kernels remain different objects.

| Diagnostic | SLP | Fixed DLP | Aligned growth |
| --- | ---: | ---: | ---: |
| Sieve positions | 184,944,003 | 233,441,013 | 100,139,772 |
| Candidate reports | 11,590,231 | 15,055,434 | 6,457,885 |
| Refinement rejections | 11,556,701 | 14,558,599 | 6,245,512 |
| Admitted atoms | 11,177 | 88,504 | 17,688 |
| LP combinations | 913 | 738 | 1,259 |
| Evictions | 445 | 78,848 | 0 |
| Split attempts | 0 | 27,592 | 11,858 |
| Rho evaluations | 0 | 15,172,034 | 6,548,004 |
| Endpoint rejections | — | 201 | 83 |
| Product rejections | — | 274,294 | 116,335 |
| Prime-residual rejections | — | 116,584 | 69,990 |
| Graph rebuilds | — | 77 | 0 |
| Final filtered rows / columns | 1,346 / 1,344 | 1,330 / 1,330 | 1,558 / 1,552 |
| Singleton removals | 106 | 116 | 33 |
| Weight-two merges | 174 | 169 | 55 |

Across all nine calls, fixed DLP has **one failed split in 249,354 attempts**;
growth has zero in 106,686. Thus the earlier zero-failure training observation
is not silently generalized to all confirmation. No duplicate, split-quota
or long-cycle rejection occurs in these calls. Fixed DLP evicts 700,416 edges
in total; growth evicts none. The instrumented training tables above supply
retained endpoint/degree/core and incidence-rank observations; accepted fresh
timing does not rebuild adjacency merely to measure it.

The final solver returns 1–4 / 1–7 / 1–14 kernel vectors for SLP / DLP /
growth. Those are available vectors, not a count of successful extractions.
There are **11 / 97 / 13 trivial dependency trials**, zero tested-cache skips,
and nine factor-producing extractions in each arm. Growth therefore does not
improve every dependency-quality statistic versus SLP. It reaches proper
factors sooner through the whole pipeline, with larger surviving matrices.

| Total stage wall seconds, nine calls | SLP | Fixed DLP | Aligned growth |
| --- | ---: | ---: | ---: |
| Collection | 1,547.104 | 2,064.531 | 842.601 |
| Preparation and verification | 116.668 | 146.388 | 157.662 |
| Filtering | 109.400 | 132.263 | 74.168 |
| Elimination | 68.934 | 81.713 | 70.100 |
| Extraction | 0.110 | 0.175 | 0.150 |

Growth increases total preparation cost, but its collection saving survives
all later work. Complete-call totals also include setup, switching, monitoring
and final output validation outside these engine timers. Filter/solve invocation
counts are 8,663/8,654 for SLP, 9,586/9,577 for fixed DLP, and 5,411/5,402 for
growth; many invocations have no useful kernel. They are not independent
useful-dependency counts. Solve-cadence optimization remains an untested future
hypothesis, not another change in this tranche.

One separate cold call on input 0 per arm completes in **184.237 / 275.436 /
131.831 seconds**, with whole-process startup-inclusive times **184.848 /
276.031 / 132.497 seconds**. Cold main-process CPU is 173.760 / 258.481 /
125.009 seconds. These single observations are diagnostics, not extra warmed
samples or a statistically accepted cold-start benefit. Profiles and charged
checkpoint/repeated-restore evidence remain in their separate sections above.

### Adopt / retain / defer

**Adopt** the explicit full-base aligned growing-DLP configuration represented
by the 50-digit `challenger` in `inputs/controls/c9_long_selected.json`, together
with the tested growing-storage API and schema-five recovery contract. The
fresh complete-factor benefit supports this measured balanced 50-digit scope
under the stated resource envelope. It is a selected bundle, not an optimized
universal digit rule or an automatic dispatch/default change.

**Retain** the historical C1 40-digit recommendation, existing 30-digit/default
behavior, legacy finite configurations and their checkpoint compatibility.
The smaller-band studies do not confirm a new growth speed gain. Preserve all
stopped and censored results with their original sources, scope and costs.

**Defer** broader size/distribution calibration, 60–100-digit capability claims,
SSS/coarse-worker growing adapters, opaque cofactors, new splitters, disk spill,
cycle-length and solve-cadence changes. C10 owns broader collector strategies;
C3 integration and E1 own the combined portfolio decision. Higher arity remains
F6 and matrix/provenance redesign retains its separate gates. Only **C9.2's
bounded source-derived yield and adaptive-storage tranche** closes here.

The final consistency check revalidates all 27 saved outputs, balanced seeds,
positive intervals, complete sample count, clean stability decision, corpus
certificates and all **209 frozen source digests**. Required control/corpus
files are versioned. Their SHA-256 values are:

- `controls/c9_long_selected.json`:
  `751d1a767f87db4e15ad4486c1c0334f60563d784dbd98756517d8414a4bd57c`.
- `corpora/c9_long_fresh_corpus.json`:
  `7a4fdc8484d9d1c7ef7c1341b4dfc9fd4ec2f2e798c22f9ba7e096b1a18737d1`.

Generated captures remain ignored under `results/c9/long-confirmation/`.
Its final capture/summary set occupies 217,862 bytes; the C9 results tree is
about 6.1 MB, excluding pre-existing C1 evidence and upstream source archives.
The capture and simultaneous QA-archive reservation checks remain below
256 MiB. No additional timing study or heavy QA follows this decision.

## Verification and reproduction

At committed revision `fcc90abe694a1514404457caa4bf7c471b7ebb08`,
`make -C v2 test` passes **594 tests**, and `make -C v2 lint` passes. A clean
`git archive` containing only committed files independently passes the same
594 tests. All **125 benchmark modules** import, and required snapshot,
A7/B4/C6/C1 and current C9 corpus/control loaders pass. Installed PyPy and
third-party dependencies are supplied externally; the archive contains no
working-tree `.venv`, ignored source or generated input dependency. Combined
final verification costs 198.327 wall / 208.271 CPU seconds, charged against
the prepaid 750-second QA reservation. The source freeze is unchanged.

The amended harness and newly committed corpus at `5c7b686` independently
pass **596 tests in both the working tree and committed-only archive**. Full
lint passes; the archive also passes **126 benchmark imports plus required
loaders**. The two additional
tests check amended cold-worker routing and restoration of the original study
identity after interruption. This verification costs 196.843 wall / 206.613
CPU seconds within the new 300-second prepaid allowance. Production source
hashes remain unchanged. Final edits after that verification change only the
seven Markdown documentation files; no implementation, runner, control or
corpus bytes change.

Graph tests compare generated cycles with an independent dense incidence
kernel, including self-loops, parallel edges and disconnected components.
Arithmetic tests reconstruct original atoms and square corrections. Growth
and explicit-FIFO cases retain referenced paths across rebuilding and repeated
restore. Fault cases cover corruption, pre-publication cap/cancellation
refusal, interrupted multi-container publication, dynamic checkpoint coexistence
and repeated restore without leaking the transient reservation.

From the repository root, with PyPy implementing Python 3.11 and an explicitly
coordinated machine lease, reproduce the original 360-second study at
**revision `c0d4924`** with the commands below. The current amended 50-digit
runner is `v2.benchmarks.qs.c9.c9_long_confirmation`; use a separate output
directory and the [amended protocol commands](c9_long_confirmation_protocol.md#reproduction).
Historical strict source checks are intentionally preserved.

Original commands:

```sh
v2/.venv/bin/python -B -m v2.benchmarks.qs.c9.c9_growth_confirmation confirm \
  --output "$PWD/v2/benchmarks/results/c9/reproduction" \
  --digits 30 --lease-seconds 600
v2/.venv/bin/python -B -m v2.benchmarks.qs.c9.c9_growth_confirmation confirm \
  --output "$PWD/v2/benchmarks/results/c9/reproduction" \
  --digits 40 --lease-seconds 1800
v2/.venv/bin/python -B -m v2.benchmarks.qs.c9.c9_growth_confirmation confirm \
  --output "$PWD/v2/benchmarks/results/c9/reproduction" \
  --digits 50 --lease-seconds 3600
v2/.venv/bin/python -B -m v2.benchmarks.qs.c9.c9_growth_confirmation summarize \
  --output "$PWD/v2/benchmarks/results/c9/reproduction"
```

Resume the same output directory after an agreed lease yield; source and
manifest checks reject retuning or corpus replacement. Each new process
performs fresh validated warmup. Do not regenerate the now-public corpus or
call it new confirmation. Use the retained source revisions for historical
training; the current confirmation loader intentionally requires the corrected
pre-confirmation freeze. Leases longer than the original default were explicitly
coordinated with C3; they change neither call limits nor the aggregate study
allowance. The frozen runner admits individual calls, so a matched group may
cross a process boundary; process IDs are retained for that limitation.

## Integration boundaries

C3 separately integrated its schema-12 implicit portfolio decoder at mainline
`d273cb8`. Integrating these branches must discriminate
`storage_policy == "grow"` before legacy DLP/SLP construction,
and test implicit allocation plus growing resume with the finite resolved child
allowance. C9's current portfolio route requires its same explicit configuration
on resume. C3/C9 also have separate current-source A7 manifests; reconcile their
pins explicitly when integrating. No merge or source retuning occurred during
these studies. No active C10 chat was present during the initial ownership
checks. C10 later started its separate worktree pinned at `129cb09` (the same frozen C9 runtime),
and was explicitly given ownership of A-product ordering and its metadata
reservations, scoring/power-plan and polynomial investigations. Growing-mode
`None` quotas and live/checkpoint/replay ownership must survive those changes.
The settled C9 decision, exact configuration and fresh-corpus exclusion hashes
were handed to C10 with explicit CPU release after confirmation. C10 proceeds
with its separate source and integrated measurement protocol.
None of C3's subsequent changes were merged into this frozen study worktree.
C9 changes no scoring, polynomial or candidate-recovery algorithm.

## Subsequent literature verification

The [11 October transfer verification](c9_literature_verification.md) checks
current source revisions and records implemented, evaluated, deferred and
previously unmapped mechanisms. Asymmetric ordered LP bounds were missing
from the earlier explicit comparison; recovered-product admission, target
feedback and cycle-weight optimization remain untested. Only the named C9.2
tranche is accepted. This addendum changes no frozen source/configuration or
historical result and does not assert exhaustive adaptation or optimality.


## Subsequent mainline integration

The [11 October integration record](c9_integration.md) reconciles C3
allocation/defaults with growing DLP and phases remaining literature transfers
into roadmap C9a. This changes no frozen control, corpus or timing verdict;
reproduce accepted timings from their historical source revisions.


## Subsequent default promotion

The user subsequently requested [50-digit default promotion](c9_default_promotion.md).
This applies the confirmed bundle and finite omitted grant without changing
this experiment's original defaults decision, observations or uncertainty.
