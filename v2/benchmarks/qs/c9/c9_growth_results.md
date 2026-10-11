# C9.2 adaptive storage and useful large-prime yield

Status: implementation, training and diagnostics complete. Fresh confirmation
is in progress; no new speed promotion is claimed at this checkpoint.

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
no overall cap increases. Research,
ordinary edits and waiting for the machine lock are excluded. Interrupted or
unpublished failed calls retain their full reserved cost. No study-cap extension
was made for growth.

All new 50-digit calls receive identical **360-second wall/CPU, 10¹³-work,
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

The 30/40-digit bands each completed 18 matched samples, with correct complete
factors from every arm. Paired signs and repeated-cell drift remained unstable
after the one permitted extension, so neither band receives a new growth
speed promotion. The 50-digit comparison is still running. Fresh inputs are
three independent balanced 50-digit, three balanced 40-digit and one 30-digit
regression integer, disjoint from versioned corpora. Their certified-prime
construction conditions p−1 on a large known factor; they are not RSA samples.
Factors and certificates are used only for validation. Seeds/repetitions do
not increase the independent-input count. The 30/40-digit challenger changes
storage on C1's existing half-base DLP bundle; the 50-digit challenger is the
selected full-base aligned growth configuration. All arms receive the same
1-GiB owned memory, 3-GiB RSS, 10¹³ work and per-band 5/30/360-second caps.

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

Graph tests compare generated cycles with an independent dense incidence
kernel, including self-loops, parallel edges and disconnected components.
Arithmetic tests reconstruct original atoms and square corrections. Growth
and explicit-FIFO cases retain referenced paths across rebuilding and repeated
restore. Fault cases cover corruption, pre-publication cap/cancellation
refusal, interrupted multi-container publication, dynamic checkpoint coexistence
and repeated restore without leaking the transient reservation.

From the repository root, with PyPy implementing Python 3.11 and an explicitly
coordinated machine lease, reproduce the frozen confirmation with:

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

C3 is developing its separate schema-12 implicit portfolio decoder. Integration
must discriminate `storage_policy == "grow"` before legacy DLP/SLP construction,
and test implicit allocation plus growing resume with the finite resolved child
allowance. C9's current portfolio route requires its same explicit configuration
on resume. C3/C9 also have separate current-source A7 manifests; reconcile their
pins explicitly when integrating. No merge or source retuning occurred during
these studies. No active C10 chat was present during ownership checks, and C9 changes no scoring,
polynomial or candidate-recovery algorithm.
