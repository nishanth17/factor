# C10a-01 v4 lower-band supplement

Status: design only, 2 November 2026. This supplements protocol_v4.md and
experiment_v4.json without changing their freezes, cells or grants. The active
plan consists of both documents. No measurement, generation or complete-factor
pilot is run by the checker. Use Astra xhigh and the parent's source/runtime,
serialization, warmup, correctness, censoring and evidence rules throughout.

## Question and order

Does the selected arithmetic bundle reduce complete factoring cost at 50/55,
and where do its setup costs or interactions with existing marking outweigh
its benefit? The existing 60 comparison remains the upper anchor. Correctness
regressions alone cannot answer this performance question.

Run the parent's clean calibration and bounded component diagnostics first;
freeze its single nominee before revealing new boundary inputs or timing this
supplement. Then run the lower screen before the later A/B/C pilots where
practical. Its results do not gate an independently valid 65/70 comparison.
Do not add a third full-bundle arm, a per-input selector or another ablation
grid. A lower-band loss may restrict future eligibility without invalidating
a larger-band hypothesis. A changed nominee requires a separate prospective
protocol, not selective repetition of the losing inputs.

## Inputs, configurations and work

The complete screen reuses a1a_integration_fresh.json cases 0, 1 and 4 at each
of 50 and 55 digits: two balanced and one uneven modulus per band. These six
moduli are already exposed despite the historical filename. Selection is by
fixed IDs and balance, not their earlier outcomes. They provide exploratory
cost/completion observations, not fresh confirmation. Hashes, N, residues and
achieved factor bits are in experiment_lower_v4.json. Validate retained prime
certificates and exact products before launch; pass only N and seed 7 to jobs.

Compare the owned 0076e25 current-mainline control with the globally selected
nominee. Freeze the lower literal configurations before timing; JSON currently
contains the original batch/windowed_pair/parity2_norm nominee template, to
be replaced only by the parent's predeclared nomination release. Use 50's
actual bounded base policy (B=30000, six A factors, pool 32, control k=1),
55's counted policy (2500 primes, seven A factors, pool 64, scored control k),
and the corresponding 51/52 and 55/56 policies. Never use the 60-digit factory
as a lower-band performance control.

Retain density+presieve+scan_fused at 50 through 55. Append parity_two in the
nominee; its initializer deliberately falls back when presieving already
supplies those levels. At 56 retain sparse control marking and parity_two for
the nominee. Existing queue, power plan, pacing, LP domains and other settings
stay literal and equal within each pair. The nominee deliberately enables
coupled automatic multiplier selection at 50/51, where control uses k=1;
report selected k, discriminant residue and fallback. This measures the whole
bundle, not a representation-only effect.

Both arms receive Tier 1: 1 GiB total owned, 3 GiB RSS. Reuse the parent's
60-band external reserve and 32 MiB owned / 256 MiB RSS parent suballocation.
Set finite row/atom/partial caps 12000/65536/32768 and split-call cap 1048576
in both arms, instead of the automatic factory's unbounded growing quotas.
These are declared operational caps; report any binding refusal. Construct
the nominal 1 GiB default policy before subtracting reserves, or use the frozen
literal mapping. Calling the factory with the reduced grant would silently
turn off calibrated policies and is forbidden.

The upper accounting envelope is the parent's 60-band bound: F <=5000,
r <=8, N <=60 digits, same row/provenance/matrix capacities and family limit.
At B=30000, the total prime count is 3245, so the 50/51 base is below that
ceiling. Smaller bases do not increase the dense-fill upper bound. Add at most
32 MiB for simultaneous presieve/scan scratch using collector_kernels.workspace
at block width 4096 and F=5000; periodic lists and scan coexistence are included.
The conditional central ownership scenario is therefore below approximately
0.801 GiB, leaving about 0.199 GiB; it is not a physical RAM or completion
claim. Worst capacity populations may still refuse, as at 60. New RSS is
unmeasured. Check actual host headroom and finite configuration admission;
never weaken coefficients or enlarge the tier after a refusal.

## Short routing-boundary checks

Run one matched uninstrumented prefix at each of 51, 52, 55 and 56 digits,
128 polynomials / 16777344 positions / 4224 blocks per call. This is eight
calls, not eight independent moduli. Reuse c10a05_plus_feas_51 and the first
55 screen input. The 55 prefix is a repeated input, not additional independent
evidence. Log actual family counts: six-factor families have only 32 Gray
polynomials, so 128 polynomials need four families at 51 and two at 52/55/56.
Early factors and censored prefixes do not count as equal-position timings.
Observers remain off; reuse the parent's observer passivity checks.

No ordinary 52/56 semiprimes are present in the retained selected corpora;
prime powers found there are unsuitable. Generate exactly one balanced
semiprime at each size, with the recorded seeds, residues 5/7, distinct ordinary
factors and no factor reuse across this supplement. Use the existing exact
certified sampler and independent certificate validation, with at most 10000
candidate pairs per slot, 1000000 random draws total, and 120 wall AND CPU
seconds including proof validation of reused inputs. Support of requested
bit-size/residue cells must be checked before generation. Accept the first
valid result; no regeneration based on factoring. Generation failure stops
that boundary cell, not unrelated screens. Freeze resulting numerical inputs
and proofs in a new execution release before any measured call. These two
pending slots are an explicit readiness gap, not invented numerical inputs.

Also perform untimed routing assertions at exact decimal thresholds around
50/51/52/55/56/60/61. Such synthetic routing values are not performance inputs.
Boundary prefixes expose setup and marking discontinuities; they do not prove
complete-factor gains at every intermediate digit size.

## Additional prospective allocation

All limits below are wall AND CPU; each measured call has a five-second margin.

| Phase | Calls | Per-call cap | Reservation |
| --- | ---: | ---: | ---: |
| 50-digit complete screen | 6 | 120 s | 750 s |
| 55-digit complete screen | 6 | 300 s | 1830 s |
| 51/52/55/56 boundary prefixes | 8 | 30 s | 280 s |
| Warmup and drift | At most 20 contexts | 25 s | 500 s |
| Input generation and proof validation | Aggregate | 120 s | 120 s |
| Analysis and validation | Aggregate | 120 s | 120 s |
| Orchestration | Aggregate | 100 s | 100 s |
| Total | 20 measured calls maximum | | 3700 s |

Retained C10a-07 combined observations completed 50-digit inputs in about
20–39 seconds and 55-digit inputs in about 73–100 seconds. Those are different
inputs/source identities; they motivate 120/300-second finite caps, not a
prediction or a scaling law. Initial total execution estimate is 15–30 minutes,
with a worst-case reservation of 61 minutes 40 seconds. Setup and later-stage
cost are included. Reassess an uninformative screen; no automatic larger retry.

For each call use the parent's exact call_work(60, seconds) numerical guard:
four times the retained 8e9 sustained rate, separately bounded setup, whole
schedule root/A charges, and dense-fill postprocessing bounds at the larger
60-band capacities. The presieve allowance is unchanged from that accounting's
supported collector work; record actual rates and refresh admission if they
exceed the parent's bootstrap. Warmup reserves 20 * 32 * call_work(60, 20).
Generation/validation and analysis each reserve 1e12 work; orchestration adds
no algorithmic work. All exact integers and cumulative sums are in JSON. Keep
the same three-second validated wall AND CPU warmup, pre/post drift sentinels,
finite invocations, lock, watchdog and interrupted-reservation rules as v4.

New output including archive is limited to 8 MiB, in addition to the parent's
64 MiB prospective output, still within the global 512 MiB ceiling. Recount
before admission. Existing bounded receipts apply; no new verbose recorder.

The parent 38600-second whole envelope plus 3700 equals 42300 seconds:
11 h 45 min wall AND CPU. With original reserved consumption 3124 instead of
its 3600 ceiling, the corresponding charged maximum is 41824 seconds.
The 60/62, 65 and 70 grants stay 9200, 7600 and 14800. No grant is transferred
from them. Confirmation is still separate and unfunded.

## Interpretation and eligibility

Report every complete and censored pair, complete CPU/wall costs, stage costs,
selected k/form, useful relations, memory and refusal. Counterbalance arm order
across inputs and reverse the first arm across sizes. Ratios are descriptive
only when both finish; capped service cost assigns the cap to failures, with
actual consumption also shown. Three independent moduli per lower band and
one 60-digit pilot pair cannot establish general gains or absence of harm.
No pooled 50–70 average may hide a lower-band regression.

Do not tune thresholds or reselect the nominee from these results. Label each
band promising, mixed, censored or unfavorable with its actual evidence;
retain current default eligibility pending confirmation. Large-size feasibility
may proceed regardless of a lower-size loss. Freeze a separately justified
size-scoped eligibility proposal before revealing confirmation inputs.

A claim of gains across 50–60 requires fresh paired complete-factor evidence
at 50, 55 AND 60, and direct evidence at changed routing bands or a narrowed
claim that does not extrapolate. At least nine independent moduli per claimed
band is a floor, not adequate power by itself. Follow the parent's pre-reveal
sampling, drift, input-level uncertainty, multiplicity and completion-harm
rules; repeated seeds or the extra 55 prefix are not new inputs. The existing
60 pilot remains exploratory and must be expanded in that separately funded
confirmation design. No default promotion follows this screen.

## Verification and remaining work

    pypy3 -m v2.benchmarks.qs.c10a01.lower_plan_v4 check

This checks both numerical plans, source pins, literal configurations and input
identities without running a factorer. The parent worker/ledger adapters remain
required, with lower SIQSConfig versus CountedSIQSConfig dispatch, composed
presieve/parity marking, shorter six-factor families, generation refusal and
global nominee substitution added to their correctness obligations. Before
execution independently validate roots, valuations, reconstruction, charged
resume and observation passivity for this lower-band composition. Configuration
construction alone does not discharge those arithmetic/integration obligations.
