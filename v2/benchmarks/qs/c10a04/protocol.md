# C10a-04 frozen experiment protocol

Stage 1 owns research, exact helpers, contracts, limited diagnosis and this
freeze. Stage 2 executes the commands in [handoff](handoff.md) without changing
sources, arms, inputs or decision rules. Neither stage implicitly promotes
production defaults or closes E1. Historical studies remain immutable.

## Configurations and population

[Research](research.md) and the [source manifest](../../inputs/controls/c10a04_sources.json)
define the source and representation assumptions. `c10a04_frozen.json` binds
runtime, all live benchmark helpers, this protocol/research, source pins,
training corpus and policy templates. Every lease refuses a source mismatch.
The frozen commit is named in the handoff; doc-only handoff/result annotations
are outside the hashed experiment sources.

| Digits | Control | Compact | Wide | Matched wall and CPU grant per call |
| --- | --- | --- | --- | ---: |
| 55 | k=1, bound 30000, M=65536, r=6, pool=32 | existing odd k, 2500 eligible primes, M=65536, r=7, pool=64 | compact with M=131072 and its coupled A target | 1200 s |
| 60 | same control | existing odd k, 5000 eligible primes, M=65536, r=8, pool=64 | compact with M=131072 and its coupled A target | 7200 s |
| 65 | same control | existing odd k, 8500 eligible primes, M=65536, r=9, pool=64 | compact with M=262144 and its coupled A target | 21600 s |

Control is the confirmed 50/51-digit collector bundle explicitly applied to
larger N, not the unrelated tiny automatic default currently used at 55+.
Every arm holds accepted powers scoring, 2 MiB plan, DLP/splitter/graph policy,
absolute LP endpoint/product bounds, collection pacing, row growth/excess 32,
100000-family allowance and 1 GiB owned memory fixed. All use existing flyer
assignments; configured Gray quota 64 means effective 32 in control and 64
in challengers. No new population or overlap policy. Count resolution scans
primes below 250000, charges exact computation/storage, and reports actual
endpoint+1, count, k, mod-8 class, singular primes and exact configuration.
Helper refusal is a refusal under its finite envelope, never silent resizing.

Three independent certified balanced products per band are versioned before
training. Revealed 55/60 pilot products are training only. Seeds 7,29,47 crossed
with these inputs yield nine matched observations per band; seeds do not add
independent inputs. The generated primes are conditioned on a large p−1
factor, not an unconditioned RSA population. This is a conditional population
and machine/backend claim only. Each group uses fixture `index%3` and seed
`(index//3)%3`, cycling exactly on the 18-observation extension. Arm order uses
`Random(1042026 + digits*100 + index).shuffle` for both phases; it is frozen,
independent of outcomes. Raw successful/unsuccessful calls are all retained.

## Timing and complete outcomes

Use PyPy implementing Python 3.11; there is no CPython fallback. Own
`/private/tmp/factor-performance.lock` for screens, timing and heavy QA.
Cold worker calls capture startup separately. Every measurement process
warms every arm for at least three validated wall **and** CPU seconds before
accepted samples, even after a short lease. Warmup is bounded by 64 calls
of 3.2 s each per arm/lease and charged to its phase. Warmup outcomes need
valid reconstruction/certainty, not a factor from each short prefix. Never
combine cold starts, warmups, profiles or Stage 1 screens with performance
observations. Inclusive cost includes count/multiplier resolution, setup,
collection, splitting, graph, filtering, matrix, extraction, classification
and independent final divisor/reconstruction validation.

Every call has 10^13 cumulative work, its arm's 1 GiB owned allowance, a 1 MiB
checkpoint allowance, a 3 GiB process RSS watchdog and the band's matched
wall/CPU limits. RSS is a process high-water mark, not isolated arm memory.
Deadlines remain cooperative; integer operations/chunks are finitely bounded
and overshoot is charged. Record exact proper divisors, unresolved cofactors,
certainty labels, work, raw/filtered rows, dependencies, stage times, owned
reservation and RSS. Search exhaustion, stalled collection, storage refusal,
work exhaustion and time exhaustion remain distinct reasons. A dependency
or an incomplete row count is never a complete factor.

## Finite resources, leases and extensions

The following independent wall and CPU caps reserve the complete worst-case
planned campaign, not a target duration or a universal size/time limit:

| Phase | Wall cap | CPU cap | Aggregate work cap |
| --- | ---: | ---: | ---: |
| Conservative prepaid Stage 1, including all QA/export checks | 3600 s | 3600 s | 10^14 |
| Training, including startup/warmup/instability repeats | 1630000 s | 1630000 s | 12×10^15 |
| Fresh generation/confirmation and repeats | 1090000 s | 1090000 s | 8×10^15 |
| Separate larger-time-grant diagnostics | 200000 s | 200000 s | 2×10^15 |
| Entire study | 2950000 s | 2950000 s | 3×10^16 |

Capture disk is at most 256 MiB in the experiment tree. No phase borrows a
reserved training/confirmation share. The training/confirmation caps cover
18 groups × three/two arms × sum(1200,7200,21600), plus setup reservations;
normal operation begins with nine groups. At most 100 process leases per
phase are permitted; repeated warmup overhead otherwise can consume the
headroom. Whole missing groups must fit both lease and phase/aggregate
allowances before starting. Default lease is 90000 s, long enough for a
three-arm 65-digit group; short leases yield with unique immutable markers.
An exhausted allowance returns an unresolved study under that allowance.

Prepaid Stage 1 covers charged computation for the limited screen, completion
pilot, corpus generation, source verification, smoke checks and heavy QA.
The screen had separate finite 1800-second/10^15-work ceilings; these were
limits, not its consumed charges. Actual screen work is separately recorded;
the conservative Stage 1 work reserve covers measured consumption and QA. No accepted sample starts in
Stage 1. Existing screen charges/failed attempts are retained locally even
if implementation fixes change the final runner's identity.

`grants --digits D` is optional **diagnosis only** after the frozen comparison:
one revealed training number, seed 7, all three arms at double that band's
wall/CPU grant (2400/14400/43200 s), same work/storage and source identity.
It has its own phase and filenames; it cannot repair a failed acceptance cell,
enter selection, or establish a resource-independent configuration speedup.
All three missing calls must fit its phase. No further extension, changed
base/LP domain, new input band, different arm or fresh retry is authorized by
this freeze. A later separately frozen tranche can exceed 65 digits or these
grants; these are not universal maxima. Preserve every prior outcome/charge.

## Interrupted execution and durable receipts

Both shared performance and per-output advisory locks serialize mutations.
Before each operation, publish/fsync a full time/work reservation. Publish
its immutable JSON receipt before settling work. Time reservations are
nonrefundable even on success, conservatively covering receipt/orchestration
overhead; actual wall/CPU (including child CPU) are also retained separately.
Cooperative overshoot increases the charged reservation. Repeated identical receipt publication is harmless;
identity/content conflicts, tampering or missing settled receipts refuse.
Restart in the **same output directory**. A crash after receipt publication
settles from that receipt once. A lost measurement before publication becomes
a validated reconstructible **censored** cell with the full reservation;
resume proceeds to other missing cells and does not selectively retry it.
Failures from correctness/contract exceptions block further measurement for
an audit. Lost diagnostic/startup/warmup attempts retain full charges; any
retry is separately charged. No clock/work or censored result is reset.

Study continuation is at receipt/group boundaries. It does not save a costly
mid-call SIQS snapshot for every block or claim process-crash recovery of
unpublished atoms. Existing SIQS checkpoints remain compatible and charged;
tests round-trip the new joint configuration and check replay. No checkpoint
schema or runtime arithmetic is deferred to Stage 2.

## Selection, uncertainty and fresh confirmation

After nine groups, extend **all arms in the band** to eighteen if any input's
paired capped-cost directions mix, early/late three-group mean effects differ
by over 0.15, or repeated-cell wall drift exceeds 15%. These are the existing
explicit instability diagnostics, not a 15% improvement floor. After eighteen,
remaining instability prevents challenger selection/acceptance; report it.
Do not hide a censored group because other calls succeeded.

Rank training arms by completed cells, then aggregate capped cost, tie to
control and then frozen arm order. Complete cells cost actual inclusive time;
unresolved cells cost max(actual time, whole cap); lost cells retain their
reserved cap/overhead. Select control if instability persists. Hash all training
receipts in a versioned selection and **commit it before generating fresh
inputs**. The generator enforces that committed identity. Exclude every integer
`n` in versioned inputs/controls/corpora, including historical exclusions and
all revealed training/screen products. Generate once with seed 202610120402;
retain certified factors and selection identity, then commit the fresh corpus
before confirmation. No inspected fresh number returns to fresh sampling.
Generation has finite 115-second / one-million-random-call / 10000-product
attempt / 4096-certificate bounds and a charged 120-second reservation; failed attempts stay charged.

Fresh confirmation compares only control and the selected challenger with
three new independent products crossed with the same three seeds, extending
to eighteen under the same rules and the **same** per-band matched envelope.
Report 4000-draw input-then-observation bootstrap 95% intervals for capped
cost and completion-rate differences. With equal full completion, acceptance
needs the entire complete-cost interval below zero and no remaining instability.
For a coverage gain, every fresh challenger call must complete, no input may
lose a completion, and the completion-gain interval's lower endpoint must
exceed zero. A reference's fast unsuccessful stop is reported separately;
a capped-cost ratio against it is never labeled a factoring speedup. Zero
correctness failures and complete planned samples are mandatory. No arbitrary
10% gain floor. Three input clusters remain a small conditional sample.

No winner, unresolved cases, work/storage refusals or unstable results retain
the control under this envelope and are publishable outcomes. Report each
band separately; do not extrapolate between/beyond bands or compound earlier
percentage gains. Promotion, if requested later, requires its own documented
scope and regression/portfolio gates. E1 remains deferred.
