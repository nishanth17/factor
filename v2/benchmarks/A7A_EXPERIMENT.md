# A7a execution specification: joint configuration search by regime

Status: committed execution design, not a record of a completed campaign or
authorization to launch timings. Do not substitute a succession of 60-digit
studies for this deliverable. Apply the current
[paired-sampling policy](../POLICY.md#cost-scaled-paired-sampling).

## Deliverable and implementation checkpoint

Deliver a table for **40, 50, 55, 60, 62, 65 and 70 decimal digits**, with
separate cold and warm rows. Each row names the complete best-found recipe,
comparison controls, independent input coverage, completion, elapsed time,
CPU, peak owned RSS, uncertainty, near ties and untested combinations. State
"best found and confirmed within this search", not global optimality.

**75 and 80 digits are separate, user-gated feasibility probes.** Do not put
them in the ordinary breadth, adaptive-search or confirmation batch. Their
only initial deliverable is the single-probe receipt and a notification.

Use the existing matching worktree
`/Users/nishanthmohan/Documents/code/factor-a7a-regime-tuning`, branch
`codex/a7a-regime-tuning`. Its implementation checkpoint is `276ebd8`, based
on mainline `e62ec18`; its recorded experimental source is `ae1cfd9`.
Documentation-only policy/plan commits may follow. Record the actual current
commit in every new freeze. Preserve all historical freezes and captures.
Do not reset another checkout, merge, push or change automatic defaults.

The existing `joint_search.py` and `regime_search.py` provide candidate
factories; **a complete automatic per-regime adaptive dispatcher is still
missing**. Implement and validate the dispatcher described below before the
ordinary search. Do not present manually generated batches as a completed
adaptive campaign. The single large-regime probe can use the existing runner.

## Exact large-regime gate: one candidate, then stop

Run 75 first when that probe is authorized. An 80-digit probe requires a
subsequent instruction; it must not be automatically submitted after 75.
After either probe, stop and notify the user regardless of its outcome.
Do not issue a second candidate, repetition, clean companion, larger grant,
alternative seed, numerical-policy change or confirmation batch.

| Field | 75-digit probe | 80-digit probe |
| --- | ---: | ---: |
| Input | `regime_calibration_75_0` | `regime_calibration_80_0` |
| Collection workers | 8 processes | 8 processes |
| Polynomials per packet | 16 | 16 |
| Actor wall allowance | 1,800 s | 3,600 s |
| Whole-owner wall ceiling, including capture/cleanup | 1,900 s | 3,700 s |
| Actor / owner aggregate CPU ceiling | 18,000 / 19,000 s | 36,000 / 37,000 s |
| Aggregate owned RAM | 22 GiB | 22 GiB |
| Worker / runtime allowance per worker | 768 / 256 MiB | 768 / 256 MiB |
| Parent QS / checkpoint allowance | 4 GiB / 256 MiB | 4 GiB / 256 MiB |
| New capture storage ceiling | 8 GiB | 8 GiB |
| Logical work ceiling | 96,000,000,000,000 | 96,000,000,000,000 |

Both use seed `104729`, nursery 4 MiB, the existing `Y0` numerical policy,
leases, persistent state, packed transport, worker verification, private
semaphore stopping, coefficient reuse, lookahead 4, state cache 64 MiB,
packet wire limit 4 MiB and atom batches of 2,048. Add the existing `wide_hits`
kernel for the large base; do not change report/product/endpoint bounds.
Record the complete resolved configuration, including residual-floor/domain
limitations. This candidate is a feasibility proposal, not an optimum.

Use **one instrumented cold call**, with full parent/all-worker CPU capture,
parent/representative-worker heap and owned RSS. Its elapsed time is diagnostic
feasibility information, not clean speedup evidence. Select no warmup and no
duplicate full factorization. Validate admission bounds before starting;
if the configuration is refused, report the refusal and stop rather than
silently shrinking it or trying another configuration.

Report the input/source identities, exact reconstruction, completion or
remaining cofactor, stop reason, charged work/wall/CPU, RSS, capture bytes,
collection progress, whether a matrix solve occurred, and CPU/heap findings.
Distinguish recorder exhaustion, candidate admission refusal and production
resource exhaustion. A timeout does not prove that every configuration or
the uninstrumented factorer is infeasible. Retain the failure and notify.

### Preparing and running that single probe

Run from the A7a worktree using PyPy implementing Python 3.11. Keep generated
files under folder-ignored `v2/benchmarks/results/a7a/`, outside Documents when
an equivalent local evidence root is available. The command below uses the
main repository's ignored evidence root.

The operator must first record a new approved grant in `grant.json` with
`identity`, `regime`, `status: "approved"`, `wall_increment`, `cpu_increment`
and `capture_increment`. The increments must cover the table above. This
document proposes those limits; committing it does not approve new timings.
Use `prior_ledger` in that file to identify the last sealed ledger. Carry its
entire prior-charge chain and captured bytes; link the unchanged M1/M2/M3
history. The old M3 allowance has only 378.775 seconds remaining and cannot
fund either large probe. Do not reset it or roll unused stage time over.

The following is an exact plan-freezing recipe; it does not run a factorer:

```python
import json
from pathlib import Path

from v2.benchmarks.qs.a7a.cell import configuration, settings
from v2.benchmarks.qs.a7a.combined import cell, corpus_identity, fixtures
from v2.benchmarks.qs.a7a.runner import source, write

root = Path("/Users/nishanthmohan/code/factor/v2/benchmarks/results/a7a")
grant = json.loads((root / "grant.json").read_text())
assert grant["status"] == "approved"
digits = grant["regime"]
assert digits in (75, 80)
actor_wall, owner_wall = (1800, 1900) if digits == 75 else (3600, 3700)
assert grant["wall_increment"] >= owner_wall
assert grant["cpu_increment"] >= 10 * owner_wall
assert grant["capture_increment"] >= 8 * 2**30
prior = json.loads(Path(grant["prior_ledger"]).read_text())
assert "capture_bytes" in prior
assert all(r["status"] != "issued" for r in prior["attempts"])
charges = list(prior["plan"].get("prior_charges", []))
charges.append(dict(identity=prior["plan"]["identity"],
                    wall=prior["wall"], cpu=prior["cpu"]))
old_storage = prior["plan"].get("prior_storage", 0) + prior["capture_bytes"]
identity = grant["identity"]
folder = root / identity
folder.mkdir(parents=True, exist_ok=False)
fixture = f"regime_calibration_{digits}_0"
config = configuration(int(fixtures("regimes")[fixture]["n"]), "Y0")
marking = config.collector.marking
if "wide_hits" not in marking.split("+"):
    marking += "+wide_hits"
features = dict(lease_budget=True, reuse_state=True, packed_transport=True,
                verification="worker", private_stop=True,
                stop_backend="semaphore", coefficient_cache=True,
                cpu_snapshot=False, lookahead=4, arrival_order=False,
                whole_family=False,
                worker_memory_bytes=768 * 2**20, state_cache_bytes=64 * 2**20,
                packet_bytes=4 * 2**20)
settings("C8", 22 * 2**30, True, 16, 2048, features)
probe = cell("C8", fixture, stage=identity, variant="single_feasibility",
             endpoint="complete", lifetime="cold", samples=1,
             warmup_seconds=3, mode="profile", observe=True,
             full_profile=True, profile_reference_seconds=120,
             tier_a=False, nursery_mb=4, seed=104729,
             wall=actor_wall, cpu=10 * actor_wall, work=96 * 10**12,
             memory=22 * 2**30, storage=8 * 2**20,
             checkpoint_bytes=256 * 2**20, qs_memory_bytes=4 * 2**30,
             polynomials_per_task=16, max_batch_atoms=2048,
             features=features,
             collector_overrides=dict(marking=marking,
                                      power_plan_bytes=2 * 2**20,
                                      sparse_queue_bytes=32 * 2**20))
plan = dict(identity=identity, source=source(), corpus_kind="regimes",
            corpus_sha256=corpus_identity("regimes"), analyze_profiles=True,
            question="Does this one candidate complete within its grant?",
            decision="Stop and notify after this one call, for every outcome.",
            wall=sum(r["wall"] for r in charges) + owner_wall,
            cpu=sum(r["cpu"] for r in charges) + 10 * owner_wall,
            storage=min(64 * 2**30, old_storage + 8 * 2**30),
            prior_charges=charges, prior_storage=old_storage,
            grant=grant,
            profile_folder=str(folder / "profiles"), cells=[probe])
assert plan["storage"] >= old_storage + 3 * 2**30
write(folder / "plan.json", plan)
print(folder / "plan.json")
print(folder / "receipts")
```

`lifetime="cold"` skips warmup; the accepted `warmup_seconds=3` metadata
default above does not execute a warmup in this call. Verify recorded warmup
work and seconds are zero.

Execute only that frozen plan, from the same worktree:

```sh
PYTHONPATH=/private/tmp/factor-a7a-vmprof-reader-20261105 \
  /opt/homebrew/bin/pypy3 - <<'PY'
import json
from pathlib import Path
from v2.benchmarks.qs.a7a.combined import run

root = Path("/Users/nishanthmohan/code/factor/v2/benchmarks/results/a7a")
grant = json.loads((root / "grant.json").read_text())
assert grant["status"] == "approved"
folder = root / grant["identity"]
plan = json.loads((folder / "plan.json").read_text())
assert len(plan["cells"]) == 1 and plan["cells"][0]["samples"] == 1
run(folder / "plan.json", folder / "receipts")
PY
```

Verify the optional vmprof reader before acquiring the window; the temporary reader path above
is the existing local installation, not a portable project dependency.
The runner acquires the shared performance lock nonblocking. If unavailable,
do not poll, run concurrently or claim a reservation. Hold it through child
cleanup, retain its file, and finish with the mandatory notification.

## Ordinary search: 40–70 digits

### Preparation before timings

1. Implement the finite adaptive dispatcher and its pure planning tests.
   It must persist per-regime coverage, pending cells, exact charges, source
   identities, censoring, paired ordering and a restart cursor. Test that a
   regime cannot be silently omitted and that a phase cannot exceed its grant.
2. Add a separately named immutable confirmation corpus and loader support:
   15 fresh certified balanced moduli for each admitted regime, using generation
   seed `202611071442`. Verify certificates, exact digit/factor sizes and
   pairwise coprimality against existing corpora. Existing generation primitives
   are `certified_prime` and `verify_certificates` in
   `v2.benchmarks.suites.build_phase_two_corpus`, with `GenerationRandom` in
   `v2.benchmarks.qs.c9a.inputs`; the previous construction is
   retained locally in `v2/audit/a7a/m_fresh_inputs_prepare.py`. Do not overwrite
   `regime_corpus.json` or train on the new confirmation inputs. The current
   three held-out moduli per regime do not supply nine independent inputs.
3. Add and validate a literal committed-mainline reference adapter. The A7a
   branch's `S0/Y0` is a same-policy serial control, not literal master. Pin the
   current mainline commit at execution; keep wrappers in benchmark code.
4. Add bounded diagnostic-window support before using shortened companions.
   Every executed candidate gets CPU and memory diagnostics, including losers.
   Calls up to 45 seconds may use a full companion. Longer repeated trials use
   a 45-second CPU/heap/RSS companion, plus one full diagnostic for each finalist
   and each concrete code change. Preserve phase/window limitations. Do not
   silently replace these with a full duplicate of every expensive repetition.
5. Verify exact arithmetic, leases, resume, output and cleanup tests; use the
   full committed-checkout suite, lint and benchmark imports for code changes.
   Heavy QA needs the shared lock. Documentation-only edits do not need
   another factoring run.

These are explicit missing implementation prerequisites, not commands claimed
to exist already. No broad timing campaign starts before they are satisfied.

### Regimes, inputs and grants

Use certified balanced semiprimes. Calibration inputs are exactly
`regime_calibration_{digits}_{0,1,2}` from `regime_corpus.json`. Known factors
validate results only. Freeze seeds `104729`, `130363`, `155921` before results;
match a seed between arms, and never choose it from outcomes.

| Digits | Per-call wall ceiling | Joint search/refinement wall ceiling |
| --- | ---: | ---: |
| 40 | 30 s | 15 min |
| 50 | 90 s | 45 min |
| 55 | 180 s | 75 min |
| 60 | 300 s | 90 min |
| 62 | 420 s | 120 min |
| 65 | 900 s | 180 min |
| 70 | 1,800 s | 300 min |

Use a machine-wide CPU ceiling of ten times wall, aggregate RAM 22 GiB,
parent QS 4 GiB, worker 512 MiB, runtime 256 MiB, checkpoint 64 MiB, packet
wire limit 4 MiB, atom batch 2,048 and logical work 96 trillion per call.
Constructor/accounting checks must admit each recipe; invalid combinations
are excluded with their reason, not silently altered. Freeze resolved base,
factor count, half-width, large-prime/report/product bounds and numerical
policies for each regime before worker comparisons. Start with current G0
and Y0 definitions; any numerical redesign gets a new identity and retuning.

Proposed ordinary-campaign ceiling: 48 hours owner wall, 480 core-hours and
64 GiB captures, including confirmation and diagnostics. This is a proposal,
not a new grant. Preserve the old expenditure separately and cumulatively.
Price the complete schedule from the first coverage observations before
committing the campaign. Report an unaffordable stage before starting it;
do not spend the allowance on one band and discover the missing coverage later.

### Search space

Use serial plus process worker counts 1–10. Start with 1, 2, 4, 6, 8 and 10;
test adjacent counts 3, 5, 7 and 9 around surviving bundles. Packet sizes are
the valid levels of `initial_design(gray_limit=config.gray_limit)` and adjacent
valid levels up to one complete family. Never submit a packet wider than the
declared family or pretend a conditionally inactive flag was exercised.

Jointly vary leases/stop backend, persistent state, packed transport,
coefficient reuse, CPU snapshots, arrival/ordered admission, lookahead
0/2/4/8, parent/worker/batch64/batch256 verification, and the collector kernels
in `joint_search.BINARY`. Include `presieve_powers` as a presieve alternative.
Include nursery 4/16/64/128 MiB, state cache 16/64/128 MiB, sparse queue
0/8/32 MiB and power plans 0/2/8 MiB where active and admitted. Wider inverse,
hit and transfer caches and pruned transfers remain eligible by regime;
an inconclusive standalone screen does not exclude their combinations.

The shared cofactor queue is unimplemented. Do not claim K=1/2 was tested.
Implement its exact ownership/lease/resume contract before those arms exist,
then compare collection/cofactor allocations with total compute workers at
most ten. Keep K=0 evidence distinct. Native/NEON work remains deferred to B5a.

### Execution order and adaptive decisions

1. **Coverage pass over 40–70:** run current production reference, same-policy
   serial and retained process incumbent on calibration input 0 in every
   regime before locally refining any one regime. Explicitly set
   `endpoint="complete"`; `regime_search.cells()` defaults to prefixes above
   55 and must not be used unchanged. Cold cells require `samples=1` and a
   fresh interpreter/pool for every repetition. Warm cells use at least three
   seconds of validated warmup and a separately declared retained lifetime.
2. Estimate A/A variability in each admitted regime and lifetime with a small
   prespecified control block. Randomize/counterbalance subsequent order and
   interleave controls; retain near ties. Record `pmset -g therm`, Foundation
   hints and host pressure outside timed intervals. Drift is uncertainty,
   not an automatic campaign-wide stop or a reason for an unbounded diagnostic
   loop. Endpoint agreement does not establish equal intermediate conditions.
3. **Joint initial search:** generate the valid `coarse()` and
   `initial_design()` pools. Force serial and one valid candidate at each
   initial worker count, then select additional joint bundles greedily by
   uncovered factor-pair levels, with deterministic tie order. Maximum total
   distinct recipes per regime/lifetime: 16 initial plus 8 refinement recipes;
   for calls over five minutes, 8 initial plus 4 refinements. Record actual
   pair coverage and untested combinations; do not claim all 32 design rows ran.
4. Race on complete-factor observations from input 0; clear losers need no
   fixed nine-repeat screen. Retain uncertain/near-tied candidates. Evaluate
   up to four survivors on calibration inputs 1 and 2. A reference timeout
   does not prevent measuring a valid candidate. If all compared arms time
   out without answering the question, diagnose and report; do not automatically
   enlarge limits. Fixed-work prefixes are mechanism diagnostics, not ranking
   evidence unless their full-run predictive validity was established first.
5. Refine by one/two option changes, adjacent workers/packets/cache levels and
   at least one joint exploration proposal outside the current neighborhood.
   Finish a coverage/refinement round in every admitted regime before another
   deep round at 60 digits. The dispatcher must enforce this barrier.
6. Review CPU/heap evidence, fix the largest concrete measured cost, revalidate
   and remeasure full factoring, then retune affected combinations. Bound each
   iteration inside the stage grant. Code/numerical/backend/lifetime changes
   invalidate dependent cost/engagement tables; do not mix source cohorts.

Prepare frozen JSON batches with `combined.cell`, the candidate factories and
`source()`/`corpus_identity()`. Always set the endpoint, lifetime, sample count,
full resource grants and source explicitly. Execute each batch only through
the existing `combined --plan ... --output ...` command above. The new
dispatcher must automate these decisions; it must not add hidden arms,
selectively retry outcomes or overwrite completed plans.

### Independent confirmation and stopping

Freeze all regime/lifetime recipes and reference identities before opening any
confirmation result. Use the first nine independent new moduli in their
committed order. For each comparison, use the following exact allocation:

| Calibrated complete-call cost tier | Initial confirmation allocation |
| --- | --- |
| Under 10 s | 9 independent inputs × 5 paired repetitions = 45 pairs |
| 10 s–5 min | 9 independent inputs × 2 paired repetitions = 18 pairs |
| Over 5 min | 9 independent inputs × 1 paired repetition = 9 pairs |

Fix the tier from calibration before confirmation; never choose repetitions
from a favorable confirmation outcome. Repetitions of one input are a cluster.
For cold measurements each repetition starts a fresh process; warm repetitions
can share only the explicitly frozen warm lifetime. Rotate arm order and
interleave controls. Compare with the best valid serial control and literal
mainline; identical controls may be deduplicated with proof of identity.

Predeclare looks after 9, 12 and 15 independent inputs, retaining the same
repetition count per input. Extend to the next three reserved moduli only
when the adjusted uncertainty interval straddles the decision and the remaining
grant funds the entire next block. Account for repeated looks and multiple
comparisons in the declared interval procedure; do not reuse unadjusted 95%
intervals at every look. Freeze and test that analysis before timing. At 15,
stop and report inconclusive if uncertainty remains. Do not retune on these
inputs or selectively retry failures. Preserve censoring; a fast refusal or
timeout is not a fast factorization. Do not insert a timeout cap as an exact
time-to-factor observation. Report completion and restricted-cost summaries
separately, and retain the protected 60/70 completion guardrail.

Confirmation has its own priced, finite plan inside the proposed campaign
ceiling. Compute its forecast from the frozen number of inputs/repetitions,
arms, measured complete costs, warmup, controls and bounded/full diagnostics.
Do not silently double it with full diagnostic copies of every long call.
Do not claim p90 precision from nine independent inputs without an adequate
uncertainty interval. An inconclusive row retains its reference and near ties.

## Final report and completion criterion

Maintain a coverage table after each round: regime/lifetime, number of distinct
recipes, independent inputs, paired observations, profiles, failures, charges,
remaining cap and next decision. Execution-cell counts are not configuration
counts or completion metrics. Preserve raw evidence locally under ignored
results/audit folders; publish concise decisions and runnable commands.

The work is complete only when the admitted regime table and its confirmation
are delivered, or a specific regime is reported infeasible/inconclusive under
its declared grant. List unimplemented mechanisms and untested interactions.
The 75/80 probes never authorize their search or confirmation automatically.
