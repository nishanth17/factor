# A1b existing fixed/list screen v1

Selection/rules/inputs/literal configurations and actual schedule were committed
in `268ca7a`, before candidate timings. Runtime base is
`4a946b8d6294d51ac5d9c36f61880cc63d5b59fd`; no C10a-01 changes.
The detailed prospective decision and reconciliation are retained locally in
`v2/audit/a1b/decision-before-screen.md`. This protocol adds executable handling;
it does not relax that selection or its advancement rule.

The sole candidate changes powers/list to fixed/list, with accepted 2 MiB power
plans and all other fields fixed. Both literal configurations request 32 MiB
queues. Only powers/list is queue-eligible. This is a **bundled scoring and
scheduling** comparison. Never infer an isolated score-arithmetic effect.

Use the pinned PyPy 3.11 binary, source manifest `a1b_source.json` and protocol
`a1b_screen.json`. The interpreter receives N and seed, never corpus factors or
proof metadata. Required inputs are committed, previously exposed data. Nine
moduli are independent; three seeds at 60 digits still represent three moduli.
No fresh corpus is generated for this screen.

Each call starts at setup and ends after its native prefix and exact result /
retained-atom validation; report native and validation-included clocks separately.
Each prefix has 2146 blocks, covering 65 complete 131073-position polynomial
intervals plus 4096 positions: 8523841 scanned positions, 66 polynomials,
with family changes. Match the actual physical cursor/polynomial/family as well
as block counts; charged work is not equal computational effort. Instrumentation
is disabled. Both arms retain multiplier policy, A stream, quota, swap_high,
interval, LP domains and collection pacing. Different conservative candidate
streams are allowed. Exact admissible-value coverage is checked independently;
all retained atoms reconstruct, and resolved/unresolved results reconstruct N.
Corpus certificates never upgrade the algorithm's certainty labels.

The original 1 GiB grant selects default kernels/plans first. Carve out 8 MiB
for runner/observation/streamed validation and use 1016 MiB native ownership in
both arms, without rerunning the memory-sensitive selector. Queues compete
inside the native allowance. Limit process RSS to 3 GiB and generated study
output/evidence to 512 MiB. No automatic refusal workaround or larger grant.

The deterministic counterbalanced schedule is committed in numerical input:
18 calls at60, six at50 and six at55, capped at30/15/20 seconds respectively,
wall AND CPU, with5 seconds reservation margin. Reservation900 seconds total.
Work2e10 per call,6e11 cumulative, derived from the current60 diagnostic's
~4.23e8 charged units/s plus setup/variation reserve. This is a cap, not a claim
that alternate configurations consume units at the same rate.

For each configuration/process, warm the first input of its band on1056 blocks
until >=3 validated wall AND CPU seconds and >=3 complete native prefixes;
last-three range/median <=20% for both clocks. A paused prefix is the declared
readiness outcome; a native refusal/limit/factor/other stop is retained and
stops this fixed-work screen for audit, never mislabeled invalid arithmetic.
Warm calls:5 seconds wall/CPU plus5 margin,1e10 work, maximum30 including end
sentinels;300-second phase cap. Starts are durable. Each band's finally block
persists attempted end observations on the same input/configuration and process.
End/start(last-three warm median) must be in[0.8,1.2] for wall AND CPU. Missing
end evidence remains missing. Resume rewarmes new processes but cannot replace
old missing endpoints or make a mixed-process band qualify.

The frozen advancement rule uses per-input paired wall/CPU ratios, taking the
median over three seeds at60. A band qualifies only if at least2/3 independent
inputs win wall time, geometric mean wall ratio<1, CPU<=1, no correctness or
coverage failure, complete identical physical prefixes and valid stability.
There is no minimum gain or training significance threshold. Ordinary variance
is reported; no selective repeats/extensions. Incomplete/censored pairs cannot
qualify. Failure keeps its charges and stops dependent work. No band qualifying
means terminal retain; supported bands alone may enter a separately frozen pilot.

Reuse `c10a04.accounting.Ledger`: reserve before calls, fsync receipts before
settlement, preserve all failures/unknown charges, recover receipt-before-ledger
crashes without rerunning, and turn lost measured cells into censored observations.
Run only unvisited cells. Full successful time reservations remain charged,
actual work settles; failed/unknown calls keep their work reservation. Source
changes require a new identity. Shared timing/heavy-QA lock:
`/private/tmp/factor-performance.lock`. No other sessions are messaged.

Stage A accounts:900 measured,300 warm,120 diagnostics (90 already charged,
including the short first warm failure),60 analysis,90 orchestration;1470 wall
AND CPU seconds. Orchestration gets a conservative90-second reservation; actual
call ledger overhead is reported separately. Expected screen100–200 seconds,
warm30–60; these estimates are recalibrated from early observations. Diagnostics
used one60-digit prefix, not a whole factorization; observer overhead12.5% on
identical work prohibits treating its attribution clocks as clean performance.

Stage B/C are NOT executable under this identity. A supported nominee needs a
separate committed pilot design, at most one matched complete pair, band caps
300/600/1800 seconds,5-second margins,4000 whole-stage ceiling. Correctness or
coverage failure stops. A completed candidate losing wall OR CPU does not
automatically advance. Unresolved/unstable observations give feasibility only.
Fresh confirmation needs a prospective complete executable design, disjoint
inputs/factors, intact readiness/drift, matched complete costs/completion and
uncertainty, within17700 seconds. No campaign is required to spend allowances.

Commands, from the task worktree:

```sh
/opt/homebrew/bin/pypy3 -B -m v2.benchmarks.qs.a1b verify
/opt/homebrew/bin/pypy3 -Bu -m v2.benchmarks.qs.a1b run
/opt/homebrew/bin/pypy3 -B -m v2.benchmarks.qs.a1b analyze
```

Do not repeat a terminal run against its ledger. Keep raw evidence ignored;
publish the concise terminal report and verify a durable local archive.
