# A1b audit against the execution request

This addendum follows the user's challenge to the completion claim at
566a461. **Actual experiments were run, but full compliance with the original
request was not established.** The 50/60 prefix evidence remains valid for its
limited purpose; the broad A1b completion claim is withdrawn pending disposition
of the gaps below. Original controls, raw evidence, failures, charges and the
verified archive remain unchanged. This audit launches no factoring calls.

## What actually ran

The native runner constructed `SIQSJob` with the frozen current-mainline
configuration and called `job.run(max_blocks=2146)`. This exercised setup,
polynomial/family transitions, score marking, scanning, recovery and admission.
These were real fixed-work factoring prefixes, not mocked timings or a
literature-only comparison. They deliberately stopped before complete factors.

| Original Stage A component | Required measured calls | Actual measured calls | Independent measured moduli |
| --- | ---: | ---: | ---: |
| 60 digits, seeds 7/29/47 | 18 | 18 | 3 |
| 50 digits, seed 7 | 6 | 6 | 3 |
| 55 digits, seed 7 | 6 | 0 | 0 |
| Total | 30 | 24 | 6 |

For example, the durable `cell-00.json` receipt records 60-digit control seed 7:
2.445351 wall /2.444273 CPU seconds, 2146 blocks, 8523841 scanned positions,
66 polynomials, and exact unresolved N. Its fixed counterpart `cell-01.json`
records 4.499584 wall seconds. Both receipts contain source-bound ledger
identity, physical cursor/polynomial, counters and retained-atom digest.
All 12 measured pairs have matching physical prefixes and retained-atom hashes.
Each measured arm passed the frozen warmup and endpoint drift interpretation.

The original output directory is `v2/benchmarks/results/a1b_screen_v1/`;
`ledger.json`, `cell-00.json` through `cell-23.json`, warm/start/end receipts and
the console transcript supply the evidence. The verified archive referenced
in [the result](a1b_results.md) preserves these files and committed source.
A separate current-control diagnostic also ran: identical-work clean/profiled
prefixes measured approximately 12.5% observer overhead. No full-factor pilot
or fresh confirmation was run.

## Requirement coverage

The numbered rows below correspond to the original user request. "Covered"
means the named work has evidence; it does not imply that later conditions
were reached or that a universal performance claim is supported.

| Section | Status | Evidence or limit |
| --- | --- | --- |
| 1. Instructions, worktree and current control | Covered | Separate task branch/worktree from exact base 4a946b8, including committed face/queue integration and fixes; no C10a-01 source, agents, inter-session messages, merge or push. Timing and heavy QA use the shared lock. |
| 2. Reconcile existing options and prior results | Covered for the focused scope | R2 frozen/baseline hashes and archived captures checked; C10 plans, A1a and C10a-07 reconciled; targeted primary source/literature review retained locally. Historical R2 loss and obsolete gate preserved. |
| 3. Establish remaining cost opportunity | Covered with stated limits | One short current 60-digit prefix, clean versus bounded observation on identical work; exclusive category clocks, inseparable costs and observer/Amdahl limitations reported. No complete-factor profile. |
| 4. Effective paths and prospective selection | Covered | Local decision note precedes runner implementation; 268ca7a commits one fixed/list candidate, literal fields, input identities, order, physical work, finite caps and advancement. Actual fixed-score queue fallback is identified as a bundled change. No queue extension or redesigned score mathematics. |
| 5. Correctness and current integration | Partial evidence mapping | Independent R2 bounds/division and current kernel/refusal/restore tests were reused; new tests check fixed/list queue fallback, presieve/scan and saved-policy rejection. Those new tests use the SLP config class, while measured arms use growing DLP. A case-by-case proof/test map for the exact selected composition was not delivered. In particular, the cited zero-norm kernel test uses powers scoring rather than fixed scoring. General suite success does not establish every requested combination. |
| 6. Stage A and band eligibility | Partial | All 18 primary 60-digit calls and all six 50-digit calls ran. No 55-digit measured calls. Warmup allocation failed for harness reasons; 55 is unmeasured, not a measured candidate loss. Per-arm native outcome/readiness and analysis/budget failure paths were not adequately tested before launch. |
| 7. Predeclared advancement | Applied to completed bands | Original input-median/geometric wall/CPU rule was frozen before timing and unchanged. Both measured bands lose all independent inputs. No conclusion about 55 throughput; its readiness evidence is incomplete. |
| 8. Complete-factor pilot | Correctly not triggered by completed screens | Neither measured band qualifies. No pilot was required merely to consume allowance. This supplies no complete-factor performance evidence. |
| 9. Fresh confirmation | Correctly not triggered | No useful pilot-supported nominee; no fresh generation, confirmation timing or adoption claim. The detailed conditional confirmation requirements were therefore not executed. |
| 10. Accounting, recovery and source freezes | Partial | Matched measurements have durable reservations/receipts, source/interpreter/configuration pins and preserved failures. However, replacement terminal analysis ran outside the already fully charged analysis reservation and was described as reporting/QA. That classification was not explicitly reconciled with the protected phase budget. Orchestration was conservatively reserved rather than fully measured. Diagnostic bounds were declared locally before execution but the numerical git freeze came later. Resource observations were below grants; RSS was checked as post-call process high-water, not continuously enforced. |
| 11. Delivery and terminal decision | Evidence delivered; completion overstated | Runner/tests/protocol/numerical inputs committed, raw evidence archived and committed-only QA passed. The retain recommendation has supporting prefix evidence. Closing the whole original request without clearly resolving the above gaps was too broad; this addendum corrects that status. |

## Specific defects and what the data still supports

1. **55-digit warmup allocation:** `a1b.py` continued reserving six endpoint
   slots after earlier endpoints had settled. It stopped the fixed arm after
   two warm calls totaling 2.767 wall /2.766 CPU seconds. Fifty warm-account
   seconds remained. This was an avoidable harness error, not a native refusal,
   evidence against the candidate, or a completed 55-digit eligibility test.
2. **Analysis preparation:** the frozen geometric-mean function applied
   `len()` to a generator and failed. The separate terminal adapter materializes
   that generator and computes the same frozen formulas, so the numerical
   correction is straightforward. Its later analysis must not be called free
   merely because it was invoked from a QA controller. The original failed
   60-second reservation remains charged; the additional analysis accounting
   is an unresolved protocol-compliance issue.
3. **Integration evidence:** existing tests cover the arithmetic components,
   but the report should distinguish inherited component evidence from direct
   coverage of the literal growing-DLP fixed/list composition. Add the missing
   case mapping/checks before claiming the entire requested correctness list.
4. **Preparation/observability:** the warm-budget and generator failure paths
   needed regression checks before launch. The resource grants and actual
   low observed peaks do not establish a stronger continuously enforced RSS
   or independently verified helper-memory contract.

Independently recomputed geometric wall/CPU ratios are 1.273336/1.273401 at
50 digits and 1.830225/1.830295 at 60. Every measured pair is slower for fixed.
Those findings support **retain current defaults for the observed existing
bundle; leave 55 unresolved**. They do not isolate score arithmetic, establish
complete-factor performance, or establish a 55-digit result. The executed
candidate bypasses queues and has more plan refusals; both interactions are
part of this bundle, not reasons to rewrite it into a new challenger.

## Remaining disposition

A complete response to the original execution request needs explicit disposition
of the missing 55 screen, exact-composition correctness evidence and analysis
accounting. A continuation, if authorized/admitted, must use a new prospective
identity, retain all old failures and charges, visit only the six unvisited
55-digit cells, rewarm in its actual new process, preserve the same candidate
and advancement rule, and show that reservations fit. The old frozen runner
must not be patched or its ledger reopened to erase the stop.

The original ledger reports 150 measured, 50 warm, 30 diagnostic and zero
analysis/orchestration seconds remaining. Those balances are not interchangeable
permission. A simple rerun cannot be promised to fit: the old readiness rule
requires at least three calls per arm plus endpoints, whereas the old warm
reservation would require at least 80 seconds in a new process. A revised
finite reservation design and an explicit analysis-budget disposition are
needed before another experiment. No confirmation funds should finance this.

Ordinary correctness/regression QA and a documentary correction remain within
the original scope. No new score mechanism, larger memory grant, favorable
retry, complete-factor campaign, merge or default change is scheduled by this
audit. The verified old archive remains immutable; this addendum is a separate
public correction.
