# Source-derived collection follow-up

The deeper source comparison found a useful missing mechanism: exact sparse
hit scheduling across working blocks. A separately implemented, frozen prototype
reduced median **fixed-prefix wall and CPU cost by 40.4%** against the existing
fused collector across nine matched pairs on three exposed 60-digit inputs.
It preserved exact output and coverage counters. This is a prefix throughput
result, not a complete-factoring speedup or adoption decision. The complete-cost
campaign remains stopped for Sol high.

The earlier recovery prepared the harness and a small collector shortlist, but
did not finish distinguishing this source gap from mechanisms already present.
That handoff was premature as a research conclusion. Its freezes and negative
probes remain unchanged; this follow-up has new source, control and output
identities.

## What the source comparison changed

Factor's bucket division stores hit primes for candidate recovery. Its scorer
still visits every factor-base prime for every block. Native sparse schedulers
instead prepare hits for the blocks they reach; the relevant paths are visible
in the pinned [FLINT source](https://github.com/flintlib/flint/blob/e473ecf2bae34bd3f171b92561c55395dedd3ada/src/qsieve/collect_relations.c)
and [Yamaquasi source](https://github.com/remyoudompheng/yamaquasi/blob/270a66177cc10ab52d49ea4f371a5e1f69664b4b/src/sieve.rs),
and official msieve 1.53's retained source. Their native thresholds and timing
claims are not transferred.

The new prototype schedules sparse prime and exact prime-power hits once per
polynomial, then replays each block's entries. Existing dense marking, thresholds,
division, certainty and reconstruction remain intact. Queues publish only when
complete, rebuild with charges on resume, and fall back exactly if unsupported
or full. It does not copy native overflow losses or the statistical small-prime
allowance described in [Bradford–Monagan–Percival §2.4](https://www.daemonology.net/papers/maplent.pdf).
Independent brute-division and score/hit tests cover exceptional roots, faces,
negative positions, tails, capacity refusal, interruption and checked restoration.

| Mechanism | Finding and owner |
| --- | --- |
| Sparse block scheduling | Missing and now implemented as an exact study-only candidate; A1a owns production integration. Distinct from C8's whole-family CRT scheduler. |
| Dense two-root interleaving | Still worth A1c investigation after scheduling; existing single-root unroll4 is different. No new timing claim. |
| Actual A-product quality and diversity | Current census already looks strong; demote an immediate C10a-02 rewrite. Unique products alone do not prove optimal useful yield. |
| Normalized odd/mod-8 polynomial forms | Real representation gap for C10a-01, requiring new root/relation/parity/checkpoint proofs. Not a safe constant substitution. |
| Small-prime refinement, presieving, Gray-root reuse | Already implemented in relevant forms. Native probabilistic loss policies do not preserve this stream's coverage requirement. |
| Base/interval/factor count and SLP/DLP policy | Coupled questions with existing owners/evidence; revisit selectively after collector cost changes, not through a broad sweep. |

Contini's cache-blocking and interval discussion, Carrier–Wagstaff's scan and
feedback sections, Alford–Pomerance's root reuse and Percival's repeated-prime
overhead discussion were revisited against executable paths. The detailed
comparison, citations, revision/file hashes and exclusions remain local in the
source-gap audit. This is targeted primary research, not exhaustive SOTA coverage.

## Retained evidence and its limits

Both arms use current e520ec9 code with the existing fused kernel explicitly
enabled at 60 digits, canonical signs, a 5000-prime base, quota64 and the same
131073-position interval. Each measured prefix executes 2146 SIQS blocks and
crosses a family boundary. Nine pairs cycle three exposed independent moduli
and seeds 7/29/47, with alternating arm order. Every arm passes process-local
wall AND CPU warmup of at least three seconds and the declared stability check.

Median candidate/control ratios are **0.596475 wall / 0.596468 CPU**. Mean CPU
ratios by input are 0.596718, 0.596300 and 0.595036. Raw CPU ratios span
0.575311–0.612580; all observations, including the slower first control, remain.
All nine pairs match retained-state hashes and scanned, candidate, refinement,
root-hit, power-hit, recovery and admission counts. No censoring or capacity
fallback occurred. The three independent inputs remain three, not nine;
the input-level two-sided sign-flip test has only eight assignments and p=.25.
Fresh population confirmation and whole-factor cost/completion are unmeasured.
No percentage is combined with the earlier fused-kernel probe.

A separate structural census found 4708–4718 of 5000 primes sparse relative to
the 4096-position block. About 76–77% of their first-polynomial block visits
were empty. No base prime was sparse relative to the full interval. This
explains why C8's earlier whole-interval eligibility result did not rule out
the present scheduling opportunity.

The first 128 A products on each input are distinct, with 127 within 1% of target.
Median errors are 0.05845%, 0.04353% and 0.03653%; each sample contains one larger
outlier. Adjacent eight-factor assignments share at most two primes and show
zero immediate flyer repetitions. The obvious spacing/flyer deficiency suspected
from source names is not present in this sample. Multiplied targets are 1,5,5
modulo 8, relevant to a future normalized-form study; no yield gain is inferred.

Native owned memory stays 1 GiB, with a **separate 32 MiB+256 KiB helper grant**
and 3 GiB RSS ceiling. Queue reserve peaked at 17,152,768 bytes. This is not a
claim at 1 GiB total ownership. A default bridge must account for that reserve
inside its intended total grant or state the larger allowance explicitly.

All 32 preparation calls—18 measured prefixes, 11 warmups and three censuses—
used 127.527 actual wall /127.500 CPU seconds, conservatively charging 1135
wall/CPU seconds and 21,837,846,918 work. Remaining: 12610 total seconds,
1780 warmup and 10830 untouched complete-screen seconds; work remaining
999,978,162,153,082. Prefix and census accounts are consumed. Output allowance
is 128 MiB; analysis and QA have separate finite accounts.

## Next action and provenance

Prioritize the exact queue's separately frozen complete-cost screen: fused
control versus sparse queue on the three exposed moduli, seed7, alternating
order, 1800 wall/CPU seconds per call. At most six calls reserve 10830 seconds
(3.01 measured hours); expected 1.5–2.5 hours. If both first-input calls remain
unresolved, stop inconclusive. Rewarm each new process, retain failures/censoring
and all charges, and report capped cost/completion separately. Do not start
the older four-arm screen simultaneously or silently replace its selected arm.

The source-gap handoff in branch `codex/c10a07-source-gap` gives exact commands.
Implementation/protocol: `42521d6`; freeze: `6c211a1`; prefix evidence:
`c58cbb1`; census/handoff: `6c3894f`. Source production runtime remains e520ec9.
The committed-only checkout passed 847 PyPy tests (9 skipped), full lint,
209 benchmark imports and frozen-loader verification in 90.670 seconds.

Freeze SHA-256:
`e2bc23abd372685bbb6c1a3163c8352b383b4858720217c828e18ebf454e4774`.
Required control/prefix/census numerical inputs are committed. Local durable
archive `factor-c10a07-source-gap/v2/audit/results/c10a07_source_gap/source-gap-evidence-v1.tar`
is 1,044,480 bytes, SHA-256
`0bd2dfd1f04c090111924698f20184bd96a8452f6d00c1c9c9057c30ac4eb780`.
It contains raw evidence, ledger, source/protocol/tests, primary references,
research and committed-only QA. Generated evidence is never force-added.
Exact ratios, accounts, identities and the census are also in the committed
[numerical summary](c10a07_source_gap_results.json).

A useful complete-cost outcome would justify a new current-automatic-default
bridge and disjoint fresh confirmation with explicit memory ownership. The
earlier independent/broader-panel design remains applicable only after updating
its identity and cost estimates. No automatic adoption, policy combination,
merge or push. E1 remains deferred.
