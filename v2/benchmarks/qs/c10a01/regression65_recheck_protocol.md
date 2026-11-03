# Host-audited 65-digit recheck

The first diagnostic stage is retained under `regression65-v1`. All four
90-second screen calls hit their native wall guard near 106–109 million
positions, before the 134218752-position target. Their before/after CPU
sentinels remained stable, but their starting rate was roughly half the
historical rate. Those censored calls cannot establish observer overhead.
The two longer v1 prefixes are **not launched**: repeating oversized fixed
work would not resolve the question.

The Mac was on AC power with `lowpowermode 1`, no reported thermal warning,
24 GiB physical RAM and 72% reported free memory. The original pilot has no
power-state receipt. A targeted unified-log query subsequently established
a transition from mode 0 to mode 1 at 15:58:20.810 on 2 November, during the
candidate and after the control ended at 15:50:49. The old pair therefore
crossed a documented host-condition change. This establishes a confound;
it does not identify a counterfactual full-run time for either arm.

The user subsequently authorized temporary disable and restoration. Before
changing AC mode, `power_probe.json` reserves exactly two 35-second calls
(same control configuration, nine warmed SIQS samples and nine native-canary
samples), then 5 seconds analysis and 5 seconds orchestration: **80 wall AND
CPU seconds**, 581000000000 total work units. Each probe has a finite
290000000000-unit guard. Mode 1 and mode 0 receipts use separate identities;
neither is overwritten. This tests a host mechanism on identical work, not
independent-input factoring performance. Restore the original AC setting
after the authorized follow-up, retaining battery settings throughout.
The command-line change required administrator authentication, so the user
disabled the setting in System Settings. A GUI-automation attempt stalled
for 2144.422 wall seconds with no factoring process active; that administrative
failure is retained separately, with tool CPU unavailable. It is not a
successful numerical call or hidden measured cost. Restoration also requires
the user's OS-authenticated action if administrator access remains unavailable.

The paired short power probes measured median SIQS CPU times of 0.940781
seconds in mode 1 and 0.546851 in mode 0 (same work and atom digest), a factor
of about 1.720. Native-canary medians were 0.031495 and 0.020032 seconds,
about 1.572 times. These are causal host-mode diagnostics on one execution
per mode, not a portable scaling law or an adjusted original factoring time.

The supported repair is to the benchmark: record active power source/mode,
requested native thread QoS, PyPy nursery and explicit GC/JIT environment,
and a fixed native OpenSSL PBKDF2 canary before and after every call.
The canary runs outside PyPy's compiled Python loops, with nine samples of
200000 iterations. Its 1800000 diagnostic work units per sentinel are
charged separately from SIQS work and included in the existing 290 billion
warm/sentinel allowance. Warmup remains at least three validated wall AND CPU
seconds, followed by nine identical-work SIQS sentinels. Cold startup remains
separate. No factorer kernel, polynomial population, multiplier or pacing is
changed merely to compensate for unexplained timing drift.

The source-pinned `regression65_screen_v2.json` first admits four clean
128-polynomial prefixes: control, candidate, candidate, control. Each has a
30-second wall AND CPU cap and 5-second reservation margin. Unlike the v1
observer comparison, these use the original uninstrumented collection path;
the repair adds only before/after observations and the existing live-budget
poller. All four calls must finish the same 16777344 physical positions with
matching within-arm work, atoms and counters before the complete pair is
admitted. Any unexpected earlier factor or refusal is retained and reviewed.

Screen ceiling: 140 measured seconds + 300 warmup + 60 analysis + 30
orchestration = **530 wall AND CPU seconds**. The expected elapsed cost is
about three minutes at the new observed rate. Exact numerical work caps use
the unchanged `plan_v4.call_work(65, 30)` analytic formula and are committed
in the JSON before execution. Full source and binary hashes are bound there.

A later separately committed `regression65_pair_v2.json` may admit one
complete matched pair, candidate then control, on the same exposed balanced
modulus. That is exploratory reproduction, not an independent input or fresh
confirmation. Set the per-call cap prospectively from the stable screen and
the selected host mode; publish the numerical work and stage sums before
launch. Both arms use the same time/memory grant and their literal v4
configurations. The control source is still identical to current committed
mainline in the operative QS/common/budget paths. No selective retries or
automatic extensions are allowed. Existing explicit live extension semantics
remain available under their new source/stage identity if the user requests
one before a live grant expires. Full serialized checkpointing remains
deferred; an exited process cannot be resumed by this runner.

For every complete call, validate product reconstruction, proper factors,
input certificates, unresolved cofactors and runtime probable-prime labels.
Report collection, preparation, filtering, elimination, work, relations and
RSS alongside total cost. Retain incomplete calls and full interrupted
reservations. A changed host state or failed 15% SIQS/native-canary drift
screen prevents a clean timing claim. Even a passing single-input pair
supports only descriptive exploratory evidence. If native execution is
stable while SIQS alone deteriorates, investigate JIT/GC or policy-specific
cost using a separately identified diagnostic; do not assert a host cause.

Matched Tier 2 remains 2 GiB owned / 5 GiB RSS, with the original explicit
parent suballocations and per-call host headroom admission. Receipts remain
at most 512 KiB, pipes 64 KiB, and each new stage at most 16 MiB within the
cumulative 512 MiB evidence cap. Refuse recorder exhaustion distinctly.
No production/default change or adoption is implied by repairing measurement.

## Admitted complete pair

All four v2 prefixes finished their exact physical work; within-arm work
and atom digests matched. Host states matched and both SIQS/native-canary
15% drift screens passed. The new pair therefore uses the original
3600-second wall AND CPU cap for **each** arm, candidate first. The screen
rate suggests roughly 60–80 minutes for both calls, with substantial
uncertainty from relation growth and downstream work. This is a planning
estimate, not a scaling law. A native limit remains an informative retained
outcome and does not authorize a retry.

| Pair reservation | Wall AND CPU seconds |
| --- | ---: |
| Two 3600-second calls plus 5-second margins | 7210 |
| Process-local warmup and before/after checks | 300 |
| Analysis | 60 |
| Orchestration | 30 |
| Total | 7600 |

`regression65_pair_v2.json` contains each numerical work allowance, inherited
setup/postprocessing bounds, exact configs, source hashes and active power
mode 0. Both arms retain Tier 2. The finite numerical follow-up allocation
is 2190 for v1 diagnosis + 80 for the power probe + 530 for the repaired
screen + 7600 for this pair = **10400 wall AND CPU seconds**. Only 660 of
the v1 allocation was used for calls; its conditional longer stage was not
launched. Unused reservations remain unused. The separately disclosed
2144.422-second administrative GUI stall is additional elapsed time with
unknown tool CPU; it is not included in the numerical-call ceiling.

The complete pair is uninstrumented inside SIQS; host/canary observations
remain in the separate before/after phases. Preserve all factors, work,
time-limit outcomes and drift failures. Inspect both returned results even
if one fails a timing screen; no result is selected away. Verify original
semantic counters/digests where the deterministic traversal reproduces.
No promotion, fresh-confirmation claim or arithmetic change follows from
one exposed modulus. Restoration of the original power setting is required
after the authorized follow-up; it must not interrupt a measured call.

Commands:

```sh
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.regression65_recheck --stage screen check
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.regression65_recheck --stage screen run
# Only after freezing and admitting the pair:
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.regression65_recheck --stage pair check
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.regression65_recheck --stage pair run
```
