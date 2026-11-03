# Benchmark decisions and limitations

## C10a-01 exploratory complete-factor pilots

On one exposed balanced modulus per size, the original frozen control and v4
candidate both completed at 60, 62 and 65 digits. Wall seconds were respectively
432.656 / 204.734, 935.933 / 708.487 and 1572.906 / 2310.505 (control/candidate).
All factors matched the input certificates; runtime primality labels remain
probable-prime labels. The 65-digit candidate's identical-work post/pre CPU
sentinel ratio was 2.504, failing the drift screen. Its regression is retained
for a separate diagnosis; these observations do not establish adoption.

The original two-hour 70-digit candidate was explicitly interrupted to restart
with four-hour grants. Its full admitted reservation remains charged; its
unmeasured final state is not a timeout or a completed factorization. A new
opt-in [live-limit runner](qs/c10a01/live_protocol.md) adds acknowledged,
journaled extensions before that restart. The frozen input, arithmetic,
policies and matched tier remain unchanged. Detailed receipts and failures
are retained under ignored results/audit directories. Defaults are unchanged.

The four-hour 70-digit candidate subsequently completed in 7603.072 wall
seconds (2 hours 6 minutes 43 seconds), with verified factors. The control was
explicitly user-interrupted after about 3 hours 31 minutes 45 seconds of
measured elapsed time. It has no returned factor or final runtime; its full
reservation remains charged. This is user censoring, not a four-hour timeout
or an identified successful-factoring speed ratio.

The 65-digit follow-up found a macOS power-daemon record switching Low Power
Mode from off to on during the original candidate, after control completion.
The old timing comparison therefore crossed a host-condition change. Short
same-work probes confirmed a substantial host-mode effect in both SIQS and
a native-code canary. The benchmark now records power source/mode, requested
thread QoS, GC/JIT environment and native-canary drift. No arithmetic repair
is justified by that evidence. A separately identified, source-pinned
[host-audited recheck](qs/c10a01/regression65_recheck_protocol.md) preserves
the original results and all censored diagnostics; defaults remain unchanged.

The 65-digit recheck completed in 1703.859 seconds for the candidate and
2007.330 for control. Both reproduced their original work, relation digest
and every non-timing SIQS statistic, and both passed SIQS drift checks with
stable recorded power mode. The candidate's elapsed time was 15.1% lower on
this exposed modulus. Control failed the additional native-canary gate
(post/pre CPU ratio 0.695), so this is qualified exploratory evidence, not
a clean accepted performance gain. Two short native-only probes showed
that three seconds of native warmup did not eliminate canary variability.
The failed gate is preserved; no additional full pair or default change follows.
