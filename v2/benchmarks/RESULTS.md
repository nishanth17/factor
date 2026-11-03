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
