# C10a-05 correctness evidence

The implementation freeze is `ea9b66cd588528d2fdfb661d25f1fa270fdcabd6`,
based on mainline `16fe8a6180b3b86553cb2d2e56e7c27c9098e2ca`.
PyPy 7.3.23 implements Python 3.11.15. The implementation freeze passed
737 tests; final integration passed **741 tests** (95.788 s test body,
96.701 s inclusive runner). `make -C v2 lint` passed Ruff checks/format and
pycodestyle (4.575 s). Thirty-five focused tests passed before the full suite.
All heavy suite runs held the machine-wide performance lock. Local transcripts
and inclusive costs are in `results/c10a05/qa/`, charged to the protocol's
1800-second QA/feasibility reservation. Committed-only verification is recorded
with the final delivery results.

The committed-only export of `f11f97dac8e401cdf0b6cac23d5802f8f229b3e3`
also passes all **741 tests** (102.419 s test body, 103.441 s runner) using
the external PyPy environment. All **176 live benchmark module/package
imports** pass, as do A7 source pins and M8/M12/M25/C1/A6 snapshot loaders.
The first broad import enumeration included a baseline source file outside
its intended package; that harness failure is retained, and the corrected
check uses the supported historical loaders. Required inputs are committed.
A separate committed-only export of measured commit `96e2aad` passes strict
runtime/helper/protocol/selection hashes and independent corpus verification.
The closing commit changes documentation only after these checks.

## Independent witnesses and state transitions

`test_c10a05_pacing.py` uses N = 10403 = 101 × 103. The atoms U = 139 and
251 have values `2*7**3*13` and `2*7*13*17**2`. Their parity rows are both
14; exhaustive dense-kernel enumeration gives {0, 3}, independently of the
filtered solver. Integer square root of their product and gcd(U1*U2−Y, N)
give 101. Both filter variants lift the same known kernel. The useful pair
is deliberately below the configured surplus of 32: the scheduler may defer
it while collecting, but direct QS and genuine terminal/storage boundaries
must test it. A singleton prefix remains available for its later partner.

The N = 77, U = 9 square gives gcd(9−2, 77) = 7 immediately despite cadence
32. N = 97, U = 49 is an independently verified trivial square and is not
retested on an unchanged store. Opposite U = ±139 form a trivial dependency;
adding U = 251 gives a proper one. Interruption after the first extractor trial
preserves its index through JSON serialization. A separate interrupted solver
has already committed a pivot before restore replays its checked prefix.

Tests exercise ordinary windows, repeated windows, actual family/stall/trivial
stopping transitions, all finite terminal reasons, full storage, refusal before
preparation/filtering/elimination/extraction, partially completed downstream
work, cancellation and work/wall/CPU refusal. Three repeated serialized pauses
assert cumulative work and prior active time; a real multi-polynomial job is
serialized after each block and compared with uninterrupted execution.

Schema 6 records the zero-row cursor and pending terminal reason. Fixed DLP,
growing DLP and ordinary collectors round-trip; incompatible policy changes
and invalid cursors are rejected. Legacy schemas retain window semantics and
omit the new default configuration field. Nested portfolio snapshots restore
with inferred or explicit configurations and reconstruct the original input.
An interrupted terminal solver also survives a monotone stop-limit extension
followed by three serialized pauses with cumulative work grants. Production
selection is compared field-for-field with the frozen 50-digit challenger and
40-digit control. Other memory allowances/backends/sizes and explicit configs
keep window policy. Legacy 50-digit standalone and portfolio snapshots retain
their concrete window configuration under the newly promoted default.
The existing suite independently covers legacy checkpoint generations,
arithmetic, certainty, storage publication, exact provenance and extensions.

One existing edge case was exposed: an empty growing-store snapshot's computed
restore cap was smaller than the public 4096-byte minimum. The temporary cap
now respects that minimum, without changing saved formats. Historical A7/C1/C9/
C10 controls remain unchanged; compatibility tests construct their legacy
configuration shapes explicitly, and a new A7 integration manifest pins the
current adapter.

`test_c10a05_study.py` verifies the single-field arm difference, two successive
yield leases, and a matched group interrupted between arms. The saved first
arm survives, the failed arm remains charged, and continuation completes the
group without overwriting evidence. These checks are independent of timing
acceptance; every measured outcome additionally validates proper divisors,
factor multiplicities, unresolved reconstruction, certainty and resource caps.

Preflight failures and fixes are retained locally. A formatter-only correction
preceded all accepted timings; both the original preflight control and final
control are versioned, with their relationship recorded in the final manifest.
No failed timing validation is silently removed or counted as a speed sample.
