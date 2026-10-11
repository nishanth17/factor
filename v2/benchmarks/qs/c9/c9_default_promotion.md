# C9 default promotion — 11 October 2026

The user explicitly requested enabling demonstrated speed improvements by
default. Promote the accepted aligned growing-DLP bundle for fresh native
50-digit SIQS/automatic calls. The source configuration comes from
`inputs/controls/c9_long_selected.json` and matches every field in its
`configurations["50"]["challenger"]` record. No new tuning/timing search or
fresh input generation is conducted. The original acceptance and its
retain-defaults verdict remain historical evidence.

## Selection and resource behavior

`v2.qs.default_siqs_config(n, backend="python-int", memory_bytes=None)` is the
public composition point. `SIQSJob(n)` uses it when config is omitted. The
bounded automatic portfolio resolves omitted SIQS once at the original input;
explicit `siqs=None`, supplied SIQSConfig objects, manual ECM plans, alternate
backends and explicit CLI QS parameter overrides preserve their choices.
The automatic/direct-SIQS CLI selects the bundle without `--qs-dlp` flags.
Children and serialized copies retain the concrete root selection.

The exact observable native band is `10**49 <= abs(n) < 10**50`. Balancedness
cannot be inferred without factors, so applying this size policy to other
50-digit classes is a user-directed extension of the measured population.
The three independent balanced inputs/nine matched groups support the reported
38.10% lower wall cost versus calibrated SLP and 53.08% versus fixed DLP.
Growth results at 30/40 digits were unstable; their recommendations remain.
No 60–100-digit or arbitrary-input benefit is inferred.

When omitted, resource grants use the confirmed 10^13 work and 1,200 wall/CPU
seconds, 1 GiB owned SIQS memory, and 16 MiB extra outer-portfolio headroom.
Explicit Budget/CLI limits remain binding, including zero and smaller grants.
This default window is not an implementation ceiling. The collector grows
within owned memory without FIFO/count quotas; finite split effort and the
256-atom cycle policy remain. The default checkpoint-byte setting matches the
confirmed configuration. Quota refusals retain their documented semantics.
RSS/JIT remains separate from owned memory; no enforced OS RSS bound is claimed.

The measured tradeoff is 326 MiB peak owned reservation versus 167/198 MiB for
SLP/fixed DLP, and more charged work despite lower wall/CPU cost. Larger input,
uneven/small-factor behavior, default portfolio handoff and combined C10 plans
retain their E1/C3/C10 acceptance boundaries. Numerical/default selection is
not a new integrated complete-run timing result.

## Compatibility and verification

Omitted CLI resume restores concrete relation settings, finite memory and
saved direct-SIQS pretest identity. Library implicit resume restores growing
collectors for allocation and ordinary SIQS snapshots. Explicit conflicting
policies still reject; saved legacy no-SIQS snapshots do not adopt the new
fallback. No serialized layout or certainty/provenance contract changes.

Independent default tests compare against the immutable accepted control,
check size/sign/backend boundaries, explicit overrides, finite parent/child
ownership, default versus explicit grants, repeated implicit resume and CLI
routing. Existing C9 arithmetic/kernel/ownership and C3 allocation tests remain.
A separate live-source A7 manifest chains to the integration manifest;
explicit benchmark arms and historical pins are unchanged. C10 can compose
its separately accepted plan through the selector after its own integration,
without modifying either study freeze. C9a remains open for untested transfers.

At runtime revision `4dcc68667eeace352710ce84899b0a8c79394a8e`,
`make -C v2 test` passes **696 tests** in a committed-files-only `git archive`
using external PyPy 3.11 dependencies. `make -C v2 lint` passes all Ruff,
format and pycodestyle checks. All **140 live non-package benchmark modules**
import, and snapshot/A7/B4/C6/C1 loaders pass. The three C9 fresh certificates
validate; the historical timing loader correctly rejects newer production
sources. The selected control and corpus SHA-256 values remain unchanged.

These final checks run under the shared performance lock, using 109.073 wall
seconds and 114.367 child CPU seconds. They are QA costs, separate from the
closed historical timing study. The temporary committed-only archive is
removed; raw receipts/logs remain ignored in `benchmarks/results/c9/`.
This receipt changes documentation only; checked runtime/tests/manifests are
unchanged. The temporary promotion worktree is removed after mainline
integration, preserving unrelated mainline edits and the C10 worktree.
