# Live limits for the 70-digit restart

The opt-in `pilots_live` runner restarts the exposed balanced 70-digit pair
with the frozen v4 nominee and owned mainline control. Both start with 14400
wall **and** CPU seconds and matched tier 3 (4 GiB owned / 8 GiB RSS).
The initial stage ceiling is 29200 wall and CPU seconds: 28810 measured,
300 warmup, 60 analysis and 30 orchestration. Exact work caps, inputs,
configurations and interrupted-run identity are in `pilot_live.json`.
The interrupted original candidate retains 7240 seconds and its full work
reservation. It is never overwritten or represented as continued work.

The runner can apply explicitly requested **increases** while the same job
continues. It reserves time and the derived exact work allowance first,
records a revision, applies it in the worker's existing Budget callback, and
acknowledges the PID, clock usage and consumed work at application. Clock
origins, work, RNG, polynomials, relations and dependencies are not reset.
The supervisor reads the same durable grant, including a pending request,
so its original watchdog cannot discard an accepted extension.

Requests use a compare-and-swap revision. The default `--arm both` updates the
active and pending arms together when both still accept changes. A completed
arm cannot have its grant changed retrospectively. A single-arm extension
must be reported as a changed resource allocation. No extension can revive a
returned timeout. Submit requests before expiry and check the acknowledgement:
`applied` means the worker has changed its Budget; `awaiting_worker` is only a
reservation, not proof of application. Extension admission closes 20 ms before
the old wall/CPU limit when observed, serializing the final worker check with
requests. Cooperative scheduling and in-flight native operations have no hard
latency guarantee. A late or closed request fails explicitly.

From the task worktree:

```sh
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.pilots_live check
/opt/homebrew/bin/pypy3 -u -m v2.benchmarks.qs.c10a01.pilots_live run
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.pilots_live status
```

For example, to raise two still-open calls from four to six hours after
checking that status reports revision 0:

```sh
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.pilots_live extend \
  --seconds 21600 --expect-revision 0 --arm both \
  --reason 'Explicitly authorized six-hour comparison grant'
```

This is an example, not an authorized or executed extension. Each accepted
request publishes the new finite stage envelope and is preserved in the
journal. Earlier requests and failures remain. There are at most 128 revisions;
each request adds a ten-second wall/CPU orchestration reservation to the new
envelope, covering controller and extra status/validation work.
State is bounded to 256 KiB, receipts to 512 KiB, child output to 64 KiB and
stage storage to 16 MiB. Source/input/policy/memory substitutions, decreases,
automatic retries, and unrestricted limits are not supported. Exact caps and
validation remain active. Both arms carry the same live polling wrapper.

The control plane uses at most 16 MiB of the existing 32 MiB harness ownership
suballocation: bounded state/JSON encoding, a bounded receipt and output buffers
coexist; no relation-store copy is created. This is reservation accounting,
not a physical PyPy RSS theorem. Full checkpoint serialization remains deferred
under its existing ownership limitation. Live extension deliberately avoids
that operation. If the worker exits or is killed, its uncheckpointed state is
lost and the latest admitted reservation stays charged. Process warmup and
nine fixed semantic sentinels precede and, where feasible, follow each call.

This is exploratory feasibility execution. Complete results, censoring,
extensions, actual costs, capped service cost and drift are all retained.
Exposed inputs and dynamically changed resources cannot establish adoption.
