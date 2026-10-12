# C10a-04 Stage 2 handoff

Stage 1 stops at this freeze. Do not merge/push or execute E1. Production
runtime/defaults, v1 and completed studies are preserved. Work is sequential;
no subagents or messages to other sessions were used.

- Managed worktree: `/Users/nishanthmohan/.codex/worktrees/c10a04-stage1/factor`.
- Branch: `codex/c10a04-stage1`.
- Starting committed mainline: `e0dafa1ebce2307f1f2ac3d787af89dd54f5ecfc`.
- Frozen experiment commit: `5dab76ec66c02049095ddf2f75e123f19d0d6627`.
  The final documentation annotation does not alter its hashed sources.
- Controls: `v2/benchmarks/inputs/controls/c10a04_frozen.json` and
  `c10a04_sources.json`; training: `inputs/corpora/c10a04_training.json`.
  The control hashes runtime, all benchmark helpers, source manifest,
  research/protocol and training inputs. `study.load()` verifies these identities.

Read [research](research.md), [protocol](protocol.md), and [Stage 1 results](results.md).
The three frozen arms are reused 50-digit control, a count-scaled compact
bundle and its wider-interval counterpart at 55/60/65 digits. Existing odd
multiplier scoring, exact A target, flyer pool and Gray reuse are coupled;
powers/plan, DLP bounds/storage mechanisms and pacing remain fixed. The source
map traces msieve 1.53, pinned YAFU/Yamaquasi/FLINT/Alpertron executable paths,
including YAFU's representation-dependent Q2 target/score. Native constants
are not PyPy optima. Even-k scoring/normalized representation remain uncovered
C10a-01 opportunities; new diversity/population, feedback, kernels and residual
algorithms stay with their named owners.

A subsequent [current-upstream review](upstream_review.md) records a material
11 October FLINT update, newly examined CUDA-MPQS calibration, and a limited
CADO-NFS SIQS screen. The freeze remains a bounded comparison of three bundles;
it does not establish a comprehensive current-source or state-of-the-art
optimum. Read that review with the original research before interpreting
Stage 2 results.

## Mechanical campaign commands

From this worktree, use the external PyPy 3.11 environment (or an equivalent
PyPy 3.11 environment with the documented dependencies). `.venv` is an ignored
local symlink, not a committed dependency. Commands acquire the shared
performance lock and refuse overlapping benchmarks/tests. Do not launch two
commands concurrently. A lease may yield at a group boundary; rerun the same
command and output directory until its finished marker exists.

```sh
PY=/Users/nishanthmohan/code/factor/v2/.venv/bin/python
$PY -B -c 'from v2.benchmarks.qs.c10a04.study import load; load()'
$PY -B -m v2.benchmarks.qs.c10a04.study train --digits 55 --lease-seconds 90000
$PY -B -m v2.benchmarks.qs.c10a04.study train --digits 60 --lease-seconds 90000
$PY -B -m v2.benchmarks.qs.c10a04.study train --digits 65 --lease-seconds 90000
$PY -B -m v2.benchmarks.qs.c10a04.study select
```

Inspect/commit the generated `c10a04_selected.json`, containing all training
receipt hashes and choices, before fresh generation. This is a mechanical
review of the frozen ranking/stability outputs, not permission to retune.

```sh
git add v2/benchmarks/inputs/controls/c10a04_selected.json
git commit -m 'Record C10a-04 frozen training selection'
$PY -B -m v2.benchmarks.qs.c10a04.study generate
git add v2/benchmarks/inputs/corpora/c10a04_fresh.json
git commit -m 'Freeze independent C10a-04 confirmation inputs'
$PY -B -m v2.benchmarks.qs.c10a04.study confirm --digits 55 --lease-seconds 90000
$PY -B -m v2.benchmarks.qs.c10a04.study confirm --digits 60 --lease-seconds 90000
$PY -B -m v2.benchmarks.qs.c10a04.study confirm --digits 65 --lease-seconds 90000
$PY -B -m v2.benchmarks.qs.c10a04.study summarize
```

The confirmation command returns without work for a band selecting control;
if all bands retain control, generate/commit the empty fresh corpus and publish
that unresolved/no-winner decision. Do not invent a winner or rerun censored
cells. Generation requires a committed selection; confirmation requires both
selection and corpus committed. Default raw output is ignored
`v2/benchmarks/results/c10a04/experiment/`. Preserve it across leases.

Optional resource-sensitivity diagnosis, **after** matched training in that band:

```sh
$PY -B -m v2.benchmarks.qs.c10a04.study grants --digits 60
```

This uses a distinct extension phase, one revealed input/seed and all arms at
the doubled time grant. It never enters selection/fresh acceptance. No other
parameter, source, input-band or resource change is permitted by this freeze.
A separately frozen future tranche can explore larger sizes/grants; current
limits are not universal maxima.

## Resources and decision procedure

No accepted training, selection or fresh confirmation was run in Stage 1.
A conservative 3600 wall/CPU seconds and 10^14 work are prepaid. Remaining
allowances are 1630000 s training, 1090000 s confirmation, 200000 s separate
extensions and 2946400 s aggregate wall **and** CPU, with phase work caps in
the protocol and 29.9×10^15 work aggregate after prepaid Stage 1. Disk is
256 MiB; each call has 10^13 work, 1 GiB owned memory, 1 MiB checkpoint and
3 GiB RSS, plus its matched 1200/7200/21600 s wall/CPU grant. These finite
worst-case reserves can span weeks; they are explicit upper allowances,
not a requirement to exhaust every phase. Inspect yields and reservations.

At least three seconds of validated PyPy wall/CPU warmup per arm/lease and
nine matched observations per band are enforced. Three independent inputs
are crossed with three seeds. Extend all arms to eighteen under the frozen
mixed-direction/15-point early-late/15% repeated-cell-drift rules. Persisting
instability retains control. Rank training by completion then capped cost,
tie to control. Fresh confirmation requires complete planned samples and zero
correctness failures. Equal full completion needs a wholly favorable 95%
input-cluster complete-cost interval; coverage gain needs all challenger calls
to complete, no input to lose completion and a strictly positive completion
interval lower bound. No arbitrary improvement floor. Preserve unresolved
outcomes and distinguish capped cost from speed against successful factors.

Receipt/group continuation is idempotent. A lost measured call becomes a
fully charged censored cell; no selective successful retry. A crash after
publication settles that exact receipt once. Correctness exceptions require
an audit. Work and time never reset; extended resource probes keep their own
charged phase. Existing SIQS checkpoint compatibility/replay is validated;
this runner deliberately does not persist every in-flight polynomial block.

## Completed Stage 1 checks and identities

The frozen commit passed in a committed-files-only Git archive at
`/private/tmp/factor-c10a04-committed-pydeqlyj`, with only an ignored external
runtime symlink added. Runtime was PyPy 7.3.23 implementing Python 3.11.15.

- `make -C v2 test PYTHON=/Users/nishanthmohan/code/factor/v2/.venv/bin/python`:
  757 tests passed, including twelve study tests; 103.742 s test execution.
- `make -C v2 lint`: Ruff check/format and pycodestyle passed (325 Python files).
- All 185 live benchmark modules imported from committed files; frozen source
  verification covered 249 files and independently certified all nine inputs.
- Frozen CLI smoke/resume: two `train --digits 55 --lease-seconds 1` commands
  with `--output v2/benchmarks/results/c10a04/smoke`; cold receipts unchanged,
  both warmups charged, two unique yields, no measured samples started.
- One-block setup preflight: all 27 input/arm cells passed, including exact
  challenger counts and reachable family streams. Receipt interruption tests
  also cover a lost extension cell and an unchanged first-nine marker.
- The main checkout remains clean at the recorded starting commit on `master`.
  No existing runtime Python, v1, baseline or historical study file changed.

SHA-256 identities:

| Versioned identity | SHA-256 |
| --- | --- |
| Frozen control | `29d9da0018727bd18db40f4bfcb4701afa37c5452bd34ea336336972debe90ca` |
| Source manifest | `2c38be3a9b0e493c5a07efbdd02a03f74af5597a030311185c7d7c77837cb107` |
| Training corpus | `31f564029964b9c8cc27a113c6b0f2b0ab551f554c12a25700ba9b2495416bc0` |

All local receipts below are relative to `v2/benchmarks/results/c10a04/`.
The source snapshots in each diagnostic retain their pre-freeze identities;
none is accepted performance evidence. These raw files are intentionally
absent from a committed checkout and are unnecessary for executing Stage 2.

| Receipt | Completed work / charge |
| --- | --- |
| `research/verified.json` | 42 inherited cached files plus newly pinned YAFU auxiliary blob verified |
| `screen-v1/ledger.json` | 15 censored calls; 1125.027 wall / 1063.013 CPU s; 146214678025 work |
| `pilot/compact-55-600.json` and `pilot/ledger.json` | Compact 55-digit factor in 390.669 s, exact reconstruction/probable labels; pilot plus warmup charged 610 wall/CPU s and 89994602634 work |
| `preflight/ledger.json` | 27 validated paused prefixes; 540 wall/CPU s and 328709181 work |
| `smoke/ledger.json` | Cold/warmup only; 150 wall/CPU s and 5906663041 work, covered by Stage 1 |
| `qa/full.log`, `qa/final-preflight.log`, `qa/preflight-2.log`, `qa/smoke.log`, `qa/committed-export.log` | Full checks, strengthened resume test, setup, CLI and committed export transcripts; earlier failed fixtures remain local |

Diagnostic ledgers total 2425.027 wall / 2363.013 CPU s and 242444652881 work.
The prepaid 3600 wall/CPU seconds and 10^14 work cover these charges and
recorded corpus/source/QA work. The smoke ledger is a diagnostic copy and
must not be used as the campaign output. The default `experiment/` directory,
selection and fresh corpus are absent. Start Stage 2 with the default commands
above; all remaining phase allowances stated earlier are intact. There is
no accepted larger-band result or default promotion at this handoff.
