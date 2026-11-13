# Checkpoints and saved resource history

Run on PyPy implementing Python 3.11 with the same committed source and numerical
configuration. A checkpoint stores published arithmetic state, cursor/RNG/search
identities and cumulative charges. Derived stores are revalidated and rebuilt
with the current invocation charged to the same root owner.

```sh
pypy3 -m v2.factor 626100403 --bounded --seed 7 --seconds 0 \
  --checkpoint /tmp/factor-progress.json
pypy3 -m v2.factor --resume /tmp/factor-progress.json --seconds 60 \
  --cpu-seconds 60 --checkpoint /tmp/factor-progress.json
```

Limits apply to recorded checkpoint charges plus the current invocation. They
are cumulative totals. Work and resources consumed after the last valid save
are not recovered when that save was skipped or a process died. A library caller
retaining the live `Budget` across retries keeps exact cumulative charges. There
is no durable ledger. Benchmark and study ledgers remain the authority for
experiment charges. Restore, replay, validation and failed encoding are charged;
paused/idle intervals between owned calls are excluded.

Resource version `observed-resources/saved-history-v3` binds this contract.
`UnsupportedResourceVersionError` names `unsupported-resource-version` for old
resource schemas; there is no migration. Integrity hashes, mathematical shape,
proper-divisor checks, provenance and exact reconstruction remain mandatory.
Changing an exhausted numerical schedule is distinct from extending root limits.

Saving is optional. `atomic_save(path, snapshot, checkpoint_bytes=None)` encodes
published state with no default size cap. CLI `--checkpoint-bytes` supplies an
optional cap. An oversized proposal records `checkpoint_size_skip` with requested
bytes, cap and scope, leaves the old file valid and continues computation.
A temporary sibling is flushed and fsynced before atomic `os.replace`, followed
by directory fsync. Before replacement, interruptions preserve the previous file. Atomic replacement
exposes a complete old or new file; a later directory-fsync error can leave the
complete new file installed. Allocation,
I/O and unsupported-encoding failures retain distinct identities and causes.
A store poisoned by interrupted publication cannot be reused or checkpointed;
recover a prior valid snapshot or the original run identity.

Wall/CPU checks are cooperative. macOS RSS sums sampled current resident bytes
for the owned tree, including parent, workers and helpers; it can overshoot and
can count shared pages more than once. See the [resource API](../README.md#resource-api).

Automatic ECM snapshots retain their checked envelope, implementation, method
seed, exact curve identities and worker policy. No live process handles or
interpreter internals are serialized. Retained contexts own startup, each active
call and close separately. Owned children drain before a final snapshot is built.

| Component | Additional continuation rules |
| --- | --- |
| [ECM](../ecm/README.md) | Finite declared tiers; exact curve/boundary replay |
| [p−1](../pm1/README.md) | Same n/base/config and charged rung verification |
| [QS/SSS](../qs/README.md) | Reverify relations, rebuild derived queues, preserve numerical quotas |
| [Execution](../execution/README.md) | Shared root owner and stage handoff |

### Automatic native execution identity

The schema-14 automatic envelope separates whole-native intent from concrete
execution state. A checkpoint saved before ECM does not load or hash native
code. Once native execution is bound, resume requires the exact saved backend,
recipe and source/build identity, even if the next stage is classification or
QS. Missing/changed native state raises rather than silently rerunning Python.
Accepted historical checkpoints with an omitted kernel retain Python;
legacy `auto` retains its old stage-one semantics. Existing resource-version
admission still applies. See the [current default scope](../README.md#automatic-ecm-execution-scope).
