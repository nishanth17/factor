# Checkpoints and limits

Bounded execution preserves checked progress and cumulative resource charges.
See [CLI options](CLI.md) or the [Python API](API.md) for starting a run.

## Save and resume from the CLI

```sh
mkdir -p v2/audit/results
pypy3 -m v2.factor 626100403 --bounded --seed 7 --work-limit 1000 \
  --checkpoint v2/audit/results/example-checkpoint.json
pypy3 -m v2.factor --resume v2/audit/results/example-checkpoint.json \
  --work-limit 2000000 --seconds 60 --cpu-seconds 60
```

Preserve the saved method and compatible settings. For explicit QS/MPQS/SIQS
jobs, repeat the method and search options:

```sh
pypy3 -m v2.factor 10002200057 --method siqs --qs-base-bound 400 \
  --seed 7 --work-limit 130000 \
  --checkpoint v2/audit/results/siqs-checkpoint.json
pypy3 -m v2.factor --method siqs --qs-base-bound 400 \
  --resume v2/audit/results/siqs-checkpoint.json --work-limit 200000000
```

Omitted allocation, tiers, relation bundle and memory settings restore saved
choices. Repeat other nondefault CLI options such as backend and Fermat
settings. Explicit incompatible overrides are rejected.

## Cumulative allowances

A resumed grant is a total allowance, not a new grant. Previously consumed
work and active-run wall/CPU time remain charged; paused time is excluded.
Verification, schedule reconstruction and disposable-cache rebuilding also
consume resources. A grant only equal to prior work can refuse reconstruction.

A deadline extension does not restart exhausted families, curves or campaigns.
ECM tiers and p−1 campaign rungs are declared up front. Direct SIQS quota
extensions require `allow_extension=True` and a compatible identity; they are
not arbitrary configuration changes. Follow the package-specific rules.

## Restore validation

Checkpoints bind the modulus, concrete configuration, backend/build and
execution policy. Restore checks reconstruction and arithmetic/provenance,
then charges replay. Checksums detect accidental corruption and are not an
authentication mechanism. Derived caches are rebuilt rather than trusted.

Supported older formats retain their saved executor, sign schedule, scoring,
primality and work identities. A legacy job is not silently upgraded to newer
defaults. Unknown or incompatible identities are rejected. Backend/build
mismatches cannot silently become a different arithmetic track.

## Resource boundaries

Time and cancellation checks are cooperative; an integer operation finishes
before the next check. Owned-workspace caps exclude interpreter/JIT memory and
caller-retained output. The bounded portfolio's input bit cap is an admission
limit, not a practical factorization-size claim.

| Component | Additional continuation rules |
| --- | --- |
| [ECM](../ecm/README.md) | Saved finite tiers, charged plans/recovery; no same-curve B1 migration |
| [p−1](../pm1/README.md) | Same n/base/config, exact rung replay and cumulative verification |
| [QS/SSS](../qs/README.md#resume-and-resource-bounds) | Original relations reverified, derived queues rebuilt, quota rules preserved |
| [Execution](../execution/README.md) | Shared child allowances and protected fallback admission |
