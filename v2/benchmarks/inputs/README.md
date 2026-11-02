# Versioned benchmark inputs

Required numerical fixtures, controls and immutable source baselines. A
committed-only checkout must be able to run tests and import benchmark
runners with these inputs; ignored research and captures are not dependencies.

| Directory | Contents |
| --- | --- |
| [corpora/](corpora/) | Independent inputs, known factorizations and certificate fixtures |
| [controls/](controls/) | Frozen configurations, source pins, ledgers and integration manifests |
| [baselines/](baselines/) | Owned historical implementations and immutable source snapshots |

Preserve historical bytes, hashes and source identities. A changed integration
needs a new manifest or control rather than repinning an old experiment.
Load a complete selected configuration bundle when reproducing a study;
engine defaults and partial parameter copies are different comparisons.

`upstream/` is folder-ignored local source material, including licenses/notices.
Regeneration that needs those external sources is separate from ordinary
fixture loading. Do not force-add sources or generated captures.

Generated output belongs under ignored `../results/<subject>/<experiment>/<run>/`.
See [benchmark guidance](../README.md) and [policy](../../POLICY.md#evidence-layout).
