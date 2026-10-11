# C3 control revisions

Protocol v1 was frozen after implementation commit df5ff28. An initial
single-call harness check rejected the result in validation because historical
corpora store repeated prime integers while C3's generator stores (prime,
exponent) pairs. No pilot, training or accepted timing capture was produced.
The original JSON freeze is preserved as c3_protocol_v1.json.

The repaired harness normalizes historical factors solely in the validation
fixture, and saves full per-call receipts separately from compact aggregate
rows to enforce the 16 MiB output cap. Settings, inputs, seeds, selection,
resource and stop rules are unchanged. The active c3_protocol.json pins this
repair before any pilot/training measurement. Its immutable baseline is the
original integrated b3b3cfb production source; actual catalog paths point to
hash-pinned committed controls. Frozen source is compiled under real file
paths for resource lookup; missing modules cannot fall through to a mutable
source package.
