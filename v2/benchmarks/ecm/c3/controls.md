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

Protocol v3 clarifies the pre-generation smaller-factor bands: a 30-digit
composite cannot have a 16-digit smaller factor, so that band is 14 digits;
the 40-digit band remains 16 digits. No confirmation inputs existed. Pilot
v2 is retained as feasibility only, with unchanged arms and allowances.
Telemetry now distinguishes selected fallback from actual setup and records
active curves at a capped stop. These reporting/generation repairs precede
training selection and all accepted measurements. Original v2 JSON and prose
remain versioned. No policy values or candidate ranking rules changed.


Revision 4 preserves revision 3's JSON/prose before any training or fresh
inputs. Review found a zero-ceiling exact-power prerequisite, eager unused
context setup, missing resume-verification byte charges and fallback telemetry
that confused selecting a stage with executing it. Fix those boundaries and
validate new progress metadata. The fresh 30-digit final uneven factor band
is 14 digits (16 would exceed half the product size). Add envelope checks
before warmup calls and reject incomplete training selection. The candidate
bundles, selection rule, seeds, total resources and acceptance thresholds stay
fixed. The earlier cold pilot remains revision-2 feasibility evidence only.


The full-suite preflight also requires a new additive A7 integration receipt:
`a7_r5_c3_sources.json` chains to the immutable B3 receipt, pins the changed
CLI/portfolio and adapter plus the new allocation module, and leaves original
A7/B3 controls untouched. It validates current arm construction, not fresh E1
performance. A genuine pre-C3 baseline test excludes the new optional field
when constructing the historical dataclass. No collector or graph changes.


Revision 5 fixes implicit SSS/SSSf configuration restoration: the nested
collector must be decoded as SieveConfig in both Python and CLI restore.
New real both-mode handoff/resume coverage exercises this contract. No SIQS,
collector, curve, parameter, seed or selection rule changes. Stop v4 training
before selection or fresh generation, retain partial rows as diagnostic,
and charge 300 seconds (258.062 elapsed plus active-call/CPU margin) against
its original 2400-second allowance. Restart training under a 2100-second cap.
Original v4 source/prose/control and the A7 receipt stay immutable; a new
chained A7 v5 integration receipt pins the repaired source. Confirmation's
original finite allowance remains unchanged.
