# C10a-01 prospective calibration and pilot protocol

Base: 4a946b8d6294d51ac5d9c36f61880cc63d5b59fd. Default arithmetic is the
literal control. The only nominee is parity2 polynomials with parity2 scoring.
All later pilots are design only in the implementation session. No defaults,
merges, pushes or E1 work are authorized by this protocol.

## Bootstrap allocation, frozen before numerical execution

`bootstrap.json` is the numerical authority. Retained 65-digit work was
432808207231 units over approximately 1800 wall seconds, about 247 million
units per CPU second. Use a deliberately larger initial sustained planning
rate of 2 billion units/s, multiplied by 4 for operational headroom, plus
10 billion units per call for setup/atomic costs omitted by that rate. This
is a finite guard, not a guarantee that time binds first. Setup count sieving,
root construction and collection charges are read from current code; later
matrix costs will be bounded independently before pilot freeze.

Prefix work: 4 * (250 + 370 + 490 + 730) billion = 7.36 trillion.
Diagnostic work: 18 * 130 + 4 * 370 + 2 * 370 billion = 4.56 trillion.
Warmup: 8 billion * 900 + 30 * 10 billion = 7.5 trillion, at most 30
compatible process/configuration warmup reservations. Generation/proof work
proxy: 100 billion; analysis: 1 trillion. Total: 20.52 trillion units.
No work limit is increased inside a call. Refusals and interrupted reservations
remain in the ledger. Generator work is explicitly a conservative bounded
proxy, not production Budget units.

Time reservations: prefixes 980; observer 360; forced policies 200; memory
100; warmup 900; generation/proof 600; analysis 120; orchestration 340.
Total 3600 wall AND CPU seconds. Each measured call includes a separate
5-second reservation margin. Ordinary development and correctness QA are
separate, with heavy QA and all performance work serialized by
`/private/tmp/factor-performance.lock`.

Generate eight exposed semiprimes once with seed 202611020100, using the
existing recursive Pocklington generator and independent verifier. Balanced
and approximately 35:65 factor-bit splits have N mod 8 respectively 1/3,
5/7, 3/1, 7/5 at 60, 62, 65, 70 digits. Reject only size, residue, equal or
already-used factors; never reject based on factoring performance. Candidate
support follows from k = 2r in p = kq + 1: consecutive r span every odd
residue mod 8 for odd q; intervals contain at least four consecutive r.
Primality and the recursive certificate condition further condition sampling;
these are not uniform RSA inputs. Maximum 10000 pair attempts per cell,
1000000 random draws, 590 wall/CPU seconds, and 100 billion proxy units.
Known factors and certificates are supplied only to validation.

Before prefixes: freeze literal configurations, exact position counts and
block order. Hold seed, capacities, LP domain, scoring/refinement, backend,
block/interval and solver policy constant except declared representation and
multiplier changes. Different A populations are intentional. Norm samples,
physical U coverage and counters describe those differences. Observer on/off
uses identical semantic work and checks passive state equality. A call that
hits its cap before its position target remains incomplete fixed work.

## Later complete-factor stages: not executed by calibration

Stage A: balanced 60 and 62, two arms each. Caps 1800 and 2400 seconds;
4 measured reservations total 8420; warmup 600, analysis 120, orchestration
60: 9200 wall/CPU seconds.
Stage B: balanced 65, two arms, 3600 seconds each; 7210 measured + 300
warmup + 60 analysis + 30 orchestration = 7600 seconds.
Stage C: balanced 70, two arms, 7200 seconds each; 14410 measured + 300
warmup + 60 analysis + 30 orchestration = 14800 seconds.
All stages plus calibration: 35200 seconds = 9 hours 46 minutes 40 seconds.
Each stage needs its own exact arithmetic, parameter reachability and resource
admission; 65-digit completion is not a prerequisite for 70-digit admission.
Numerical post-calibration work limits and memory tiers must be frozen before
any pilot command can execute. No automatic additional pairs or larger retry.

Complete cost includes multiplier/base/A setup, collection, preparation,
filtering, solving, extraction and final classification/reconstruction.
Cold import/startup is separate; warm each arm/configuration/process for at
least 3 validated wall AND CPU seconds. Preserve every failed/censored call.
Successful/successful pairs permit descriptive elapsed-cost ratios only;
otherwise compare finite-budget completion and capped service cost (cap for
unfinished jobs), never successful-factor speed against a failure. One modulus
per pilot size cannot establish general speedup or input-level uncertainty.

A nominee requires a new confirmation identity, selection committed before
revealing disjoint inputs, at least nine independent moduli per claimed band,
matched resources, input-level paired uncertainty and fresh repeatability.
Repeated seeds are not independent moduli. No universal percentage threshold.
Further redesign uses exploratory evidence and cannot reuse confirmation as
untouched data. Extension is not authorized automatically.
