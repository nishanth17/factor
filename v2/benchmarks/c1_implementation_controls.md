# C1 implementation controls and accounting notes

Before training, preserve owned production source from `ea22953` in the
hash-checked `inputs/baselines/c1_slp_baseline.json`. This is the calibrated
SLP execution control, including its original relation layout and dispatch.
The five frozen arms and allowances in `c1_implementation_protocol.md` stay
unchanged. Training and confirmation load the control into a separate package;
its exception/budget classes stay together. No upstream licensed implementation
is loaded into an experiment.

The first training input uses seed7 and the second seed29. Each arm/band gets
at least three seconds of validated warmup before its one run per input.
Confirmation rotates seeds7/29/47, pairs sample indices, and extends both arms
of an unstable input together. Cold processes are separate captures. Timer
end follows output reconstruction, expected-factor and certainty validation.
Failure cost is the assigned wall cap, with actual wall/CPU, reason, factors
and unresolved cofactors retained alongside it. Corpus certification remains
outside factoring time and is never supplied to the algorithm.

The graph reservation is 32KiB plus4KiB per retained forest edge. It covers at
most two vertices per edge, two coexisting map sets during staged rebuilding,
parent/component/depth/size maps, adjacency dictionaries and traversal/update
scratch. Labels are at most40bits; atomic IDs are bounded shared strings.
Atoms and combined rows have their own reservations; checkpoint JSON has the
existing separate encoding allowance. Evicted atoms remain counted until the
transaction commits. Reservation bytes are an owned-storage bound, not RSS.
Long paths exceeding256atoms are reported yield losses. Split calls remain
charged even if cancellation prevents publishing their atom. Quota exhaustion
rejects further composite residuals but still permits full/SLP collection.

A cycle contains one new closing atom and its unique forest path. Pinning
those path atoms means later FIFO eviction cannot break an emitted cycle.
Every admitted closing atom is unique; therefore the fundamental cycles are
independent in the retained LP incidence kernel. This applies in all components
and to self-loops/parallel edges. The independent exhaustive boolean oracle
checks the complete span of small retained multigraphs, including eviction.
The production cycle-length cap deliberately limits that completeness claim.

Opt-in support is for serial `SieveCollector`/`SIQSJob`, and portfolio dispatch
through the same explicit SIQS configuration. It uses SIQS checkpoint version4;
SLP remains version3. No shared portfolio or stage-job schema changes are made.
SSS and worker exporters reject the DLP configuration before collection.
