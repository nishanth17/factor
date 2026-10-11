# C10a-05 fixed collection pacing

Implementation/protocol freeze: `ea9b66cd588528d2fdfb661d25f1fa270fdcabd6`.
Training selection: `51f98fe`; fresh corpus: `96e2aad`.
Starting mainline: `master`, `16fe8a6180b3b86553cb2d2e56e7c27c9098e2ca`.
Work was isolated on managed branch `codex/c10a-05-pacing`; unrelated original
checkout edits and v1 were preserved.

See [research](research.md), [state invariants](design.md),
[correctness](correctness.md), and the immutable [protocol](protocol.md).
The one challenger retains the existing cadence/excess across shared SIQS
polynomial windows, preserves immediate zero-row checks and forces charged
processing at genuine collection stops. All other configuration fields match.
This does not repeat R3's cadence-32 experiment or implement C9a feedback.

**Decision:** adopt collection pacing only for the existing native 50-digit
automatic bundle at exactly 1 GiB owned memory. Retain window pacing at
40 digits and the reference/small/30-digit configurations. Keep explicit
configurations and concrete legacy checkpoints. The fresh 50-digit complete
cost gate passes; the 40-digit gate fails its uncertainty and stability tests.
This closes only C10a-05, not the C10a umbrella, C9a or E1. No earlier C1/C9/C10
percentage is compounded into this result.

## Training, before fresh inputs

Cost change is candidate/control − 1. Intervals are the frozen 4000-draw
hierarchical bootstrap, resampling inputs then their matched observations.
Each band has three distinct input numbers; seeds and repeat cells are not
additional independent numbers. All timed calls completed and validated.
The small regression products are distinct, but two share a prime factor;
they are arithmetic regression fixtures, not an independently sampled population.

| Digits | Pairs | Control total s | Collection total s | Cost change (95% interval) | Training decision |
| --- | ---: | ---: | ---: | --- | --- |
| 8 | 18 | 0.151056 | 0.160445 | +6.22% (−15.96%, +40.07%) | Retain; unstable |
| 30 | 18 | 7.118154 | 7.190944 | +1.02% (−6.66%, +14.18%) | Retain; unstable |
| 40 | 18 | 122.799365 | 120.232865 | −2.09% (−3.17%, +1.09%) | Advance by frozen cost ranking; unstable |
| 50 | 9 | 1150.598235 | 840.591088 | −26.94% (−27.72%, −26.13%) | Advance; stable |

The 40-digit first input was 1.31% slower in aggregate while the other two
were 2.68%/2.77% faster. The 30-digit third input was 15.46% slower despite
the other two improving. These are reasons against broad scheduling promotion.
Tiny 8-digit durations are particularly noisy; the declared extension was
retained and did not establish a speed benefit. No completion regression was
observed in these controls. They do not establish tail coverage.

The 50-digit per-input training improvements were 26.25%, 27.61% and 26.82%.
Its initial nine pairs passed stability; the other bands retained all eighteen
pairs after their mandatory extensions. Selection was committed before six
new confirmation inputs were generated. The generator excludes versioned
corpora and historical exclusions and uses independent prime certificates;
its balanced primes are conditioned on a large factor of p−1, not an
RSA-distributed population.

## Fresh confirmation and complete cost

| Digits | Pairs | Complete, each arm | Control total s | Collection total s | Cost change (95% interval) | Decision |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| 40 | 18 | 18/18 | 119.064695 | 117.252903 | −1.52% (−3.49%, +2.10%) | Retain window; unstable |
| 50 | 9 | 9/9 | 1005.917189 | 744.523235 | −25.99% (−28.24%, −22.87%) | Adopt scoped collection |

The 50-digit per-input improvements were 24.39%, 22.79% and 28.30%, with no
mixed paired signs or extension trigger. The 40-digit effects were −2.96%,
−1.03% and +2.21%; the first two inputs still had mixed paired signs after the
predeclared extension. Fewer solves do not rescue its failed complete-cost
gate. All 180 timed calls across 90 pairs completed correctly; no accepted
call was censored. Completion was equal, not improved. Zero correctness
failures remained; no samples or failed runs were discarded.

Nine 50-digit fresh pairs means three independent numbers with three seeds
each, not nine independent semiprimes. Training used 63 pairs and fresh
confirmation used 27. Intervals are conditional on these small, conditioned
input populations and this machine, not a universal size-band guarantee.
The declared five-percentage-point class completion regression limit was not
crossed. Other sizes, residual domains, backends and owned allowances were
not promoted.

Fresh stage sums below include unsuccessful attempts. Setup is factor-base,
assignment and root work; the remainder includes orchestration, classification,
validation and time outside the measured stage scopes. Inclusive parent
`collection_to_extraction` is deliberately not added a second time.

| Seconds, aggregate | 40 control | 40 collection | 50 control | 50 collection |
| --- | ---: | ---: | ---: | ---: |
| Collection | 105.816 | 109.588 | 647.864 | 653.941 |
| Preparation | 3.321 | 1.427 | 156.316 | 23.061 |
| Filtering | 3.101 | 0.517 | 73.733 | 8.162 |
| Elimination | 1.067 | 0.168 | 71.980 | 1.464 |
| Extraction | 0.069 | 0.038 | 0.160 | 0.103 |
| Setup | 1.098 | 1.076 | 10.280 | 11.885 |
| Remainder | 4.593 | 4.439 | 45.584 | 45.908 |
| Complete wall | 119.065 | 117.253 | 1005.917 | 744.523 |
| Complete CPU | 113.479 | 111.863 | 953.367 | 705.988 |

| Fresh aggregate or maximum | 40 control | 40 collection | 50 control | 50 collection |
| --- | ---: | ---: | ---: | ---: |
| Filter / solve calls | 1808 / 1790 | 226 / 36 | 5436 / 5427 | 479 / 9 |
| Solves without a proper factor | 1772 | 18 | 5418 | 0 |
| Trivial extraction trials | 246 | 30 | 18 | 8 |
| Final available kernel vectors | 72 | 780 | 40 | 504 |
| Proper factor outcomes | 18 | 18 | 9 | 9 |
| Polynomials | 2792 | 2904 | 7717 | 7821 |
| Collected positions | 365955816 | 379767622 | 1011490341 | 1024405124 |
| Retained rows, sum / peak | 5950 / 343 | 6664 / 389 | 14585 / 1642 | 15061 / 1697 |
| Peak partials | 5468 | 5634 | 13998 | 14066 |
| Peak owned MiB | 116.217 | 122.151 | 332.009 | 337.179 |
| Charged work | 192221130118 | 198494246200 | 826027713373 | 589678562032 |

Kernel vectors are available dependencies, not a count of useful factors;
extraction stops after a proper divisor. The 50-digit gain costs 104 extra
polynomials, 476 extra retained rows in aggregate and 5.17 MiB extra maximum
owned reservation. Additional collection is measured, not assumed free.
The shared-process RSS high-water maxima were 255.22 MiB at 40 digits and
251.69 MiB at 50 digits; these are not isolated per-arm memory measurements.
Owned reservations are conservative accounting bounds, not measured RSS.

Training resource totals provide the regression controls:

| Digits | CPU s, control → collection | Work, control → collection | Peak owned MiB, control → collection |
| --- | --- | --- | --- |
| 8 | 0.151 → 0.160 | 7883250 → 7884936 | 13.610 → 13.610 |
| 30 | 7.112 → 7.185 | 49404265752 → 53808402526 | 46.739 → 48.874 |
| 40 | 116.843 → 114.451 | 195277775302 → 200389389056 | 119.349 → 122.401 |
| 50 | 1091.473 → 797.938 | 893343509013 → 614581180134 | 329.437 → 332.909 |

## Attribution, replay and resource costs

The source diagnosis and complete feasibility traces are in [research](research.md).
On the inspected 50-digit control, all 640 changed-store polynomial-final
solves bypassed growth 32; the first observed nonempty dependency attempt
found the factor. Collection pacing used one solve after more collection.
The analogous 40-digit trace had 52 forced solves. This establishes repeated
work, not the earliest mathematically useful untested prefix. Independent
fixtures demonstrate useful below-surplus and immediate square dependencies;
they are the correctness oracle for bounded terminal opportunity.

Two separate 30-second profiles both stopped at `wall_limit` after 18
polynomials, retaining unresolved reconstruction. Their six-row prefixes
are too early for representative downstream attribution; they support no
speed claim. Complete boundary traces and the accepted stage counters supply
that attribution. The profiles, including censored outcomes, remain preserved.

The checkpoint diagnostic made two JSON round trips per arm after two-block
pauses, with a common 16 MiB encoding cap and cumulative 60-second grant.
Control payloads were 6559/6714 bytes; collection payloads 6656/6805 bytes.
Control encoding took 0.01904 s total and restore 0.10435 s; collection took
0.00470 s and 0.06892 s. Each arm charged 6,038,215 reconstruction work units;
total job work was 15,698,895 versus 15,698,912. Polynomial identity survived.
Encoding records zero abstract work but consumes wall/CPU time; rebuilding is
charged. These short, order-sensitive diagnostics are not performance claims
about checkpoint speed. Independent tests cover interruptions later inside
preparation, filtering, solver and extractor state.

Twelve separate cold subprocess probes consumed 47.111 s including startup.
Six completed (8/30-digit training and the first 40-digit training number);
six were censored at five seconds (50-digit and fresh 40-digit probes).
Startup-inclusive durations ranged from 0.424 to 5.431 s. They do not enter
the warmed comparison. Each of twelve band/arm warmups validated at least
three seconds of wall and CPU execution; their 601 calls consumed 85.529 s
including wrapper overhead. Many were tiny-input repetitions, not independent
observations. The 180 sample wrappers cost 4249.958 s, profiles 60.535 s,
checkpoint diagnostics 0.733 s and corpus generation 0.316 s.

The durable aggregate ledger charges **6244.182 s** against the frozen
14,400-second ceiling, including the conservative 1800-second prepayment for
all feasibility and QA, including retained failed preflight checks. No
allowance was reset. Sum-of-stage timings and accepted call timings are not
substituted for this inclusive resource ledger. Final QA below remains within
that prepaid reservation. Machine: Apple M4, 24 GiB, macOS 26.6.2 arm64;
accepted measurements and heavy QA held the machine-wide performance lock.

## Delivery validation

`make -C v2 test` passes **741 tests** in the worktree and in a committed-only
export of integration commit `f11f97d`. `make -C v2 lint` passes all three
checks. The committed export imports **176 live benchmark modules/packages**
and verifies A7 and the M8/M12/M25/C1/A6 historical loaders. Measured commit
`96e2aad` separately passes its strict frozen source/selection/corpus checks.
See [correctness](correctness.md) for independent witnesses and compatibility.
All subsequent closing changes are documentation only.

## Reproduction and evidence

Use PyPy implementing Python 3.11 (this study used PyPy 7.3.23 / Python 3.11.15).
Checkout the frozen study/corpus commit `96e2aad` to reproduce the measured
runtime and strict source pins. From the repository root, using an empty
ignored output directory and the versioned controls/corpus:

```sh
pypy3 -m v2.benchmarks.qs.c10a05.study train --digits 50 --lease-seconds 5400
pypy3 -m v2.benchmarks.qs.c10a05.study train --digits 40 --lease-seconds 5400
pypy3 -m v2.benchmarks.qs.c10a05.study train --digits 30 --lease-seconds 5400
pypy3 -m v2.benchmarks.qs.c10a05.study train --digits 8 --lease-seconds 5400
pypy3 -m v2.benchmarks.qs.c10a05.study confirm --digits 50 --lease-seconds 7200
pypy3 -m v2.benchmarks.qs.c10a05.study confirm --digits 40 --lease-seconds 7200
pypy3 -m v2.benchmarks.qs.c10a05.study diagnose
pypy3 -m v2.benchmarks.qs.c10a05.study summarize
```

The recorded training and confirmation each used one process across their
bands, in the order shown, with a fresh validated warmup for every band/arm.
Separate CLI leases are supported and repeat warmup; compare new results as
a new execution, not byte-identical timings. To exactly match process grouping,
call `study.measure(study.OUTPUT, phase, digits, lease_seconds)` for those bands
inside one `study.machine_window()` context after `study.c9.runtime()`.
The runner resumes partial groups from files and never resets its expense
ledger. Use a fresh output namespace for independent reproduction; never
remove inconvenient rows or charges.

Training per-call wall/CPU caps were 5/10/60/360 seconds for 8/30/40/50 digits;
fresh 40/50 caps are 60/480 seconds. These were chosen from feasibility, not a
universal 120- or 1200-second ceiling. Every call also has 10^13 work units,
matched owned storage and a 3 GiB RSS watchdog. A timeout is a censored result,
not proof that an input cannot be factored. The protocol fixes a 14,400-second
aggregate charge including prepaid QA/feasibility, with persisted reservations
for interruption. Cold probes, profiles and checkpoint rebuild diagnostics
are charged and reported separately from accepted complete-call timings.

Raw receipts, profiles, checkpoints and QA transcripts remain under ignored
`v2/benchmarks/results/c10a05/`. Required fixtures, certificates, configurations,
source pins and selection hashes are versioned under `benchmarks/inputs/`.
Historical C10 sources, captures and verdicts are unchanged. The new runner's
unique lease markers repair the prerequisite repeated-lease problem without
editing the frozen old runner. This does not close unrelated C10a tranches.

The 180-receipt SHA-256 index is preserved locally as
`experiment/accepted-receipt-sha256.txt`, digest
`c3816d722702b6210d2e1ba8e10bfc6a03208944e22c1ece9e0ca83815e7044f`.
The frozen summary digest is
`6ab8feab57a8ad8f11603f1a504103d53c3b9fff292d710132ed8e6ad8d399f5`;
the expense ledger digest is
`35688488f40288cd1c491ec5967dcc47c6a71c50892d2a8cb545824979640058`.
Production integration changes the default selector and A7 integration pin;
the strict historical study must be rerun at `96e2aad`, not silently repinned
to the delivery checkout.
