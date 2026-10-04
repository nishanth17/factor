# Benchmarks

Run from the repository root on PyPy implementing Python 3.11:

```sh
make -C v2 benchmark WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-two WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-three-reference WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-three-collector WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-three-pipeline WARMUP_SECONDS=3 REPETITIONS=9
```

Use unique output names when supplying BENCHMARK_OUTPUT. Outputs remain
local and Git-ignored: raw JSON/gzip samples, stdout transcripts, profiles,
generated freeze files and validation/verification dumps. They are useful
working evidence, rather than source files to publish for every run.

Committed inputs include the independent Phase 2 corpora, competitor metadata,
exact old-source baselines in audit/, and the frozen P3.3 configuration used
by `phase_three_ownership.py`. Those snapshots are hash-checked before use.
The preserved v1 comparison is emulated Python 2 through `lib2to3`, not a
native Python 2 measurement.

## Recorded P3.3 result

A matched warmed comparison on 16 small held-out balanced inputs reduced
complete QS cohort median time from **27.220 to 21.446 ms**, about **21.2%**.
Both arms used the same frozen bucket recovery and filtering settings,
including setup through verified factor reconstruction.

The measurement used PyPy 7.3.23 / Python 3.11.15 on macOS arm64, default JIT,
at least three seconds of validated warmup, and nine samples. Training seed
329 preceded a configuration freeze; the final held-out seed was 335.
The fixed-cohort bootstrap median-ratio interval was 0.726–0.863.

Reproduce the final comparison with:

```sh
pypy3 -m v2.benchmarks.phase_three_ownership \
  --warmup-seconds 3 --repetitions 9 \
  --output v2/benchmarks/qs_comparison_LOCAL.json
```

These are 23–26-bit inputs, not evidence of large SIQS scalability. Cold
startup has no demonstrated gain. Tiny collector windows still trail
exhaustive enumeration. Resieving and tighter candidate scoring did not
improve complete-run time and remain optional.

Archive full raw results locally when needed. Publish selected summaries
with environment, commands, correctness checks and scope; avoid committing
a transcript or profiler dump for every development attempt.
