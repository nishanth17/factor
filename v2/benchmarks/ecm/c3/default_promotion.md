# C3 default decision following the initial source screen

On 11 October 2026 the user explicitly requested bringing the observed
improvements into defaults. This supersedes the earlier defer decision in
[source-screen results](source_preset_results.md), whose frozen assessment and
raw observations remain unchanged. No new search, 18/27-sample extension or
fresh confirmation is initiated for this change.

Fresh bounded calls with omitted ECM tiers use compact64, `(2000, 50000, 64)`,
for exactly 30 decimal digits, and Yamaquasi's escalating ECM prefix,
`(200, 7700, 10), (2000, 81000, 30), (10000, 554000, 100)`, for exactly 40 digits.
Other sizes retain fixed32. This scope avoids extrapolating a threshold into
unmeasured size bands. Dispatch sees only the absolute original integer, never
its expected factors. It resolves once; recursive children retain that plan.

Compact64 has the lowest observed 30-digit CPU mean, 0.1367 seconds, versus
0.2245 for fixed32 and 0.1678 for escalating140. At 40 digits escalating140's
mean is 2.1747 seconds versus 3.9918 for fixed32 and 3.4840 for fitted64.
Across the complete screen its paired CPU reduction is 44.44%
[14.37%, 66.47%] versus fixed32 and 35.30% [1.52%, 59.49%] versus fitted64.
These intervals are conditional on six revealed subjects and do not correct
for selection among candidates. The hybrid itself was not measured as a
separate arm, and this document claims no fresh or broad hybrid speedup.

Alpertron115 has a slightly lower 40-digit point estimate, but its comparison
against fitted64 crosses zero. Its larger workspace also exceeds the default
CLI headroom with the existing executor settings. Escalating140 fits that
headroom, has positive lower bounds against both screen controls, and remains
available alongside Alpertron115 as an explicit preset. Native ECM families
and stage-two methods differ; upstream counts are numerical transfer inputs,
not equal-probability promises.

The implementation preserves every other configuration field and the total
work/wall/CPU grants. It does not create or enlarge fallback reserves. Existing
allocation floors still protect SIQS admission; a larger curve ceiling never
authorizes extra work. The larger schedule must coexist with unchanged SIQS
storage under the existing outer cap; otherwise automatic selection retains
fixed32. Explicit paired executors keep their prior plan. Explicit tiers,
including empty schedules and fixed32, override selection.

Only concrete numerical tiers enter checkpoints. Old schemas retain their
saved tiers, and new resumes do not rerun dispatch. Explicit conflicting
configurations still reject. The unbounded single-bound `factorize()` API is
unchanged; bounded API and CLI defaults use this decision.

The screen completed 54 matched groups / 324 valid complete calls with warmed
PyPy 3.11, nine seeds, fixed SIQS settings and 10^13-work / 5-or-30-second / 288
MiB grants. It did not meet its declared interval-width gate. Fresh inputs,
other structures/backends, default-grant economics and joint policy calibration
remain unconfirmed. C3/G1/E1 acceptance remains open. Historical manifests are
not repinned; reproduce the timing screen at `e73dbac`.

## Verification

The committed-only archive at `42075b7` passes all 663 tests, full lint, 154
benchmark module imports, 35 certified corpora and all 57 training fixtures
under PyPy 7.3.23 / Python 3.11.15. The existing A7 loader passes its separate
additive source-compatibility manifest; its explicit engine settings and old
manifests remain unchanged. Verification took 107.072 seconds under the shared
machine lock. A legacy test fixture now copies serialized config fields rather
than private routing metadata into the immutable mainline constructor; its
original resume/reconstruction/work assertions remain intact.

The focused checks cover exact decimal boundaries and sign, explicit/empty
schedules, the CLI's original 80 MiB / 64 MiB SIQS headroom, tight-memory
retention, protected fallback refusal, root-plan freezing, new/legacy numerical
resume and direct-SIQS CLI precedence. Reproduce the required checks from a
checkout containing only committed files with the project PyPy environment:

    make -C v2 test PYTHON=.venv/bin/python
    make -C v2 lint

Raw QA receipts and logs remain ignored under
`results/c3/default-promotion/`. No timing screen or fresh-input generation was
run during this promotion.
