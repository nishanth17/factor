# ECM-local two-curve NEON feasibility prototype

This optional prototype is separate from the scalar C/GMP and OpenMP executors.
It is never selected by factorize_bounded, the production packet, or QS. It
supports exactly 200/266/333-bit unsigned operands and two independent curve lanes.
Necessary local packing/carries are restricted to this screen. Reduction is the
existing ordinary GMP mpz_mod operation; no new reduction context is provided.

Build explicitly with PyPy 3.11 on ARM64 macOS:

```sh
pypy3 -B -m v2.ecm.native.neon_build --gmp-prefix /opt/homebrew/opt/gmp
```

Build products and disassembly live in ignored native/build/neon. GMP/CFFI and
the compiler remain explicit dependencies. Source, binary, dependency and flag
identities are recorded. No LTO and no automatic loop/SLP vectorization are used
in the wrapper/scalar control. Admission inspects the explicit NEON symbol for
vector widening products and verifies that its scalar counterpart lacks them.
The prebuilt GMP library's own implementation remains unchanged.

The radix is 2^26, with 8/11/13 digits respectively. Each NEON uint64 lane sums a
Comba product diagonal using vmull_u32. At most 13 products plus a carry<13*2^26
fit below13*2^52<2^56. Carry propagation produces canonical 26-bit digits and the
final carry fits one digit. Highest partial input digits and single-curve tails
are checked against full GMP products. The squaring candidate uses the same
kernel with equal operands; it has no specialized symmetry optimization. The
GMP square control uses its specialized mpn_sqr.

| Endpoint / arm | Timed work |
| --- | --- |
| arithmetic-gmp-mpn | Two mpn_mul_n or mpn_sqr calls on prepacked 64-bit limbs and disjoint preallocated outputs |
| arithmetic-scalar26 | Two scalar products with exactly the NEON layout, diagonal order and carries; auto-vectorization disabled |
| arithmetic-neon26 | One two-lane NEON product/square on prepacked radix26 input, preallocated output |
| complete-production-mpz | Two reusable-object mpz_mul followed by mpz_mod, matching the actual B4a production operation pattern; equal operands select squaring |
| complete-hybrid-scalar26 | Export, interleave, scalar26 product, deinterleave, leading-zero normalization, import and GMP reduction |
| complete-hybrid-neon26 | The identical hybrid path with the explicit NEON product |
| layout-only-diagnostic | Identical packing/import/normalization/reduction around a precomputed correct product; not a candidate multiplication or a subtractable exact overhead |

The scalar26-versus-NEON26 comparison isolates SIMD under the same representation
and surrounding work. Comparing either to GMP also changes radix/algorithm and
bookkeeping. A hybrid scalar win is not automatically a SIMD win. Arithmetic-only
and complete-hybrid summaries are separate; a multiplication-only win cannot
justify executor integration. All inputs are nonnegative and width-bounded;
future integration would need to handle actual point intermediates and signs.
No arbitrary-width arithmetic or whole-factor speed claim follows.

Timing runs only as the separately admitted Stage A screen. Its fixed 64-pair
panel uses varied uniformly drawn canonical residues modulo 2^width-159. The
correctness lane separately tests zeros/ones, maximal values, carries, partial
limbs, random products/squares, malformed canonical inputs and inactive tails.
All objects/buffers are retained across useful batches. Every pair enters the
same opaque noinline sink in a separate translation unit, with a checked low26
checksum update; no LTO permits cross-unit elimination/hoisting. Full products or
residues for every distinct panel entry and the complete batch checksum are
verified against GMP outside the timed region. Repeated panel offsets rotate.

Each endpoint group uses identical operation counts/repetitions and input order.
Calibration targets at least 10 ms for the fastest matched batch and at most 100 ms
for the slowest; max 65,536 repetitions, native cooperative cap 250 ms. Interleave
arms and retain calibrated and warmup rows, then 9 measured batches/arm. Warmup is
at least 3 seconds per arm of the validated PyPy driver path, including its checks;
report that wall time honestly, rather than labeling it 3 seconds of C-only work.
Raw JSONL records relative order/time, every batch's wall/process CPU/checksum,
operation count and QoS result. Summaries report median/min/max/MAD and paired
complete-path uncertainty. These are exploratory batch observations, not
independent factoring moduli or promotion evidence.

Before each batch, on the single arithmetic thread, call exactly
pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0). Record its result,
pthread_get_qos_class_np readback, relative priority and thread ID. No privileged
settings or affinity changes occur. Successful QoS is not a P-core/E-core pin.
Failed/inconsistent requests remain visible and prevent a clean follow-up signal.

Untimed admission uses neon_qa, which builds, checks exactness/code generation and
runs lint under its own finite shared-lock grant. It never invokes timed batches
or benchmark warmup. Ordinary scalar/OpenMP QA does not require this optional
prototype; future Stage Q --with-neon explicitly enables its additional admission
and committed-source build checks. Only a frozen v3 Stage A cell may run timing.
A negative complete-path result closes this prototype without blocking other
B4a arms or ruling out every SIMD design.

The authorized 4 November screen is complete on source c485f2b. All six complete
hybrid cases lose to the B4a scalar GMP path by 2.27–3.14x. This prototype is
closed without integration; scalar/OpenMP evaluation remains independent. See
[recorded results](../../../benchmarks/RESULTS.md#b4a-optional-ecm-neon-screen-4-november-2026).

Arithmetic alone also loses to GMP mpn at every width/operation (1.53–3.54x
cost). ECM SIMD is parked; scalar B4a/native threads remain the ECM focus and
NEON work remains under B5a QS scans/root updates.
