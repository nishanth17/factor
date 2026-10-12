# CADO-NFS SIQS calibration trace for C10a

Reviewed at revision `692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b`.
This completes the source trace relevant to C10a-04's calibration questions,
following the initial shallow screen. It covers the active factorization
route, not every optional CADO-NFS algorithm. Native code was not executed;
no upstream performance or optimality claim was reproduced.

The [source manifest](../../inputs/controls/c10a04_upstream_review.json)
pins fetched files and tree-verified blobs. Source links below resolve at the
same revision. Captures and independent checks remain in ignored
`v2/benchmarks/results/c10a04/upstream-audit/`.

## Active caller and parameter path

| Stage | Executable path and result |
| --- | --- |
| Profile lookup | `scripts/cadofactor/toplevel.py::find_default_parameter_file` searches `params.qs.c<size>` for QS factorization, with nearby sizes at distance at most three below 200 digits. The shipped profiles are c40/c60/c70/c80/c90. Consequently there is no automatic 55- or 65-digit QS profile within that search radius; an explicit profile/configuration is required. These files are selected presets, not an executable autotuner. |
| Task composition | `cadotask.py::CompleteFactorization` selects `PolyselQSTask`, one-sided `FactorBaseTask`, `QuadraticSievingTask`, duplicate/Galois filtering, purge, merge, linear algebra and square root. QS forces Galois action `_y`. |
| Input polynomial | `PolyselQSTask.run` writes `x²-N` for factorization. This path applies no searched SIQS multiplier: its D is N. The discriminant adjustment for class groups is a distinct branch. |
| CLI and build | `cadoprograms.py::Siqs` maps q0/q1 to `qidx0/qidx1`; `sieve/CMakeLists.txt` builds `siqs` with `SIQS_SIEVE`. `las.cpp` binds `ALGO=SIQS`, then `do_one_special_q` selects the area, prepares roots/norms and sieves it. |

Sources: [profile lookup](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/scripts/cadofactor/toplevel.py),
[task wiring](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/scripts/cadofactor/cadotask.py).

## A products, Gray enumeration and normalized units

`siqs_todo_list::special_q_from_index` decodes a combinatorial index into
r distinct prime indices and multiplies them into q, corresponding to v2's A.
`get_qfac_prime` generates primes with Legendre(D,p)=1 in the prescribed
factor range. The active upper check is **inclusive**; a parameter-file comment
calling it exclusive is not the executable rule. Ordinary scheduling walks
indices; random sampling is optional. qrange counts A products per work unit,
not consecutive integer A values. The prescribed pool and factor count choose
A magnitude; this route has no corrective-factor search around an executable
sqrt(2N)/M target. Documentation's discriminant heuristic is not such a caller.

`siqs_special_q_data` constructs centered CRT terms. Its Gray root for j is
B_j; the coordinate conversion gives U=q i+B_j and b=1.
`lognorm_base::norm<siqs_special_q_data>` evaluates the exact norm
F_j=((q i+B_j)²-N)/q, including the j-dependent correction.
Thus U²-N=q F_j, the classic representation used by v2, with D=N here.
CADO-NFS centers individual CRT terms; v2 recenters each resulting B modulo A.
Those choices can shift the sampled integer window even with the same residue.

`siever_config::parse_default<SIQS>` sets logA=I+(r-1).
The SIQS `choose_sieve_area_impl` recovers logI=I and sets J=2^(r-1).
This visits 2^I positions for each of all J Gray variants, for
2^(I+r-1) positions per A. **logA here describes area, not log₂ of coefficient A.**
The equivalent half-width is 2^(I-1), with one endpoint convention differing
from v2's 2M+1 positions.

Sources: [A enumeration](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/las-todo-list.cpp),
[CRT/Gray](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/las-qlattice.hpp),
[area](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/las-choose-sieve-area.cpp),
[norm](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/las-norms.cpp).

## Presets, residuals and usable relations

| Profile | Prime endpoint | Half-width | r | Gray variants | Positions/A | LP bits / cofactor bits |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| c40 | 7500 | 4096 | 6 | 32 | 262144 | 15 / 30 |
| c60 | 50000 | 8192 | 8 | 128 | 2097152 | 20 / 42 |
| c70 | 150000 | 16384 | 10 | 512 | 16777216 | 21 / 42 |
| c80 | 100000 | 2048 | 12 | 2048 | 8388608 | 19 / 38 |

The [presets](https://github.com/cado-nfs/cado-nfs/tree/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/parameters/factor)
are visibly nonmonotone. They supply hypotheses, not demonstrated PyPy scaling.

`makefb_with_powers` visits primes p<=lim and stores polynomial roots and
lifted prime powers. `FactorBaseTask` supplies maxbits from I; `fb_factorbase::read`
applies lim/powlim. The endpoint is neither eligible-prime cardinality nor matrix
column count. Independently counting distinct eligible primes, including 2
and excluding power/root multiplicity, gives 2567/2544/2533 for the three
frozen 60-digit inputs at D=N and lim=50000. Their actual endpoints are
49957/49999/49993. This is an interpreted count, not a native matrix dump.

`check_leftover_norm` applies mfb/lpb bit limits and impossible-product/primality
checks before facul cofactorization. In the ordinary resieve route,
A-dividing primes omitted from sieving are divided explicitly from F;
`cofac_standalone::get_relation` then reattaches q's known prime factors.
Duplicate/Galois filtering and purge determine usable excess; purge shortage
propagates back through task notifications to request more relations. Raw
relations and the nominal base size therefore do not specify a usable matrix.
Bucket allocations, threaded scratch and external relation files also differ
from v2's owned-memory/receipt contract; none are portable capacity constants.

## Exact checks and implications

Independent PyPy 3.11 checks used N=77 and exhaustive small CRT families:
35 A products, 1540 positions, 1120 nonsingular root-transform cases and
84 combinatorial indices. Reconstruction, modulo-sign uniqueness and the
index-ranking identity passed. These validate this interpretation, not the
native implementation. The ignored check script and JSON preserve reproduction.

The v2 caller takes indices 0 through gray_limit-1 from the ordinary Gray walk.
Quota 64 is a **prefix**, not a demonstrated diverse selection. With r=8/9,
some free CRT signs do not vary in that prefix. This corrects the terminology
in the hashed Stage 1 research without editing that frozen document.

The source trace reinforces a concrete missing question: smaller intervals
with fuller within-A reuse may amortize differently from the frozen compact/
wide bundles. It does not prove they win. C10a-07 now owns source reconciliation,
quota attribution and conditional traversal-policy experiments. C10a-02 retains
selection of different A products; A1a owns transition kernels and C8 the sorted
CRT scheduler. C10a-04's three arms and acceptance identities remain unchanged.
