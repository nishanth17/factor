# Roadmap, research and reference tools

[TODOs](TODOS.md) records phased work and acceptance/experiment gates.
Research notes cover [QS/SIQS](quadratic_sieve_research.md),
[matrix algorithms](gf2_matrix_research.md),
[portfolio optimization](phase_two_optimization_research.md) and
[C sieve ideas](sieve_port_review.md).

The directory retains small source/citation manifests and exact source
snapshots needed by tests or benchmark comparisons. Milestone filenames on
those inputs identify their provenance; they are not disposable run output.

Generated validation captures, comparison samples, profiles, HTML/text audit
exports and detailed verification dumps stay local and are excluded from Git.
Historical references to them in the roadmap/research are shown as plain
labels. Public results and rerun commands are in
[the benchmark guide](../benchmarks/README.md); accepted behavior is summarized
in [the changelog](../../CHANGELOG.md).

## Reference checks

From the repository root, with PyPy implementing Python 3.11:

```sh
pypy3 v2/audit/test_prac_reference.py
pypy3 v2/audit/reduction_bench.py
```

PRAC remains experimental; production uses the checked ladder. The prototype
does not establish correctness for every exceptional projective pair.

Other historical audit scripts inspect the preserved v1 baseline through
`lib2to3` and can expose known v1 failures. Profile/analysis helpers may also
need local run captures. They are diagnostic tools, not the native v2
acceptance suite; use `make -C v2 test` for that suite.
