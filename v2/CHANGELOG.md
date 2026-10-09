# v2 development changes

## B2 paired ECM continuation — bounded opt-in tranche, 9 October 2026

- Add opt-in bounded reusable +/- stage-two execution through
  `PortfolioConfig.ecm_pair_distance`. Retain Python integers, the Montgomery
  ladder, unpaired streamed execution and existing production defaults.
- Keep immutable coverage records in the accepted finite program store and
  baby/giant points, term products and recovery in each curve's private job.
  Charge construction, decoding, arithmetic and recovery; preserve exact
  eligible primes, segmented tails and direct exceptions.
- Recover saturated paired terms by replaying both certified prime scalars.
  Validate every returned divisor. Add version-7 paired checkpoints, canonical
  int/GMP resume and predeclared finite campaign continuation with cumulative
  allowances. Preserve existing unpaired schemas and accounting.
- Freeze committed `94caf40` streamed/reusable controls, independent certified
  training/held-out inputs, D choices, budgets and selection criteria. Complete
  exclusive-window PyPy comparisons: 30 training arms, 15 held-out arms and
  135 separate cold starts. Extend four unstable training arms to 31 samples.
- Retain defaults after negative matched evidence. Selected paired execution
  is 81.3%, 154.2% and 45.0% slower than streamed complete factoring on the
  held-out small, medium and uneven cohorts; the nonsplitting campaign is
  15.0% slower. Completion is unchanged. Reusable unpaired programs also win.
  D choices 24/64/768/2048 are experimental selections, not recommendations.
  Segment boundaries eliminate pairing at the largest campaign choices;
  smaller D's product savings do not offset total overhead.
- Pass 367 system-PyPy tests (two optional skips), 370 PyPy/GMP tests with no
  skips and full `make -C v2 lint`, in the coordinated correctness window.
  Repeat both suites and import all 59 benchmark modules from committed files.
- Defer increased-B1 continuation because A6's exact schedule-ratio contract
  is absent. Adding curves or changing bounds on an exhausted checkpoint,
  wheel/common-Z, production PRAC, kernels and allocation remain separate.

This scoped changelog is under `v2/` to honor the B2 work boundary; the root
changelog and `v1/` are unchanged.
