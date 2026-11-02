# V1 baseline

The [v1 directory](../../v1/) preserves the original Python 2 implementation.
Its files, including README.md, are protected by the baseline hash manifest.
They remain unchanged. Use [Factor v2](../README.md) for current usage,
repairs, bounded execution and supported runtime instructions.

## Source map

| Original file | Component |
| --- | --- |
| [factor.py](../../v1/factor.py) | Interactive factoring dispatcher |
| [pollardPm1.py](../../v1/pollardPm1.py) | Two-stage Pollard p−1 |
| [pollardRho.py](../../v1/pollardRho.py) | Brent-style Pollard rho |
| [ecm.py](../../v1/ecm.py) | Montgomery ECM, Suyama setup and original PRAC helper |
| [primeSieve.py](../../v1/primeSieve.py) | Prime-sieve implementations |
| [utils.py](../../v1/utils.py) | Arithmetic and primality helpers |

## Historical usage and evidence

The original entry point is `factor.py`. Its syntax, arithmetic behavior and
API belong to the historical implementation. The [original README](../../v1/README.md)
retains usage transcripts and references. Its “50–60 digits within a minute”
statement is not a supported performance guarantee for current Factor.

Compatibility comparisons use a separately labelled adapter and validate
outputs. See [benchmark results](../benchmarks/RESULTS.md#historical-v1-comparisons)
for the measured scope. Documentation cleanup does not rewrite the baseline
or weaken its hash check.
