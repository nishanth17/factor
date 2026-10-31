# A1a frozen execution results — 31 October 2026

Subsequent [user-directed adoption and cleanup](a1a_default_promotion.md)
enables the integration at 50–55 digits. The decisions below remain the
original frozen results. Execution-time paths identify the original worktrees;
their ignored receipts and journals were archived locally before removal.

The selected integration passes its frozen 50-digit gate. Its 55-digit result remains inconclusive because repeated-cell drift persists after the one permitted full extension. The combined selected-power-period/four-stride extension loses at both sizes; retain its prior integrated control. These decisions authorize no automatic defaults or portfolio/E1 claim.

| Comparison | Digits | Mean complete wall seconds, control → selected | Change and 99.5% interval | Exact one-sided p | Decision |
| --- | ---: | ---: | --- | ---: | --- |
| Integration | 50 | 46.538 → 43.742 | -6.008% [-7.795% to -4.373%] | 0.001953125 | Pass, scoped |
| Integration | 55 | 227.993 → 213.650 | -6.291% [-7.575% to -5.310%] | 0.001953125 | Inconclusive; retain control (drift) |
| Extension | 50 | 81.748 → 86.068 | +5.285% [+4.008% to +6.245%] | 1.000000000 | Retain control (higher cost) |
| Extension | 55 | 270.654 → 278.758 | +2.994% [+2.432% to +3.501%] | 1.000000000 | Retain control (higher cost) |

## Endpoint and interpretation

Each band has nine independent certified fresh inputs, seed 7, and two opposite-order repetitions per arm. Integration at 55 digits retains four repetitions per arm after the required full stability extension. Repetitions are averaged within input; the table reports the equally weighted input means. Absolute seconds and differences are observed costs. Uncertainty is the frozen 50,000-draw input bootstrap interval for the total cost ratio change; no additional absolute-cost interval or inferential look was added.

All 180 measured calls completed and reconstructed exactly. Proper divisors, matching deterministic paths, runtime certainty labels, work and memory bounds were checked; fixture certificates did not upgrade probable-prime API labels. There are no measured timeouts, censored results, interrupted stages, failed reservations or overhead overruns. All process warmups passed the validated three wall-and-CPU-second minimum and stable semantically matching prefix criteria. Base repeated-cell wall and CPU drift had to be at most 15%; only integration at 55 digits required an extension. All-data drift still failed there, so its negative cost estimate is descriptive and is not an accepted speed improvement.

These are warmed complete-factor costs, including configuration/job/count setup, factoring, proper-divisor and reconstruction checks, certainty classification and explicit teardown. Cold imports/startup and instrumented profiles are outside the endpoint. Native exact 50/55-digit SIQS uses a 1 GiB owned cap, 3 GiB cooperative RSS cap, half-width 65536, block 4096, list/powers/bucket collection and a 2 MiB power plan. The 50-digit bundle has bound 30000, multiplier 1, six A factors, pool 32 and 32 actual Gray variants; 55 digits counts exactly 2500 eligible entries under ceiling 250000, with selected multiplier, seven A factors, pool 64 and 64 actual variants. Per-call work is capped at 10^13, wall and CPU at 480/1159 seconds, respectively.

The conditioned p−1 prime populations do not represent arbitrary semiprimes. Inputs, controls and charges remain separate. Do not pool these studies, chain their percentages, compare absolute costs across their different panels/time periods, infer component effects, or infer extended-versus-pre-integration-mainline performance. Forty digits retains control; historical clean 60-digit evidence remains inconclusive and 65-digit completion remains unproved. A direct combined comparison and any default adoption require separate authorization and evidence.

## Experiment 1: integration

Control is owned pre-integration mainline `92543797f4ed3fae347ecf5a3e938361d7ec7047`; selected is `density+presieve+scan_fused`. The historical isolated 5.945%/3.214% results are context only.

### 50 digits

36/36 calls completed; 18 matched pairs, nine independent inputs. Observed mean wall difference is -2.796 seconds per factorization; mean CPU is 45.294 → 42.502 seconds. The base schedule passed stability; no extension was required.

Per-input costs average all retained repetitions. IDs below are suffixes of `a1a_integration_50_`. Drift is range/median in percent; C/S means control/selected. Every listed repetition completed.

| Input | Wall seconds C/S | CPU seconds C/S | Wall change | Wall drift % C/S | CPU drift % C/S |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 30.260/27.872 | 29.346/27.075 | -7.894% | 9.00/8.48 | 7.80/7.81 |
| 1 | 45.416/42.931 | 44.158/41.730 | -5.472% | 5.10/1.15 | 4.32/0.69 |
| 2 | 37.462/35.459 | 36.523/34.551 | -5.346% | 3.66/2.54 | 3.23/2.25 |
| 3 | 45.487/41.964 | 44.354/40.816 | -7.745% | 2.00/5.57 | 1.90/4.75 |
| 4 | 39.509/35.914 | 38.449/34.993 | -9.097% | 3.35/1.29 | 2.77/1.22 |
| 5 | 48.598/45.814 | 47.390/44.531 | -5.728% | 1.85/0.51 | 1.51/0.32 |
| 6 | 49.658/46.446 | 48.170/45.078 | -6.469% | 1.57/1.61 | 1.48/1.24 |
| 7 | 86.312/83.410 | 84.013/80.764 | -3.362% | 4.18/0.46 | 3.76/0.42 |
| 8 | 36.137/33.863 | 35.242/32.980 | -6.292% | 0.05/2.03 | 0.05/1.62 |

| Arm | Mean work/call | Maximum owned bytes | Maximum RSS bytes |
| --- | ---: | ---: | ---: |
| control | 69312216774.667 | 355723176 | 280297472 |
| selected | 69619733633.111 | 367650728 | 280281088 |

Owned workspace is the runtime accounting measure; RSS is the cooperative process high-water observation, not a fresh per-call allocation peak.

### 55 digits

72/72 calls completed; 36 matched pairs, nine independent inputs. Observed mean wall difference is -14.342 seconds per factorization; mean CPU is 218.341 → 204.631 seconds. One full extension retained every base observation. Stability failed in all 18 input/arm cells on both wall and CPU; no further extension is permitted.

Per-input costs average all retained repetitions. IDs below are suffixes of `a1a_integration_55_`. Drift is range/median in percent; C/S means control/selected. Every listed repetition completed.

| Input | Wall seconds C/S | CPU seconds C/S | Wall change | Wall drift % C/S | CPU drift % C/S |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 258.664/235.594 | 248.359/226.292 | -8.919% | 66.05/71.41 | 64.51/69.63 |
| 1 | 201.077/186.939 | 192.415/178.992 | -7.031% | 51.79/53.47 | 50.86/52.28 |
| 2 | 247.361/233.773 | 236.840/223.741 | -5.493% | 52.97/52.46 | 51.89/51.48 |
| 3 | 252.910/236.370 | 242.171/226.301 | -6.540% | 52.81/52.70 | 51.81/51.76 |
| 4 | 208.811/198.598 | 199.859/190.167 | -4.891% | 52.19/53.59 | 51.28/52.46 |
| 5 | 202.415/191.842 | 193.872/183.750 | -5.224% | 53.21/53.06 | 52.07/52.05 |
| 6 | 210.229/196.425 | 201.318/188.091 | -6.566% | 53.14/52.92 | 52.04/51.85 |
| 7 | 226.493/211.879 | 216.737/202.764 | -6.452% | 51.38/51.54 | 50.52/50.66 |
| 8 | 243.975/231.435 | 233.493/221.583 | -5.140% | 52.13/52.89 | 51.15/51.86 |

| Arm | Mean work/call | Maximum owned bytes | Maximum RSS bytes |
| --- | ---: | ---: | ---: |
| control | 79720474601.556 | 429532064 | 280707072 |
| selected | 80325494501.889 | 445129632 | 280707072 |

Owned workspace is the runtime accounting measure; RSS is the cooperative process high-water observation, not a fresh per-call allocation peak.

### Exact identities and retained accounts

Prepared source: `fee6a41c9af2702c34f964536e71ad3011d4cb22` (261 frozen source bindings).

- Control SHA-256: `707a22d6b08500efdf250b51e04ead9f859b2f134f642b24002c00735848e1e5`
- Protocol SHA-256: `a13f5c2d5f8c1ab07ec2a7b7a97f030464b19d44de6ba8394adc5a5a1c195c2e`
- Corpus SHA-256: `c688817f8020c3645c0ec1a5f757f7d26ca3aaa06c6c060955dca706479b9f71`
- Baseline SHA-256: `e0ce2f9fc837ff2d80c118fb94b3fbf617420985bdf6efbdcfef5d1168a361d4`

Checked 140 hashed ledger receipts; output including the owned replay package occupies 1906251 of 134217728 allowed bytes. Source/control/corpus pins and both prepared checkouts remain unchanged. Full successful time reservations remain charged; measured work is charged as recorded. Remaining amounts are protected ceilings, not permission to repeat completed stages or transfer between phases.

| Protected phase | Charged wall / CPU seconds | Remaining wall / CPU seconds | Charged work | Remaining work |
| --- | ---: | ---: | ---: | ---: |
| analysis | 600.000 / 600.000 | 0.000 / 0.000 | 0 | 0 |
| base50 | 17460.000 / 17460.000 | 0.000 / 0.000 | 2500775107340 | 357499224892660 |
| base55 | 41904.000 / 41904.000 | 0.000 / 0.000 | 2880827443862 | 357119172556138 |
| extension50 | 0.000 / 0.000 | 17460.000 / 17460.000 | 0 | 360000000000000 |
| extension55 | 41904.000 / 41904.000 | 0.000 / 0.000 | 2880827443862 | 357119172556138 |
| generation | 120.000 / 120.000 | 0.000 / 0.000 | 0 | 0 |
| orchestration | 180.000 / 180.000 | 420.000 / 420.000 | 0 | 0 |
| total | 103078.000 / 103078.000 | 22570.000 / 22570.000 | 8296435038238 | 3031703564961762 |
| warm50 | 280.000 / 280.000 | 2520.000 / 2520.000 | 19496871904 | 799980503128096 |
| warm55 | 630.000 / 630.000 | 2170.000 / 2170.000 | 14508171270 | 799985491828730 |

## Experiment 2: extension

Control is owned integration `fee6a41c9af2702c34f964536e71ad3011d4cb22` with `density+presieve+scan_fused`; selected is the combined `density+presieve_powers+unroll4+scan_fused`. Both tested bands retain that control. This does not reject either component separately.

### 50 digits

36/36 calls completed; 18 matched pairs, nine independent inputs. Observed mean wall difference is +4.320 seconds per factorization; mean CPU is 78.082 → 82.200 seconds. The base schedule passed stability; no extension was required.

Per-input costs average all retained repetitions. IDs below are suffixes of `a1a_extension_50_`. Drift is range/median in percent; C/S means control/selected. Every listed repetition completed.

| Input | Wall seconds C/S | CPU seconds C/S | Wall change | Wall drift % C/S | CPU drift % C/S |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 62.412/65.250 | 59.601/62.340 | +4.547% | 1.75/0.63 | 1.71/0.69 |
| 1 | 136.033/144.873 | 129.983/138.352 | +6.499% | 0.39/0.53 | 0.25/0.41 |
| 2 | 62.838/66.489 | 60.001/63.447 | +5.810% | 1.12/0.34 | 1.05/0.30 |
| 3 | 83.115/87.992 | 79.447/84.088 | +5.868% | 0.45/0.15 | 0.47/0.13 |
| 4 | 73.555/75.872 | 70.139/72.444 | +3.151% | 0.05/1.27 | 0.20/1.09 |
| 5 | 103.446/110.247 | 98.862/105.316 | +6.574% | 0.36/0.31 | 0.44/0.23 |
| 6 | 64.199/67.176 | 61.286/64.207 | +4.637% | 1.51/0.71 | 1.48/0.72 |
| 7 | 53.592/56.682 | 51.222/54.128 | +5.767% | 0.53/0.11 | 0.51/0.07 |
| 8 | 96.540/100.033 | 92.195/95.477 | +3.618% | 0.97/0.41 | 0.97/0.24 |

| Arm | Mean work/call | Maximum owned bytes | Maximum RSS bytes |
| --- | ---: | ---: | ---: |
| control | 68044601151.333 | 363495744 | 273547264 |
| selected | 68045597304.000 | 367755584 | 273301504 |

Owned workspace is the runtime accounting measure; RSS is the cooperative process high-water observation, not a fresh per-call allocation peak.

### 55 digits

36/36 calls completed; 18 matched pairs, nine independent inputs. Observed mean wall difference is +8.104 seconds per factorization; mean CPU is 258.416 → 266.125 seconds. The base schedule passed stability; no extension was required.

Per-input costs average all retained repetitions. IDs below are suffixes of `a1a_extension_55_`. Drift is range/median in percent; C/S means control/selected. Every listed repetition completed.

| Input | Wall seconds C/S | CPU seconds C/S | Wall change | Wall drift % C/S | CPU drift % C/S |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 220.289/227.268 | 210.299/216.984 | +3.168% | 0.24/0.67 | 0.23/0.66 |
| 1 | 342.044/353.374 | 326.540/337.411 | +3.312% | 0.23/0.04 | 0.20/0.04 |
| 2 | 259.101/268.763 | 247.405/256.563 | +3.729% | 0.43/0.18 | 0.40/0.09 |
| 3 | 252.842/262.576 | 241.435/250.586 | +3.850% | 0.01/0.53 | 0.02/0.48 |
| 4 | 275.999/283.897 | 263.535/271.094 | +2.862% | 0.26/0.19 | 0.26/0.21 |
| 5 | 214.774/219.117 | 205.077/209.165 | +2.022% | 0.21/0.27 | 0.07/0.27 |
| 6 | 349.170/356.704 | 333.381/340.682 | +2.158% | 0.03/0.12 | 0.02/0.12 |
| 7 | 306.925/316.477 | 293.102/302.058 | +3.112% | 0.00/0.38 | 0.04/0.32 |
| 8 | 214.746/220.645 | 204.968/210.580 | +2.747% | 0.16/0.34 | 0.17/0.39 |

| Arm | Mean work/call | Maximum owned bytes | Maximum RSS bytes |
| --- | ---: | ---: | ---: |
| control | 83576972311.778 | 452801616 | 293011456 |
| selected | 83588426994.222 | 457061456 | 292995072 |

Owned workspace is the runtime accounting measure; RSS is the cooperative process high-water observation, not a fresh per-call allocation peak.

### Exact identities and retained accounts

Prepared source: `80d8cda1b11f2ca2724d428078c206930999e1ac` (268 frozen source bindings).

- Control SHA-256: `ac461c7178e6f83b1987ffd0d52a2deac36c45dbf996f8a686818a8bce582bb4`
- Protocol SHA-256: `ac505f5ce8df777d177dc53557d574e47971d3a5502c2b223ee12441c466f560`
- Corpus SHA-256: `fc66410532754e31ecd138cbffa9f5d4fa3df6b6fd3ef30163ac8673094140c2`
- Baseline SHA-256: `81adbf563b3a049947c21fdb2a1f1a17128bb327e091c331e1c62b2b6be49d82`

Checked 103 hashed ledger receipts; output including the owned replay package occupies 1614718 of 134217728 allowed bytes. Source/control/corpus pins and both prepared checkouts remain unchanged. Full successful time reservations remain charged; measured work is charged as recorded. Remaining amounts are protected ceilings, not permission to repeat completed stages or transfer between phases.

| Protected phase | Charged wall / CPU seconds | Remaining wall / CPU seconds | Charged work | Remaining work |
| --- | ---: | ---: | ---: | ---: |
| analysis | 600.000 / 600.000 | 0.000 / 0.000 | 0 | 0 |
| base50 | 17460.000 / 17460.000 | 0.000 / 0.000 | 2449623572196 | 357550376427804 |
| base55 | 41904.000 / 41904.000 | 0.000 / 0.000 | 3008977187508 | 356991022812492 |
| extension50 | 0.000 / 0.000 | 17460.000 / 17460.000 | 0 | 360000000000000 |
| extension55 | 0.000 / 0.000 | 41904.000 / 41904.000 | 0 | 360000000000000 |
| generation | 120.000 / 120.000 | 0.000 / 0.000 | 0 | 0 |
| orchestration | 120.000 / 120.000 | 480.000 / 480.000 | 0 | 0 |
| total | 61114.000 / 61114.000 | 64534.000 / 64534.000 | 5509030195848 | 3034490969804152 |
| warm50 | 630.000 / 630.000 | 2170.000 / 2170.000 | 43946600472 | 799956053399528 |
| warm55 | 280.000 / 280.000 | 2520.000 / 2520.000 | 6482835672 | 799993517164328 |

## Execution record and verification

Preserved worktrees are `/Users/nishanthmohan/Documents/code/factor-a1a-integration` and `/Users/nishanthmohan/Documents/code/factor-a1a-extension`. In each, durable receipts are under `v2/benchmarks/results/a1a_NAME/frozen/`: `ledger.json`, `environment.json`, per-call JSON, `base-finished-{50,55}.json`, `analysis-{50,55}.json` and `decision-{50,55}.json`. Integration also has `extension-finished-55.json` and its full extension receipts. Detailed per-input results, identity checks and account reconciliation remain ignored under each `v2/audit/a1a_NAME/execution_summary.json`; CLI logs and the sequential controller journal are in the corresponding audit execution folders.

Execution used the user-approved recorded `/opt/homebrew/bin/pypy3`: PyPy 7.3.23, Python 3.11.15, macOS 26.6.2 arm64. The project venv and recorded paths resolve to the same binary, SHA-256 `3c45fc49fa01571d19bac1594495d7300c30589c79c0cad2c5bee8d59f9e5ef2`. The initial venv-path invocation was rejected by the immutable executable receipt before any calibration or reservation. Its evidence and the explicit interpreter-path authorization are preserved locally; no receipt was overwritten to bypass the guard.

Executed sequence, with an explicit working directory for every command:

```sh
# cwd: /Users/nishanthmohan/Documents/code/factor-a1a-integration
/opt/homebrew/bin/pypy3 -B -m v2.benchmarks.qs.a1a_integration.study verify
/opt/homebrew/bin/pypy3 -Bu -m v2.benchmarks.qs.a1a_integration.study run --digits 50
/opt/homebrew/bin/pypy3 -B -m v2.benchmarks.qs.a1a_integration.study analyze --digits 50
/opt/homebrew/bin/pypy3 -Bu -m v2.benchmarks.qs.a1a_integration.study run --digits 55
/opt/homebrew/bin/pypy3 -Bu -m v2.benchmarks.qs.a1a_integration.study extend --digits 55
/opt/homebrew/bin/pypy3 -B -m v2.benchmarks.qs.a1a_integration.study analyze --digits 55

# cwd: /Users/nishanthmohan/Documents/code/factor-a1a-extension
/opt/homebrew/bin/pypy3 -B -m v2.benchmarks.qs.a1a_extension.study verify
/opt/homebrew/bin/pypy3 -Bu -m v2.benchmarks.qs.a1a_extension.study run --digits 50
/opt/homebrew/bin/pypy3 -B -m v2.benchmarks.qs.a1a_extension.study analyze --digits 50
/opt/homebrew/bin/pypy3 -Bu -m v2.benchmarks.qs.a1a_extension.study run --digits 55
/opt/homebrew/bin/pypy3 -B -m v2.benchmarks.qs.a1a_extension.study analyze --digits 55
```

These are execution records, not restart instructions. All stages are terminal; do not rerun generation, freeze, measurements or final inference against these ledgers. Analysis was performed once per band. Acceptance required every call complete, stable validated warmup, unchanged deterministic paths, all-data repeated-cell wall/CPU drift ≤15%, a 99.5% interval upper endpoint below zero and exact one-sided p≤0.005. No minimum percentage floor was added.

After all measurement processes ended, receipt audits and committed-files-only tests, lint, benchmark imports and frozen bindings were verified serially under the shared performance lock and the separate cumulative 3600-second wall/CPU QA allowance. Detailed commands, successful and earlier failed charges, clone identities and final remaining QA balances are preserved in the local audit ledger. The report branch changes documentation only; frozen runtime, runners, protocols and inputs are unchanged.
