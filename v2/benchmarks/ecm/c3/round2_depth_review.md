# C3 curve-count boundary and upstream effort

Reviewed 11 October 2026 after the user questioned the inherited 32-curve
schedule. PortfolioConfig limits the number of tiers to 64, but accepts any
nonnegative integer curve count in a tier. Thirty-two is a numerical default,
not an algorithmic maximum. Explicit total work/time limits remain binding.
The committed round-two candidate already doubles the count to 64.

## Actual upstream policies

Curve counts need their bounds, target factor size and caller. An ECM-only
campaign and a cheap pre-SIQS test optimize different costs.

| Source | Observed effort | Interpretation |
| --- | --- | --- |
| [GMP-ECM README](https://github.com/sethtroisi/gmp-ecm/blob/main/README), table 1 | Expected counts: 74 at B1/B2 11000/1.9e6 for a 20-digit factor; 214 at 50000/1.3e7 for 25 digits; 430 at 250000/1.3e8 for 30 digits; 7553 at 43000000/2.4e11 for 50 digits, using its stated continuation. | These are factor digits, not modulus digits or automatic SIQS cutoffs. One expected count still has a model miss probability about exp(-1). The native continuation and curve-family assumptions differ from v2. |
| [Zimmermann's parameter note](https://members.loria.fr/PZimmermann/records/ecm/params.html), preferred June 2019 table | A different time optimization for GMP-ECM 7 on a 512-bit modulus gives 41 expected curves at B1=24433 for about 20-digit factors, 339 at 445657 for 30 digits, and 8177 at 46919468 for 50 digits. | Different bound/continuation/cost assumptions produce different counts. Do not mix the older table on that page with this optimization or treat any table as a universal cap. |
| [YAFU scheduler](https://raw.githubusercontent.com/bbuhrow/yafu/master/factor/autofactor.c), schedule_work | Normal target factor depth is 4/13 of current modulus digits; light is 2/9, deep 1/3. get_curves_required credits completed effort and derives more work across tiers. A target below 15 digits on inputs <=45 digits can skip ECM. | No universal 32-curve cutoff. A depth-target rule and its native QS/NFS crossover cannot be transferred as a PyPy timing rule. The previous pinned research source shows the same mechanism. |
| [Alpertron ecmCurve](https://raw.githubusercontent.com/alpertron/calculators/master/ecm.c) | The bounds ladder allocates 25 at B1=2000, 90 at 11000, 300 at 50000 and deeper tiers. Separately, the automatic SIQS handoff table starts with threshold 10 and grows to 350. | The handoff check precedes execution of the threshold curve. Its NumberLength*9 proxy, NextEC mode and input range govern whether that check applies; the bounds ladder alone is not the pretest policy. |
| [Pinned Yamaquasi ecm_auto / ecm_only](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/ecm.rs) | Auto uses eight 200/7700 curves at 65..160 input bits, increasing to counts 20,25,40,100,150,100,250 in successive larger bands. ECM-only escalates through counts 10,30,100,100,200,600,2000,10000,15000 with much larger bounds. | The same implementation deliberately invests little where SIQS is inexpensive and much more for ECM-only searches. The inspected native Edwards/Suyama-11 backend differs from v2. |

The online GMP, YAFU and Alpertron sources were rechecked for this question.
The exact Yamaquasi revision was unavailable through the web fetch on this
recheck; its downloaded source and content hash from the original research
remain the evidence. No upstream code was copied or executed.

## What the present calibration can support

Both fitted counts reached the maximum offered, 64, so neither is an
established optimum. Inspection of the completed immutable captures gives:

| Band / ladder | Root entries | Surviving seed assignments after the full ladder | Independent surviving subjects |
| --- | ---: | ---: | ---: |
| 30 / compact64 | 171 | 3 | 3 |
| 30 / wide64 | 171 | 2 | 2 |
| 30 / deep20 | 171 | 17 | 9 |
| 40 / compact64 | 171 | 82 | 12 |
| 40 / wide64 | 171 | 74 | 12 |
| 40 / deep20 | 171 | 88 | 13 |

All 513 paired wide/deep captures agree on the first sixteen root curve
seeds, bounds and outcomes, including subjects solved before ECM. This is a
trajectory check, not a timing result. Root entries comprise nineteen
independent inputs x nine seeds per band; seeds are not independent subjects.

The 40-digit survivors provide the clearest training support for a separately
frozen 128/256-curve continuation and a stronger-bound tier after 64 failures.
At 30 digits only three compact64 survivors remain, below the already frozen
twelve-subject model floor. More seeds on those same inputs cannot repair
that population limitation. A larger independently generated training cohort
would be needed to estimate that tail under the existing support rule.

The immediate complete-call comparison still tests the committed 64-curve
candidate against both controls. Any deeper search must preserve this result,
use a new finite source/assignment freeze, and commit its resulting candidate
before new confirmation. No observation beyond 64 curves exists yet. Neither
this review nor the instrumented fit establishes a production default change.
