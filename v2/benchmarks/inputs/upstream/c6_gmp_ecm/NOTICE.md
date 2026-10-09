# C6 upstream reference inputs

Exact files from GMP-ECM commit
`8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e`, fetched 9 October 2026:
https://github.com/sethtroisi/gmp-ecm/tree/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e

`LucasChainGen.c`, `LucasChainGen.h`, `LCG_macros.h` and `README.upstream`
come from LucasChainGenerator (README renamed only); `ecm.c` is the upstream
decoder and arithmetic reference. Copyright notices are preserved. These
sources are LGPL-3.0-or-later; COPYING and COPYING.LIB accompany them.
They are independent build/verification inputs, not linked into production.

The C6 build script copies these files to ignored scratch storage, reduces
MAX_CODE_OR_PRIME_COUNT from 6,000,000 to 4,096, and uses one generator
thread with B1=2,000. This changes allocation capacity, not search logic.
The decoder's declarations and functions are extracted verbatim, with
ASSERT enabled and a small standalone I/O wrapper. Generated records are
independently interpreted with exact Python integers before use.

No CADO source is adapted. Its pinned bytecode implementation at
`692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b` and LGPL-2.1 COPYING were
inspected as research references only. New Python register allocation and
integer verification code are written independently of the C decoder.
