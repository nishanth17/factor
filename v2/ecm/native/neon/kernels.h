#ifndef FACTOR_ECM_NEON_KERNELS_H
#define FACTOR_ECM_NEON_KERNELS_H
#include <stddef.h>
#include <stdint.h>

#define B4N_LIMBS 13
#define B4N_RADIX 26
#define B4N_MASK ((UINT32_C(1) << B4N_RADIX) - 1)

/* [digit][curve lane], with disjoint preallocated output; exactly two lanes. */
void b4n_neon_product(uint32_t out[2 * B4N_LIMBS][2],
                      const uint32_t a[B4N_LIMBS][2],
                      const uint32_t b[B4N_LIMBS][2], unsigned limbs);
void b4n_scalar_product(uint32_t out[2 * B4N_LIMBS][2],
                        const uint32_t a[B4N_LIMBS][2],
                        const uint32_t b[B4N_LIMBS][2], unsigned limbs);
uint64_t b4n_consume(uint64_t checksum, uint32_t low0, uint32_t low1,
                      const void *output, size_t bytes);
#endif
