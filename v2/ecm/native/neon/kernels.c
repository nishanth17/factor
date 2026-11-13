/* ECM-local, two-curve feasibility kernel; no modular arithmetic here. */
#include "kernels.h"
#if !defined(__aarch64__)
#error "The optional ECM NEON screen requires AArch64"
#endif
#include <arm_neon.h>

/* With B=2^26 and L<=13, carry<13B and each diagonal accumulator is
 * <=13(B-1)^2+(13B-1)<13B^2<2^56. No lane can overflow uint64.
 * The final carry fits one radix digit because each input has <=L digits.
 */
__attribute__((noinline))
void b4n_neon_product(uint32_t out[2 * B4N_LIMBS][2],
                      const uint32_t a[B4N_LIMBS][2],
                      const uint32_t b[B4N_LIMBS][2], unsigned limbs) {
    uint64x2_t carry = vdupq_n_u64(0);
    const uint64x2_t mask = vdupq_n_u64(B4N_MASK);
    for (unsigned k = 0; k < 2 * limbs - 1; ++k) {
        uint64x2_t total = carry;
        unsigned first = k < limbs ? 0 : k - limbs + 1;
        unsigned last = k < limbs ? k : limbs - 1;
        for (unsigned i = first; i <= last; ++i) {
            uint32x2_t left = vld1_u32(a[i]);
            uint32x2_t right = vld1_u32(b[k - i]);
            total = vaddq_u64(total, vmull_u32(left, right));
        }
        vst1_u32(out[k], vmovn_u64(vandq_u64(total, mask)));
        carry = vshrq_n_u64(total, B4N_RADIX);
    }
    vst1_u32(out[2 * limbs - 1], vmovn_u64(carry));
}

/* Same representation, diagonal order and carries. The build disables
 * automatic loop/SLP vectorization; inspect both emitted symbols in admission.
 */
__attribute__((noinline))
void b4n_scalar_product(uint32_t out[2 * B4N_LIMBS][2],
                        const uint32_t a[B4N_LIMBS][2],
                        const uint32_t b[B4N_LIMBS][2], unsigned limbs) {
    uint64_t carry0 = 0, carry1 = 0;
    for (unsigned k = 0; k < 2 * limbs - 1; ++k) {
        uint64_t total0 = carry0, total1 = carry1;
        unsigned first = k < limbs ? 0 : k - limbs + 1;
        unsigned last = k < limbs ? k : limbs - 1;
        for (unsigned i = first; i <= last; ++i) {
            total0 += (uint64_t)a[i][0] * b[k - i][0];
            total1 += (uint64_t)a[i][1] * b[k - i][1];
        }
        out[k][0] = (uint32_t)(total0 & B4N_MASK);
        out[k][1] = (uint32_t)(total1 & B4N_MASK);
        carry0 = total0 >> B4N_RADIX;
        carry1 = total1 >> B4N_RADIX;
    }
    out[2 * limbs - 1][0] = (uint32_t)carry0;
    out[2 * limbs - 1][1] = (uint32_t)carry1;
}
