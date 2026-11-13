#include "kernels.h"

/* Separate translation unit, no LTO: the caller must materialize every full
 * output before this opaque call. The same low26 checksum update is used by
 * every arm; the full products/residues are checked outside the timed region.
 */
__attribute__((noinline))
uint64_t b4n_consume(uint64_t checksum, uint32_t low0, uint32_t low1,
                      const void *output, size_t bytes) {
    __asm__ __volatile__("" : : "r"(output), "r"(bytes) : "memory");
    return ((checksum << 7) | (checksum >> 57)) ^
           ((uint64_t)low0 + UINT64_C(0x9e3779b97f4a7c15) * low1);
}
