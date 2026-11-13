/* Independent canonical-residue feasibility probe, never a factoring backend.
 * GMP owns carries and division. Inputs are nonnegative residues below n.
 * Distinct product, quotient and remainder buffers avoid all partial overlap.
 */
#define _POSIX_C_SOURCE 200809L
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

_Static_assert(GMP_NUMB_BITS == 64 && GMP_NAIL_BITS == 0,
               "probe requires ordinary 64-bit GMP limbs");

static uint64_t now(void) {
    struct timespec value;
    if (clock_gettime(CLOCK_MONOTONIC, &value)) abort();
    return (uint64_t)value.tv_sec * 1000000000u + value.tv_nsec;
}

static void export_limbs(mp_limb_t *out, const mpz_t value, size_t width) {
    size_t count;
    memset(out, 0, width * sizeof(*out));
    mpz_export(out, &count, -1, sizeof(*out), 0, 0, value);
    if (count > width) abort();
}

static void reduce(mp_limb_t *out, mp_limb_t *product,
                   const mp_limb_t *modulus, size_t width) {
    mp_limb_t quotient[7], remainder[6] = {0};
    size_t used = 2 * width;
    while (used && !product[used - 1]) --used;
    /* If used<width, the product is below B^(width-1)<=n. Otherwise GMP's
     * tdiv contract applies: nn>=dn, top divisor limb nonzero, qxn=0.
     * Quotient needs at most (2*width-width+1) limbs; remainder needs width.
     */
    if (used < width) memcpy(remainder, product, used * sizeof(*product));
    else mpn_tdiv_qr(quotient, remainder, 0, product, used, modulus, width);
    memcpy(out, remainder, width * sizeof(*out));
}

static uint64_t pipeline(mpz_t output, const mpz_t n, const mpz_t initial,
                          const mpz_t multiplier, unsigned loops, int limb) {
    size_t width = mpz_size(n);
    if (width < 4 || width > 6 || mpz_sgn(initial) < 0 ||
        mpz_cmp(initial, n) >= 0 || mpz_sgn(multiplier) < 0 ||
        mpz_cmp(multiplier, n) >= 0) abort();
    mp_limb_t x[6], a[6], modulus[6], product[12];
    export_limbs(x, initial, width);
    export_limbs(a, multiplier, width);
    export_limbs(modulus, n, width);
    mpz_set(output, initial);
    mpz_realloc2(output, 2 * width * GMP_NUMB_BITS);
    uint64_t start = now();
    for (unsigned index = 0; index < loops; ++index) {
        if (limb) {
            mpn_mul_n(product, x, a, width);
            reduce(x, product, modulus, width);
            mpn_sqr(product, x, width);
            reduce(x, product, modulus, width);
        } else {
            mpz_mul(output, output, multiplier);
            mpz_mod(output, output, n);
            mpz_mul(output, output, output);
            mpz_mod(output, output, n);
        }
    }
    uint64_t elapsed = now() - start;
    if (limb) mpz_import(output, width, -1, sizeof(*x), 0, 0, x);
    return elapsed;
}

static void verify(void) {
    mpz_t n, x, a, ordinary, limbs;
    mpz_inits(n, x, a, ordinary, limbs, NULL);
    unsigned cases = 0;
    for (unsigned width = 4; width <= 6; ++width) {
        for (unsigned shape = 0; shape < 3; ++shape) {
            mpz_set_ui(n, 1);
            mpz_mul_2exp(n, n, shape == 2 ? 64 * width : 64 * (width - 1));
            if (shape == 2) mpz_sub_ui(n, n, 1);
            else mpz_add_ui(n, n, shape ? 0xffff : 3);
            for (unsigned left = 0; left < 4; ++left) {
                for (unsigned right = 0; right < 4; ++right) {
                    if (left < 2) mpz_set_ui(x, left);
                    else mpz_sub_ui(x, n, left - 1);
                    if (right < 2) mpz_set_ui(a, right);
                    else mpz_sub_ui(a, n, right - 1);
                    pipeline(ordinary, n, x, a, 127, 0);
                    pipeline(limbs, n, x, a, 127, 1);
                    if (mpz_cmp(ordinary, limbs) || mpz_cmp(limbs, n) >= 0)
                        abort();
                    ++cases;
                }
            }
        }
    }
    printf("{\"verified_cases\":%u,\"steps_per_case\":127}\n", cases);
    mpz_clears(n, x, a, ordinary, limbs, NULL);
}

int main(int argc, char **argv) {
    if (argc == 2 && !strcmp(argv[1], "verify")) { verify(); return 0; }
    if (argc != 6) return 2;
    mpz_t n, x, a, output;
    mpz_inits(n, x, a, output, NULL);
    if (mpz_set_str(n, argv[1], 16) || mpz_set_str(x, argv[2], 16) ||
        mpz_set_str(a, argv[3], 16)) return 2;
    unsigned loops = (unsigned)strtoul(argv[4], NULL, 10);
    int mode = atoi(argv[5]);
    if (!loops || loops > 4000000 || (mode != 0 && mode != 1)) return 2;
    uint64_t start = now();
    uint64_t elapsed = pipeline(output, n, x, a, loops, mode);
    uint64_t complete = now() - start;
    char *value = mpz_get_str(NULL, 16, output);
    printf("{\"ns\":%llu,\"complete_ns\":%llu,\"loops\":%u,\"mode\":%d,\"result\":\"%s\"}\n",
           (unsigned long long)elapsed, (unsigned long long)complete,
           loops, mode, value);
    free(value);
    mpz_clears(n, x, a, output, NULL);
    return 0;
}
