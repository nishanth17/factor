/* Full canonical-residue SME feasibility, never selected by Factor.
 * The scalar controls and ordinary-division proof are shared with the owned
 * limb probe. Transposed limbs persist through the entire arithmetic batch.
 */
#define main scalar_pipeline_main
#include "limb_pipeline_probe.c"
#undef main
#include <arm_sve.h>

__arm_locally_streaming __attribute__((noinline))
static unsigned lanes(void) { return (unsigned)svcntd(); }

__arm_locally_streaming __attribute__((noinline))
static void vector_product(mp_limb_t product[12][32],
                            mp_limb_t left[6][32], mp_limb_t right[6][32],
                            unsigned width, unsigned count) {
    svbool_t pg = svwhilelt_b64_u64(0, count);
    svuint64_t zero = svdup_u64(0);
    for (unsigned k = 0; k < 2 * width; ++k) svst1_u64(pg, product[k], zero);
    for (unsigned i = 0; i < width; ++i) {
        svuint64_t a = svld1_u64(pg, left[i]), carry = zero;
        for (unsigned j = 0; j < width; ++j) {
            svuint64_t b = svld1_u64(pg, right[j]);
            svuint64_t lo = svmul_u64_x(pg, a, b);
            svuint64_t hi = svmulh_u64_x(pg, a, b);
            svuint64_t old = svld1_u64(pg, product[i + j]);
            svuint64_t first = svadd_u64_x(pg, lo, old);
            svbool_t c1 = svcmplt_u64(pg, first, lo);
            svuint64_t second = svadd_u64_x(pg, first, carry);
            svbool_t c2 = svcmplt_u64(pg, second, first);
            carry = svadd_u64_x(pg, hi, svdup_u64_z(c1, 1));
            carry = svadd_u64_x(pg, carry, svdup_u64_z(c2, 1));
            /* a*b+old+carry_in <= (B-1)^2+2(B-1)=B^2-1.
             * Thus the final carry fits one limb even when both adds carry. */
            svst1_u64(pg, product[i + j], second);
        }
        svst1_u64(pg, product[i + width], carry);
    }
}

static void vector_reduce(mp_limb_t x[6][32], mp_limb_t p[12][32],
                           const mp_limb_t *n, unsigned width, unsigned count) {
    for (unsigned lane = 0; lane < count; ++lane) {
        mp_limb_t product[12], remainder[6];
        for (unsigned j = 0; j < 2 * width; ++j) product[j] = p[j][lane];
        reduce(remainder, product, n, width);
        for (unsigned j = 0; j < width; ++j) x[j][lane] = remainder[j];
    }
}

static uint64_t batch(mpz_t *output, const mpz_t n, mpz_t *x, mpz_t *a,
                       unsigned loops, unsigned mode, unsigned count) {
    uint64_t start = now();
    if (mode < 2) {
        for (unsigned i = 0; i < count; ++i)
            pipeline(output[i], n, x[i], a[i], loops, mode);
    } else {
        unsigned width = (unsigned)mpz_size(n);
        if (width < 4 || width > 6 || count > lanes() || count > 8) abort();
        mp_limb_t xp[6][32] = {{0}}, ap[6][32] = {{0}}, product[12][32];
        mp_limb_t modulus[6], temporary[6];
        export_limbs(modulus, n, width);
        for (unsigned i = 0; i < count; ++i) {
            export_limbs(temporary, x[i], width);
            for (unsigned j = 0; j < width; ++j) xp[j][i] = temporary[j];
            export_limbs(temporary, a[i], width);
            for (unsigned j = 0; j < width; ++j) ap[j][i] = temporary[j];
        }
        for (unsigned i = 0; i < loops; ++i) {
            vector_product(product, xp, ap, width, count);
            vector_reduce(xp, product, modulus, width, count);
            vector_product(product, xp, xp, width, count);
            vector_reduce(xp, product, modulus, width, count);
        }
        for (unsigned i = 0; i < count; ++i) {
            for (unsigned j = 0; j < width; ++j) temporary[j] = xp[j][i];
            mpz_import(output[i], width, -1, sizeof(*temporary), 0, 0, temporary);
        }
    }
    return now() - start;
}

int main(int argc, char **argv) {
    if (argc != 7) return 2;
    mpz_t n, initial, multiplier, x[8], a[8], output[8];
    mpz_inits(n, initial, multiplier, NULL);
    if (mpz_set_str(n, argv[1], 16) || mpz_set_str(initial, argv[2], 16) ||
        mpz_set_str(multiplier, argv[3], 16)) return 2;
    unsigned loops = (unsigned)strtoul(argv[4], NULL, 10);
    unsigned mode = (unsigned)atoi(argv[5]), count = (unsigned)atoi(argv[6]);
    if (!loops || loops > 1000000 || mode > 2 || !count || count > 8) return 2;
    for (unsigned i = 0; i < count; ++i) {
        mpz_inits(x[i], a[i], output[i], NULL);
        mpz_add_ui(x[i], initial, 101 * i); mpz_mod(x[i], x[i], n);
        mpz_add_ui(a[i], multiplier, 103 * i); mpz_mod(a[i], a[i], n);
    }
    uint64_t elapsed = batch(output, n, x, a, loops, mode, count);
    printf("{\"ns\":%llu,\"loops\":%u,\"mode\":%u,\"curves\":%u,\"lanes\":%u,\"results\":[",
           (unsigned long long)elapsed, loops, mode, count, lanes());
    for (unsigned i = 0; i < count; ++i) {
        char *value = mpz_get_str(NULL, 16, output[i]);
        printf("%s\"%s\"", i ? "," : "", value);free(value);
        mpz_clears(x[i], a[i], output[i], NULL);
    }
    puts("]}");mpz_clears(n, initial, multiplier, NULL);
    return 0;
}
