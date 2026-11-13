/* Fixed-width ECM feasibility harness. Reduction is ordinary GMP mpz_mod. */
#define _DARWIN_C_SOURCE 1
#include "screen.h"
#include "kernels.h"
#include <gmp.h>
#include <pthread.h>
#include <pthread/qos.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if GMP_NUMB_BITS != 64 || GMP_NAIL_BITS != 0
#error "This optional M4 screen requires 64-bit GMP limbs without nails"
#endif
#define ENTRIES 64
#define WORDS 6
#define INITIAL_CHECKSUM UINT64_C(0x6a09e667f3bcc909)

typedef struct {
    mpz_t inputs[4], product[2][2], residue[2][2], result[2];
    mp_limb_t words[4][WORDS];
    uint32_t a[B4N_LIMBS][2], b[B4N_LIMBS][2];
    uint32_t work_a[B4N_LIMBS][2], work_b[B4N_LIMBS][2];
    uint32_t packed[4][B4N_LIMBS], normalized[2][2 * B4N_LIMBS];
    uint32_t known[2][2 * B4N_LIMBS][2];
    union {
        mp_limb_t mpn[2][2 * WORDS];
        uint32_t digits[2 * B4N_LIMBS][2];
    } raw;
    unsigned active;
} entry;

struct b4n_context {
    mpz_t modulus, exported[2];
    entry entries[ENTRIES];
    uint64_t ready;
    unsigned width, limbs, words, last_arm, last_operation, last_index;
    unsigned last_repeats, last_offset, has_result, has_batch;
    b4n_observation last_observation;
};

unsigned b4n_abi(void) { return 1; }
const char *b4n_gmp_version(void) { return gmp_version; }

static double seconds(clockid_t clock) {
    struct timespec value;
    if (clock_gettime(clock, &value)) return -1;
    return (double)value.tv_sec + (double)value.tv_nsec / 1e9;
}

static int canonical(b4n_value value, unsigned maximum) {
    return value.data && value.size && value.size <= maximum &&
           (value.size == 1 || value.data[0] != 0);
}

static void read_value(mpz_t out, b4n_value value) {
    mpz_import(out, value.size, 1, 1, 1, 0, value.data);
}

b4n_context *b4n_create(unsigned width, b4n_value modulus) {
    if ((width != 200 && width != 266 && width != 333) ||
        !canonical(modulus, (width + 7) / 8)) return NULL;
    b4n_context *c = calloc(1, sizeof(*c));
    if (!c) return NULL;
    c->width = width;
    c->limbs = (width + B4N_RADIX - 1) / B4N_RADIX;
    c->words = (width + 63) / 64;
    mpz_init2(c->modulus, width);
    mpz_init2(c->exported[0], 2 * width + 64);
    mpz_init2(c->exported[1], 2 * width + 64);
    read_value(c->modulus, modulus);
    if (mpz_sizeinbase(c->modulus, 2) != width || mpz_even_p(c->modulus)) {
        mpz_clears(c->modulus, c->exported[0], c->exported[1], NULL);
        free(c);
        return NULL;
    }
    for (unsigned i = 0; i < ENTRIES; ++i) {
        entry *e = &c->entries[i];
        for (unsigned j = 0; j < 4; ++j) mpz_init2(e->inputs[j], width);
        for (unsigned op = 0; op < 2; ++op)
            for (unsigned lane = 0; lane < 2; ++lane) {
                mpz_init2(e->product[op][lane], 2 * width + 64);
                mpz_init2(e->residue[op][lane], 2 * width + 64);
            }
        for (unsigned lane = 0; lane < 2; ++lane)
            mpz_init2(e->result[lane], 2 * width + 64);
    }
    return c;
}

void b4n_destroy(b4n_context *c) {
    if (!c) return;
    for (unsigned i = 0; i < ENTRIES; ++i) {
        entry *e = &c->entries[i];
        for (unsigned j = 0; j < 4; ++j) mpz_clear(e->inputs[j]);
        for (unsigned op = 0; op < 2; ++op)
            for (unsigned lane = 0; lane < 2; ++lane) {
                mpz_clear(e->product[op][lane]);
                mpz_clear(e->residue[op][lane]);
            }
        for (unsigned lane = 0; lane < 2; ++lane) mpz_clear(e->result[lane]);
    }
    mpz_clears(c->modulus, c->exported[0], c->exported[1], NULL);
    free(c);
}

static void pack(entry *e, unsigned limbs, unsigned square) {
    /* mpz_export writes nothing for zero: clear the complete bounded buffer.
     * Six nail bits give canonical radix26 words using the documented API.
     */
    memset(e->packed, 0, sizeof(e->packed));
    for (unsigned lane = 0; lane < 2; ++lane) {
        mpz_export(e->packed[2 * lane], NULL, -1, sizeof(uint32_t), 0, 6,
                   e->inputs[2 * lane]);
        if (!square)
            mpz_export(e->packed[2 * lane + 1], NULL, -1, sizeof(uint32_t), 0, 6,
                       e->inputs[2 * lane + 1]);
    }
    for (unsigned i = 0; i < limbs; ++i)
        for (unsigned lane = 0; lane < 2; ++lane) {
            e->work_a[i][lane] = e->packed[2 * lane][i];
            e->work_b[i][lane] = square ? e->work_a[i][lane] :
                                          e->packed[2 * lane + 1][i];
        }
}

int b4n_set(b4n_context *c, unsigned index, const b4n_value values[4],
             unsigned active) {
    if (!c || !values || index >= ENTRIES || (active != 1 && active != 2)) return -1;
    c->ready &= ~(UINT64_C(1) << index);
    c->has_result = c->has_batch = 0;
    entry *e = &c->entries[index];
    for (unsigned j = 0; j < 2 * active; ++j) {
        if (!canonical(values[j], (c->width + 7) / 8)) return -1;
        read_value(e->inputs[j], values[j]);
        if (mpz_sizeinbase(e->inputs[j], 2) > c->width) return -1;
    }
    /* A one-curve tail never reads the absent lane's input descriptors. */
    if (active == 1) mpz_set_ui(e->inputs[2], 0), mpz_set_ui(e->inputs[3], 0);
    e->active = active;
    memset(e->words, 0, sizeof(e->words));
    for (unsigned j = 0; j < 4; ++j)
        mpz_export(e->words[j], NULL, -1, sizeof(mp_limb_t), 0, 0, e->inputs[j]);
    pack(e, c->limbs, 0);
    memcpy(e->a, e->work_a, sizeof(e->a));
    memcpy(e->b, e->work_b, sizeof(e->b));
    for (unsigned op = 0; op < 2; ++op) {
        memset(e->known[op], 0, sizeof(e->known[op]));
        for (unsigned lane = 0; lane < 2; ++lane) {
            unsigned right = 2 * lane + (op == 0);
            mpz_mul(e->product[op][lane], e->inputs[2 * lane], e->inputs[right]);
            mpz_mod(e->residue[op][lane], e->product[op][lane], c->modulus);
            memset(e->normalized[lane], 0, sizeof(e->normalized[lane]));
            mpz_export(e->normalized[lane], NULL, -1, sizeof(uint32_t), 0, 6,
                       e->product[op][lane]);
            for (unsigned i = 0; i < 2 * c->limbs; ++i)
                e->known[op][i][lane] = e->normalized[lane][i];
        }
    }
    c->ready |= UINT64_C(1) << index;
    return 0;
}

static void reduce_digits(b4n_context *c, entry *e) {
    for (unsigned lane = 0; lane < 2; ++lane) {
        for (unsigned i = 0; i < 2 * c->limbs; ++i)
            e->normalized[lane][i] = e->raw.digits[i][lane];
        unsigned count = 2 * c->limbs;
        while (count && !e->normalized[lane][count - 1]) --count;
        mpz_import(e->result[lane], count, -1, sizeof(uint32_t), 0, 6,
                   e->normalized[lane]);
        mpz_mod(e->result[lane], e->result[lane], c->modulus);
    }
}

static void evaluate(b4n_context *c, entry *e, unsigned arm, unsigned op) {
    if (arm == 0) {
        for (unsigned lane = 0; lane < 2; ++lane) {
            if (op) mpn_sqr(e->raw.mpn[lane], e->words[2 * lane], c->words);
            else mpn_mul_n(e->raw.mpn[lane], e->words[2 * lane],
                            e->words[2 * lane + 1], c->words);
        }
    } else if (arm == 1 || arm == 2) {
        if (arm == 1)
            b4n_scalar_product(e->raw.digits, e->a, op ? e->a : e->b, c->limbs);
        else b4n_neon_product(e->raw.digits, e->a, op ? e->a : e->b, c->limbs);
    } else if (arm == 3) {
        /* Exactly the production packet's ordinary-GMP multiply/square then
         * reduction pattern, with reusable objects and no Python boundary.
         */
        for (unsigned lane = 0; lane < 2; ++lane) {
            unsigned right = 2 * lane + (op == 0);
            mpz_mul(e->result[lane], e->inputs[2 * lane], e->inputs[right]);
            mpz_mod(e->result[lane], e->result[lane], c->modulus);
        }
    } else {
        pack(e, c->limbs, op);
        if (arm == 4)
            b4n_scalar_product(e->raw.digits, e->work_a, e->work_b, c->limbs);
        else if (arm == 5)
            b4n_neon_product(e->raw.digits, e->work_a, e->work_b, c->limbs);
        else memcpy(e->raw.digits, e->known[op],
                     2 * c->limbs * 2 * sizeof(uint32_t));
        reduce_digits(c, e);
    }
}

static void materialize(b4n_context *c, entry *e, unsigned arm) {
    for (unsigned lane = 0; lane < 2; ++lane) {
        if (arm == 0)
            mpz_import(c->exported[lane], 2 * c->words, -1,
                       sizeof(mp_limb_t), 0, 0, e->raw.mpn[lane]);
        else if (arm <= 2) {
            for (unsigned i = 0; i < 2 * c->limbs; ++i)
                e->normalized[lane][i] = e->raw.digits[i][lane];
            mpz_import(c->exported[lane], 2 * c->limbs, -1,
                       sizeof(uint32_t), 0, 6, e->normalized[lane]);
        } else mpz_set(c->exported[lane], e->result[lane]);
    }
}

int b4n_once(b4n_context *c, unsigned arm, unsigned op, unsigned index) {
    if (!c || arm > 6 || op > 1 || index >= ENTRIES ||
        !(c->ready & (UINT64_C(1) << index))) return -1;
    c->has_result = c->has_batch = 0;
    evaluate(c, &c->entries[index], arm, op);
    materialize(c, &c->entries[index], arm);
    for (unsigned lane = 0; lane < 2; ++lane)
        if (mpz_cmp(c->exported[lane], arm <= 2 ?
                    c->entries[index].product[op][lane] :
                    c->entries[index].residue[op][lane])) return -2;
    c->last_arm = arm; c->last_operation = op; c->last_index = index;
    c->has_result = 1;
    return 0;
}

int b4n_export(b4n_context *c, b4n_output *out) {
    if (!c || !out || !c->has_result) return -1;
    materialize(c, &c->entries[c->last_index], c->last_arm);
    memset(out, 0, sizeof(*out));
    for (unsigned lane = 0; lane < 2; ++lane) {
        size_t count = 0;
        if (mpz_sizeinbase(c->exported[lane], 2) > 2 * c->width) return -1;
        mpz_export(out->bytes[lane], &count, 1, 1, 1, 0, c->exported[lane]);
        out->size[lane] = count ? (unsigned)count : 1;
    }
    return 0;
}

static uint64_t mix(uint64_t checksum, uint32_t low0, uint32_t low1) {
    return ((checksum << 7) | (checksum >> 57)) ^
           ((uint64_t)low0 + UINT64_C(0x9e3779b97f4a7c15) * low1);
}

void b4n_qos(b4n_observation *out) {
    if (!out) return;
    out->requested_qos = QOS_CLASS_USER_INTERACTIVE;
    out->qos_set = pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
    qos_class_t observed = QOS_CLASS_UNSPECIFIED;
    out->qos_get = pthread_get_qos_class_np(pthread_self(), &observed,
                                           &out->relative_priority);
    out->observed_qos = observed;
    uint64_t thread_id = 0;
    out->thread_id_status = pthread_threadid_np(NULL, &thread_id);
    out->thread_id = thread_id;
}

int b4n_batch(b4n_context *c, unsigned arm, unsigned op,
               unsigned repeats, unsigned offset, b4n_observation *out) {
    if (!c || !out || arm > 6 || op > 1 || !repeats || repeats > 65536 ||
        offset >= ENTRIES || c->ready != UINT64_MAX) return -1;
    for (unsigned i = 0; i < ENTRIES; ++i)
        if (c->entries[i].active != 2) return -1;
    c->has_result = c->has_batch = 0;
    memset(out, 0, sizeof(*out));
    b4n_qos(out);
    double wall = seconds(CLOCK_MONOTONIC), cpu = seconds(CLOCK_PROCESS_CPUTIME_ID);
    if (wall < 0 || cpu < 0) return -1;
    uint64_t checksum = INITIAL_CHECKSUM;
    unsigned r;
    for (r = 0; r < repeats; ++r) {
        /* Same cooperative watchdog check in every arm, outside the pair loop.
         * It bounds an accidental oversized calibration without per-op clocks.
         */
        if (!(r & 15)) {
            double now = seconds(CLOCK_MONOTONIC);
            if (now < 0 || now - wall >= 0.25) break;
        }
        for (unsigned j = 0; j < ENTRIES; ++j) {
            entry *e = &c->entries[(j + r + offset) & (ENTRIES - 1)];
            evaluate(c, e, arm, op);
            uint32_t low0, low1;
            if (!arm) low0 = e->raw.mpn[0][0] & B4N_MASK,
                      low1 = e->raw.mpn[1][0] & B4N_MASK;
            else if (arm <= 2) low0 = e->raw.digits[0][0] & B4N_MASK,
                              low1 = e->raw.digits[0][1] & B4N_MASK;
            else low0 = mpz_get_ui(e->result[0]) & B4N_MASK,
                 low1 = mpz_get_ui(e->result[1]) & B4N_MASK;
            checksum = b4n_consume(checksum, low0, low1, e, sizeof(*e));
        }
    }
    double wall_end = seconds(CLOCK_MONOTONIC);
    double cpu_end = seconds(CLOCK_PROCESS_CPUTIME_ID);
    if (wall_end < 0 || cpu_end < 0) return -1;
    out->wall = wall_end - wall;
    out->cpu = cpu_end - cpu;
    out->checksum = checksum;
    out->pairs = (uint64_t)r * ENTRIES;
    c->last_arm = arm; c->last_operation = op; c->last_index = 0;
    c->last_repeats = r; c->last_offset = offset;
    c->has_batch = c->has_result = r != 0;
    c->last_observation = *out;
    return r == repeats ? 0 : 1;
}

int b4n_validate(b4n_context *c, const b4n_observation *out) {
    if (!c || !out || !c->has_batch || out->pairs != c->last_observation.pairs ||
        out->checksum != c->last_observation.checksum) return -1;
    unsigned arm = c->last_arm, op = c->last_operation;
    /* Full-width comparison is deliberately outside the timed region. */
    for (unsigned index = 0; index < ENTRIES; ++index) {
        entry *e = &c->entries[index];
        materialize(c, e, arm);
        for (unsigned lane = 0; lane < 2; ++lane)
            if (mpz_cmp(c->exported[lane], arm <= 2 ? e->product[op][lane] :
                                                      e->residue[op][lane]))
                return -1;
    }
    uint64_t expected = INITIAL_CHECKSUM;
    for (unsigned r = 0; r < c->last_repeats; ++r)
        for (unsigned j = 0; j < ENTRIES; ++j) {
            entry *e = &c->entries[(j + r + c->last_offset) & (ENTRIES - 1)];
            mpz_srcptr value0 = arm <= 2 ? e->product[op][0] : e->residue[op][0];
            mpz_srcptr value1 = arm <= 2 ? e->product[op][1] : e->residue[op][1];
            expected = mix(expected, mpz_get_ui(value0) & B4N_MASK,
                                      mpz_get_ui(value1) & B4N_MASK);
        }
    return expected == out->checksum ? 0 : -1;
}
