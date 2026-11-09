/* Independent implementation of Factor's existing point/program formulas.
 * Only documented GMP mpz APIs are used. No allocator hooks or global state.
 */
#define _POSIX_C_SOURCE 200809L
#include "packet.h"
#include <gmp.h>
#include <stdint.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef B4_OPENMP
#include <omp.h>
#endif

#define ABI 1
#define POINTS 16

struct b4_context {
    mpz_t n, a24, scalar, x[POINTS], z[POINTS], scratch[7], guard;
    unsigned char modulus_bytes[128], a24_bytes[128];
    size_t modulus_size, a24_size;
    const b4_record *records;
    size_t count, record, opcode;
    long bit;
    unsigned mode, done, completed, starting;
    /* Per-context deterministic clock scripts; inert in normal execution. */
    uint64_t clock_script[2][16];
    unsigned clock_count[2], clock_index[2];
};

unsigned b4_abi(void) { return ABI; }
const char *b4_gmp_version(void) { return gmp_version; }

b4_context *b4_create(void) {
    b4_context *c = calloc(1, sizeof(*c));
    if (!c) return NULL;
    mpz_inits(c->n, c->a24, c->scalar, c->guard, NULL);
    for (unsigned i = 0; i < POINTS; ++i) mpz_inits(c->x[i], c->z[i], NULL);
    for (unsigned i = 0; i < 7; ++i) mpz_init(c->scratch[i]);
    return c;
}

void b4_destroy(b4_context *c) {
    if (!c) return;
    mpz_clears(c->n, c->a24, c->scalar, c->guard, NULL);
    for (unsigned i = 0; i < POINTS; ++i) mpz_clears(c->x[i], c->z[i], NULL);
    for (unsigned i = 0; i < 7; ++i) mpz_clear(c->scratch[i]);
    free(c);
}

void b4_abort(b4_context *c) {
    if (!c) return;
    c->mode = c->done = 0;
    c->records = NULL;
    c->count = 0;
}

static int canonical(b4_value v, size_t cap) {
    return v.data && v.size && v.size <= cap &&
           (v.size == 1 || v.data[0] != 0);
}

static void import_value(mpz_t dest, b4_value v) {
    mpz_import(dest, v.size, 1, 1, 1, 0, v.data);
}

static int begin(b4_context *c, unsigned abi, b4_value n, b4_value x,
                 b4_value z, b4_value a24) {
    if (c) c->mode = c->done = 0;
    if (!c || abi != ABI || !canonical(n, 128) || !canonical(x, 128) ||
        !canonical(z, 128) || !canonical(a24, 128)) return -1;
    c->mode = c->done = 0;
    /* Rebind by canonical bytes, never by a caller's implicit context claim.
     * Constant modulus/curve data normally survive many useful packets.
     */
    if (c->modulus_size != n.size || memcmp(c->modulus_bytes, n.data, n.size)) {
        import_value(c->n, n);
        memcpy(c->modulus_bytes, n.data, n.size);
        c->modulus_size = n.size;
    }
    import_value(c->x[0], x);
    import_value(c->z[0], z);
    if (c->a24_size != a24.size || memcmp(c->a24_bytes, a24.data, a24.size)) {
        import_value(c->a24, a24);
        memcpy(c->a24_bytes, a24.data, a24.size);
        c->a24_size = a24.size;
    }
    if (mpz_cmp_ui(c->n, 3) < 0 || mpz_even_p(c->n) ||
        mpz_cmp(c->x[0], c->n) >= 0 || mpz_cmp(c->z[0], c->n) >= 0 ||
        mpz_cmp(c->a24, c->n) >= 0 ||
        (!mpz_sgn(c->x[0]) && !mpz_sgn(c->z[0]))) return -1;
    mpz_set_ui(c->guard, 1);
    c->completed = 0;
    return 0;
}

/* Stage the two results before writing: destination may alias any input. */
static void doubling(b4_context *c, unsigned dest, unsigned src) {
    mpz_t *s = c->scratch;
    mpz_add(s[0], c->x[src], c->z[src]);
    mpz_sub(s[1], c->x[src], c->z[src]);
    mpz_mul(s[0], s[0], s[0]); mpz_mod(s[0], s[0], c->n);
    mpz_mul(s[1], s[1], s[1]); mpz_mod(s[1], s[1], c->n);
    mpz_sub(s[2], s[0], s[1]);
    mpz_mul(s[3], c->a24, s[2]);
    /* a24=(A+2)/4 requires BB (s[1]) in this bracket. */
    mpz_add(s[3], s[3], s[1]);
    mpz_mul(s[4], s[0], s[1]); mpz_mod(s[4], s[4], c->n);
    mpz_mul(s[5], s[2], s[3]); mpz_mod(s[5], s[5], c->n);
    mpz_set(c->x[dest], s[4]); mpz_set(c->z[dest], s[5]);
}

static void addition(b4_context *c, unsigned dest, unsigned left,
                     unsigned right, unsigned difference) {
    mpz_t *s = c->scratch;
    mpz_sub(s[0], c->x[left], c->z[left]);
    mpz_add(s[1], c->x[right], c->z[right]);
    mpz_mul(s[0], s[0], s[1]); mpz_mod(s[0], s[0], c->n);
    mpz_add(s[1], c->x[left], c->z[left]);
    mpz_sub(s[2], c->x[right], c->z[right]);
    mpz_mul(s[1], s[1], s[2]); mpz_mod(s[1], s[1], c->n);
    mpz_add(s[2], s[0], s[1]); mpz_sub(s[3], s[0], s[1]);
    mpz_mul(s[2], s[2], s[2]); mpz_mul(s[3], s[3], s[3]);
    mpz_mul(s[4], c->z[difference], s[2]);
    mpz_mul(s[5], c->x[difference], s[3]);
    mpz_mod(c->x[dest], s[4], c->n); mpz_mod(c->z[dest], s[5], c->n);
}

/* Match Factor's accepted fused ladder: the doubled adjacent point shares
 * its sum/difference with differential addition. No representation changes.
 */
static void ladder_bit(b4_context *c, unsigned q, unsigned r) {
    mpz_t *s = c->scratch;
    mpz_add(s[0], c->x[q], c->z[q]);
    mpz_sub(s[1], c->x[q], c->z[q]);
    mpz_mul(s[2], s[0], s[0]); mpz_mod(s[2], s[2], c->n);
    mpz_mul(s[3], s[1], s[1]); mpz_mod(s[3], s[3], c->n);
    mpz_sub(s[4], s[2], s[3]);
    mpz_add(s[5], c->x[r], c->z[r]);
    mpz_mul(s[5], s[1], s[5]); mpz_mod(s[5], s[5], c->n);
    mpz_sub(s[6], c->x[r], c->z[r]);
    mpz_mul(s[6], s[0], s[6]); mpz_mod(s[6], s[6], c->n);
    mpz_add(s[0], s[5], s[6]); mpz_sub(s[1], s[5], s[6]);
    mpz_mul(s[5], s[0], s[0]); mpz_mul(s[5], s[5], c->z[0]);
    mpz_mod(s[5], s[5], c->n);
    mpz_mul(s[6], s[1], s[1]); mpz_mul(s[6], s[6], c->x[0]);
    mpz_mod(s[6], s[6], c->n);
    mpz_mul(s[0], s[2], s[3]); mpz_mod(s[0], s[0], c->n);
    mpz_mul(s[1], c->a24, s[4]); mpz_add(s[1], s[1], s[3]);
    mpz_mul(s[1], s[1], s[4]); mpz_mod(s[1], s[1], c->n);
    mpz_set(c->x[q], s[0]); mpz_set(c->z[q], s[1]);
    mpz_set(c->x[r], s[5]); mpz_set(c->z[r], s[6]);
}

int b4_ladder(b4_context *c, unsigned abi, b4_value n, b4_value x,
              b4_value z, b4_value a24, b4_value scalar) {
    if (c) c->mode = c->done = 0;
    if (!canonical(scalar, 512) || begin(c, abi, n, x, z, a24)) return -1;
    import_value(c->scalar, scalar);
    c->mode = 1;
    mpz_set(c->x[1], c->x[0]); mpz_set(c->z[1], c->z[0]);
    if (!mpz_sgn(c->scalar) || !mpz_sgn(c->z[0])) {
        mpz_set_ui(c->x[1], 1); mpz_set_ui(c->z[1], 0); c->done = 1;
    } else if (mpz_cmp_ui(c->scalar, 1) == 0) {
        c->done = 1;
    } else {
        doubling(c, 2, 0);
        c->bit = (long)mpz_sizeinbase(c->scalar, 2) - 2;
    }
    return 0;
}

int b4_records(b4_context *c, unsigned abi, b4_value n, b4_value x,
               b4_value z, b4_value a24, const b4_record *records, size_t count) {
    if (c) c->mode = c->done = 0;
    if (!records || !count || count > 16) return -1;
    /* Structural safety is checked independently of Python's scalar/frontier proof. */
    for (size_t i = 0; i < count; ++i) {
        const b4_record *r = &records[i];
        if (!r->code || !r->masks || r->code_size % 4 || r->code_size > 2048 ||
            r->masks_size != r->code_size / 4 + 1 || !r->slots ||
            r->slots > POINTS || r->output >= r->slots) return -1;
        unsigned initialized = 1;
        if (r->masks[0] > 3) return -1;
        for (size_t j = 0; j < r->code_size; j += 4) {
            const unsigned char *op = r->code + j;
            if (op[0] >= r->slots || op[1] >= r->slots || op[2] >= r->slots ||
                !(initialized & (1u << op[1])) || !(initialized & (1u << op[2])) ||
                r->masks[j / 4 + 1] > 3) return -1;
            if (op[3] == 255) {
                if (op[1] != op[2]) return -1;
            } else if (op[3] >= r->slots || !(initialized & (1u << op[3]))) return -1;
            initialized |= 1u << op[0];
        }
        if (!(initialized & (1u << r->output))) return -1;
    }
    if (begin(c, abi, n, x, z, a24)) return -1;
    c->records = records; c->count = count;
    c->record = c->opcode = 0; c->starting = 1; c->mode = 2;
    return 0;
}

static void guard_point(b4_context *c, unsigned slot, unsigned mask) {
    /* A zero mask contributes no coordinate; an already-zero guard remains
     * zero. Continue executing the record to preserve its canonical output.
     */
    if (!mask || !mpz_sgn(c->guard)) return;
    if (mask & 1) mpz_mul(c->guard, c->guard, c->x[slot]);
    if (mask & 2) mpz_mul(c->guard, c->guard, c->z[slot]);
    mpz_mod(c->guard, c->guard, c->n);
}

static uint64_t nanoseconds(b4_context *c, clockid_t clock) {
    unsigned kind = clock == CLOCK_MONOTONIC ? 0 : 1;
    if (c && c->clock_count[kind]) {
        unsigned index = c->clock_index[kind]++;
        if (index >= c->clock_count[kind]) index = c->clock_count[kind] - 1;
        return c->clock_script[kind][index];
    }
    struct timespec value;
    if (clock_gettime(clock, &value)) return 0;
    return (uint64_t)value.tv_sec * 1000000000u + (uint64_t)value.tv_nsec;
}

int b4_test_clocks(b4_context *c, const unsigned long long *wall, size_t nw,
                    const unsigned long long *cpu, size_t nc) {
    if (!c || nw > 16 || nc > 16 || (nw && !wall) || (nc && !cpu)) return -1;
    for (unsigned k = 0; k < 2; ++k) {
        size_t count = k ? nc : nw;
        const unsigned long long *values = k ? cpu : wall;
        for (size_t i = 0; i < count; ++i) c->clock_script[k][i] = values[i];
        c->clock_count[k] = (unsigned)count;
        c->clock_index[k] = 0;
    }
    return 0;
}

unsigned b4_test_clock_reads(const b4_context *c, unsigned kind) {
    return c && kind < 2 ? c->clock_index[kind] : 0;
}

static int step_until(b4_context *c, uint64_t wall_end, uint64_t cpu_end,
                      unsigned actions, const atomic_uint *cancel) {
    uint64_t cpu_wall = 0;
    for (unsigned i = 0; i < actions && !c->done; ++i) {
        if (cancel && atomic_load(cancel)) return -3;
        uint64_t wall = nanoseconds(c, CLOCK_MONOTONIC);
        if (!wall) return -1;
        if (wall >= wall_end) break;
        /* Keep the cheap wall check at every action boundary: the first
         * boundary after an exceptionally slow GMP operation must return.
         * Process CPU is aggregate across the team and relatively expensive;
         * sample at entry and after one ms of wall progress. Root resource
         * polls still run before and after every bounded foreign slice. */
        if (!cpu_wall || wall < cpu_wall || wall - cpu_wall >= 1000000u) {
            uint64_t cpu = nanoseconds(c, CLOCK_PROCESS_CPUTIME_ID);
            if (!cpu) return -1;
            if (cpu >= cpu_end) break;
            cpu_wall = wall;
        }
        if (c->mode == 1) {
            if (mpz_tstbit(c->scalar, (mp_bitcnt_t)c->bit)) {
                ladder_bit(c, 2, 1);
            } else {
                ladder_bit(c, 1, 2);
            }
            if (--c->bit < 0) c->done = 1;
        } else {
            const b4_record *r = &c->records[c->record];
            if (c->starting) {
                guard_point(c, 0, r->masks[0]); c->starting = 0;
            }
            if (c->opcode < r->code_size) {
                const unsigned char *op = r->code + c->opcode;
                if (op[3] == 255) doubling(c, op[0], op[1]);
                else addition(c, op[0], op[1], op[2], op[3]);
                c->opcode += 4;
                guard_point(c, op[0], r->masks[c->opcode / 4]);
            }
            if (c->opcode == r->code_size) {
                mpz_set(c->x[0], c->x[r->output]);
                mpz_set(c->z[0], c->z[r->output]);
                ++c->record; ++c->completed;
                c->opcode = 0; c->starting = 1;
                if (c->record == c->count || !mpz_sgn(c->guard)) c->done = 1;
            }
        }
    }
    return c->done ? 1 : 0;
}

int b4_ready(const b4_context *c) { return c && c->mode; }

int b4_slice(b4_context *c, unsigned long long wall_end,
             unsigned long long cpu_end, unsigned actions, const void *cancel) {
    if (!b4_ready(c) || !actions || actions > 128) return -1;
    return step_until(c, wall_end, cpu_end, actions, cancel);
}

int b4_step(b4_context *c, unsigned milliseconds, unsigned actions) {
    if (!c || !c->mode || !milliseconds || milliseconds > 5 ||
        !actions || actions > 128) return -1;
    uint64_t wall = nanoseconds(c, CLOCK_MONOTONIC);
    uint64_t cpu = nanoseconds(c, CLOCK_PROCESS_CPUTIME_ID);
    if (!wall || !cpu) return -1;
    uint64_t interval = (uint64_t)milliseconds * 1000000u;
    return step_until(c, wall + interval, cpu + interval, actions, NULL);
}

unsigned b4_openmp(void) {
#ifdef B4_OPENMP
    return 1;
#else
    return 0;
#endif
}

int b4_team_step(b4_context *const *contexts, size_t count, unsigned workers,
                 unsigned schedule, unsigned milliseconds, unsigned actions,
                 int *statuses, unsigned *actual_workers) {
    if (!contexts || !statuses || !actual_workers || !count || count > 32 ||
        (workers != 1 && workers != 2 && workers != 4) || schedule > 2 ||
        !milliseconds || milliseconds > 5 || !actions || actions > 128)
        return -1;
    *actual_workers = 0;
    for (size_t i = 0; i < count; ++i) {
        if (!contexts[i] || !contexts[i]->mode) return -1;
        for (size_t j = 0; j < i; ++j)
            if (contexts[i] == contexts[j]) return -1;
        statuses[i] = -1;
    }
#ifndef B4_OPENMP
    return -2;
#else
    /* No nesting, process x team hybrid, global schedule mutation, callbacks
     * or shared mutable GMP objects. A size mismatch refuses before arithmetic.
     * Cold runtime creation and OS scheduling are outside the cooperative
     * five-ms work target; the external owner bounds the whole process.
     */
    if (omp_in_parallel()) return -2;
    uint64_t wall = nanoseconds(NULL, CLOCK_MONOTONIC);
    uint64_t cpu = nanoseconds(NULL, CLOCK_PROCESS_CPUTIME_ID);
    if (!wall || !cpu) return -1;
    uint64_t interval = (uint64_t)milliseconds * 1000000u;
    uint64_t wall_end = wall + interval, cpu_end = cpu + workers * interval;
    unsigned actual = 0;
#pragma omp parallel num_threads(workers) shared(actual)
    {
#pragma omp single
        actual = (unsigned)omp_get_num_threads();
        /* single's barrier makes actual visible before any work dispatch. */
        if (actual == workers) {
            if (schedule == 0) {
#pragma omp for schedule(dynamic, 1)
                for (size_t i = 0; i < count; ++i)
                    statuses[i] = step_until(contexts[i], wall_end, cpu_end, actions, NULL);
            } else if (schedule == 1) {
#pragma omp for schedule(static, 1)
                for (size_t i = 0; i < count; ++i)
                    statuses[i] = step_until(contexts[i], wall_end, cpu_end, actions, NULL);
            } else {
#pragma omp for schedule(guided, 1)
                for (size_t i = 0; i < count; ++i)
                    statuses[i] = step_until(contexts[i], wall_end, cpu_end, actions, NULL);
            }
        }
    }
    *actual_workers = actual;
    return actual == workers ? 0 : -2;
#endif
}

static void export_value(unsigned char *data, size_t *size, const mpz_t value) {
    *size = 0;
    mpz_export(data, size, 1, 1, 1, 0, value);
    if (!*size) { data[0] = 0; *size = 1; }
}

int b4_result(const b4_context *c, b4_output *out) {
    if (!c || !out || !c->mode || !c->done) return -1;
    memset(out, 0, sizeof(*out));
    unsigned slot = c->mode == 1 ? 1 : 0;
    export_value(out->x, &out->x_size, c->x[slot]);
    export_value(out->z, &out->z_size, c->z[slot]);
    export_value(out->guard, &out->guard_size, c->guard);
    out->completed = c->completed;
    return 0;
}
