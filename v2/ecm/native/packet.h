/* Factor's optional scalar ECM packet. All integers use canonical unsigned BE. */
#ifndef FACTOR_ECM_PACKET_H
#define FACTOR_ECM_PACKET_H
#include <stddef.h>

typedef struct b4_context b4_context;
typedef struct { const unsigned char *data; size_t size; } b4_value;
typedef struct {
    const unsigned char *code;
    const unsigned char *masks;
    size_t code_size, masks_size;
    unsigned slots, output;
} b4_record;
typedef struct {
    unsigned char x[128], z[128], guard[128];
    size_t x_size, z_size, guard_size;
    unsigned completed;
} b4_output;

unsigned b4_abi(void);
unsigned b4_limb_division(void);
const char *b4_gmp_version(void);
typedef struct {
    unsigned long long birth, rss, user, system;
} b4_root_observation;
/* 0=observed childless root, 1=use full tree provider, negative=errno. */
int b4_observe_root(int pid, b4_root_observation *output);
b4_context *b4_create(void);
void b4_destroy(b4_context *ctx);
void b4_abort(b4_context *ctx);
int b4_ladder(b4_context *ctx, unsigned abi, b4_value n, b4_value x,
              b4_value z, b4_value a24, b4_value scalar);
int b4_records(b4_context *ctx, unsigned abi, b4_value n, b4_value x,
               b4_value z, b4_value a24, const b4_record *records, size_t count);
/* 0=yield, 1=done, -1=invalid; clocks include this process's native threads. */
int b4_step(b4_context *ctx, unsigned milliseconds, unsigned actions);
/* Internal C-team boundary; cancel is an optional C11 atomic unsigned. */
int b4_slice(b4_context *ctx, unsigned long long wall_end,
             unsigned long long cpu_end, unsigned actions, const void *cancel);
int b4_ready(const b4_context *ctx);
/* Per-context deterministic test clocks, copied and bounded; zero counts
 * restore real clocks. No Python callbacks or global clock replacement. */
int b4_test_clocks(b4_context *ctx, const unsigned long long *wall, size_t nw,
                    const unsigned long long *cpu, size_t nc);
unsigned b4_test_clock_reads(const b4_context *ctx, unsigned kind);
typedef struct b4_pthread_team b4_pthread_team;
unsigned b4_pthread(void);
b4_pthread_team *b4_pthread_create(void);
int b4_pthread_shutdown(b4_pthread_team *team);
void b4_pthread_destroy(b4_pthread_team *team);
/* Creates, joins, dispatched generations and worker generation entries. */
unsigned long long b4_pthread_counter(b4_pthread_team *team, unsigned kind);
void b4_pthread_cancel(b4_pthread_team *team);
int b4_pthread_reset(b4_pthread_team *team);
/* Dynamic=0, static cyclic chunks=1. Each generation drains before return;
 * fixed workers join on shutdown/close, including partial creation failure.
 * -2=create failure/unavailable, -3=cancelled, -4=worker failure. */
int b4_pthread_step(b4_pthread_team *team, b4_context *const *contexts,
                    size_t count, unsigned workers, unsigned schedule,
                    unsigned chunk, unsigned milliseconds, unsigned actions,
                    int *statuses, unsigned *actual_workers);
/* Deterministic per-owner fault/gate controls; -1 disables each failure.
 * No process-global hooks and no Python callback inside a C worker. */
int b4_pthread_inject(b4_pthread_team *team, int create_at, int job_at,
                      unsigned pause);
void b4_pthread_wait(b4_pthread_team *team, unsigned waiting);
void b4_pthread_release(b4_pthread_team *team);
/* schedule: 0=dynamic,1=static,2=guided, all chunk1. Scalar builds refuse.
 * All contexts must be distinct. Return only after the OpenMP region drains.
 * -2=team unavailable/size mismatch, -1=invalid, 0=checked statuses written.
 */
unsigned b4_openmp(void);
int b4_team_step(b4_context *const *contexts, size_t count, unsigned workers,
                 unsigned schedule, unsigned milliseconds, unsigned actions,
                 int *statuses, unsigned *actual_workers);
int b4_result(const b4_context *ctx, b4_output *output);
/* Whole-curve ordinary-mpz owner. Metadata is immutable and copied once. */
typedef struct b4_prime_block b4_prime_block;
typedef struct {
    unsigned long long prime, power, reserve;
    b4_record fast, unit;
} b4_curve_unit;
typedef struct {
    unsigned phase, done, setup_done, chain_chunks, chain_replays, boundary;
    unsigned need, chain, strict_state, operation;
    unsigned long long b1, b2, center, distance, replay_power;
    unsigned long long left, next, hi, work, setup_ops, stage_one_ops, stage_two_ops;
    size_t index, powers_count, baby_count, terms_count, replay_index, recovery;
    unsigned long long phase_ns[3];
} b4_curve_info;
typedef struct { unsigned char data[128]; size_t size; } b4_atom;
/* Exactness hook: 1=limb path, 0=ordinary fallback, -1=invalid input. */
int b4_test_product(b4_context *ctx, b4_value a, b4_value b,
                     unsigned negative_a, unsigned negative_b,
                     unsigned alias, b4_atom *output);
unsigned b4_whole_ecm(void);
b4_prime_block *b4_block_create(unsigned long long lo, unsigned long long hi,
                               unsigned long long bound, unsigned chain,
                               const b4_curve_unit *units, size_t count);
int b4_block_destroy(b4_prime_block *block);
int b4_block_prepare(b4_prime_block *block, unsigned long long b1,
                      unsigned long long b2);
int b4_curve_begin(b4_context *ctx, b4_value n, b4_value sigma,
                    unsigned long long b1, unsigned long long b2,
                    unsigned chunk, unsigned gcd_batch, unsigned segment);
int b4_curve_feed(b4_context *ctx, b4_prime_block *block);
int b4_curve_info_get(const b4_context *ctx, b4_curve_info *info);
/* Canonical fields: a24, value X/Z, replay X/Z, product, factor,
 * baby X/Z, terms, giant X/Z, previous X/Z. No partial arithmetic exports. */
int b4_curve_get(const b4_context *ctx, unsigned field, size_t index, b4_atom *out);
int b4_curve_power(const b4_context *ctx, size_t index, unsigned long long *pair);
/* Restore complete canonical actions; speculative scratch is discarded. */
int b4_curve_restore(b4_context *ctx, const b4_curve_info *info,
                      const b4_value *fields, const b4_value *baby,
                      const b4_value *terms, b4_prime_block *powers,
                      b4_prime_block *block);
int b4_curve_options(b4_context *ctx, unsigned atomic_actions, unsigned phase_clocks);
/* Inert per-owner deterministic failure/cancellation hooks. Zero disables. */
int b4_curve_faults(b4_context *ctx, unsigned preparation_failure,
                     unsigned cancel_after_actions);
void b4_curve_cancel(b4_context *ctx);
#endif
