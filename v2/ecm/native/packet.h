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
const char *b4_gmp_version(void);
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
typedef struct b4_pthread_team b4_pthread_team;
unsigned b4_pthread(void);
b4_pthread_team *b4_pthread_create(void);
void b4_pthread_destroy(b4_pthread_team *team);
void b4_pthread_cancel(b4_pthread_team *team);
int b4_pthread_reset(b4_pthread_team *team);
/* Dynamic=0, static cyclic chunks=1. All threads join before return.
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
#endif
