/* Scoped whole-curve OpenMP owner. Runtime workers never call Python. */
#define _POSIX_C_SOURCE 200809L
#include "packet.h"
#include <stdlib.h>

#ifdef B4_OPENMP
#include <omp.h>
#include <stdatomic.h>
#include <stdint.h>
#include <time.h>

struct b4_omp_team {
    atomic_uint active, cancelled, failed, entered, pause, participants;
    atomic_ullong counters[4];
    int job_at;
    /* Each rank writes only its own slot; region exit drains all writers. */
    struct { uint64_t tasks, curves; unsigned char spacing[112]; }
        activity[B4_NATIVE_MAX_WORKERS];
};

unsigned b4_omp_whole_kernel(void) { return 1; }

b4_omp_team *b4_omp_create(void) {
    b4_omp_team *t = calloc(1, sizeof(*t));
    if (!t) return NULL;
    atomic_init(&t->active, 0); atomic_init(&t->cancelled, 0);
    atomic_init(&t->failed, 0); atomic_init(&t->entered, 0);
    atomic_init(&t->pause, 0); atomic_init(&t->participants, 0);
    for (unsigned i = 0; i < 4; ++i) atomic_init(&t->counters[i], 0);
    t->job_at = -1;
    return t;
}

int b4_omp_shutdown(b4_omp_team *t) {
    return !t || atomic_load(&t->active) ? -1 : 0;
}

void b4_omp_destroy(b4_omp_team *t) {
    if (!t) return;
    if (b4_omp_shutdown(t)) abort();
    /* No implicit task can reference the owner after the region barrier.
     * Runtime-owned idle workers and the DSO have process lifetime. */
    free(t);
}

unsigned long long b4_omp_counter(b4_omp_team *t, unsigned kind) {
    return !t || kind >= 4 ? 0 : atomic_load(&t->counters[kind]);
}

unsigned b4_omp_participants(b4_omp_team *t) {
    return t ? atomic_load(&t->participants) : 0;
}

unsigned long long b4_omp_activity(b4_omp_team *t, unsigned rank,
                                    unsigned field) {
    if (!t || rank >= B4_NATIVE_MAX_WORKERS || field > 1 ||
        atomic_load(&t->active)) return UINT64_MAX;
    return field ? t->activity[rank].curves : t->activity[rank].tasks;
}

void b4_omp_cancel(b4_omp_team *t) {
    if (t) atomic_store(&t->cancelled, 1);
}

int b4_omp_reset(b4_omp_team *t) {
    if (!t || atomic_load(&t->active)) return -1;
    atomic_store(&t->cancelled, 0);
    return 0;
}

int b4_omp_inject(b4_omp_team *t, int create_at, int job_at, unsigned pause) {
    /* libomp owns thread creation: actual constrained teams test refusal.
     * A pthread creation-rank injection cannot impersonate that mechanism. */
    if (!t || atomic_load(&t->active) || create_at != -1 ||
        job_at < -1 || job_at > 31 || pause > 1) return -1;
    t->job_at = job_at;
    atomic_store(&t->pause, pause);
    atomic_store(&t->entered, 0);
    return 0;
}

static void pause_briefly(void) {
    struct timespec delay = {0, 50000};
    nanosleep(&delay, NULL);
}

void b4_omp_wait(b4_omp_team *t, unsigned waiting) {
    if (!t || !waiting || waiting > B4_NATIVE_MAX_WORKERS) return;
    while (atomic_load(&t->entered) < waiting && !atomic_load(&t->cancelled))
        pause_briefly();
}

void b4_omp_release(b4_omp_team *t) {
    if (t) atomic_store(&t->pause, 0);
}

static uint64_t clock_ns(clockid_t id) {
    struct timespec ts;
    if (clock_gettime(id, &ts)) return 0;
    return (uint64_t)ts.tv_sec * 1000000000u + (uint64_t)ts.tv_nsec;
}

static void execute_one(b4_omp_team *t, b4_context *context, size_t index,
                        int *status, uint64_t wall_end, uint64_t cpu_end,
                        unsigned actions) {
    if (atomic_load(&t->cancelled)) return;
    int result = (int)index == t->job_at ? -4 :
        b4_slice(context, wall_end, cpu_end, actions, &t->cancelled);
    unsigned rank = (unsigned)omp_get_thread_num();
    ++t->activity[rank].tasks;
    if (result == 1) ++t->activity[rank].curves;
    *status = result;
    if (result < 0 && result != -3) {
        atomic_store(&t->failed, 1);
        atomic_store(&t->cancelled, 1);
    }
}

int b4_omp_step(b4_omp_team *t, b4_context *const *contexts, size_t count,
                unsigned workers, unsigned schedule, unsigned chunk,
                unsigned milliseconds, unsigned actions, int *statuses,
                unsigned *actual_workers) {
    if (actual_workers) *actual_workers = 0;
    if (!t || !actual_workers || count > 32 ||
        (count && (!contexts || !statuses)) || !workers || workers > B4_NATIVE_MAX_WORKERS ||
        schedule > 1 || !chunk || chunk > 32 || !milliseconds || milliseconds > 5 ||
        !actions || actions > b4_slice_actions() || omp_get_level()) return -1;
    unsigned inactive = 0;
    if (!atomic_compare_exchange_strong(&t->active, &inactive, 1)) return -1;
    for (size_t i = 0; i < count; ++i) {
        if (!b4_ready(contexts[i])) { atomic_store(&t->active, 0); return -1; }
        for (size_t j = 0; j < i; ++j)
            if (contexts[i] == contexts[j]) { atomic_store(&t->active, 0); return -1; }
        statuses[i] = -3;
    }
    if (!count || atomic_load(&t->cancelled)) {
        int result = atomic_load(&t->cancelled) ? -3 : 0;
        atomic_store(&t->active, 0);
        return result;
    }
    atomic_store(&t->failed, 0); atomic_store(&t->entered, 0);
    atomic_store(&t->participants, 0);
    unsigned actual = 0;
    uint64_t wall_end = 0, cpu_end = 0;
#pragma omp parallel num_threads(workers) default(none) \
    shared(t, actual, wall_end, cpu_end, contexts, statuses, count, workers, \
           schedule, chunk, milliseconds, actions)
    {
        unsigned rank = (unsigned)omp_get_thread_num();
        atomic_fetch_or(&t->participants, 1u << rank);
        atomic_fetch_add(&t->entered, 1);
        while (atomic_load(&t->pause) && !atomic_load(&t->cancelled)) pause_briefly();
#pragma omp single
        {
            actual = (unsigned)omp_get_num_threads();
            atomic_fetch_add(&t->counters[3], actual);
            if (actual == workers) {
                uint64_t wall = clock_ns(CLOCK_MONOTONIC);
                uint64_t cpu = clock_ns(CLOCK_PROCESS_CPUTIME_ID);
                if (!wall || !cpu) {
                    atomic_store(&t->failed, 1); atomic_store(&t->cancelled, 1);
                } else {
                    wall_end = wall + (uint64_t)milliseconds * 1000000u;
                    cpu_end = cpu + (uint64_t)workers * milliseconds * 1000000u;
                    atomic_fetch_add(&t->counters[2], 1);
                }
            }
        }
        /* The single barrier makes actual team size and deadlines visible.
         * All threads take the same branch; mismatch executes no kernel. */
        if (actual == workers && wall_end && cpu_end) {
            if (!schedule) {
#pragma omp for schedule(dynamic, chunk) nowait
                for (size_t i = 0; i < count; ++i)
                    execute_one(t, contexts[i], i, &statuses[i], wall_end, cpu_end, actions);
            } else {
#pragma omp for schedule(static, chunk) nowait
                for (size_t i = 0; i < count; ++i)
                    execute_one(t, contexts[i], i, &statuses[i], wall_end, cpu_end, actions);
            }
        }
        /* The parallel region's implicit barrier drains every borrower. */
    }
    *actual_workers = actual;
    int result = actual != workers ? -2 : atomic_load(&t->failed) ? -4 :
        atomic_load(&t->cancelled) ? -3 : 0;
    atomic_store(&t->active, 0);
    return result;
}
#else
unsigned b4_omp_whole_kernel(void) { return 0; }
b4_omp_team *b4_omp_create(void) { return NULL; }
int b4_omp_shutdown(b4_omp_team *t) { (void)t; return -2; }
void b4_omp_destroy(b4_omp_team *t) { (void)t; }
unsigned long long b4_omp_counter(b4_omp_team *t, unsigned k) {
    (void)t; (void)k; return 0;
}
unsigned b4_omp_participants(b4_omp_team *t) { (void)t; return 0; }
unsigned long long b4_omp_activity(b4_omp_team *t, unsigned r, unsigned f) {
    (void)t; (void)r; (void)f; return ~0ull;
}
void b4_omp_cancel(b4_omp_team *t) { (void)t; }
int b4_omp_reset(b4_omp_team *t) { (void)t; return -2; }
int b4_omp_inject(b4_omp_team *t, int c, int j, unsigned p) {
    (void)t; (void)c; (void)j; (void)p; return -2;
}
void b4_omp_wait(b4_omp_team *t, unsigned n) { (void)t; (void)n; }
void b4_omp_release(b4_omp_team *t) { (void)t; }
int b4_omp_step(b4_omp_team *t, b4_context *const *c, size_t n, unsigned w,
                unsigned s, unsigned k, unsigned m, unsigned a, int *r, unsigned *v) {
    (void)t; (void)c; (void)n; (void)w; (void)s; (void)k; (void)m; (void)a; (void)r;
    if (v) *v = 0;
    return -2;
}
#endif
