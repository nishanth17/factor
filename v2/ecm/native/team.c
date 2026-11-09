/* Owned pthread dispatch. Threads borrow distinct GMP contexts until joined. */
#define _POSIX_C_SOURCE 200809L
#include "packet.h"
#include <stdlib.h>
#ifdef B4_PTHREAD
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <time.h>

struct b4_pthread_team {
    pthread_mutex_t mutex;
    pthread_cond_t condition;
    atomic_uint cancelled, failed;
    atomic_size_t next;
    unsigned active, started, waiting, pause;
    int create_at, job_at;
    b4_context *const *contexts;
    int *statuses;
    size_t count;
    unsigned workers, schedule, chunk, actions;
    uint64_t wall_end, cpu_end;
};

typedef struct { b4_pthread_team *team; unsigned rank; } worker_arg;

unsigned b4_pthread(void) { return 1; }

b4_pthread_team *b4_pthread_create(void) {
    b4_pthread_team *t = calloc(1, sizeof(*t));
    if (!t) return NULL;
    if (pthread_mutex_init(&t->mutex, NULL)) { free(t); return NULL; }
    if (pthread_cond_init(&t->condition, NULL)) {
        pthread_mutex_destroy(&t->mutex); free(t); return NULL;
    }
    atomic_init(&t->cancelled, 0);
    atomic_init(&t->failed, 0);
    atomic_init(&t->next, 0);
    t->create_at = t->job_at = -1;
    return t;
}

void b4_pthread_destroy(b4_pthread_team *t) {
    if (!t) return;
    /* The Python owner forbids close during dispatch. All workers have joined. */
    pthread_cond_destroy(&t->condition);
    pthread_mutex_destroy(&t->mutex);
    free(t);
}

void b4_pthread_cancel(b4_pthread_team *t) {
    if (!t) return;
    atomic_store(&t->cancelled, 1);
    pthread_mutex_lock(&t->mutex);
    pthread_cond_broadcast(&t->condition);
    pthread_mutex_unlock(&t->mutex);
}

int b4_pthread_reset(b4_pthread_team *t) {
    if (!t) return -1;
    pthread_mutex_lock(&t->mutex);
    int result = t->active ? -1 : 0;
    if (!result) atomic_store(&t->cancelled, 0);
    pthread_mutex_unlock(&t->mutex);
    return result;
}

int b4_pthread_inject(b4_pthread_team *t, int create_at, int job_at,
                      unsigned pause) {
    if (!t || create_at < -1 || create_at > 3 || job_at < -1 ||
        job_at > 31 || pause > 1) return -1;
    pthread_mutex_lock(&t->mutex);
    int result = t->active ? -1 : 0;
    if (!result) {
        t->create_at = create_at; t->job_at = job_at; t->pause = pause;
    }
    pthread_mutex_unlock(&t->mutex);
    return result;
}

void b4_pthread_wait(b4_pthread_team *t, unsigned waiting) {
    if (!t || !waiting || waiting > 4) return;
    pthread_mutex_lock(&t->mutex);
    while (t->waiting < waiting && !atomic_load(&t->cancelled))
        pthread_cond_wait(&t->condition, &t->mutex);
    pthread_mutex_unlock(&t->mutex);
}

void b4_pthread_release(b4_pthread_team *t) {
    if (!t) return;
    pthread_mutex_lock(&t->mutex);
    t->pause = 0;
    pthread_cond_broadcast(&t->condition);
    pthread_mutex_unlock(&t->mutex);
}

static void *worker(void *opaque) {
    worker_arg *arg = opaque;
    b4_pthread_team *t = arg->team;
    pthread_mutex_lock(&t->mutex);
    ++t->waiting;
    pthread_cond_broadcast(&t->condition);
    while ((!t->started || t->pause) && !atomic_load(&t->cancelled))
        pthread_cond_wait(&t->condition, &t->mutex);
    pthread_mutex_unlock(&t->mutex);

    size_t start = (size_t)arg->rank * t->chunk;
    while (!atomic_load(&t->cancelled)) {
        if (!t->schedule) start = atomic_fetch_add(&t->next, t->chunk);
        if (start >= t->count) break;
        size_t end = start + t->chunk;
        if (end > t->count) end = t->count;
        for (size_t i = start; i < end; ++i) {
            if (atomic_load(&t->cancelled)) break;
            int result = (int)i == t->job_at ? -4 :
                b4_slice(t->contexts[i], t->wall_end, t->cpu_end,
                         t->actions, &t->cancelled);
            t->statuses[i] = result;
            if (result < 0 && result != -3) {
                atomic_store(&t->failed, 1);
                b4_pthread_cancel(t);
                break;
            }
        }
        if (t->schedule) start += (size_t)t->workers * t->chunk;
    }
    return NULL;
}

static uint64_t clock_ns(clockid_t id) {
    struct timespec ts;
    if (clock_gettime(id, &ts)) return 0;
    return (uint64_t)ts.tv_sec * 1000000000u + (uint64_t)ts.tv_nsec;
}

int b4_pthread_step(b4_pthread_team *t, b4_context *const *contexts,
                    size_t count, unsigned workers, unsigned schedule,
                    unsigned chunk, unsigned milliseconds, unsigned actions,
                    int *statuses, unsigned *actual) {
    if (actual) *actual = 0;
    if (!t || !actual || count > 32 ||
        (count && (!contexts || !statuses)) ||
        (workers != 1 && workers != 2 && workers != 4) || schedule > 1 ||
        !chunk || chunk > 32 || !milliseconds || milliseconds > 5 ||
        !actions || actions > 128) return -1;
    for (size_t i = 0; i < count; ++i) {
        if (!b4_ready(contexts[i])) return -1;
        for (size_t j = 0; j < i; ++j)
            if (contexts[i] == contexts[j]) return -1;
        statuses[i] = -3;
    }
    pthread_mutex_lock(&t->mutex);
    if (t->active) { pthread_mutex_unlock(&t->mutex); return -1; }
    t->active = 1; t->started = t->waiting = 0;
    t->contexts = contexts; t->statuses = statuses; t->count = count;
    t->workers = workers; t->schedule = schedule; t->chunk = chunk;
    t->actions = actions;
    atomic_store(&t->next, 0); atomic_store(&t->failed, 0);
    pthread_mutex_unlock(&t->mutex);

    pthread_t threads[4];
    worker_arg args[4];
    unsigned created = 0;
    int result = 0;
    if (!count || atomic_load(&t->cancelled)) goto finish;
    for (unsigned i = 0; i < workers; ++i) {
        args[i].team = t; args[i].rank = i;
        if ((int)i == t->create_at ||
            pthread_create(&threads[i], NULL, worker, &args[i])) {
            result = -2;
            b4_pthread_cancel(t);
            break;
        }
        ++created;
    }
    /* No arithmetic before the complete team exists. Set shared deadlines
     * after creation so short work is not starved by thread startup. CPU is
     * the process total, already including all native threads exactly once. */
    uint64_t wall = clock_ns(CLOCK_MONOTONIC);
    uint64_t cpu = clock_ns(CLOCK_PROCESS_CPUTIME_ID);
    if (!wall || !cpu) { result = -4; b4_pthread_cancel(t); }
    t->wall_end = wall + (uint64_t)milliseconds * 1000000u;
    t->cpu_end = cpu + (uint64_t)workers * milliseconds * 1000000u;
    pthread_mutex_lock(&t->mutex);
    t->started = 1;
    pthread_cond_broadcast(&t->condition);
    pthread_mutex_unlock(&t->mutex);
    for (unsigned i = 0; i < created; ++i) {
        /* Valid, joinable, exclusively owned handles cannot return EINVAL,
         * ESRCH or EDEADLK. If that invariant fails, returning borrowed data
         * to Python would be unsafe; fail the process, never free live state. */
        if (pthread_join(threads[i], NULL)) abort();
    }
finish:
    *actual = created;
    if (!result && atomic_load(&t->failed)) result = -4;
    if (!result && atomic_load(&t->cancelled)) result = -3;
    pthread_mutex_lock(&t->mutex);
    t->active = 0; t->contexts = NULL; t->statuses = NULL;
    pthread_mutex_unlock(&t->mutex);
    return result;
}
#else
unsigned b4_pthread(void) { return 0; }
b4_pthread_team *b4_pthread_create(void) { return NULL; }
void b4_pthread_destroy(b4_pthread_team *t) { (void)t; }
void b4_pthread_cancel(b4_pthread_team *t) { (void)t; }
int b4_pthread_reset(b4_pthread_team *t) { (void)t; return -2; }
int b4_pthread_inject(b4_pthread_team *t, int c, int j, unsigned p) {
    (void)t; (void)c; (void)j; (void)p; return -2;
}
void b4_pthread_wait(b4_pthread_team *t, unsigned n) { (void)t; (void)n; }
void b4_pthread_release(b4_pthread_team *t) { (void)t; }
int b4_pthread_step(b4_pthread_team *t, b4_context *const *c, size_t n,
                    unsigned w, unsigned s, unsigned k, unsigned m, unsigned a,
                    int *r, unsigned *actual) {
    (void)t; (void)c; (void)n; (void)w; (void)s; (void)k; (void)m;
    (void)a; (void)r;
    if (actual) *actual = 0;
    return -2;
}
#endif
