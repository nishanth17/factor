/* One owned persistent pthread team; borrowed contexts drain each generation. */
#define _POSIX_C_SOURCE 200809L
#include "packet.h"
#include <stdlib.h>
#ifdef B4_PTHREAD
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <time.h>

#define MAX_WORKERS B4_NATIVE_MAX_WORKERS

typedef struct {
    uint64_t cpu, wall, one, two, slices, errors, start, finish;
    uint64_t one_cpu, two_cpu, mixed_cpu, unsupported;
    uint64_t tasks, curves, identity;
    unsigned char spacing[136];
} b4_worker_metric;
_Static_assert(sizeof(b4_worker_metric) == 256, "separate rank writes");
_Static_assert(MAX_WORKERS <= 16, "participant mask capacity");

static uint64_t clock_ns(clockid_t id);

typedef struct {
    b4_pthread_team *team;
    unsigned rank;
    unsigned long long generation;
} worker_arg;

struct b4_pthread_team {
    pthread_mutex_t mutex;
    pthread_cond_t condition;
    pthread_t threads[MAX_WORKERS];
    worker_arg args[MAX_WORKERS];
    atomic_uint cancelled, failed;
    atomic_size_t next;
    unsigned active, live, completed, waiting, pause, stop, observers;
    unsigned long long generation, counters[4];
    int create_at, job_at;
    b4_context *const *contexts;
    int *statuses;
    size_t count;
    unsigned workers, schedule, chunk, actions;
    uint64_t wall_end, cpu_end;
    b4_worker_metric metrics[MAX_WORKERS];
    uint64_t metric_generations, metric_span, metric_gap, metric_errors;
    uint64_t metric_small_generations, metric_small_span, metric_full_span;
    unsigned diagnostic_enabled, participants;
};

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

/* Only idle owners may shut down. Join before changing handles/arguments or
 * releasing the DSO; parked workers still execute code from that DSO. */
int b4_pthread_shutdown(b4_pthread_team *t) {
    if (!t) return -1;
    pthread_mutex_lock(&t->mutex);
    if (t->active) { pthread_mutex_unlock(&t->mutex); return -1; }
    t->stop = 1;
    pthread_cond_broadcast(&t->condition);
    unsigned live = t->live;
    pthread_mutex_unlock(&t->mutex);
    for (unsigned i = 0; i < live; ++i) {
        if (pthread_join(t->threads[i], NULL)) abort();
    }
    pthread_mutex_lock(&t->mutex);
    t->counters[1] += live;
    t->live = t->workers = t->waiting = t->stop = 0;
    pthread_mutex_unlock(&t->mutex);
    return 0;
}

void b4_pthread_destroy(b4_pthread_team *t) {
    if (!t) return;
    if (b4_pthread_shutdown(t)) abort();
    pthread_cond_destroy(&t->condition);
    pthread_mutex_destroy(&t->mutex);
    free(t);
}

unsigned long long b4_pthread_counter(b4_pthread_team *t, unsigned kind) {
    if (!t || kind >= 4) return 0;
    pthread_mutex_lock(&t->mutex);
    unsigned long long value = t->counters[kind];
    pthread_mutex_unlock(&t->mutex);
    return value;
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
    if (!t || create_at < -1 || create_at >= MAX_WORKERS || job_at < -1 ||
        job_at > 31 || pause > 1) return -1;
    /* A creation injection explicitly rearms startup on this idle owner.
     * Normal dispatch never destroys/recreates an already-live team. */
    if (create_at >= 0 && b4_pthread_shutdown(t)) return -1;
    pthread_mutex_lock(&t->mutex);
    int result = t->active ? -1 : 0;
    if (!result) {
        t->create_at = create_at; t->job_at = job_at; t->pause = pause;
        t->waiting = 0;
    }
    pthread_mutex_unlock(&t->mutex);
    return result;
}

void b4_pthread_wait(b4_pthread_team *t, unsigned waiting) {
    if (!t || !waiting || waiting > MAX_WORKERS) return;
    pthread_mutex_lock(&t->mutex);
    ++t->observers;
    while (t->waiting < waiting && !atomic_load(&t->cancelled))
        pthread_cond_wait(&t->condition, &t->mutex);
    --t->observers;
    pthread_mutex_unlock(&t->mutex);
}

void b4_pthread_release(b4_pthread_team *t) {
    if (!t) return;
    pthread_mutex_lock(&t->mutex);
    t->pause = 0;
    pthread_cond_broadcast(&t->condition);
    pthread_mutex_unlock(&t->mutex);
}

static int diagnostic_slice(b4_pthread_team *t, unsigned rank,
                            b4_context *context) {
    b4_curve_info before, after;
    b4_worker_metric *m = &t->metrics[rank];
    int valid = !b4_curve_info_get(context, &before);
    uint64_t started = clock_ns(CLOCK_THREAD_CPUTIME_ID);
    int result = b4_slice(context, t->wall_end, t->cpu_end,
                          t->actions, &t->cancelled);
    uint64_t finished = clock_ns(CLOCK_THREAD_CPUTIME_ID);
    ++m->slices;
    if (!started || finished < started) ++m->errors;
    if (valid && !b4_curve_info_get(context, &after)) {
        if (after.stage_one_ops < before.stage_one_ops ||
            after.stage_two_ops < before.stage_two_ops) ++m->errors;
        else {
            uint64_t one = after.stage_one_ops - before.stage_one_ops;
            uint64_t two = after.stage_two_ops - before.stage_two_ops;
            m->one += one; m->two += two;
            uint64_t elapsed = finished >= started ? finished - started : 0;
            if (one && !two) m->one_cpu += elapsed;
            else if (two && !one) m->two_cpu += elapsed;
            else m->mixed_cpu += elapsed;
        }
    } else ++m->unsupported;
    return result;
}

static void work(b4_pthread_team *t, unsigned rank) {
    size_t start = (size_t)rank * t->chunk;
    while (!atomic_load(&t->cancelled)) {
        if (!t->schedule) start = atomic_fetch_add(&t->next, t->chunk);
        if (start >= t->count) break;
        size_t end = start + t->chunk;
        if (end > t->count) end = t->count;
        for (size_t i = start; i < end; ++i) {
            if (atomic_load(&t->cancelled)) break;
            int result = (int)i == t->job_at ? -4 :
                (t->diagnostic_enabled && !(t->generation & 7)) ? diagnostic_slice(t, rank, t->contexts[i]) :
                b4_slice(t->contexts[i], t->wall_end, t->cpu_end,
                         t->actions, &t->cancelled);
            ++t->metrics[rank].tasks;
            if (result == 1) ++t->metrics[rank].curves;
            t->statuses[i] = result;
            if (result < 0 && result != -3) {
                atomic_store(&t->failed, 1);
                b4_pthread_cancel(t);
                break;
            }
        }
        if (t->schedule) start += (size_t)t->workers * t->chunk;
    }
}

static void *worker(void *opaque) {
    worker_arg *arg = opaque;
    b4_pthread_team *t = arg->team;
    unsigned long long generation = arg->generation;
    pthread_mutex_lock(&t->mutex);
    t->metrics[arg->rank].identity = (uint64_t)(uintptr_t)pthread_self();
    for (;;) {
        while (!t->stop && generation == t->generation)
            pthread_cond_wait(&t->condition, &t->mutex);
        if (t->stop) break;
        generation = t->generation;
        ++t->waiting; ++t->counters[3];
        t->participants |= 1u << arg->rank;
        /* Only explicit entry observers need this intermediate notification.
         * The normal caller waits for every borrower to drain. */
        if (t->observers) pthread_cond_broadcast(&t->condition);
        while (t->pause && !atomic_load(&t->cancelled))
            pthread_cond_wait(&t->condition, &t->mutex);
        pthread_mutex_unlock(&t->mutex);
        unsigned metered = t->diagnostic_enabled && !(generation & 7);
        uint64_t cpu_at = metered ? clock_ns(CLOCK_THREAD_CPUTIME_ID) : 0;
        uint64_t wall_at = metered ? clock_ns(CLOCK_MONOTONIC) : 0;
        work(t, arg->rank);
        if (metered) {
            b4_worker_metric *m = &t->metrics[arg->rank];
            uint64_t wall_now = clock_ns(CLOCK_MONOTONIC);
            uint64_t cpu_now = clock_ns(CLOCK_THREAD_CPUTIME_ID);
            if (!wall_at || !cpu_at || wall_now < wall_at || cpu_now < cpu_at)
                ++m->errors;
            else { m->cpu += cpu_now-cpu_at; m->wall += wall_now-wall_at; }
            m->start=wall_at; m->finish=wall_now;
        }
        pthread_mutex_lock(&t->mutex);
        ++t->completed;
        if (t->completed == t->workers)
            pthread_cond_broadcast(&t->condition);
    }
    pthread_mutex_unlock(&t->mutex);
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
        (!workers || workers > MAX_WORKERS) || schedule > 1 ||
        !chunk || chunk > 32 || !milliseconds || milliseconds > 5 ||
        !actions || actions > b4_slice_actions()) return -1;
    for (size_t i = 0; i < count; ++i) {
        if (!b4_ready(contexts[i])) return -1;
        for (size_t j = 0; j < i; ++j)
            if (contexts[i] == contexts[j]) return -1;
        statuses[i] = -3;
    }
    pthread_mutex_lock(&t->mutex);
    if (t->active || (t->live && t->workers != workers)) {
        pthread_mutex_unlock(&t->mutex); return -1;
    }
    if (!count || atomic_load(&t->cancelled)) {
        int result = atomic_load(&t->cancelled) ? -3 : 0;
        pthread_mutex_unlock(&t->mutex); return result;
    }
    t->active = 1;
    /* Workers cannot see a new generation until all requested handles exist.
     * Partial creation failure publishes no arithmetic and joins its subset. */
    if (!t->live) {
        t->workers = workers;
        for (unsigned i = 0; i < workers; ++i) {
            t->args[i] = (worker_arg){t, i, t->generation};
            if ((int)i == t->create_at ||
                pthread_create(&t->threads[i], NULL, worker, &t->args[i])) {
                *actual = t->live;
                atomic_store(&t->cancelled, 1);
                t->active = 0;
                pthread_mutex_unlock(&t->mutex);
                if (b4_pthread_shutdown(t)) abort();
                return -2;
            }
            ++t->live; ++t->counters[0];
        }
    }
    *actual = t->live;
    uint64_t wall = clock_ns(CLOCK_MONOTONIC);
    uint64_t cpu = clock_ns(CLOCK_PROCESS_CPUTIME_ID);
    if (!wall || !cpu) {
        t->active = 0;
        pthread_mutex_unlock(&t->mutex);
        return -4;
    }
    t->contexts = contexts; t->statuses = statuses; t->count = count;
    t->schedule = schedule; t->chunk = chunk; t->actions = actions;
    t->wall_end = wall + (uint64_t)milliseconds * 1000000u;
    t->cpu_end = cpu + (uint64_t)workers * milliseconds * 1000000u;
    t->completed = t->waiting = t->participants = 0;
    atomic_store(&t->next, 0); atomic_store(&t->failed, 0);
    ++t->generation; ++t->counters[2];
    pthread_cond_broadcast(&t->condition);
    while (t->completed < t->live)
        pthread_cond_wait(&t->condition, &t->mutex);
    if (t->diagnostic_enabled && !(t->generation & 7)) {
        uint64_t first = UINT64_MAX, earliest_end = UINT64_MAX, last = 0;
        unsigned valid = 1;
        for (unsigned rank=0; rank<t->workers; ++rank) {
            b4_worker_metric *m = &t->metrics[rank];
            if (!m->start || m->finish < m->start) valid = 0;
            if (m->start < first) first = m->start;
            if (m->finish < earliest_end) earliest_end = m->finish;
            if (m->finish > last) last = m->finish;
        }
        if (valid) {
            ++t->metric_generations; t->metric_span += last-first;
            t->metric_gap += last-earliest_end;
            if (t->count < t->workers) {
                ++t->metric_small_generations;
                t->metric_small_span += last-first;
            } else t->metric_full_span += last-first;
        } else ++t->metric_errors;
    }
    /* Every worker has relinquished this generation's borrowed buffers.
     * They may remain parked, but cannot use contexts after this barrier. */
    int result = atomic_load(&t->failed) ? -4 :
        atomic_load(&t->cancelled) ? -3 : 0;
    t->active = 0; t->contexts = NULL; t->statuses = NULL;
    pthread_mutex_unlock(&t->mutex);
    return result;
}
#else
unsigned b4_pthread(void) { return 0; }
b4_pthread_team *b4_pthread_create(void) { return NULL; }
int b4_pthread_shutdown(b4_pthread_team *t) { (void)t; return -2; }
void b4_pthread_destroy(b4_pthread_team *t) { (void)t; }
unsigned long long b4_pthread_counter(b4_pthread_team *t, unsigned k) {
    (void)t; (void)k; return 0;
}
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

#ifdef B4_PTHREAD
unsigned long long b4_pthread_diagnostic(b4_pthread_team *t,
                                         unsigned rank, unsigned field) {
    if (!t || rank > MAX_WORKERS || field > 14) return UINT64_MAX;
    pthread_mutex_lock(&t->mutex);
    uint64_t result = UINT64_MAX;
    if (!t->active) {
        if (rank == MAX_WORKERS) {
            switch (field) {
            case 0: result = t->metric_generations; break;
            case 1: result = t->metric_span; break;
            case 2: result = t->metric_gap; break;
            case 3: result = t->metric_errors; break;
            case 4: result = t->metric_small_generations; break;
            case 5: result = t->metric_small_span; break;
            case 6: result = t->metric_full_span; break;
            default: break;
            }
        } else {
            b4_worker_metric *m = &t->metrics[rank];
            switch (field) {
            case 0: result=m->cpu; break; case 1: result=m->wall; break;
            case 2: result=m->one; break; case 3: result=m->two; break;
            case 4: result=m->slices; break; case 5: result=m->errors; break;
            case 6: result=m->start; break; case 7: result=m->finish; break;
            case 8: result=m->one_cpu; break; case 9: result=m->two_cpu; break;
            case 10: result=m->mixed_cpu; break; case 11: result=m->unsupported; break;
            case 12: result=m->tasks; break; case 13: result=m->curves; break;
            case 14: result=m->identity; break;
            }
        }
    }
    pthread_mutex_unlock(&t->mutex);
    return result;
}
#endif

#ifdef B4_PTHREAD
int b4_pthread_measure(b4_pthread_team *t, unsigned enabled) {
    if (!t || enabled > 1) return -1;
    pthread_mutex_lock(&t->mutex);
    int status = t->active ? -1 : 0;
    if (!status) t->diagnostic_enabled = enabled;
    pthread_mutex_unlock(&t->mutex);
    return status;
}
unsigned b4_pthread_participants(b4_pthread_team *t) {
    if (!t) return 0;
    pthread_mutex_lock(&t->mutex);
    unsigned result = t->active ? 0 : t->participants;
    pthread_mutex_unlock(&t->mutex);
    return result;
}
#else
unsigned long long b4_pthread_diagnostic(b4_pthread_team *t,
                                         unsigned r, unsigned f) {
    (void)t; (void)r; (void)f; return ~0ull;
}
int b4_pthread_measure(b4_pthread_team *t, unsigned e) {
    (void)t; (void)e; return -2;
}
unsigned b4_pthread_participants(b4_pthread_team *t) { (void)t; return 0; }
#endif
unsigned b4_team_max_workers(void) {
#if defined(B4_PTHREAD) || defined(B4_OPENMP)
    return B4_NATIVE_MAX_WORKERS;
#else
    return 0;
#endif
}
