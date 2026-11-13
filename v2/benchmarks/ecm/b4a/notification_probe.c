/* A pre-registered entry observer must progress without paused workers. */
#define _POSIX_C_SOURCE 200809L
#define _DARWIN_C_SOURCE 1
#include <assert.h>
#include <pthread.h>
#include <sched.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include B4_PACKET_SOURCE
#include B4_TEAM_SOURCE

static atomic_uint observed;

static void *wait_for_entry(void *team) {
    b4_pthread_wait(team, 2);
    atomic_store(&observed, 1);
    return NULL;
}

static uint64_t now_ns(void) {
    struct timespec timestamp;
    assert(!clock_gettime(CLOCK_MONOTONIC, &timestamp));
    return (uint64_t)timestamp.tv_sec * 1000000000u + timestamp.tv_nsec;
}

int main(void) {
    unsigned char modulus[] = {3, 245}, x[] = {2}, z[] = {1}, a[] = {2};
    unsigned char scalar[512];
    b4_context *contexts[8];
    int statuses[8];
    unsigned actual = 0;
    memset(scalar, 255, sizeof(scalar));
    for (unsigned i = 0; i < 8; ++i) {
        contexts[i] = b4_create();
        assert(contexts[i]);
        assert(!b4_ladder(contexts[i], 1, (b4_value){modulus, 2},
            (b4_value){x, 1}, (b4_value){z, 1}, (b4_value){a, 1},
            (b4_value){scalar, sizeof(scalar)}));
    }
    b4_pthread_team *team = b4_pthread_create();
    assert(team);
    assert(!b4_pthread_step(team, contexts, 8, 4, 0, 1, 5, 16,
                           statuses, &actual));
    assert(actual == 4);

    /* Reset the entry predicate after the first generation. Inspecting the
     * private counter under its mutex makes registration deterministic. */
    assert(!b4_pthread_inject(team, -1, -1, 0));
    atomic_init(&observed, 0);
    pthread_t waiter;
    assert(!pthread_create(&waiter, NULL, wait_for_entry, team));
    uint64_t end = now_ns() + 1000000000u;
    for (;;) {
        pthread_mutex_lock(&team->mutex);
        unsigned active = team->observers;
        pthread_mutex_unlock(&team->mutex);
        if (active) break;
        assert(now_ns() < end);
        sched_yield();
    }
    assert(!b4_pthread_step(team, contexts, 8, 4, 0, 1, 5, 16,
                           statuses, &actual));
    assert(!pthread_join(waiter, NULL));
    assert(atomic_load(&observed));

    assert(!b4_pthread_shutdown(team));
    b4_pthread_destroy(team);
    for (unsigned i = 0; i < 8; ++i) b4_destroy(contexts[i]);
    puts("pthread entry observer progressed");
    return 0;
}
