/* Deterministic libproc failure proof; compile only in an admitted QA window. */
#define proc_pid_rusage b4_fake_pid_rusage
#define proc_listchildpids b4_fake_child_pids
#ifndef B4_PACKET_SOURCE
#define B4_PACKET_SOURCE "../../../ecm/native/packet.c"
#endif
#include B4_PACKET_SOURCE
#include <assert.h>

static int child_mode, usage_error, child_calls;

int b4_fake_pid_rusage(int pid, int flavor, rusage_info_t *buffer) {
    assert(pid == 42 && flavor == RUSAGE_INFO_V0);
    if (usage_error) { errno = usage_error; return -1; }
    struct rusage_info_v0 *usage = (struct rusage_info_v0 *)buffer;
    memset(usage, 0, sizeof(*usage));
    usage->ri_proc_start_abstime = 7;
    usage->ri_resident_size = 900;
    usage->ri_user_time = 100;
    usage->ri_system_time = 200;
    return 0;
}

int b4_fake_child_pids(pid_t pid, void *buffer, int size) {
    assert(pid == 42 && buffer && size == sizeof(pid_t));
    ++child_calls;
    if (child_mode == 1) { *(pid_t *)buffer = 123; return sizeof(pid_t); }
    if (child_mode == 2) { errno = EPERM; return 0; }
    if (child_mode == 3) { errno = EACCES; return -1; }
    if (child_mode == 4) { errno = EIO; return 0; }
    if (child_mode == 5) return -1;
    if (child_mode == 6) { errno = ESRCH; return 0; }
    return 0;
}

int main(void) {
    b4_root_observation output;
    const int expected[] = {0, 1, -EPERM, -EACCES, -EIO, -EIO, -ESRCH};
    for (child_mode = 0; child_mode < 7; ++child_mode) {
        child_calls = 0;
        errno = EPERM; /* A successful childless query must ignore stale errno. */
        assert(b4_observe_root(42, &output) == expected[child_mode]);
        assert(child_calls == 1);
        if (!child_mode) {
            assert(output.birth == 7 && output.rss == 900);
            assert(output.user == 100 && output.system == 200);
        }
    }
    for (unsigned i = 0; i < 3; ++i) {
        const int errors[] = {EPERM, EACCES, ESRCH};
        usage_error = errors[i]; child_calls = 0;
        assert(b4_observe_root(42, &output) == -usage_error);
        assert(child_calls == 0);
    }
    assert(b4_observe_root(0, &output) == -EINVAL);
    assert(b4_observe_root(42, NULL) == -EINVAL);
    return 0;
}
