"""Fresh-job diagnostic of parent CPU guarding, never a production policy."""

import time


def execute_packet(task, state=None):
    """Keep packet-private limits and CPU publication; parent guards the sum.

    This benchmark-only adapter reuses the unchanged arithmetic packet. The
    parent still reserves private grants, polls aggregate CPU while waiting
    and validating, signals stop, charges final snapshots and drains workers.
    Global stop latency differs from peer polling; no production checkpoint
    may be saved or resumed with this unversioned diagnostic adapter.
    """
    from ....qs import parallel
    from ....qs.worker_packet import collect_packet

    process = state is None
    slot = None
    if process:
        state, slot = parallel._PROCESS_STATE
    last_poll = 0.0

    def publish(budget):
        nonlocal last_poll
        now = time.monotonic()
        if process and now - last_poll >= 0.001:
            state["cpu"][slot] = time.process_time()
            state["rss"][slot] = parallel._rss()
            last_poll = now

    try:
        # collect_packet retains immediate work refusal, private CPU/wall
        # limits, <=64-action stop checks and forced publication checks.
        return collect_packet(
            task,
            cancelled=state["stop"].is_set,
            publish_cpu=publish,
        )
    finally:
        if process:
            state["cpu"][slot] = time.process_time()
            state["rss"][slot] = parallel._rss()


def main():
    from ....qs import current_parallel
    from .cell import main as cell_main

    # Only this fresh diagnostic subprocess selects the adapter. Ordinary
    # imports and production entry points retain the incumbent peer guard.
    current_parallel.execute_packet = execute_packet
    cell_main()


if __name__ == "__main__":
    main()
