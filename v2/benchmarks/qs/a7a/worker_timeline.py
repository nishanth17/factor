"""Optional worker service/queue spans, aggregated once per minute."""

import json
import os
import time
from collections import defaultdict
from functools import wraps
from pathlib import Path


def initialize_worker(state):
    """Run the production initializer before installing benchmark observers."""
    from multiprocessing import queues, util

    from ....qs import worker_packet
    from ....qs.parallel import _initialize_process

    _initialize_process(state)
    folder = Path(os.environ["FACTOR_A7A_TIMELINE_DIR"])
    path = folder / f"worker-{os.getpid()}.json"
    totals = defaultdict(lambda: [0, 0.0, 0.0])
    last = [time.monotonic() - 60]
    failure = [None]

    def save(force=False):
        now = time.monotonic()
        if force or now - last[0] >= 60:
            try:
                temporary = path.with_suffix(".tmp")
                temporary.write_text(
                    json.dumps(
                        dict(
                            pid=os.getpid(),
                            monotonic=now,
                            spans=dict(totals),
                            error=failure[0],
                        )
                    )
                )
                temporary.replace(path)
                last[0] = now
            except Exception as error:
                failure[0] = repr(error)

    def wrap(function, label):
        @wraps(function)
        def measured(*args, **kwargs):
            wall, cpu = time.monotonic(), time.process_time()
            try:
                return function(*args, **kwargs)
            finally:
                row = totals[label]
                row[0] += 1
                row[1] += time.monotonic() - wall
                row[2] += time.process_time() - cpu
                save()

        return measured

    worker_packet.collect_packet = wrap(worker_packet.collect_packet, "active")
    queues.SimpleQueue.put = wrap(queues.SimpleQueue.put, "publication")
    queues.Queue.get = wrap(queues.Queue.get, "waiting_for_task")
    util.Finalize(None, save, args=(True,), exitpriority=1)
    save(True)
