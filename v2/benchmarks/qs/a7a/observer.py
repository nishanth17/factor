"""Optional two-Hz RSS samples of this call's owned processes only."""

import os
import subprocess
import threading
import time


class OwnedRSS:
    """Bounded outer samples; child HWM sums are never called aggregate RSS."""

    def __init__(self, pool, *, limit, samples=16000):
        self.pool, self.limit, self.maximum = pool, limit, samples
        self.rows, self.error, self.exceeded = [], None, False
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        while not self.stop.is_set():
            processes = []
            if self.pool is not None and self.pool.mode == "process":
                processes = list(
                    (self.pool.executor._processes or {}).values()
                )
            pids = [os.getpid()] + [p.pid for p in processes if p.pid]
            try:
                from multiprocessing import resource_tracker

                helper = resource_tracker._resource_tracker._pid
                if helper is not None:
                    pids.append(helper)
                result = subprocess.run(
                    ["ps", "-o", "pid=,rss=", "-p", ",".join(map(str, pids))],
                    capture_output=True,
                    text=True,
                    timeout=1,
                )
                values = [line.split() for line in result.stdout.splitlines()]
                rss = {int(pid): 1024 * int(value) for pid, value in values}
                total = sum(rss.values())
                self.rows.append(
                    dict(
                        time=time.monotonic(),
                        bytes=total,
                        processes=rss,
                        helper=helper,
                    )
                )
                if total > self.limit or len(self.rows) >= self.maximum:
                    self.exceeded = True
                    return
            except Exception as error:
                self.error = str(error)
                return
            self.stop.wait(0.5)

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *args):
        self.stop.set()
        self.thread.join(timeout=2)
