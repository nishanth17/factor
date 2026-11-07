"""Bounded progress snapshots; observation failures never cancel arithmetic."""

import json
import resource
import threading
import time
from collections import deque
from pathlib import Path

from .runner import thermal_status


class Progress:
    """One durable row per minute plus the final observation."""

    def __init__(self, path, counters, observer, *, interval=60):
        if interval != 60:
            raise ValueError("progress interval must be 60 seconds")
        self.path = Path(path)
        self.counters, self.observer = counters, observer
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self._run, daemon=True)
        self.history = deque(maxlen=32)
        self.rows = 0
        self.error = None
        self.job = None

    def sample(self):
        now = time.monotonic()
        row = self.counters.snapshot()
        if row.get("job") != self.job:
            self.history.clear()
            self.job = row.get("job")
        useful = row.get("usable_relations")
        if useful is not None:
            self.history.append((now, useful))
            rates = {}
            for seconds in (300, 900):
                points = [p for p in self.history if now - p[0] <= seconds]
                if len(points) > 1 and useful >= points[0][1]:
                    rates[str(seconds)] = (useful - points[0][1]) / (
                        now - points[0][0]
                    )
            row["rolling_relations_per_second"] = rates
        row.update(
            monotonic=now,
            cpu=resource.getrusage(resource.RUSAGE_SELF).ru_utime
            + resource.getrusage(resource.RUSAGE_SELF).ru_stime,
            rss=None
            if not self.observer or not self.observer.rows
            else self.observer.rows[-1],
            thermal=thermal_status(),
        )
        if self.rows >= 1000:
            raise RuntimeError("progress observation allowance exhausted")
        with self.path.open("a") as stream:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
        self.rows += 1

    def _run(self):
        while not self.stop.is_set():
            try:
                self.sample()
            except Exception as error:
                self.error = repr(error)
                return
            self.stop.wait(60)

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *args):
        self.stop.set()
        self.thread.join(timeout=5)
        if self.thread.is_alive():
            self.error = "progress recorder did not stop"
        elif self.error is None:
            try:
                self.sample()
            except Exception as error:
                self.error = repr(error)
