"""Timeline capture must preserve the requested worker CPU/heap observer."""

import json
import os
import tempfile
import unittest
from multiprocessing import queues, util
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.qs.a7a.worker_timeline import initialize_worker
from v2.qs import worker_packet


class TimelineCompositionTests(unittest.TestCase):
    def test_requested_diagnostics_are_initialized_before_queue_observation(
        self,
    ):
        for diagnostic in (False, True):
            with tempfile.TemporaryDirectory() as temporary:
                calls = []
                with (
                    patch.dict(
                        os.environ,
                        {
                            "FACTOR_A7A_TIMELINE_DIR": temporary,
                            "FACTOR_A7A_DIAGNOSTIC_DIR": "capture"
                            if diagnostic
                            else "",
                        },
                    ),
                    patch(
                        "v2.qs.parallel._initialize_process",
                        side_effect=lambda state: calls.append("base"),
                    ),
                    patch(
                        "v2.benchmarks.qs.a7a.diagnostic_hook.initialize_worker",
                        side_effect=lambda state: calls.append("diagnostics"),
                    ),
                    patch.object(worker_packet, "collect_packet"),
                    patch.object(queues.SimpleQueue, "put"),
                    patch.object(queues.Queue, "get"),
                    patch.object(util, "Finalize"),
                ):
                    initialize_worker({})
                self.assertEqual(
                    calls, ["diagnostics" if diagnostic else "base"]
                )
                record = json.loads(
                    next(Path(temporary).glob("worker-*.json")).read_text()
                )
                self.assertIsNone(record["error"])
                self.assertEqual(record["spans"], {})
