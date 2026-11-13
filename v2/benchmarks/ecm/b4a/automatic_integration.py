"""Source-pinned automatic ECM-to-QS comparison entry point.

Only benchmark adapters come from adapter_root. Production imports must come
from the declared committed source_root, including a pre-activation control.
"""

import json
import sys
from pathlib import Path


def main():
    request = json.loads(sys.stdin.readline())
    root = Path(request["source_root"]).resolve()
    sys.path.insert(0, str(root))
    import v2

    adapters = str(Path(request["adapter_root"]).resolve() / "v2")
    if adapters not in v2.__path__:
        v2.__path__.append(adapters)
    from v2.benchmarks.ecm.b4a.confirmation_runtime import close_spawn_helper
    from v2.benchmarks.ecm.b4a.integration_call import run_request

    try:
        row = run_request(request)
        modules = (
            "v2.portfolio",
            "v2.execution.resources",
            "v2.ecm.auto",
            "v2.ecm.workers",
            "v2.qs.current_parallel",
        )
        paths = {}
        for name in modules:
            module = sys.modules.get(name)
            if module is not None:
                path = Path(module.__file__).resolve()
                if not path.is_relative_to(root):
                    raise ValueError(
                        "integration production source escaped its pin"
                    )
                paths[name] = str(path)
        row["production_modules"] = paths
        row["source"] = request.get("source")
        if row["ecm_execution"] is None:
            jobs = []
            snapshots = (row.get("observation") or {}).get("snapshots", ())
            counts = [
                len(s["native"].get("ecm_worker_pids", ())) for s in snapshots
            ]
            processes = max(counts, default=0)
            for event in row["events"]:
                if event["stage"] in ("ecm", "ecm_workers"):
                    jobs.append(
                        dict(
                            backend="python",
                            implementation="python",
                            topology="processes" if processes else "scalar",
                            workers=processes or 1,
                            native_entries=0,
                            event=event,
                        )
                    )
            row["ecm_execution"] = dict(
                schema="historical-python-events-v1",
                requested="omitted",
                backend="python" if jobs else "not-engaged",
                jobs=jobs,
                native_entries=0,
                reason="historical_committed_default",
            )
        print(json.dumps(row), flush=True)
        return 1 if row.get("native_assertion_failure") else 0
    finally:
        close_spawn_helper()


if __name__ == "__main__":
    raise SystemExit(main())
