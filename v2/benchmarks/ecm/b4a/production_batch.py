"""Complete portfolio calls through retained production ECM owners."""

import os
import time


class ProductionSession:
    """Retain one declared recipe and settle final owner cleanup explicitly."""

    def __init__(self):
        self.owner = None
        self.key = None
        self.closed = False

    @staticmethod
    def make_owner(request):
        from ....ecm.workers import ECMExecutor
        from ....execution.worker_contract import WorkerSettings

        role = request["role"]
        native = (
            {}
            if request.get("metadata_cache") is None
            else dict(metadata_cache=request["metadata_cache"])
        )
        if role in ("python", "b12", "scalar"):
            values = dict(request.get("worker_settings", {}))
            values["workers"] = request.get("workers", 1)
            return ECMExecutor(
                seed=request.get("seed", 73),
                settings=WorkerSettings(**values),
                force_workers=role == "b12",
                **native,
            )
        if role == "pthread":
            from ....execution.ecm_native import NativeECMExecutor

            owner_type = NativeECMExecutor
        elif role == "openmp":
            from ....execution.ecm_openmp import OpenMPExecutor

            owner_type = OpenMPExecutor
        else:
            raise ValueError("unknown production ECM role")
        values = dict(request.get("worker_settings", {}))
        values["workers"] = request.get("workers", 4)
        return owner_type(
            seed=request.get("seed", 73),
            settings=WorkerSettings(**values),
            schedule=request.get("schedule", "dynamic"),
            chunk_size=request.get("chunk_size", 1),
            capacity=request.get("capacity", 32),
            **native,
        )

    def state(self):
        owner = self.owner
        team = getattr(owner, "team", None)
        if team is not None:
            packets = team.packets
            return dict(
                contexts=len(packets),
                live_contexts=sum(p.context is not None for p in packets),
                actual_workers=0
                if team.actual is None
                else int(team.actual[0]),
                slices=team.slices,
                returns=getattr(team, "foreign_returns", team.slices),
                return_reasons=dict(getattr(team, "return_reasons", {})),
                native_segment=self.segment_state(team),
                retained_bytes=team.blocks.used,
                blocks=len(team.blocks.blocks),
                backend=team.variant,
                metadata=self.metadata_state(team.blocks),
                participants=(
                    int(team.lib.b4_omp_participants(team.owner))
                    if team.variant == "openmp" and team.owner is not None
                    else None
                ),
            )
        engine = getattr(owner, "serial_engine", None)
        source = None if engine is None else engine.programs or engine.context
        packet = getattr(source, "native_packet", None)
        store = getattr(packet, "blocks", None)
        return dict(
            contexts=int(packet is not None),
            live_contexts=int(
                packet is not None and packet.context is not None
            ),
            returns=getattr(packet, "calls", 0),
            return_reasons=dict(getattr(packet, "return_reasons", {})),
            native_segment=self.segment_state(packet),
            retained_bytes=getattr(store, "used", 0),
            blocks=len(getattr(store, "blocks", {})),
            metadata=self.metadata_state(store),
            pool_started=bool(
                getattr(owner, "slots", ()) and owner.pool_started()
            ),
            processes=len(getattr(owner, "slots", ())),
        )

    @staticmethod
    def metadata_state(store):
        metadata = getattr(store, "metadata", None)
        index = getattr(store, "index", None)
        return (
            None
            if metadata is None
            else dict(
                stats=dict(metadata.stats),
                indexed_native_blocks=0
                if index is None
                else int(store.lib.b4_index_stat(index, 3)),
                indexed_owned_bytes=0
                if index is None
                else int(store.lib.b4_index_stat(index, 1)),
                indexed_metadata_work=0
                if index is None
                else int(store.lib.b4_index_stat(index, 2)),
                catalogue_bytes=metadata.used,
                cache_directory=None
                if metadata.cache_directory is None
                else str(metadata.cache_directory),
            )
        )

    @staticmethod
    def segment_state(owner):
        segment = getattr(owner, "_segment", None)
        return (
            None
            if segment is None
            else dict(
                returns=dict(segment.returns),
                waves=segment.waves,
                refills=segment.refills,
                observations=segment.observations,
            )
        )

    def call(self, request):
        from ....execution.budget import Budget
        from ....execution.portfolio_config import PortfolioConfig
        from ....portfolio import factorize_bounded

        if self.closed:
            raise RuntimeError("production session is closed")
        key = tuple(
            request.get(name)
            for name in (
                "role",
                "workers",
                "schedule",
                "chunk_size",
                "capacity",
                "seed",
                "metadata_cache",
            )
        )
        key += (tuple(sorted(request.get("worker_settings", {}).items())),)
        if self.owner is None:
            self.owner, self.key = self.make_owner(request), key
        elif self.key != key:
            raise ValueError("retained production recipe changed")
        role = request["role"]
        config = PortfolioConfig(
            trial_bound=request.get("trial_bound", 2),
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=((request["b1"], request["b2"], request["curves"]),),
            ecm_chain_mode="auto",
            ecm_kernel="python" if role in ("python", "b12") else "c-gmp-full",
            siqs=None,
            max_input_bits=1024,
        )
        grant = Budget(
            seconds=request.get("seconds", 30),
            cpu_seconds=request.get("cpu_seconds", 120),
            rss_bytes=8 * 2**30,
        )
        started = time.monotonic()
        before = self.state()
        try:
            run = factorize_bounded(
                request["n"],
                seed=request.get("seed", 73),
                config=config,
                budget=grant,
                ecm_executor=self.owner,
                qs_workers=None,
            )
            if run.result.reconstruct() != request["n"]:
                raise ValueError("production result failed reconstruction")
            events = [e for e in run.events if e["stage"] == "ecm_workers"]
            stats = [e["stats"] for e in events]
            state = self.state()
            dispatched = sum(s.get("worker_dispatches", 0) for s in stats)
            completed = sum(
                s.get("completed", s.get("completed_after", 0)) for s in stats
            )
            if role == "b12" and events and not dispatched:
                raise ValueError(
                    "B12 control failed to dispatch actual workers"
                )
            if role in ("pthread", "openmp") and events:
                workers = request.get("workers", 4)
                if state["actual_workers"] != workers:
                    raise ValueError("production compute-team mismatch")
                if (
                    role == "openmp"
                    and state["participants"] != (1 << workers) - 1
                ):
                    raise ValueError(
                        "OpenMP caller/team participation mismatch"
                    )
            result = dict(
                request,
                reason=run.reason,
                complete=run.result.complete,
                factors=[
                    [int(p.value), p.exponent, p.certainty.value]
                    for p in run.result.factors
                ],
                remaining=[int(n) for n in run.result.remaining],
                work=run.work_used,
                observed_rss=grant.observed_rss,
                ecm_events=events,
                completed_curves=completed,
                state=state,
                before_state=before,
                worker_dispatches=dispatched,
                manifests=[s.get("manifest") for s in stats],
                resource_stop=run.resource_stop,
            )
        finally:
            grant.close()
        result.update(api_wall=time.monotonic() - started, cpu=grant.cpu_used)
        return result

    @staticmethod
    def cpu_marker():
        """Read quiescent lifetime counters already converted by H3."""
        from ....execution.resources import MacResources

        samples = MacResources().sample({os.getpid(): None})
        return dict(
            cpu=sum(s.cpu + s.reaped_cpu for s in samples),
            rss=sum(s.rss for s in samples),
            pids=[s.pid for s in samples],
        )

    def close(self):
        if self.closed:
            return
        if self.owner is not None:
            self.owner.close()
        self.closed = True
